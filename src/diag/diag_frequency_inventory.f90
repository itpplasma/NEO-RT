module diag_frequency_inventory
    !! Surface frequency and resonance inventory for the configured NEO-RT model.
    !! Model 0 is evaluated through native Boozer canonical frequencies; model
    !! 2 is evaluated through the direct real-space return map.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_is_finite
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_freq, only: Om_th, Om_ph, Om_tB
    use neort_profiles, only: vth, Om_tE, mi, qi
    use driftorbit, only: etatp, sign_vpar, frequency_model, mph, &
        FREQUENCY_MODEL_LEGACY, FREQUENCY_MODEL_GC_FULL
    use do_magfie_mod, only: s, bfac, q
    use neort, only: mth_max_abs, harmonic_bounds, set_to_passing_region, &
        set_to_trapped_region
    use neort_gc_orbit_integrator, only: GC_ORBIT_PASSING, GC_ORBIT_TRAPPED
    use neort_gc_frequency_provider, only: &
        GC_FREQUENCY_SUCCESS, gc_frequency_context_t, &
        gc_full_orbit_frequency_result_t, initialize_gc_frequency_context, &
        evaluate_gc_full_orbit_frequency
    use neort_gc_full_resonance, only: find_gc_resonances
    use util, only: files_exist
    implicit none
    private

    integer, parameter :: inventory_points = 128
    integer, parameter :: resonance_scan_intervals = 256
    integer, parameter :: max_inventory_roots = 64
    real(dp), parameter :: residual_tolerance = 1.0e-10_dp
    real(dp), parameter :: eta_tolerance = 1.0e-9_dp

    public :: run_frequency_inventory_diag
    public :: frequency_inventory_model_label
    public :: frequency_inventory_class_label

contains

    subroutine run_frequency_inventory_diag(runname)
        character(len=*), intent(in) :: runname
        character(len=512) :: config_file
        integer :: unit_frequency, unit_resonance, model, setup_status

        config_file = trim(runname)//'.in'
        call neort_init(config_file, 'in_file', 'in_file_pert')
        model = frequency_model
        if (model /= FREQUENCY_MODEL_LEGACY .and. &
            model /= FREQUENCY_MODEL_GC_FULL) then
            error stop 'frequency_inventory supports frequency_model=0 or 2'
        end if
        if (files_exist('plasma.in', 'profile.in')) then
            call neort_prepare_splines('plasma.in', 'profile.in')
        end if
        call neort_setup_at_s(s)

        open (newunit=unit_frequency, file=trim(runname)//'_frequency_inventory.out', &
            status='replace', action='write', form='formatted')
        open (newunit=unit_resonance, file=trim(runname)//'_resonance_inventory.out', &
            status='replace', action='write', form='formatted')
        write (unit_frequency, '(A)') '# frequency_inventory_v1'
        write (unit_frequency, '(A)') &
            '# columns model class sign eta etatp eta_over_etatp omega_b omega_phi omega_magnetic omega_electric status orbit_status'
        write (unit_frequency, '(A)') &
            '# real2 component columns omega_magnetic and omega_electric are NaN: direct provider returns canonical frequencies only'
        write (unit_frequency, '(A)') '# status=0 success; nonzero rows retain invalid frequencies'
        write (unit_resonance, '(A)') '# resonance_inventory_v1'
        write (unit_resonance, '(A)') &
            '# columns model class sign mph mth eta etatp eta_over_etatp residual dres_deta status orbit_status'
        write (unit_resonance, '(A)') '# status=0 success; 3 partial scan; nonzero rows retain status'

        setup_status = GC_FREQUENCY_SUCCESS
        if (model == FREQUENCY_MODEL_GC_FULL) then
            call write_full_inventory(s, unit_frequency, unit_resonance, setup_status)
        else
            call write_legacy_inventory(unit_frequency, unit_resonance)
        end if
        close (unit_frequency)
        close (unit_resonance)
        if (setup_status /= GC_FREQUENCY_SUCCESS) then
            write (0, '(A,I0)') 'frequency_inventory full-orbit setup status: ', setup_status
        end if
    end subroutine run_frequency_inventory_diag

    subroutine write_legacy_inventory(unit_frequency, unit_resonance)
        integer, intent(in) :: unit_frequency, unit_resonance
        integer :: sign_index, orbit_class, k, mth_min, mth_max, status
        real(dp) :: eta_lo, eta_hi, eta, ratio, omega_b, omega_phi
        real(dp) :: omega_magnetic, omega_electric
        real(dp) :: roots(max_inventory_roots), derivatives(max_inventory_roots)

        call harmonic_bounds(mph, q, mth_max_abs, mth_min, mth_max)
        do orbit_class = GC_ORBIT_TRAPPED, GC_ORBIT_PASSING
            call region_bounds(orbit_class, eta_lo, eta_hi)
            do sign_index = -1, 1, 2
                sign_vpar = real(sign_index, dp)
                do k = 0, inventory_points
                    eta = eta_lo + (eta_hi - eta_lo)*real(k, dp) &
                        /real(inventory_points, dp)
                    call legacy_frequency(eta, orbit_class, omega_b, omega_phi, &
                        omega_magnetic, omega_electric, status)
                    ratio = eta/etatp
                    call write_frequency_row(unit_frequency, FREQUENCY_MODEL_LEGACY, &
                        orbit_class, sign_index, eta, etatp, ratio, omega_b, &
                        omega_phi, omega_magnetic, omega_electric, status, 0)
                end do
                do k = mth_min, mth_max
                    call find_inventory_roots(eta_lo, eta_hi, orbit_class, sign_index, &
                        k, roots, derivatives, status)
                    call write_roots(unit_resonance, FREQUENCY_MODEL_LEGACY, &
                        orbit_class, sign_index, k, roots, derivatives, status)
                end do
            end do
        end do
    end subroutine write_legacy_inventory

    subroutine write_full_inventory(surface, unit_frequency, unit_resonance, setup_status)
        real(dp), intent(in) :: surface
        integer, intent(in) :: unit_frequency, unit_resonance
        integer, intent(out) :: setup_status
        type(gc_frequency_context_t) :: context
        type(gc_full_orbit_frequency_result_t) :: result
        integer :: sign_index, orbit_class, k, mth_min, mth_max, status
        real(dp) :: eta_lo, eta_hi, eta, ratio
        real(dp) :: roots(max_inventory_roots), derivatives(max_inventory_roots)

        call initialize_gc_frequency_context(surface, 0.0_dp, bfac, Om_tE, mi, qi, &
            vth, context, setup_status)
        if (setup_status /= GC_FREQUENCY_SUCCESS) then
            call harmonic_bounds(mph, q, mth_max_abs, mth_min, mth_max)
            do orbit_class = GC_ORBIT_TRAPPED, GC_ORBIT_PASSING
                call region_bounds(orbit_class, eta_lo, eta_hi)
                do sign_index = -1, 1, 2
                    do k = 0, inventory_points
                        eta = eta_lo + (eta_hi - eta_lo)*real(k, dp) &
                            /real(inventory_points, dp)
                        call write_frequency_row(unit_frequency, FREQUENCY_MODEL_GC_FULL, &
                            orbit_class, sign_index, eta, etatp, eta/etatp, 0.0_dp, &
                            nan_value(), nan_value(), nan_value(), setup_status, 0)
                    end do
                    do k = mth_min, mth_max
                        call write_resonance_row(unit_resonance, &
                            FREQUENCY_MODEL_GC_FULL, orbit_class, sign_index, k, &
                            nan_value(), nan_value(), nan_value(), nan_value(), &
                            setup_status, 0)
                    end do
                end do
            end do
            return
        end if

        call harmonic_bounds(mph, q, mth_max_abs, mth_min, mth_max)
        do orbit_class = GC_ORBIT_TRAPPED, GC_ORBIT_PASSING
            call region_bounds(orbit_class, eta_lo, eta_hi)
            do sign_index = -1, 1, 2
                sign_vpar = real(sign_index, dp)
                do k = 0, inventory_points
                    eta = eta_lo + (eta_hi - eta_lo)*real(k, dp) &
                        /real(inventory_points, dp)
                    call full_frequency(context, eta, orbit_class, sign_index, &
                        result, status)
                    ratio = eta/etatp
                    call write_frequency_row(unit_frequency, FREQUENCY_MODEL_GC_FULL, &
                        orbit_class, sign_index, eta, etatp, ratio, result%omega_b, &
                        result%omega_phi, nan_value(), nan_value(), status, &
                        result%orbit_status)
                end do
                do k = mth_min, mth_max
                    call find_inventory_roots(eta_lo, eta_hi, orbit_class, sign_index, &
                        k, roots, derivatives, status, context)
                    call write_roots(unit_resonance, FREQUENCY_MODEL_GC_FULL, &
                        orbit_class, sign_index, k, roots, derivatives, status, context)
                end do
            end do
        end do
    end subroutine write_full_inventory

    subroutine region_bounds(orbit_class, eta_lo, eta_hi)
        integer, intent(in) :: orbit_class
        real(dp), intent(out) :: eta_lo, eta_hi

        if (orbit_class == GC_ORBIT_PASSING) then
            call set_to_passing_region(eta_lo, eta_hi)
        else
            call set_to_trapped_region(eta_lo, eta_hi)
        end if
    end subroutine region_bounds

    subroutine legacy_frequency(eta, orbit_class, omega_b, omega_phi, &
            omega_magnetic, omega_electric, status)
        real(dp), intent(in) :: eta
        integer, intent(in) :: orbit_class
        real(dp), intent(out) :: omega_b, omega_phi, omega_magnetic, omega_electric
        integer, intent(out) :: status
        real(dp) :: d1, d2

        status = 0
        omega_b = 0.0_dp
        omega_phi = 0.0_dp
        omega_magnetic = 0.0_dp
        omega_electric = Om_tE
        if (eta <= 0.0_dp) then
            status = 1
            return
        end if
        call Om_th(vth, eta, omega_b, d1, d2)
        call Om_ph(vth, eta, omega_phi, d1, d2)
        call Om_tB(vth, eta, omega_magnetic, d1, d2)
        if (.not. ieee_is_finite(omega_b) .or. .not. ieee_is_finite(omega_phi)) then
            status = 1
        end if
        if (orbit_class /= merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, eta > etatp)) then
            status = 1
        end if
    end subroutine legacy_frequency

    subroutine full_frequency(context, eta, orbit_class, sign_index, result, status)
        type(gc_frequency_context_t), intent(in) :: context
        real(dp), intent(in) :: eta
        integer, intent(in) :: orbit_class, sign_index
        type(gc_full_orbit_frequency_result_t), intent(out) :: result
        integer, intent(out) :: status
        real(dp) :: period_estimate, omega_b, d1, d2

        period_estimate = 1.0_dp
        call Om_th(vth, eta, omega_b, d1, d2)
        if (ieee_is_finite(omega_b) .and. abs(omega_b) > tiny(omega_b)) then
            period_estimate = 2.0_dp*acos(-1.0_dp)/abs(omega_b)
        end if
        call evaluate_gc_full_orbit_frequency(context, eta, sign_index, orbit_class, &
            period_estimate, result, status, vth)
    end subroutine full_frequency

    subroutine write_frequency_row(unit, model, orbit_class, sign_index, eta, etatp_value, &
            ratio, omega_b, omega_phi, omega_magnetic, omega_electric, status, orbit_status)
        integer, intent(in) :: unit, model, orbit_class, sign_index, status, orbit_status
        real(dp), intent(in) :: eta, etatp_value, ratio, omega_b, omega_phi
        real(dp), intent(in) :: omega_magnetic, omega_electric

        write (unit, '(I0,1X,A,1X,I0,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,I0,1X,I0)') &
            model, frequency_inventory_class_label(orbit_class), sign_index, eta, &
            etatp_value, ratio, omega_b, omega_phi, omega_magnetic, omega_electric, &
            status, orbit_status
    end subroutine write_frequency_row

    subroutine write_roots(unit, model, orbit_class, sign_index, mth_value, roots, &
            derivatives, search_status, context)
        integer, intent(in) :: unit, model, orbit_class, sign_index, mth_value
        real(dp), intent(in) :: roots(:), derivatives(:)
        integer, intent(in) :: search_status
        type(gc_frequency_context_t), intent(in), optional :: context
        integer :: k, nroots, status, orbit_status, row_status
        real(dp) :: residual, omega_b, omega_phi, d1, d2
        type(gc_full_orbit_frequency_result_t) :: result

        nroots = count(roots /= 0.0_dp)
        if (nroots == 0) then
            if (search_status /= 0) then
                call write_resonance_row(unit, model, orbit_class, sign_index, &
                    mth_value, nan_value(), nan_value(), nan_value(), nan_value(), &
                    search_status, 0)
            end if
            return
        end if
        do k = 1, nroots
            orbit_status = 0
            row_status = search_status
            if (model == FREQUENCY_MODEL_GC_FULL) then
                call full_frequency(context, roots(k), orbit_class, sign_index, result, status)
                orbit_status = result%orbit_status
                if (status == GC_FREQUENCY_SUCCESS) then
                    residual = real(mth_value, dp)*result%omega_b &
                        + real(mph, dp)*result%omega_phi
                else
                    residual = nan_value()
                end if
            else
                call legacy_frequency(roots(k), orbit_class, omega_b, omega_phi, d1, d2, status)
                if (status == 0) then
                    residual = real(mth_value, dp)*omega_b &
                        + real(mph, dp)*omega_phi
                else
                    residual = nan_value()
                end if
            end if
            if (status /= 0) row_status = status
            call write_resonance_row(unit, model, orbit_class, sign_index, mth_value, &
                roots(k), etatp, roots(k)/etatp, residual, row_status, orbit_status, &
                derivatives(k))
        end do
    end subroutine write_roots

    subroutine write_resonance_row(unit, model, orbit_class, sign_index, mth_value, &
            eta, etatp_value, ratio, residual, status, orbit_status, derivative)
        integer, intent(in) :: unit, model, orbit_class, sign_index, mth_value
        real(dp), intent(in) :: eta, etatp_value, ratio, residual
        integer, intent(in) :: status, orbit_status
        real(dp), intent(in), optional :: derivative
        real(dp) :: derivative_value

        derivative_value = nan_value()
        if (present(derivative)) derivative_value = derivative
        write (unit, '(I0,1X,A,1X,I0,1X,I0,1X,I0,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,ES24.16E3,1X,I0,1X,I0)') &
            model, frequency_inventory_class_label(orbit_class), sign_index, mph, &
            mth_value, eta, etatp_value, ratio, residual, derivative_value, status, &
            orbit_status
    end subroutine write_resonance_row

    subroutine find_inventory_roots(eta_lo, eta_hi, orbit_class, sign_index, mth_value, &
            roots, derivatives, search_status, context)
        real(dp), intent(in) :: eta_lo, eta_hi
        integer, intent(in) :: orbit_class, sign_index, mth_value
        real(dp), intent(out) :: roots(:), derivatives(:)
        integer, intent(out) :: search_status
        type(gc_frequency_context_t), intent(in), optional :: context
        integer :: nroots

        roots = 0.0_dp
        derivatives = 0.0_dp
        call find_gc_resonances(residual_at, eta_lo, eta_hi, &
            resonance_scan_intervals, residual_tolerance, eta_tolerance, roots, &
            derivatives, nroots, search_status)

    contains

        subroutine residual_at(eta_value, residual, status)
            real(dp), intent(in) :: eta_value
            real(dp), intent(out) :: residual
            integer, intent(out) :: status
            real(dp) :: omega_b, omega_phi, d1, d2
            type(gc_full_orbit_frequency_result_t) :: value

            if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
                call full_frequency(context, eta_value, orbit_class, sign_index, value, status)
                omega_b = value%omega_b
                omega_phi = value%omega_phi
            else
                call legacy_frequency(eta_value, orbit_class, omega_b, omega_phi, d1, d2, status)
            end if
            residual = real(mth_value, dp)*omega_b + real(mph, dp)*omega_phi
        end subroutine residual_at
    end subroutine find_inventory_roots

    pure function nan_value() result(value)
        real(dp) :: value
        value = ieee_value(0.0_dp, ieee_quiet_nan)
    end function nan_value

    pure function frequency_inventory_model_label(model) result(label)
        integer, intent(in) :: model
        character(len=8) :: label

        select case (model)
        case (FREQUENCY_MODEL_LEGACY)
            label = 'boozer0'
        case (FREQUENCY_MODEL_GC_FULL)
            label = 'real2'
        case default
            label = 'unknown'
        end select
    end function frequency_inventory_model_label

    pure function frequency_inventory_class_label(orbit_class) result(label)
        integer, intent(in) :: orbit_class
        character(len=8) :: label

        select case (orbit_class)
        case (GC_ORBIT_PASSING)
            label = 'passing'
        case (GC_ORBIT_TRAPPED)
            label = 'trapped'
        case default
            label = 'unknown'
        end select
    end function frequency_inventory_class_label

end module diag_frequency_inventory
