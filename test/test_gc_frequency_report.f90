program test_gc_frequency_report
    !! Write a directly comparable legacy/GC frequency and resonance table.
    !! The report is diagnostic output; the production frequency path remains
    !! selected by frequency_model in the normal NEO-RT executable.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s
    use driftorbit, only: FREQUENCY_MODEL_LEGACY, FREQUENCY_MODEL_GC_THIN, &
        frequency_model, magdrift, magdrift_passing, sign_vpar, &
        etatp, etadt, mth, mph, nlev
    use neort_freq, only: init_canon_freq_trapped_spline, &
        init_canon_freq_passing_spline, Om_th, Om_ph, Om_tB
    use neort_gc_frequency_splines, only: GC_SPLINE_SUCCESS, &
        evaluate_gc_spline
    use neort_gc_orbit_integrator, only: GC_ORBIT_PASSING, GC_ORBIT_TRAPPED
    use neort_magfie, only: init_flux_surface_average
    use neort_orbit, only: th0
    use neort_profiles, only: vth, Om_tE
    use neort_resonance, only: driftorbit_coarse, driftorbit_root
    use util, only: qe, mu, qi, mi

    implicit none

    integer, parameter :: nfrequency = 120
    integer, parameter :: nmodels = 2
    integer, parameter :: model_ids(nmodels) = [ &
        FREQUENCY_MODEL_LEGACY, FREQUENCY_MODEL_GC_THIN]
    real(dp), parameter :: surface = 0.25_dp
    real(dp), parameter :: velocity = 1.0e8_dp
    real(dp), parameter :: mass = 2.014_dp*mu
    real(dp), parameter :: electric_frequency = 4.0e3_dp
    character(len=1024) :: eqdsk_file, output_prefix
    real(dp) :: eta, eta_min, eta_max, eta_fraction
    real(dp) :: omega_b, omega_phi, omega_magnetic
    real(dp) :: omega_electric
    real(dp) :: gc_omega_b, gc_domega_bdv, gc_domega_bdeta
    real(dp) :: gc_omega_magnetic, gc_domega_magnetic_dv
    real(dp) :: gc_domega_magnetic_deta, gc_domega_electric_deta
    integer :: gc_status, gc_orbit_class
    real(dp) :: domega_bdv, domega_bdeta, domega_phdv, domega_phdeta
    real(dp) :: domega_bmagdv, domega_bmagdeta
    real(dp) :: roots(nlev, 3), root_value(2), residual
    integer :: model_index, model, class_index, k, nroots, root_index
    integer :: unit_frequency, unit_resonance
    character(len=16) :: model_name, class_name

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    if (len_trim(eqdsk_file) == 0) error stop 'EQDSK_FILE is required'
    call get_environment_variable('GC_FREQUENCY_REPORT_PREFIX', output_prefix)
    if (len_trim(output_prefix) == 0) output_prefix = 'gc_frequency'

    open (newunit=unit_frequency, file=trim(output_prefix)//'_scan.out', &
        status='replace', action='write')
    open (newunit=unit_resonance, file=trim(output_prefix)//'_resonances.out', &
        status='replace', action='write')
    write (unit_frequency, '(a)') &
        '# model class eta_over_etatp omega_b omega_phi omega_magnetic omega_electric'
    write (unit_resonance, '(a)') &
        '# model class mth eta_over_etatp resonance_residual dres_deta'

    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))
    call set_s(surface)
    call init_magfie_at_s()
    call init_flux_surface_average(surface)
    qi = qe
    mi = mass
    vth = velocity
    Om_tE = electric_frequency
    sign_vpar = 1.0_dp
    mph = 3
    magdrift = .true.
    magdrift_passing = 1

    do model_index = 1, nmodels
        model = model_ids(model_index)
        frequency_model = model
        call init_canon_freq_trapped_spline()
        call init_canon_freq_passing_spline()
        model_name = model_label(model)

        do class_index = 1, 2
            if (class_index == 1) then
                class_name = 'passing'
                eta_min = etatp*1.0e-3_dp
                eta_max = etatp*(1.0_dp - 1.0e-3_dp)
            else
                class_name = 'trapped'
                eta_min = etatp*(1.0_dp + 1.0e-3_dp)
                eta_max = etatp + (etadt - etatp)*(1.0_dp - 1.0e-3_dp)
            end if
            do k = 0, nfrequency - 1
                eta_fraction = real(k, dp)/real(nfrequency - 1, dp)
                eta = eta_min + eta_fraction*(eta_max - eta_min)
                call Om_th(velocity, eta, omega_b, domega_bdv, domega_bdeta)
                call Om_ph(velocity, eta, omega_phi, domega_phdv, domega_phdeta)
                call Om_tB(velocity, eta, omega_magnetic, domega_bmagdv, &
                    domega_bmagdeta)
                if (model == FREQUENCY_MODEL_GC_THIN) then
                    gc_orbit_class = merge(GC_ORBIT_PASSING, GC_ORBIT_TRAPPED, &
                        class_index == 1)
                    call evaluate_gc_spline(velocity, eta, int(sign_vpar), &
                        gc_orbit_class, gc_omega_b, gc_domega_bdv, &
                        gc_domega_bdeta, gc_omega_magnetic, &
                        gc_domega_magnetic_dv, gc_domega_magnetic_deta, &
                        omega_electric, gc_domega_electric_deta, gc_status)
                    if (gc_status /= GC_SPLINE_SUCCESS) then
                        error stop 'GC frequency report spline evaluation failed'
                    end if
                else
                    omega_electric = Om_tE
                end if
                write (unit_frequency, '(a,1x,a,1x,5(es20.12,1x))') &
                    trim(model_name), trim(class_name), eta/etatp, omega_b, &
                    omega_phi, omega_magnetic, omega_electric
            end do
        end do

        do class_index = 1, 2
            if (class_index == 1) then
                class_name = 'passing'
                eta_min = etatp*1.0e-3_dp
                eta_max = etatp*(1.0_dp - 1.0e-3_dp)
            else
                class_name = 'trapped'
                eta_min = etatp*(1.0_dp + 1.0e-3_dp)
                eta_max = etatp + (etadt - etatp)*(1.0_dp - 1.0e-3_dp)
            end if
            do mth = -9, 9
                call driftorbit_coarse(velocity, eta_min, eta_max, roots, nroots)
                do root_index = 1, nroots
                    root_value = driftorbit_root(velocity, &
                        1.0e-8_dp*abs(Om_tE), roots(root_index, 1), &
                        roots(root_index, 2))
                    if (root_value(1) < 0.0_dp) cycle
                    call Om_th(velocity, root_value(1), omega_b, domega_bdv, &
                        domega_bdeta)
                    call Om_ph(velocity, root_value(1), omega_phi, domega_phdv, &
                        domega_phdeta)
                    residual = mth*omega_b + mph*omega_phi
                    write (unit_resonance, '(a,1x,a,1x,i4,1x,3(es20.12,1x))') &
                        trim(model_name), trim(class_name), mth, &
                        root_value(1)/etatp, residual, &
                        real(mth, dp)*domega_bdeta + real(mph, dp)*domega_phdeta
                end do
            end do
        end do
    end do

    close (unit_frequency)
    close (unit_resonance)

contains

    function model_label(value) result(label)
        integer, intent(in) :: value
        character(len=16) :: label

        select case (value)
        case (FREQUENCY_MODEL_LEGACY)
            label = 'legacy'
        case (FREQUENCY_MODEL_GC_THIN)
            label = 'gc_thin'
        case default
            label = 'unknown'
        end select
    end function model_label

end program test_gc_frequency_report
