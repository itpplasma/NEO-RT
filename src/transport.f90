module neort_transport
    use iso_fortran_env, only: dp => real64
    use util, only: imun, pi, c, qi
    use logger, only: trace, debug, warning, error
    use do_magfie_mod, only: do_magfie, s, a, R0, iota, q, psi_pr, eps, inp_swi, &
        bphcov, dbthcovds, dbphcovds, q, dqds, sign_theta, Bthcov
    use do_magfie_pert_mod, only: do_magfie_pert_amp
    use neort_magfie, only: dVds, B0
    use neort_profiles, only: ni1, Om_tE
    use neort_nonlin, only: nonlinear_attenuation
    use neort_freq, only: Om_th, Om_ph
    use neort_gc_frequency_splines, only: &
        evaluate_gc_phase_average_surface, gc_spline_q, &
        evaluate_gc_full_orbit_frequency_surface, &
        evaluate_gc_full_orbit_phase_average_surface
    use neort_gc_frequency_provider, only: gc_full_orbit_frequency_result_t, &
        GC_FREQUENCY_SUCCESS, GC_FREQUENCY_ORBIT_ERROR
    use neort_gc_full_resonance, only: GC_RESONANCE_SUCCESS, GC_RESONANCE_PARTIAL, &
        GC_RESONANCE_INVALID_INPUT, GC_RESONANCE_BOUNDARY_INVALID, &
        GC_RESONANCE_SAMPLE_VALID, GC_RESONANCE_SAMPLE_BOUNDARY, &
        GC_RESONANCE_SAMPLE_UNCONFINED, GC_RESONANCE_SAMPLE_WALL, &
        GC_RESONANCE_SAMPLE_RADIAL_DOMAIN, GC_RESONANCE_SAMPLE_INVALID, &
        gc_resonance_diagnostics_t, find_gc_resonances
    use neort_gc_orbit_integrator, only: GC_ORBIT_SUCCESS, GC_ORBIT_TRAPPED, &
        GC_ORBIT_PASSING, GC_ORBIT_FIELD_ERROR, GC_ORBIT_STATE_ERROR, &
        GC_ORBIT_START_ROOT_ERROR, GC_ORBIT_INTEGRATOR_ERROR, GC_ORBIT_NO_RETURN, &
        GC_ORBIT_PERTURBATION_ERROR, GC_ORBIT_UNCONFINED, GC_ORBIT_WALL_LOSS, &
        GC_ORBIT_RADIAL_DOMAIN, gc_orbit_average_t, &
        gc_orbit_status_is_physical_loss
    use neort_orbit, only: bounce_fast, nvar, noshear, poloidal_velocity
    use neort_resonance, only: driftorbit_coarse, driftorbit_root
    use driftorbit, only: vth, mth, mph, mi, B0, Bmin, Bmax, comptorque, epsmn, &
        etamin, etamax, A1, A2, nlev, pertfile, nonlin, m0, etatp, etadt, &
        sign_vpar_htheta, sign_vpar, frequency_model, &
        FREQUENCY_MODEL_GC_THIN, FREQUENCY_MODEL_GC_FULL

    implicit none

    integer, parameter, public :: GC_TRANSPORT_SUCCESS = 0
    integer, parameter, public :: GC_TRANSPORT_FULL_ORBIT_FAILURE = 71

    ! Residual callback statuses use a separate namespace from both resonance
    ! and orbit statuses.  In particular, GC_RESONANCE_BOUNDARY_INVALID has
    ! the same integer value as GC_ORBIT_INTEGRATOR_ERROR.
    integer, parameter :: GC_FULL_RESIDUAL_ORBIT_BASE = 100
    integer, parameter :: GC_FULL_RESIDUAL_SETUP_ERROR = 200
    integer, parameter :: GC_FULL_RESIDUAL_ZERO_THIN_FREQUENCY = 201
    integer, parameter :: GC_FULL_RESIDUAL_ORBIT_FIELD_ERROR = &
        GC_FULL_RESIDUAL_ORBIT_BASE + GC_ORBIT_FIELD_ERROR
    integer, parameter :: GC_FULL_RESIDUAL_ORBIT_STATE_ERROR = &
        GC_FULL_RESIDUAL_ORBIT_BASE + GC_ORBIT_STATE_ERROR
    integer, parameter :: GC_FULL_RESIDUAL_ORBIT_START_ROOT_ERROR = &
        GC_FULL_RESIDUAL_ORBIT_BASE + GC_ORBIT_START_ROOT_ERROR
    integer, parameter :: GC_FULL_RESIDUAL_ORBIT_INTEGRATOR_ERROR = &
        GC_FULL_RESIDUAL_ORBIT_BASE + GC_ORBIT_INTEGRATOR_ERROR
    integer, parameter :: GC_FULL_RESIDUAL_ORBIT_NO_RETURN = &
        GC_FULL_RESIDUAL_ORBIT_BASE + GC_ORBIT_NO_RETURN
    integer, parameter :: GC_FULL_RESIDUAL_ORBIT_PERTURBATION_ERROR = &
        GC_FULL_RESIDUAL_ORBIT_BASE + GC_ORBIT_PERTURBATION_ERROR
    integer, parameter :: GC_FULL_RESIDUAL_ORBIT_UNCONFINED = &
        GC_FULL_RESIDUAL_ORBIT_BASE + GC_ORBIT_UNCONFINED
    integer, parameter :: GC_FULL_RESIDUAL_ORBIT_WALL_LOSS = &
        GC_FULL_RESIDUAL_ORBIT_BASE + GC_ORBIT_WALL_LOSS
    integer, parameter :: GC_FULL_RESIDUAL_ORBIT_RADIAL_DOMAIN = &
        GC_FULL_RESIDUAL_ORBIT_BASE + GC_ORBIT_RADIAL_DOMAIN

    type, public :: gc_transport_failure_t
        integer :: resonance_partial = 0
        integer :: resonance_failures = 0
        integer :: frequency_failures = 0
        integer :: frequency_setup_failures = 0
        integer :: phase_failures = 0
        integer :: orbit_failures = 0
        integer :: lost_orbits = 0
        integer :: unconfined_samples = 0
        integer :: wall_orbits = 0
        integer :: radial_domain_orbits = 0
        integer :: numerical_samples = 0
        integer :: measure_failures = 0
        integer :: unknown_measure_cells = 0
        integer :: component_count = 0
        integer :: field_error_samples = 0
        integer :: state_error_samples = 0
        integer :: start_root_error_samples = 0
        integer :: integrator_error_samples = 0
        integer :: no_return_samples = 0
        integer :: perturbation_error_samples = 0
        integer :: unknown_status_samples = 0
        integer :: scanned_samples = 0
        integer :: confined_samples = 0
        integer :: boundary_samples = 0
        real(dp) :: canonical_scan_measure = 0.0_dp
        real(dp) :: canonical_confined_measure = 0.0_dp
        real(dp) :: canonical_physical_measure = 0.0_dp
        real(dp) :: canonical_boundary_measure = 0.0_dp
        real(dp) :: canonical_unresolved_measure = 0.0_dp
        real(dp) :: unknown_measure_coordinate_span = 0.0_dp
        real(dp) :: thermal_scan_measure = 0.0_dp
        real(dp) :: thermal_confined_measure = 0.0_dp
        real(dp) :: thermal_physical_measure = 0.0_dp
        real(dp) :: thermal_boundary_measure = 0.0_dp
        real(dp) :: thermal_unresolved_measure = 0.0_dp
        logical :: canonical_measure_certified = .false.
        logical :: component_identity_certified = .false.
        real(dp) :: canonical_coverage_fraction = 0.0_dp
        real(dp) :: physical_canonical_coverage_fraction = 0.0_dp
        real(dp) :: canonical_unresolved_fraction = 0.0_dp
        real(dp) :: thermal_coverage_fraction = 0.0_dp
        real(dp) :: physical_thermal_coverage_fraction = 0.0_dp
        real(dp) :: thermal_unresolved_fraction = 0.0_dp
        real(dp) :: confined_coverage_fraction = 0.0_dp
        real(dp) :: physical_coverage_fraction = 0.0_dp
    end type gc_transport_failure_t

    real(dp) :: Omth, dOmthdv, dOmthdeta

contains

    pure integer function gc_transport_failure_code(failures)
        type(gc_transport_failure_t), intent(in) :: failures

        gc_transport_failure_code = GC_TRANSPORT_SUCCESS
        if (failures%resonance_partial > 0 .or. failures%resonance_failures > 0 &
            .or. failures%frequency_failures > 0 &
            .or. failures%frequency_setup_failures > 0 &
            .or. failures%phase_failures > 0 &
            .or. failures%orbit_failures > 0 .or. failures%unconfined_samples > 0 &
            .or. failures%numerical_samples > 0 &
            .or. failures%radial_domain_orbits > 0 &
            .or. failures%unknown_status_samples > 0 .or. &
            failures%measure_failures > 0 .or. &
            (failures%scanned_samples > 0 &
            .and. .not. failures%canonical_measure_certified) .or. &
            (failures%scanned_samples > 0 &
            .and. .not. failures%component_identity_certified)) then
            gc_transport_failure_code = GC_TRANSPORT_FULL_ORBIT_FAILURE
        end if
    end function gc_transport_failure_code

    pure function fmt_dbg(msg1, v1, msg2, v2, msg3, v3, msg4, v4) result(s)
        ! Helper to compose a short debug line
        character(*), intent(in) :: msg1, msg2
        character(*), intent(in), optional :: msg3, msg4
        real(dp), intent(in) :: v1, v2
        real(dp), intent(in), optional :: v3, v4
        character(len=256) :: s
        character(len=64) :: a1, a2, a3, a4
        a3 = ''; a4 = ''
        write(a1,'(ES12.5)') v1
        write(a2,'(ES12.5)') v2
        if (present(v3)) write(a3,'(ES12.5)') v3
        if (present(v4)) write(a4,'(ES12.5)') v4
        if (present(msg4)) then
            s = trim(msg1)//trim(a1)//' '//trim(msg2)//trim(a2)//' '//trim(msg3)//trim(a3)//' '//trim(msg4)//trim(a4)
        else if (present(msg3)) then
            s = trim(msg1)//trim(a1)//' '//trim(msg2)//trim(a2)//' '//trim(msg3)//trim(a3)
        else
            s = trim(msg1)//trim(a1)//' '//trim(msg2)//trim(a2)
        end if
    end function fmt_dbg

    ! original contains follows

    pure function D11int(ux, taub, Hmn2)
        real(dp) :: D11int
        real(dp), intent(in) :: ux, taub, Hmn2

        D11int = pi**(3.0_dp / 2.0_dp) * mph**2 * c**2 * q * vth &
            / (qi**2 * dVds * abs(psi_pr)) * ux**3 * exp(-ux**2) * taub * Hmn2
    end function D11int

    pure function D12int(ux, taub, Hmn2)
        real(dp) :: D12int
        real(dp), intent(in) :: ux, taub, Hmn2

        D12int = D11int(ux, taub, Hmn2) * ux**2
    end function D12int

    pure function Tphi_int(ux, taub, Hmn2)
        real(dp) :: Tphi_int
        real(dp), intent(in) :: ux, taub, Hmn2

        Tphi_int = sign(1.0_dp, psi_pr * q * sign_theta) * pi**(3.0_dp / 2.0_dp) * mph**2 * ni1 * &
            c * vth / qi &
            * ux**3 * exp(-ux**2) * taub * Hmn2 * (A1 + A2 * ux**2)
    end function Tphi_int

    subroutine compute_transport_integral(vmin, vmax, vsteps, D, T)
        ! compute transport integral via midpoint rule
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(dp), intent(out) :: D(2), T ! Transport coefficients D and torque density T
        real(dp) :: D_plateau, dsdreff ! Plateau diffusion coefficient and ds/dreff=<|grad s|>
        real(dp) :: ux, du, dD11, dD12, dT, v, eta
        real(dp) :: eta_res(2)
        real(dp) :: taub, bounceavg(nvar), Omph, dOmphdv, dOmphdeta
        real(dp) :: q_fieldline
        integer :: istate_dv
        integer :: direct_status, orbit_class
        type(gc_orbit_average_t) :: direct_average
        type(gc_full_orbit_frequency_result_t) :: full_frequency
        real(dp) :: Hmn2, attenuation_factor
        real(dp) :: roots(nlev, 3)
        integer :: nroots, kr, ku
        real(dp) :: full_root_values(nlev), full_root_derivatives(nlev)
        type(gc_transport_failure_t) :: full_failures
        real(dp) :: resonance_scan_min, resonance_scan_max

        call debug(fmt_dbg('compute_transport_integral: vmin=', vmin, ' vmax=', vmax, ' vsteps=', dble(vsteps)))

        D = 0.0_dp
        T = 0.0_dp
        full_failures = gc_transport_failure_t()
        du = (vmax - vmin) / (vsteps * vth)
        ux = vmin / vth + du / 2.0_dp

        do ku = 1, vsteps
            v = ux * vth
            if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
                call collect_full_orbit_roots(v, full_root_values, &
                    full_root_derivatives, nroots, direct_status)
                if (direct_status /= GC_RESONANCE_SUCCESS .and. &
                    direct_status /= GC_RESONANCE_PARTIAL) then
                    call warning(fmt_dbg('full GC resonance search failed: v=', &
                        v, ' status=', dble(direct_status)))
                    nroots = 0
                    full_failures%resonance_failures = &
                        full_failures%resonance_failures + 1
                else if (direct_status == GC_RESONANCE_PARTIAL) then
                    call warning(fmt_dbg('full GC resonance search retained partial roots: v=', &
                        v, ' nroots=', dble(nroots)))
                    full_failures%resonance_partial = &
                        full_failures%resonance_partial + 1
                end if
            else
                call driftorbit_coarse(v, etamin, etamax, roots, nroots)
            end if
            ! No explicit nroots==0 guard: the do-loop below is empty when
            ! nroots==0, and a bare `cycle` here would skip the trailing
            ! `ux = ux + du` velocity-grid increment and stall the sweep.
            do kr = 1, nroots
                if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
                    eta_res = [full_root_values(kr), full_root_derivatives(kr)]
                else
                    eta_res = driftorbit_root(v, 1.0e-8_dp * abs(Om_tE), &
                        roots(kr, 1), roots(kr, 2))
                end if
                if (eta_res(1) < 0.0_dp) cycle
                if (abs(eta_res(2)) <= sqrt(epsilon(eta_res(2))) &
                    *max(1.0_dp, abs(Om_tE)) &
                    /max(etamax - etamin, tiny(eta_res(2)))) then
                    if (frequency_model == FREQUENCY_MODEL_GC_FULL) &
                        full_failures%resonance_failures = &
                        full_failures%resonance_failures + 1
                    call warning(fmt_dbg('ill-conditioned resonance skipped: eta=', &
                        eta_res(1), ' derivative=', eta_res(2)))
                    cycle
                end if
                eta = eta_res(1)

                call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)

                taub = 2.0_dp * pi / abs(Omth)
                if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
                    orbit_class = merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, &
                        eta > etatp)
                    call evaluate_gc_full_orbit_frequency_surface(v, eta, &
                        int(sign_vpar), orbit_class, taub, full_frequency, &
                        direct_status)
                    if (direct_status /= GC_FREQUENCY_SUCCESS) then
                        call record_frequency_outcome(full_failures, &
                            full_frequency%orbit_status, direct_status)
                        cycle
                    end if
                    Omth = full_frequency%omega_b
                    Omph = full_frequency%omega_phi
                    taub = full_frequency%period
                    call evaluate_gc_full_orbit_phase_average_surface(v, eta, &
                        int(sign_vpar), orbit_class, taub, Omth, Omph, mth, mph, &
                        evaluate_direct_perturbation, direct_average, direct_status)
                    if (direct_status /= GC_FREQUENCY_SUCCESS) then
                        call record_phase_outcome(full_failures, direct_average%status)
                        cycle
                    end if
                    bounceavg = 0.0_dp
                    bounceavg(3) = real(direct_average%perturbation_average)
                    bounceavg(4) = aimag(direct_average%perturbation_average)
                    bounceavg(5) = direct_average%inverse_b_average
                    bounceavg(6) = direct_average%b_average
                    istate_dv = 2
                else if (frequency_model == FREQUENCY_MODEL_GC_THIN) then
                    orbit_class = merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, &
                        eta > etatp)
                    call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
                    q_fieldline = gc_spline_q()
                    call evaluate_gc_phase_average_surface(v, eta, &
                        int(sign_vpar), orbit_class, taub, Omth, Omph, &
                        q_fieldline, mth, mph, evaluate_direct_perturbation, &
                        direct_average, direct_status)
                    if (direct_status /= GC_FREQUENCY_SUCCESS) then
                        call error(fmt_dbg('direct GC orbit average failed: mth=', &
                            dble(mth), ' eta=', eta))
                    end if
                    bounceavg = 0.0_dp
                    bounceavg(3) = real(direct_average%perturbation_average)
                    bounceavg(4) = aimag(direct_average%perturbation_average)
                    bounceavg(5) = direct_average%inverse_b_average
                    bounceavg(6) = direct_average%b_average
                    istate_dv = 2
                else
                    call bounce_fast(v, eta, taub, bounceavg, &
                        timestep_transport, istate_dv)
                end if
                if (istate_dv == -1) then
                    call error(fmt_dbg('VODE MXSTEP: mth=', dble(mth), ' ux=', ux, ' eta=', eta, ' taub=', taub))
                else if (istate_dv /= 2) then
                    call warning(fmt_dbg('dvode istate=', dble(istate_dv), ' at mth=', dble(mth), ' ux=', ux, ' eta=', eta))
                else
                    if (abs(eta - etatp) < 1.0e-8_dp*etatp) then
                        call trace(fmt_dbg('near etatp: mth=', dble(mth), ' ux=', ux, ' eta=', eta, ' taub=', taub))
                    end if
                end if
                Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2) * (mi * (ux * vth)**2 / 2.0_dp)**2
                attenuation_factor = nonlinear_attenuation(ux, eta, bounceavg, Omth, &
                    dOmthdv, dOmthdeta, Hmn2)

                dD11 = du * D11int(ux, taub, Hmn2) / abs(eta_res(2))
                dD12 = du * D12int(ux, taub, Hmn2) / abs(eta_res(2))
                D(1) = D(1) + dD11 * attenuation_factor
                D(2) = D(2) + dD12 * attenuation_factor

                if (comptorque) then
                    dT = du * Tphi_int(ux, taub, Hmn2) / abs(eta_res(2))
                    T = T + dT * attenuation_factor
                end if
            end do
            ux = ux + du
        end do

        if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
            call update_transport_coverage(full_failures)
            call emit_gc_transport_failure(full_failures)
            if (gc_transport_failure_code(full_failures) /= GC_TRANSPORT_SUCCESS) then
                D = 0.0_dp
                T = 0.0_dp
                error stop GC_TRANSPORT_FULL_ORBIT_FAILURE
            end if
        end if

        D_plateau = pi * vth**3 / (16.0_dp * R0 * iota * (qi * B0 / (mi * c))**2)
        dsdreff = 2.0_dp / a * sqrt(s) ! TODO: Use exact value instead of this approximation
        D = dsdreff**(-2) * D / D_plateau

        call debug("compute_transport_integral complete")

    contains

        subroutine collect_full_orbit_roots(velocity, root_values, &
                root_derivatives, root_count, root_status)
            real(dp), intent(in) :: velocity
            real(dp), intent(out) :: root_values(:), root_derivatives(:)
            integer, intent(out) :: root_count, root_status

            real(dp) :: region_roots(size(root_values))
            real(dp) :: region_derivatives(size(root_values))
            real(dp) :: lower, upper
            integer :: region_count, region_status
            real(dp) :: thermal_weight
            type(gc_resonance_diagnostics_t) :: region_diagnostics

            root_values = 0.0_dp
            root_derivatives = 0.0_dp
            root_count = 0
            root_status = GC_RESONANCE_SUCCESS

            lower = etamin
            upper = min(etamax, etatp*(1.0_dp - 1.0e-8_dp))
            if (upper > lower) then
                resonance_scan_min = lower
                resonance_scan_max = upper
                call find_gc_resonances(full_residual, lower, upper, nlev, &
                    max(1.0e-8_dp*abs(Om_tE), 1.0e-6_dp), 1.0e-10_dp, &
                    region_roots, region_derivatives, region_count, region_status, &
                    region_diagnostics, classify_full_residual, &
                    full_nonlocal_measure_unavailable)
                thermal_weight = thermal_measure_weight(velocity)
                call accumulate_resonance_diagnostics(full_failures, &
                    region_diagnostics, thermal_weight)
                root_values(1:region_count) = region_roots(1:region_count)
                root_derivatives(1:region_count) = &
                    region_derivatives(1:region_count)
                root_count = region_count
                if (region_status /= GC_RESONANCE_SUCCESS) root_status = region_status
            end if

            lower = max(etamin, etatp*(1.0_dp + 1.0e-8_dp))
            upper = etamax
            if (upper > lower .and. root_count < size(root_values)) then
                resonance_scan_min = lower
                resonance_scan_max = upper
                call find_gc_resonances(full_residual, lower, upper, nlev, &
                    max(1.0e-8_dp*abs(Om_tE), 1.0e-6_dp), 1.0e-10_dp, &
                    region_roots, region_derivatives, region_count, region_status, &
                    region_diagnostics, classify_full_residual, &
                    full_nonlocal_measure_unavailable)
                thermal_weight = thermal_measure_weight(velocity)
                call accumulate_resonance_diagnostics(full_failures, &
                    region_diagnostics, thermal_weight)
                region_count = min(region_count, size(root_values) - root_count)
                root_values(root_count + 1:root_count + region_count) = &
                    region_roots(1:region_count)
                root_derivatives(root_count + 1:root_count + region_count) = &
                    region_derivatives(1:region_count)
                root_count = root_count + region_count
                if (region_status /= GC_RESONANCE_SUCCESS) root_status = region_status
            else if (upper > lower) then
                root_status = GC_RESONANCE_PARTIAL
            end if
        end subroutine collect_full_orbit_roots

        subroutine full_residual(pitch, residual, residual_status)
            real(dp), intent(in) :: pitch
            real(dp), intent(out) :: residual
            integer, intent(out) :: residual_status

            real(dp) :: thin_omega_b, thin_dv, thin_deta, period_estimate
            integer :: local_class
            type(gc_full_orbit_frequency_result_t) :: local_frequency

            residual = 0.0_dp
            call Om_th(v, pitch, thin_omega_b, thin_dv, thin_deta)
            if (thin_omega_b == 0.0_dp) then
                if (is_open_scan_endpoint(pitch)) then
                    residual_status = GC_RESONANCE_BOUNDARY_INVALID
                else
                    full_failures%frequency_failures = full_failures%frequency_failures + 1
                    full_failures%frequency_setup_failures = &
                        full_failures%frequency_setup_failures + 1
                    residual_status = GC_FULL_RESIDUAL_ZERO_THIN_FREQUENCY
                end if
                return
            end if
            period_estimate = 2.0_dp*pi/abs(thin_omega_b)
            local_class = merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, &
                pitch > etatp)
            call evaluate_gc_full_orbit_frequency_surface(v, pitch, &
                int(sign_vpar), local_class, period_estimate, local_frequency, &
                residual_status)
            if (residual_status /= GC_FREQUENCY_SUCCESS) then
                call record_frequency_outcome(full_failures, local_frequency%orbit_status, &
                    residual_status)
                if (local_frequency%orbit_status == GC_ORBIT_SUCCESS) then
                    residual_status = GC_FULL_RESIDUAL_SETUP_ERROR
                else
                    residual_status = encode_full_residual_orbit_status( &
                        local_frequency%orbit_status)
                end if
                return
            end if
            residual = real(mth, dp)*local_frequency%omega_b &
                +real(mph, dp)*local_frequency%omega_phi
        end subroutine full_residual

        logical function is_open_scan_endpoint(pitch)
            real(dp), intent(in) :: pitch

            is_open_scan_endpoint = pitch == resonance_scan_min &
                .or. pitch == resonance_scan_max
        end function is_open_scan_endpoint

        integer function classify_full_residual(sample_status)
            integer, intent(in) :: sample_status

            select case (sample_status)
            case (GC_RESONANCE_SUCCESS)
                classify_full_residual = GC_RESONANCE_SAMPLE_VALID
            case (GC_RESONANCE_BOUNDARY_INVALID)
                classify_full_residual = GC_RESONANCE_SAMPLE_BOUNDARY
            case (GC_FULL_RESIDUAL_ORBIT_UNCONFINED)
                classify_full_residual = GC_RESONANCE_SAMPLE_UNCONFINED
            case (GC_FULL_RESIDUAL_ORBIT_WALL_LOSS)
                classify_full_residual = GC_RESONANCE_SAMPLE_WALL
            case (GC_FULL_RESIDUAL_ORBIT_RADIAL_DOMAIN)
                classify_full_residual = GC_RESONANCE_SAMPLE_RADIAL_DOMAIN
            case (GC_FULL_RESIDUAL_SETUP_ERROR, GC_FULL_RESIDUAL_ZERO_THIN_FREQUENCY, &
                    GC_FULL_RESIDUAL_ORBIT_FIELD_ERROR, &
                    GC_FULL_RESIDUAL_ORBIT_STATE_ERROR, &
                    GC_FULL_RESIDUAL_ORBIT_START_ROOT_ERROR, &
                    GC_FULL_RESIDUAL_ORBIT_INTEGRATOR_ERROR, &
                    GC_FULL_RESIDUAL_ORBIT_NO_RETURN, &
                    GC_FULL_RESIDUAL_ORBIT_PERTURBATION_ERROR)
                classify_full_residual = GC_RESONANCE_SAMPLE_INVALID
            case default
                classify_full_residual = GC_RESONANCE_SAMPLE_INVALID
            end select
        end function classify_full_residual

        integer function encode_full_residual_orbit_status(orbit_status)
            integer, intent(in) :: orbit_status

            encode_full_residual_orbit_status = GC_FULL_RESIDUAL_ORBIT_BASE &
                +orbit_status
        end function encode_full_residual_orbit_status

        subroutine full_nonlocal_measure_unavailable(pitch, density, measure_status)
            !! Model 2 deliberately does not enter the old local eta torque
            !! integral.  Buchholz et al. Eq. 17 needs a nonlocal fixed-H,
            !! fixed-J_perp class coordinate and |d psi_star/dx| at its
            !! resonance root; a callback over the local eta scan cannot
            !! supply that transform without changing the transport integral.
            real(dp), intent(in) :: pitch
            real(dp), intent(out) :: density
            integer, intent(out) :: measure_status

            associate (unused_pitch => pitch)
            end associate
            density = 0.0_dp
            measure_status = GC_RESONANCE_INVALID_INPUT
        end subroutine full_nonlocal_measure_unavailable

        real(dp) function thermal_measure_weight(velocity)
            real(dp), intent(in) :: velocity
            real(dp) :: speed_ratio

            thermal_measure_weight = 0.0_dp
            if (vth <= 0.0_dp .or. velocity <= 0.0_dp) return
            speed_ratio = velocity/vth
            thermal_measure_weight = du*speed_ratio**3*exp(-speed_ratio**2)
        end function thermal_measure_weight

        subroutine accumulate_resonance_diagnostics(failures, diagnostics, thermal_weight)
            type(gc_transport_failure_t), intent(inout) :: failures
            type(gc_resonance_diagnostics_t), intent(in) :: diagnostics
            real(dp), intent(in) :: thermal_weight
            logical :: first_scan

            first_scan = failures%scanned_samples == 0
            failures%scanned_samples = failures%scanned_samples + diagnostics%scan_samples
            failures%confined_samples = failures%confined_samples + diagnostics%confined_samples
            failures%boundary_samples = failures%boundary_samples + diagnostics%boundary_samples
            ! full_residual records every non-success orbit status while it
            ! evaluates a scan sample.  Do not add the same status again from
            ! the scan diagnostics; transport counters also include root and
            ! phase evaluations outside the retained scan cells.
            failures%numerical_samples = failures%numerical_samples &
                +diagnostics%numerical_samples
            failures%measure_failures = failures%measure_failures &
                +diagnostics%measure_failures
            failures%unknown_measure_cells = failures%unknown_measure_cells &
                +diagnostics%unknown_measure_cells
            failures%component_count = failures%component_count &
                +diagnostics%component_count
            if (first_scan) then
                failures%canonical_measure_certified = &
                    diagnostics%canonical_measure_certified
                failures%component_identity_certified = &
                    diagnostics%component_identity_certified
            else
                failures%canonical_measure_certified = &
                    failures%canonical_measure_certified .and. &
                    diagnostics%canonical_measure_certified
                failures%component_identity_certified = &
                    failures%component_identity_certified .and. &
                    diagnostics%component_identity_certified
            end if
            failures%canonical_scan_measure = failures%canonical_scan_measure &
                +diagnostics%canonical_scan_measure
            failures%canonical_confined_measure = failures%canonical_confined_measure &
                +diagnostics%canonical_confined_measure
            failures%canonical_physical_measure = failures%canonical_physical_measure &
                +diagnostics%canonical_physical_measure
            failures%canonical_boundary_measure = failures%canonical_boundary_measure &
                +diagnostics%canonical_boundary_measure
            failures%canonical_unresolved_measure = failures%canonical_unresolved_measure &
                +diagnostics%canonical_unresolved_measure
            failures%unknown_measure_coordinate_span = &
                failures%unknown_measure_coordinate_span &
                +diagnostics%unknown_measure_coordinate_span
            failures%thermal_scan_measure = failures%thermal_scan_measure &
                +thermal_weight*diagnostics%canonical_scan_measure
            failures%thermal_confined_measure = failures%thermal_confined_measure &
                +thermal_weight*diagnostics%canonical_confined_measure
            failures%thermal_physical_measure = failures%thermal_physical_measure &
                +thermal_weight*diagnostics%canonical_physical_measure
            failures%thermal_boundary_measure = failures%thermal_boundary_measure &
                +thermal_weight*diagnostics%canonical_boundary_measure
            failures%thermal_unresolved_measure = failures%thermal_unresolved_measure &
                +thermal_weight*diagnostics%canonical_unresolved_measure
        end subroutine accumulate_resonance_diagnostics

        subroutine record_frequency_outcome(failures, orbit_status, frequency_status)
            type(gc_transport_failure_t), intent(inout) :: failures
            integer, intent(in) :: orbit_status
            integer, intent(in) :: frequency_status

            call record_orbit_status(failures, orbit_status)
            if (gc_orbit_status_is_physical_loss(orbit_status)) then
                call record_physical_orbit(failures, orbit_status)
                return
            end if
            failures%frequency_failures = failures%frequency_failures + 1
            if (frequency_status /= GC_FREQUENCY_ORBIT_ERROR) then
                failures%frequency_setup_failures = &
                    failures%frequency_setup_failures + 1
            end if
            if (orbit_status /= GC_ORBIT_SUCCESS) then
                failures%orbit_failures = failures%orbit_failures + 1
            end if
        end subroutine record_frequency_outcome

        subroutine record_phase_outcome(failures, orbit_status)
            type(gc_transport_failure_t), intent(inout) :: failures
            integer, intent(in) :: orbit_status

            call record_orbit_status(failures, orbit_status)
            if (gc_orbit_status_is_physical_loss(orbit_status)) then
                call record_physical_orbit(failures, orbit_status)
                return
            end if
            failures%phase_failures = failures%phase_failures + 1
            if (orbit_status /= GC_ORBIT_SUCCESS) then
                failures%orbit_failures = failures%orbit_failures + 1
            end if
        end subroutine record_phase_outcome

        subroutine record_orbit_status(failures, orbit_status)
            type(gc_transport_failure_t), intent(inout) :: failures
            integer, intent(in) :: orbit_status

            select case (orbit_status)
            case (GC_ORBIT_SUCCESS)
                continue
            case (GC_ORBIT_FIELD_ERROR)
                failures%field_error_samples = failures%field_error_samples + 1
            case (GC_ORBIT_STATE_ERROR)
                failures%state_error_samples = failures%state_error_samples + 1
            case (GC_ORBIT_START_ROOT_ERROR)
                failures%start_root_error_samples = failures%start_root_error_samples + 1
            case (GC_ORBIT_INTEGRATOR_ERROR)
                failures%integrator_error_samples = failures%integrator_error_samples + 1
            case (GC_ORBIT_NO_RETURN)
                failures%no_return_samples = failures%no_return_samples + 1
            case (GC_ORBIT_PERTURBATION_ERROR)
                failures%perturbation_error_samples = &
                    failures%perturbation_error_samples + 1
            case (GC_ORBIT_UNCONFINED)
                failures%unconfined_samples = failures%unconfined_samples + 1
            case (GC_ORBIT_WALL_LOSS)
                failures%wall_orbits = failures%wall_orbits + 1
            case (GC_ORBIT_RADIAL_DOMAIN)
                failures%radial_domain_orbits = failures%radial_domain_orbits + 1
            case default
                failures%unknown_status_samples = failures%unknown_status_samples + 1
            end select
        end subroutine record_orbit_status

        subroutine record_physical_orbit(failures, orbit_status)
            type(gc_transport_failure_t), intent(inout) :: failures
            integer, intent(in) :: orbit_status

            if (gc_orbit_status_is_physical_loss(orbit_status)) then
                failures%lost_orbits = failures%lost_orbits + 1
            end if
        end subroutine record_physical_orbit

    end subroutine compute_transport_integral

    subroutine update_transport_coverage(failures)
        type(gc_transport_failure_t), intent(inout) :: failures

        real(dp) :: canonical_eligible, thermal_eligible

        if (.not. failures%canonical_measure_certified) then
            failures%canonical_coverage_fraction = 0.0_dp
            failures%physical_canonical_coverage_fraction = 0.0_dp
            failures%canonical_unresolved_fraction = 0.0_dp
            failures%thermal_coverage_fraction = 0.0_dp
            failures%physical_thermal_coverage_fraction = 0.0_dp
            failures%thermal_unresolved_fraction = 0.0_dp
            failures%confined_coverage_fraction = 0.0_dp
            failures%physical_coverage_fraction = 0.0_dp
            return
        end if

        canonical_eligible = failures%canonical_scan_measure &
            -failures%canonical_boundary_measure
        if (canonical_eligible > 0.0_dp) then
            failures%canonical_coverage_fraction = failures%canonical_confined_measure &
                /canonical_eligible
            failures%physical_canonical_coverage_fraction = &
                (failures%canonical_confined_measure &
                +failures%canonical_physical_measure)/canonical_eligible
            failures%canonical_unresolved_fraction = &
                failures%canonical_unresolved_measure/canonical_eligible
        end if

        thermal_eligible = failures%thermal_scan_measure &
            -failures%thermal_boundary_measure
        if (thermal_eligible > 0.0_dp) then
            failures%thermal_coverage_fraction = failures%thermal_confined_measure &
                /thermal_eligible
            failures%physical_thermal_coverage_fraction = &
                (failures%thermal_confined_measure &
                +failures%thermal_physical_measure)/thermal_eligible
            failures%thermal_unresolved_fraction = &
                failures%thermal_unresolved_measure/thermal_eligible
        end if
        failures%confined_coverage_fraction = failures%canonical_coverage_fraction
        failures%physical_coverage_fraction = &
            failures%physical_canonical_coverage_fraction
    end subroutine update_transport_coverage

    subroutine emit_gc_transport_failure(failures)
        type(gc_transport_failure_t), intent(in) :: failures

        write (0, '(A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0)') &
            'GC_TRANSPORT_STATUS code=', gc_transport_failure_code(failures), &
            'resonance_partial=', failures%resonance_partial, &
            'resonance_failures=', failures%resonance_failures, &
            'frequency_failures=', failures%frequency_failures, &
            'phase_failures=', failures%phase_failures, &
            'orbit_failures=', failures%orbit_failures
        write (0, '(A,I0)') 'GC_TRANSPORT_FREQUENCY frequency_setup_failures=', &
            failures%frequency_setup_failures
        write (0, '(A,I0,1X,A,ES16.8)') &
            'GC_TRANSPORT_MEASURE unknown_measure_cells=', failures%unknown_measure_cells, &
            'unknown_measure_coordinate_span=', failures%unknown_measure_coordinate_span
        write (0, '(A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0)') &
            'GC_TRANSPORT_LOSS lost_orbits=', failures%lost_orbits, &
            'wall_orbits=', failures%wall_orbits, &
            'radial_domain_orbits=', failures%radial_domain_orbits, &
            'unconfined_samples=', failures%unconfined_samples, &
            'numerical_samples=', failures%numerical_samples, &
            'unknown_status_samples=', failures%unknown_status_samples, &
            'measure_failures=', failures%measure_failures
        write (0, '(A,I0,1X,A,I0,1X,A,I0)') &
            'GC_TRANSPORT_TOPOLOGY component_count=', failures%component_count, &
            'canonical_measure_certified=', merge(1, 0, failures%canonical_measure_certified), &
            'component_identity_certified=', merge(1, 0, failures%component_identity_certified)
        write (0, '(A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0,1X,A,I0)') &
            'GC_TRANSPORT_ORBIT_STATUS field_errors=', failures%field_error_samples, &
            'state_errors=', failures%state_error_samples, &
            'start_root_errors=', failures%start_root_error_samples, &
            'integrator_errors=', failures%integrator_error_samples, &
            'no_return=', failures%no_return_samples, &
            'perturbation_errors=', failures%perturbation_error_samples
        write (0, '(A,ES16.8,1X,A,ES16.8,1X,A,ES16.8,1X,A,ES16.8,1X,A,ES16.8,1X,A,ES16.8)') &
            'GC_TRANSPORT_COVERAGE canonical_coverage_fraction=', &
            failures%canonical_coverage_fraction, &
            'physical_canonical_coverage_fraction=', &
            failures%physical_canonical_coverage_fraction, &
            'canonical_unresolved_fraction=', failures%canonical_unresolved_fraction, &
            'thermal_coverage_fraction=', failures%thermal_coverage_fraction, &
            'physical_thermal_coverage_fraction=', &
            failures%physical_thermal_coverage_fraction, &
            'thermal_unresolved_fraction=', failures%thermal_unresolved_fraction
    end subroutine emit_gc_transport_failure

    subroutine evaluate_direct_perturbation(position, bmod, amplitude, status)
        real(dp), intent(in) :: position(3), bmod
        complex(dp), intent(out) :: amplitude
        integer, intent(out) :: status
        complex(dp) :: field_amplitude

        amplitude = (0.0_dp, 0.0_dp)
        status = 1
        if (bmod <= 0.0_dp) return
        if (pertfile) then
            call do_magfie_pert_amp(position, field_amplitude)
            amplitude = epsmn*field_amplitude/bmod
        else
            amplitude = epsmn*exp(imun*real(m0, dp)*position(3))
        end if
        status = 0
    end subroutine evaluate_direct_perturbation

    subroutine timestep_transport(v, eta, neq, t, y, ydot)
        !
        !  Timestep function for orbit integration.
        !  Includes poloidal angle theta and parallel velocity.
        !  More integrands may be added starting from y(3)
        !

        real(dp), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq)
        real(dp), intent(out) :: ydot(neq)

        ! BEGIN TODO: remove all of this after refactoring and re-use routine in orbit
        ! for y(1:3)
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3), Om_tB_v
        real(dp) :: t0, parallel_phase
        complex(dp) :: epsn, Hn ! relative amplitude of perturbation field epsn=Bn/B0
        ! and Hamiltonian Hn = (H - H0)_n

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        call poloidal_velocity(v, eta, bmod, hctrvr(3), hder(3), y(2), ydot)
        if (inp_swi == 11) then
            ydot(7) = y(2)*hctrvr(2)
            parallel_phase = mph*y(7)
        else
            ydot(7) = 0.0_dp
            parallel_phase = q*mph*y(1)
        end if

        ! evaluate orbit averages of Hamiltonian perturbation
        if (pertfile) then
            call do_magfie_pert_amp(x, epsn)
            epsn = epsmn * epsn / bmod
        else
            epsn = epsmn * exp(imun * m0 * y(1))
        end if

        if (eta > etatp) then
            !t0 = 0.25*2*pi/Omth ! Different starting position in orbit
            t0 = 0.0_dp
            Hn = (2.0_dp - eta * bmod) * epsn * exp(imun * (parallel_phase - mth * (t - &
                t0) * Omth))
        else
            Hn = (2.0_dp - eta * bmod) * epsn * exp(imun * (parallel_phase - (mth + q * mph) &
                * t * Omth))
        end if
        ydot(3) = real(Hn)
        ydot(4) = aimag(Hn)

        ! evaluate orbit averages for nonlinear attenuation
        if (nonlin) then
            ydot(5) = 1.0_dp / bmod
            ydot(6) = bmod
        else
            ydot(5:6) = 0.0_dp
        end if

        ! Zero any trailing integrands this routine does not compute (the unused
        ! abs(B) slot, y(7)). The DOP853 solver carries every component through
        ! its stage combinations, so an uninitialised derivative leaves a
        ! denormal in the state that trips the underflow FPE trap.
        if (inp_swi /= 11) ydot(7:) = 0.0_dp
    end subroutine timestep_transport

end module neort_transport
