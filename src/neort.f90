module neort
    use iso_fortran_env, only: dp => real64
    use logger, only: debug, set_log_level, get_log_level, log_result, LOG_INFO
    use neort_datatypes, only: magfie_data_t, transport_data_t, transport_harmonic_t
    use neort_profiles, only: init_thermodynamic_forces, init_profiles, vth, dvthds, ni1, &
                              dni1ds, Ti1, dTi1ds, qi, mi, mu, qe
    use neort_magfie, only: init_flux_surface_average
    use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline
    use neort_transport, only: compute_transport_integral
    use do_magfie_mod, only: s
    use neort_gc_eqdsk_nonlocal_transport, only: gc_eqdsk_nonlocal_factory_t
    use neort_gc_model2_transport_dispatch, only: &
        GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE, GC_MODEL2_DISPATCH_NOT_CERTIFIED, &
        GC_MODEL2_DISPATCH_SUCCESS, gc_model2_backend_evidence_t, &
        gc_model2_observed_evidence_t, gc_model2_transport_execution_t, &
        gc_model2_transport_options_t, execute_gc_model2_transport, &
        finalize_gc_model2_transport_execution

    implicit none

    ! Number of integration steps in v, set 0 for adaptive integration by quadpack
    integer :: vsteps = 256

    ! Maximum absolute poloidal/orbit resonance harmonic. A negative value
    ! preserves the historical q-dependent automatic range.
    integer :: mth_max_abs = -1

    ! Upper velocity-space cutoff in units of the thermal velocity. Default 4.0
    ! captures the far-tail resonance; the historical hard-coded 3.0 bound
    ! truncated it. Set vmax_over_vth = 3.0 to reproduce pre-2026-07-20 results.
    real(dp) :: vmax_over_vth = 4.0_dp

    type(gc_eqdsk_nonlocal_factory_t), pointer, save :: model2_factory => null()
    type(gc_model2_transport_options_t), save :: model2_options
    type(gc_model2_backend_evidence_t), save :: model2_backend
    type(gc_model2_transport_execution_t), save :: model2_execution
    integer, save :: model2_harmonic_m = 0
    integer, save :: model2_harmonic_n = 0
    integer, save :: model2_last_status = GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE
    logical, save :: model2_configured = .false.

contains

    subroutine configure_model2_transport(factory, harmonic_m, harmonic_n, backend, &
            options, status)
        type(gc_eqdsk_nonlocal_factory_t), target, intent(inout) :: factory
        integer, intent(in) :: harmonic_m, harmonic_n
        type(gc_model2_backend_evidence_t), intent(in) :: backend
        type(gc_model2_transport_options_t), intent(in) :: options
        integer, intent(out) :: status

        model2_configured = .false.
        model2_last_status = GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE
        model2_execution = gc_model2_transport_execution_t()
        if (.not. factory%initialized) then
            status = GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE
            return
        end if
        if (.not. factory%field_ready .or. .not. factory%profile_ready) then
            status = GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE
            return
        end if
        if (.not. factory%perturbation_ready .or. .not. factory%wall_ready) then
            status = GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE
            return
        end if
        model2_factory => factory
        model2_harmonic_m = harmonic_m
        model2_harmonic_n = harmonic_n
        model2_backend = backend
        model2_options = options
        model2_configured = .true.
        model2_last_status = GC_MODEL2_DISPATCH_SUCCESS
        status = GC_MODEL2_DISPATCH_SUCCESS
    end subroutine configure_model2_transport

    subroutine clear_model2_transport()
        nullify(model2_factory)
        model2_options = gc_model2_transport_options_t()
        model2_backend = gc_model2_backend_evidence_t()
        model2_execution = gc_model2_transport_execution_t()
        model2_harmonic_m = 0
        model2_harmonic_n = 0
        model2_last_status = GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE
        model2_configured = .false.
    end subroutine clear_model2_transport

    subroutine get_model2_transport_execution(execution, status)
        type(gc_model2_transport_execution_t), intent(out) :: execution
        integer, intent(out) :: status

        execution = model2_execution
        status = model2_last_status
    end subroutine get_model2_transport_execution

    pure subroutine harmonic_bounds(mph_value, q_value, max_abs, mth_min, mth_max)
        integer, intent(in) :: mph_value, max_abs
        real(dp), intent(in) :: q_value
        integer, intent(out) :: mth_min, mth_max

        if (max_abs >= 0) then
            mth_min = -max_abs
            mth_max = max_abs
        else
            mth_min = -ceiling(2*abs(mph_value*q_value))
            mth_max = ceiling(2*abs(mph_value*q_value))
        end if
    end subroutine harmonic_bounds

    subroutine init
        use do_magfie_mod, only: psi_pr, q
        use driftorbit, only: nopassing, sign_vpar, etamin, etamax, comptorque

        call debug('init')
        call init_flux_surface_average(s)
        call init_canon_freq_trapped_spline
        if (.not. nopassing) call init_canon_freq_passing_spline
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax)
        if (comptorque) call init_thermodynamic_forces(psi_pr, q)
        call debug('init complete')
    end subroutine init

    pure subroutine set_to_trapped_region(eta_min, eta_max)
        use driftorbit, only: epst, etatp, etadt

        real(dp), intent(out) :: eta_min, eta_max

        eta_min = (1 + epst)*etatp
        eta_max = (1 - epst)*etadt
    end subroutine set_to_trapped_region

    pure subroutine set_to_passing_region(eta_min, eta_max)
        use driftorbit, only: epsp, etatp

        real(dp), intent(out) :: eta_min, eta_max

        eta_min = epsp*etatp
        eta_max = (1 - epsp)*etatp
    end subroutine set_to_passing_region

    subroutine check_magfie(data)
        use util, only: pi, c, imun, eV
        use do_magfie_mod, only: do_magfie, sign_theta, R0, a, psi_pr, iota, q, Bthcov, Bphcov, &
            dBthcovds, dBphcovds, eps
        use do_magfie_pert_mod, only: do_magfie_pert_amp, mph
        use driftorbit, only: B0, etatp, etadt, M_t, Om_tE, m0, pertfile, nonlin, &
            epsmn, pertfile_scale, dVds

        type(magfie_data_t), intent(out) :: data

        integer, parameter :: nth = 50
        integer :: k
        real(dp) :: thmin, thmax
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: D_plateau, Drp
        real(dp) :: dpp, dhh, fpeff
        complex(dp) :: bn

        D_plateau = pi * vth**3 / (16.0_dp * R0 * iota * (qi * B0 / (mi * c))**2)
        Drp = 4 * mph * q / (eps**2 * sqrt(pi))  ! actually this is Drp/Dp

        thmin = -pi
        thmax = pi
        x(1) = s
        x(2) = 0.0_dp
        x(3) = 0.0_dp

        data%params%s = s
        data%params%R0 = R0
        data%params%a = a
        data%params%eps = eps
        data%params%psi_pr = psi_pr
        data%params%B0 = B0
        data%params%Bthcov = Bthcov
        data%params%Bphcov = Bphcov
        data%params%dBthcovds = dBthcovds
        data%params%dBphcovds = dBphcovds
        data%params%q = q
        data%params%iota = iota
        data%params%dVds = dVds
        data%params%M_t = M_t
        data%params%Om_tE = Om_tE
        data%params%Om_tBref = c * mi * vth**2 / (2.0_dp * qi * sign_theta * psi_pr)
        data%params%vth = vth
        data%params%T_in_eV = mi / 2.0_dp * vth**2 / eV
        data%params%m0 = 1.0_dp * m0
        data%params%n0 = 1.0_dp * mph
        data%params%Dp = D_plateau
        data%params%Drp_over_Dp = Drp
        data%params%etatp = etatp
        data%params%etadt = etadt
        data%params%pertfile = pertfile
        data%params%nonlin = nonlin
        dpp = 0.0_dp
        dhh = 0.0_dp
        fpeff = 0.0_dp
        data%params%dpp = dpp
        data%params%dhh = dhh
        data%params%fpeff = fpeff

        data%n_points = nth
        if (allocated(data%tensors)) deallocate(data%tensors)
        allocate(data%tensors(nth))

        do k = 1, nth
            x(3) = thmin + (k - 1)*(thmax - thmin)/(nth - 1)
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            if (pertfile) then
                call do_magfie_pert_amp(x, bn)
                bn = pertfile_scale*bn/bmod
            else
                bn = epsmn*exp(imun*m0*x(3))
            end if
            data%tensors(k)%theta = x(3)
            data%tensors(k)%bmod = bmod
            data%tensors(k)%sqrtg = sqrtg
            data%tensors(k)%hder = hder
            data%tensors(k)%hcovar = hcovar
            data%tensors(k)%hctrvr = hctrvr
            data%tensors(k)%hcurl = hcurl
            data%tensors(k)%bn = bn
            data%tensors(k)%eps_exp = epsmn*exp(imun*m0*x(3))
        end do
    end subroutine check_magfie

    subroutine compute_transport(result_)
        use driftorbit, only: mth, M_t, R0, etamin, etamax, sign_vpar, nopassing, comptorque, &
            dVds, mph, dOm_tEds, dM_tds, supban, frequency_model, &
            FREQUENCY_MODEL_GC_FULL
        use neort_profiles, only: Om_tE
        use do_magfie_mod, only: q

        type(transport_data_t), intent(out) :: result_

        integer :: j, idx

        ! Transport coefficients D11, D12 in approximate effective radius r=2.0_dp/a*sqrt(s)
        real(dp) :: Dco(2), Dctr(2), Dt(2)

        ! Torque density dTphi_int/ds for integration over normalized toroidal flux s
        real(dp) :: Tco, Tctr, Tt
        integer :: mthmin, mthmax

        call debug('compute_transport')

        if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
            call compute_model2_transport(result_)
            return
        end if

        if (allocated(result_%harmonics)) deallocate(result_%harmonics)

        Dco = 0.0_dp
        Dctr = 0.0_dp
        Dt = 0.0_dp
        Tco = 0.0_dp
        Tctr = 0.0_dp
        Tt = 0.0_dp

        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0 + M_t*dvthds/R0

        call harmonic_bounds(mph, q, mth_max_abs, mthmin, mthmax)

        if (supban) then
            ! Superbanana plateau is the trapped ell=0 resonance only.
            mthmin = 0
            mthmax = 0
        end if

        if (mthmax < mthmin) then
            ! Edge case: no valid harmonics (upper bound below lower bound).
            ! In this situation we intentionally allocate an empty array.
            allocate(result_%harmonics(0))
        else
            allocate(result_%harmonics(mthmax - mthmin + 1))
        end if

        idx = 0
        do j = mthmin, mthmax
            idx = idx + 1
            call compute_transport_harmonic(j, Dco, Dctr, Dt, Tco, Tctr, Tt, result_%harmonics(idx))
        end do

        result_%summary%M_t = M_t
        result_%summary%Dco = Dco
        result_%summary%Dctr = Dctr
        result_%summary%Dt = Dt

        result_%torque%has_torque = comptorque
        if (comptorque) then
            result_%torque%s = s
            result_%torque%dVds = dVds
            result_%torque%M_t = M_t
            result_%torque%Tco = Tco
            result_%torque%Tctr = Tctr
            result_%torque%Tt = Tt
        else
            result_%torque%s = 0.0_dp
            result_%torque%dVds = 0.0_dp
            result_%torque%M_t = 0.0_dp
            result_%torque%Tco = 0.0_dp
            result_%torque%Tctr = 0.0_dp
            result_%torque%Tt = 0.0_dp
        end if
        call debug('compute_transport complete')
    end subroutine compute_transport

    subroutine compute_model2_transport(result_)
        use driftorbit, only: M_t, comptorque, dVds

        type(transport_data_t), intent(out) :: result_

        type(gc_model2_observed_evidence_t) :: observed
        type(gc_model2_transport_execution_t) :: execution
        integer :: status, i

        if (allocated(result_%harmonics)) deallocate(result_%harmonics)
        result_%summary%M_t = 0.0_dp
        result_%summary%Dco = 0.0_dp
        result_%summary%Dctr = 0.0_dp
        result_%summary%Dt = 0.0_dp
        result_%torque%has_torque = .false.
        result_%torque%s = 0.0_dp
        result_%torque%dVds = 0.0_dp
        result_%torque%M_t = 0.0_dp
        result_%torque%Tco = 0.0_dp
        result_%torque%Tctr = 0.0_dp
        result_%torque%Tt = 0.0_dp
        execution = gc_model2_transport_execution_t()
        if (.not. model2_configured .or. .not. associated(model2_factory)) then
            model2_last_status = GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE
            model2_execution = execution
            return
        end if

        call execute_gc_model2_transport(model2_factory%provider, model2_harmonic_m, &
            model2_harmonic_n, model2_backend, model2_options, execution, status)
        if (status /= GC_MODEL2_DISPATCH_SUCCESS) then
            model2_last_status = status
            model2_execution = execution
            return
        end if

        observed = execution%observed
        call update_model2_backend_execution_certificates(observed)
        call finalize_gc_model2_transport_execution(execution, model2_backend, observed, &
            status)
        model2_last_status = status
        model2_execution = execution
        if (status /= GC_MODEL2_DISPATCH_SUCCESS) return

        allocate(result_%harmonics(size(execution%harmonics)))
        do i = 1, size(execution%harmonics)
            result_%harmonics(i)%Dresco = [execution%harmonics(i)%d11(1), &
                execution%harmonics(i)%d12(1)]
            result_%harmonics(i)%Dresctr = [execution%harmonics(i)%d11(2), &
                execution%harmonics(i)%d12(2)]
            result_%harmonics(i)%Drest = [execution%harmonics(i)%d11(3), &
                execution%harmonics(i)%d12(3)]
            result_%harmonics(i)%Tresco = execution%harmonics(i)%torque(1)
            result_%harmonics(i)%Tresctr = execution%harmonics(i)%torque(2)
            result_%harmonics(i)%Trest = execution%harmonics(i)%torque(3)
            result_%harmonics(i)%vminp_over_vth = 0.0_dp
            result_%harmonics(i)%vmaxp_over_vth = 0.0_dp
            result_%harmonics(i)%vmint_over_vth = 0.0_dp
            result_%harmonics(i)%vmaxt_over_vth = 0.0_dp
            result_%harmonics(i)%mth = &
                execution%harmonics(i)%poloidal_harmonic
        end do
        result_%summary%M_t = M_t
        result_%summary%Dco = [execution%d11_class(1), execution%d12_class(1)]
        result_%summary%Dctr = [execution%d11_class(2), execution%d12_class(2)]
        result_%summary%Dt = [execution%d11_class(3), execution%d12_class(3)]
        result_%torque%has_torque = comptorque
        if (comptorque) then
            result_%torque%s = s
            result_%torque%dVds = dVds
            result_%torque%M_t = M_t
            result_%torque%Tco = execution%torque_class(1)
            result_%torque%Tctr = execution%torque_class(2)
            result_%torque%Tt = execution%torque_class(3)
        else
            result_%torque%s = 0.0_dp
            result_%torque%dVds = 0.0_dp
            result_%torque%M_t = 0.0_dp
            result_%torque%Tco = 0.0_dp
            result_%torque%Tctr = 0.0_dp
            result_%torque%Tt = 0.0_dp
        end if
    end subroutine compute_model2_transport

    subroutine update_model2_backend_execution_certificates(observed)
        type(gc_model2_observed_evidence_t), intent(in) :: observed
        logical :: component_certificate

        model2_backend%topology_certified = .false.
        if (observed%topology_certification_attempts > 0) then
            model2_backend%topology_certified = &
                observed%topology_certification_successes == &
                observed%topology_certification_attempts
        end if

        model2_backend%canonical_measure_certified = .false.
        if (model2_factory%provider%node_ready) then
            if (model2_factory%provider%node_context%charge_c_locked) then
                if (.not. associated(model2_factory%provider%node_context% &
                    canonical_conversion_provider)) then
                    model2_backend%canonical_measure_certified = .true.
                end if
            end if
        end if

        component_certificate = .false.
        if (model2_factory%provider%node_ready) then
            if (model2_factory%provider%node_context%components_enumerated) then
                if (model2_factory%provider%node_class_result%class_complete) then
                    if (allocated(model2_factory%provider%node_components)) then
                        component_certificate = size(model2_factory%provider%node_components) == &
                            model2_factory%provider%node_class_result%nclasses
                    end if
                end if
            end if
        end if
        model2_backend%component_identity_certified = component_certificate
    end subroutine update_model2_backend_execution_certificates

    subroutine compute_transport_harmonic(j, Dco, Dctr, Dt, Tco, Tctr, Tt, harmonic)
        use driftorbit, only: mth, M_t, etamin, etamax, sign_vpar, nopassing, supban

        integer, intent(in) :: j
        real(dp), intent(inout) :: Dco(2), Dctr(2), Dt(2), Tco, Tctr, Tt
        type(transport_harmonic_t), intent(out) :: harmonic

        real(dp) :: Dresco(2), Dresctr(2), Drest(2), Tresco, Tresctr, Trest
        real(dp) :: vminp, vmaxp, vmint, vmaxt
        character(len=256) :: buffer

        write(buffer, '(A,ES12.5,A,I0)') "compute_transport_harmonic: M_t = ", M_t, ", mth = ", j
        call debug(buffer)

        mth = j

        Drest = 0.0_dp
        Dresco = 0.0_dp
        Dresctr = 0.0_dp
        Tresco = 0.0_dp
        Tresctr = 0.0_dp
        Trest = 0.0_dp

        vminp = 1.0e-6_dp * vth
        vmaxp = vmax_over_vth * vth
        vmint = vminp
        vmaxt = vmaxp

        if (supban) then
            ! Shaing superbanana plateau: trapped ell=0 resonance only, with a
            ! wide velocity bracket so the analytic Omph = 0 resonance is
            ! captured across the thermal range.
            vmint = 0.01_dp * vth
            vmaxt = 5.0_dp * vth
            sign_vpar = 1
            call set_to_trapped_region(etamin, etamax)
            call compute_transport_integral(vmint, vmaxt, vsteps, Drest, Trest)
            Dt = Dt + Drest
            Tt = Tt + Trest

            harmonic%mth = mth
            harmonic%Dresco = Dresco
            harmonic%Dresctr = Dresctr
            harmonic%Drest = Drest
            harmonic%Tresco = Tresco
            harmonic%Tresctr = Tresctr
            harmonic%Trest = Trest
            harmonic%vminp_over_vth = vminp/vth
            harmonic%vmaxp_over_vth = vmaxp/vth
            harmonic%vmint_over_vth = vmint/vth
            harmonic%vmaxt_over_vth = vmaxt/vth

            call debug('compute_transport_harmonic complete (supban)')
            return
        end if

        ! Passing resonance (co-passing)
        if (.not. nopassing) then
            sign_vpar = 1
            call set_to_passing_region(etamin, etamax)
            call compute_transport_integral(vminp, vmaxp, vsteps, Dresco, Tresco)
            Dco = Dco + Dresco
            Tco = Tco + Tresco
        end if

        ! Passing resonance (counter-passing)
        if (.not. nopassing) then
            sign_vpar = -1
            call set_to_passing_region(etamin, etamax)
            call compute_transport_integral(vminp, vmaxp, vsteps, Dresctr, Tresctr)
            Dctr = Dctr + Dresctr
            Tctr = Tctr + Tresctr
        end if

        ! Trapped resonance (trapped)
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax)
        call compute_transport_integral(vmint, vmaxt, vsteps, Drest, Trest)
        Dt = Dt + Drest
        Tt = Tt + Trest

        call log_result("")
        write(buffer, '(A,ES12.5,A,I0)') "test_flux: Mt = ", M_t, ", mth = ", mth
        call log_result(buffer)
        write(buffer, '(4ES12.2,2F12.2)') Dresco(1), Dresctr(1), &
            Drest(1), Dresco(1) + Dresctr(1) + Drest(1), &
            vminp/vth, vmint/vth
        call log_result(buffer)
        write(buffer, '(4ES12.2,2F12.2)') Dresco(2), Dresctr(2), &
            Drest(2), Dresco(2) + Dresctr(2) + Drest(2), &
            vmaxp/vth, vmaxt/vth
        call log_result(buffer)

        call log_result("")
        write(buffer, '(A,ES12.5,A,I0)') "compute_torque: Mt = ", M_t, ", mth = ", mth
        call log_result(buffer)
        write(buffer, '(4ES12.2)') Tresco, Tresctr, Trest, Tresco + Tresctr + Trest
        call log_result(buffer)

        harmonic%mth = mth
        harmonic%Dresco = Dresco
        harmonic%Dresctr = Dresctr
        harmonic%Drest = Drest
        harmonic%Tresco = Tresco
        harmonic%Tresctr = Tresctr
        harmonic%Trest = Trest
        harmonic%vminp_over_vth = vminp/vth
        harmonic%vmaxp_over_vth = vmaxp/vth
        harmonic%vmint_over_vth = vmint/vth
        harmonic%vmaxt_over_vth = vmaxt/vth

        call debug('compute_transport_harmonic complete')
    end subroutine compute_transport_harmonic

    subroutine write_magfie_data_to_files(data, base_path)
        type(magfie_data_t), intent(in) :: data
        character(len=*), intent(in) :: base_path

        integer :: k
        integer, parameter :: unit = 9

        open (unit=unit, file=trim(adjustl(base_path))//"_magfie_param.out", recl=1024)
        write (unit, *) "-------------------------"
        write (unit, *) "check_magfie: s         = ", data%params%s
        write (unit, *) "check_magfie: R0        = ", data%params%R0
        write (unit, *) "check_magfie: a         = ", data%params%a
        write (unit, *) "check_magfie: eps       = ", data%params%eps
        write (unit, *) "check_magfie: A         = ", 1.0_dp / data%params%eps
        write (unit, *) "check_magfie: psi_pr    = ", data%params%psi_pr
        write (unit, *) "check_magfie: B0        = ", data%params%B0
        write (unit, *) "check_magfie: Bthcov    = ", data%params%Bthcov
        write (unit, *) "check_magfie: Bphcov    = ", data%params%Bphcov
        write (unit, *) "check_magfie: dBthcovds = ", data%params%dBthcovds
        write (unit, *) "check_magfie: dBphcovds = ", data%params%dBphcovds
        write (unit, *) "check_magfie: q         = ", data%params%q
        write (unit, *) "check_magfie: iota      = ", data%params%iota
        write (unit, *) "check_magfie: dVds      = ", data%params%dVds
        write (unit, *) "check_magfie: M_t       = ", data%params%M_t
        write (unit, *) "check_magfie: Om_tE     = ", data%params%Om_tE
        write (unit, *) "check_magfie: Om_tBref  = ", data%params%Om_tBref
        write (unit, *) "check_magfie: vth       = ", data%params%vth
        write (unit, *) "check_magfie: T [eV]    = ", data%params%T_in_eV
        write (unit, *) "check_magfie: m0        = ", data%params%m0
        write (unit, *) "check_magfie: n0        = ", data%params%n0
        write (unit, *) "check_magfie: Dp        = ", data%params%Dp
        write (unit, *) "check_magfie: Drp/Dp    = ", data%params%Drp_over_Dp
        write (unit, *) "check_magfie: etatp     = ", data%params%etatp
        write (unit, *) "check_magfie: etadt     = ", data%params%etadt
        write (unit, *) "-------------------------"
        write (unit, *) "check_magfie: pertfile  = ", data%params%pertfile
        if (data%params%nonlin) then
            write (unit, *) "-------------------------"
            write (unit, *) "check_magfie: dpp       = ", data%params%dpp
            write (unit, *) "check_magfie: dhh       = ", data%params%dhh
            write (unit, *) "check_magfie: dfpeff    = ", data%params%fpeff
            write (unit, *) "-------------------------"
        end if
        close (unit=unit)

        open (unit=unit, file=trim(adjustl(base_path))//"_magfie.out", recl=1024)
        do k = 1, data%n_points
            write (unit, *) data%tensors(k)%theta, data%tensors(k)%bmod, data%tensors(k)%sqrtg, &
                            data%tensors(k)%hder(1), data%tensors(k)%hder(2), &
                            data%tensors(k)%hder(3), data%tensors(k)%hcovar(1), &
                            data%tensors(k)%hcovar(2), data%tensors(k)%hcovar(3), &
                            data%tensors(k)%hctrvr(1), data%tensors(k)%hctrvr(2), &
                            data%tensors(k)%hctrvr(3), data%tensors(k)%hcurl(1), &
                            data%tensors(k)%hcurl(2), data%tensors(k)%hcurl(3), &
                            real(data%tensors(k)%bn), aimag(data%tensors(k)%bn), &
                            real(data%tensors(k)%eps_exp), aimag(data%tensors(k)%eps_exp)
        end do
        close (unit=unit)
    end subroutine write_magfie_data_to_files

    subroutine write_transport_data_to_files(data, base_path)
        type(transport_data_t), intent(in) :: data
        character(len=*), intent(in) :: base_path

        integer :: k
        real(dp) :: total_D1, total_D2
        integer, parameter :: unit1 = 9
        integer, parameter :: unit2 = 10

        open (unit=unit1, file=trim(adjustl(base_path))//".out", recl=1024)
        write (unit1, *) "# M_t D11co D11ctr D11t D11 D12co D12ctr D12t D12"
        total_D1 = data%summary%Dco(1) + data%summary%Dctr(1) + data%summary%Dt(1)
        total_D2 = data%summary%Dco(2) + data%summary%Dctr(2) + data%summary%Dt(2)
        write (unit1, *) data%summary%M_t, data%summary%Dco(1), data%summary%Dctr(1), &
                        data%summary%Dt(1), total_D1, data%summary%Dco(2), data%summary%Dctr(2), &
                        data%summary%Dt(2), total_D2
        close (unit=unit1)

        if (data%torque%has_torque) then
            open (unit=unit1, file=trim(adjustl(base_path))//"_torque.out", recl=1024)
            write (unit1, *) "# s dVds M_t Tco Tctr Tt"
            write (unit1, *) data%torque%s, data%torque%dVds, data%torque%M_t, data%torque%Tco, &
                             data%torque%Tctr, data%torque%Tt
            close (unit=unit1)
        end if

        open (unit=unit1, file=trim(adjustl(base_path))//"_integral.out", recl=1024)
        open (unit=unit2, file=trim(adjustl(base_path))//"_torque_integral.out", recl=1024)
        do k = 1, size(data%harmonics)
            total_D1 = data%harmonics(k)%Dresco(1) + data%harmonics(k)%Dresctr(1) + &
                       data%harmonics(k)%Drest(1)
            total_D2 = data%harmonics(k)%Dresco(2) + data%harmonics(k)%Dresctr(2) + &
                       data%harmonics(k)%Drest(2)
            write (unit1, *) data%summary%M_t, data%harmonics(k)%mth, data%harmonics(k)%Dresco(1), &
                             data%harmonics(k)%Dresctr(1), data%harmonics(k)%Drest(1), &
                             total_D1, data%harmonics(k)%Dresco(2), data%harmonics(k)%Dresctr(2), &
                             data%harmonics(k)%Drest(2), total_D2, &
                             data%harmonics(k)%vminp_over_vth, data%harmonics(k)%vmaxp_over_vth, &
                             data%harmonics(k)%vmint_over_vth, data%harmonics(k)%vmaxt_over_vth

            write (unit2, *) data%harmonics(k)%mth, data%harmonics(k)%Tresco, &
                             data%harmonics(k)%Tresctr, data%harmonics(k)%Trest
        end do
        close (unit=unit1)
        close (unit=unit2)
    end subroutine write_transport_data_to_files

end module neort
