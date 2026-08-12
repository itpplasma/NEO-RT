module neort_freq
    use iso_fortran_env, only: dp => real64
    use logger, only: debug, trace, get_log_level, LOG_TRACE
    use util, only: pi
    use spline, only: spline_coeff, spline_val_0
    use neort_orbit, only: nvar, bounce_fast, bounce_time, timestep
    use neort_separatrix, only: separatrix_model_t, fit_separatrix_model, &
        evaluate_separatrix_model
    use neort_profiles, only: vth, Om_tE, dOm_tEds
    use driftorbit, only: etamin, etamax, etatp, etadt, epsst_spl, epst_spl, magdrift, &
        magdrift_passing, &
        epssp_spl, epsp_spl, sign_vpar, sign_vpar_htheta, mph, nonlin, supban
    use shaing, only: omph_shaing
    use do_magfie_mod, only: iota, s, Bthcov, Bphcov, q
    implicit none

    ! For splining in the trapped eta region
    integer, parameter :: netaspl = 100
    real(dp), allocatable :: OmtB_spl_coeff(:, :)
    real(dp), allocatable :: Omth_spl_coeff(:, :)
    real(dp), allocatable :: vres_spl_coeff(:, :)

    ! For splining in the passing eta region
    integer, parameter :: netaspl_pass = 100
    real(dp), allocatable :: OmtB_pass_spl_coeff(:, :)
    real(dp), allocatable :: Omth_pass_spl_coeff(:, :)
    real(dp), allocatable :: vres_pass_spl_coeff(:, :)

    type(separatrix_model_t) :: separatrix_trapped, separatrix_passing

    ! Initialization flags for threadprivate allocatable arrays
    logical :: freq_trapped_initialized = .false.
    logical :: freq_passing_initialized = .false.

    ! Magic sentinel for auto-initializing threadprivate state in worker threads
    integer, parameter :: FREQ_INIT_SENTINEL = 161803398
    integer, save :: freq_thread_init_state = 0

    !$omp threadprivate (freq_trapped_initialized, freq_passing_initialized)
    !$omp threadprivate (freq_thread_init_state)
    !$omp threadprivate (OmtB_spl_coeff, Omth_spl_coeff, vres_spl_coeff)
    !$omp threadprivate (OmtB_pass_spl_coeff, Omth_pass_spl_coeff, vres_pass_spl_coeff)
    !$omp threadprivate (separatrix_trapped, separatrix_passing)

contains

    subroutine freq_thread_init()
        ! Initialize threadprivate variables for this thread
        ! Must be called once per thread before using frequency routines
        freq_trapped_initialized = .false.
        freq_passing_initialized = .false.
        separatrix_trapped = separatrix_model_t()
        separatrix_passing = separatrix_model_t()
    end subroutine freq_thread_init

    subroutine init_canon_freq_trapped_spline
        ! Initialise splines for canonical frequencies of trapped orbits

        real(dp) :: etarange(netaspl), Om_tB_v(netaspl), Omth_v(netaspl)
        real(dp) :: taub_v(netaspl)
        integer :: k
        real(dp) :: aa, b
        real(dp) :: v, eta, taub, taub_est, bounceavg(nvar)
        logical :: fit_ok

        ! Auto-initialize threadprivate state if not yet done for this thread
        if (freq_thread_init_state /= FREQ_INIT_SENTINEL) then
            call freq_thread_init()
            freq_thread_init_state = FREQ_INIT_SENTINEL
        end if

        call trace('init_canon_freq_trapped_spline')

        v = vth
        ! Start the spline knots where the spline is actually evaluated
        ! (eta > etatp*(1+epst_spl)); knots closer to the trapped-passing
        ! boundary sit in the logarithmic taub singularity and make the global
        ! cubic spline ring. The region below is covered by the log extrapolation.
        etamin = (1.0_dp + epst_spl) * etatp
        etamax = etatp + (etadt - etatp) * (1.0_dp - epsst_spl)
        Om_tB_v = 0.0_dp
        taub_v = 0.0_dp
        ! Allocate coefficient arrays for trapped region splines (safe for
        ! undefined allocation status).
        if (.not. freq_trapped_initialized) then
            if (allocated(Omth_spl_coeff)) deallocate(Omth_spl_coeff)
            if (allocated(OmtB_spl_coeff)) deallocate(OmtB_spl_coeff)
            if (allocated(vres_spl_coeff)) deallocate(vres_spl_coeff)
            allocate(Omth_spl_coeff(netaspl - 1, 5))
            allocate(OmtB_spl_coeff(netaspl - 1, 5))
            allocate(vres_spl_coeff(netaspl - 1, 5))
            freq_trapped_initialized = .true.
        end if
        if (get_log_level() >= LOG_TRACE) then
            write(*,'(A)') '[TRACE] init_canon_freq_trapped_spline state:'
            write(*,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  v =', v, 'Om_tE =', Om_tE
            write(*,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  etatp =', etatp, 'etadt =', etadt
            write(*,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  etamin =', etamin, 'etamax =', etamax
            write(*,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  epsst_spl =', epsst_spl, 'epst_spl =', epst_spl
            write(*,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') '  s =', s, 'q =', q, 'iota =', iota
            write(*,'(A,1X,I0,2X,A,1X,I0,2X,A,1X,L1)') '  mph =', mph, 'sign_vpar =', int(sign_vpar), 'nonlin =', nonlin
        end if

        ! logspace for eta
        b = log(epst_spl)
        aa = 1.0_dp / (netaspl - 1.0_dp) * (log(etamax / etamin - 1.0_dp) - b)

        do k = netaspl - 1, 0, -1
            eta = etamin * (1.0_dp + exp(aa * k + b))
            etarange(k + 1) = eta
            if (get_log_level() >= LOG_TRACE) then
                write(*,'(A,I4,A,ES12.5)') '[TRACE] init_canon_freq_trapped_spline k=', k, ' eta=', eta
            end if
            if (k == netaspl - 1) then
                taub_est = bounce_time(v, eta)
            else
                taub_est = bounce_time(v, eta, taub_estimate=taub_est)
            end if
            taub = taub_est
            call bounce_fast(v, eta, taub, bounceavg, timestep)
            if (get_log_level() >= LOG_TRACE) then
                write(*,'(A,I4,A,ES12.5,A,ES12.5)') '[TRACE] init_canon_freq_trapped_spline k=', k, ' eta=', eta, ' taub=', taub
            end if
            if (magdrift) Om_tB_v(k + 1) = bounceavg(3)
            taub_v(k + 1) = v*taub
            Omth_v(k + 1) = 2*pi/(v*taub)
        end do

        call fit_separatrix_model((etarange - etatp)/etatp, taub_v, &
            Om_tB_v, separatrix_trapped, fit_ok)
        if (.not. fit_ok) error stop "trapped separatrix fit failed"
        Omth_spl_coeff = spline_coeff(etarange, Omth_v)

        if (magdrift) then
            OmtB_spl_coeff = spline_coeff(etarange, Om_tB_v)
        end if

        call trace('init_canon_freq_trapped_spline complete')

    end subroutine init_canon_freq_trapped_spline

    subroutine init_canon_freq_passing_spline
        ! Initialise splines for canonical frequencies of passing orbits

        real(dp) :: etarange(netaspl_pass), Om_tB_v(netaspl_pass), Omth_v(netaspl_pass)
        real(dp) :: taub_v(netaspl_pass)
        real(dp) :: aa, b
        integer :: k
        real(dp) :: v, eta, taub, taub_est, bounceavg(nvar)
        logical :: fit_ok

        ! Auto-initialize threadprivate state if not yet done for this thread
        if (freq_thread_init_state /= FREQ_INIT_SENTINEL) then
            call freq_thread_init()
            freq_thread_init_state = FREQ_INIT_SENTINEL
        end if

        call trace('init_canon_freq_passing_spline')

        v = vth
        etamin = etatp*epssp_spl
        etamax = etatp
        Om_tB_v = 0.0_dp
        taub_v = 0.0_dp
        ! Allocate coefficient arrays for passing region splines (safe for undefined allocation status)
        if (.not. freq_passing_initialized) then
            if (allocated(Omth_pass_spl_coeff)) deallocate(Omth_pass_spl_coeff)
            if (allocated(OmtB_pass_spl_coeff)) deallocate(OmtB_pass_spl_coeff)
            if (allocated(vres_pass_spl_coeff)) deallocate(vres_pass_spl_coeff)
            allocate(Omth_pass_spl_coeff(netaspl_pass - 1, 5))
            allocate(OmtB_pass_spl_coeff(netaspl_pass - 1, 5))
            allocate(vres_pass_spl_coeff(netaspl_pass - 1, 5))
            freq_passing_initialized = .true.
        end if
        if (get_log_level() >= LOG_TRACE) then
            write(*,'(A)') '[TRACE] init_canon_freq_passing_spline state:'
            write(*,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  v =', v, 'Om_tE =', Om_tE
            write(*,'(A,1X,ES12.5)') '  etatp =', etatp
            write(*,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  etamin =', etamin, 'etamax =', etamax
            write(*,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  epssp_spl =', epssp_spl, 'epsp_spl =', epsp_spl
            write(*,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') '  s =', s, 'q =', q, 'iota =', iota
            write(*,'(A,1X,I0,2X,A,1X,I0,2X,A,1X,L1)') '  mph =', mph, 'sign_vpar =', int(sign_vpar), 'nonlin =', nonlin
        end if

        b = log((etamax - etamin) / etamax)
        aa = 1.0_dp / (netaspl_pass - 1.0_dp) * (log(epsp_spl) - b)

        do k = netaspl_pass - 1, 0, -1
            eta = etamax * (1.0_dp - exp(aa * k + b))
            etarange(k + 1) = eta
            if (get_log_level() >= LOG_TRACE) then
                write(*,'(A,I4,A,ES12.5)') '[TRACE] init_canon_freq_passing_spline k=', k, ' eta=', eta
            end if
            if (k == netaspl_pass - 1) then
                taub_est = bounce_time(v, eta)
            else
                taub_est = bounce_time(v, eta, taub_estimate=taub_est)
            end if
            taub = taub_est
            call bounce_fast(v, eta, taub, bounceavg, timestep)
            if (get_log_level() >= LOG_TRACE) then
                write(*,'(A,I4,A,ES12.5,A,ES12.5)') '[TRACE] init_canon_freq_passing_spline k=', k, ' eta=', eta, ' taub=', taub
            end if
            if (magdrift_passing > 0) Om_tB_v(k + 1) = bounceavg(3)
            taub_v(k + 1) = v*taub
            Omth_v(k + 1) = 2*pi/(v*taub)
        end do

        call fit_separatrix_model((etatp - etarange)/etatp, taub_v, &
            Om_tB_v, separatrix_passing, fit_ok)
        if (.not. fit_ok) error stop "passing separatrix fit failed"
        Omth_pass_spl_coeff = spline_coeff(etarange, Omth_v)

        if (magdrift_passing > 0) then
            OmtB_pass_spl_coeff = spline_coeff(etarange, Om_tB_v)
        end if
        call trace('init_canon_freq_passing_spline complete')
    end subroutine init_canon_freq_passing_spline

    subroutine Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
        ! returns bounce averaged toroidal magnetic drift frequency
        ! and derivatives w.r.t. v and eta
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: OmtB, dOmtBdv, dOmtBdeta
        real(dp) :: splineval(3)
        real(dp) :: x, tau_model, dtau_dx, omega_model, domega_dx
        real(dp) :: drift_model, ddrift_dx
        logical :: fit_ok
        OmtB = 0.0_dp
        dOmtBdv = 0.0_dp
        dOmtBdeta = 0.0_dp
        if (eta > etatp) then
            if (.not. magdrift) return
            x = (eta - etatp)/etatp
            if (x <= separatrix_trapped%xscale) then
                call evaluate_separatrix_model(separatrix_trapped, x, tau_model, &
                    dtau_dx, omega_model, domega_dx, drift_model, ddrift_dx, fit_ok)
                if (.not. fit_ok) error stop "trapped separatrix drift evaluation failed"
                OmtB = v**2*drift_model
                dOmtBdv = 2.0_dp*v*drift_model
                dOmtBdeta = v**2*ddrift_dx/etatp
                return
            else
                splineval = spline_val_0(OmtB_spl_coeff, eta)
            end if
        else if (magdrift_passing <= 0) then
            ! MARS-K has no passing magnetic-drift path; this reproduces that
            ! scope exactly rather than approximately.
            splineval = 0.0_dp
        else
            x = (etatp - eta)/etatp
            if (x <= separatrix_passing%xscale) then
                call evaluate_separatrix_model(separatrix_passing, x, tau_model, &
                    dtau_dx, omega_model, domega_dx, drift_model, ddrift_dx, fit_ok)
                if (.not. fit_ok) error stop "passing separatrix drift evaluation failed"
                OmtB = v**2*drift_model
                dOmtBdv = 2.0_dp*v*drift_model
                dOmtBdeta = -v**2*ddrift_dx/etatp
                return
            else
                splineval = spline_val_0(OmtB_pass_spl_coeff, eta)
            end if
        end if
        OmtB = splineval(1) * v**2
        dOmtBdv = 2.0_dp * splineval(1) * v
        dOmtBdeta = splineval(2) * v**2
    end subroutine Om_tB

    subroutine Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        ! returns canonical toroidal frequency
        ! and derivatives w.r.t. v and eta
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: Omph, dOmphdv, dOmphdeta
        real(dp) :: Omth, dOmthdv, dOmthdeta
        real(dp) :: OmtB, dOmtBdv, dOmtBdeta
        real(dp) :: deta, dv

        if (supban .and. (eta > etatp)) then
            ! Shaing superbanana-plateau: the analytic bounce-averaged toroidal
            ! drift frequency replaces the numeric magnetic-drift precession.
            ! Derivatives are taken by centred finite differences because the
            ! resonance root finder needs d(res)/d(eta) for its Jacobian.
            Omph = omph_shaing(v, eta)
            deta = 1.0e-6_dp*eta
            dv = 1.0e-6_dp*v
            dOmphdeta = (omph_shaing(v, eta + deta) - omph_shaing(v, eta - deta)) &
                        /(2.0_dp*deta)
            dOmphdv = (omph_shaing(v + dv, eta) - omph_shaing(v - dv, eta)) &
                      /(2.0_dp*dv)
            return
        end if

        if (eta > etatp) then
            Omph = Om_tE
            dOmphdv = 0.0_dp
            dOmphdeta = 0.0_dp
            if (magdrift) then
                call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
                Omph = Omph + OmtB
                dOmphdv = dOmphdv + dOmtBdv
                dOmphdeta = dOmphdeta + dOmtBdeta
            end if
        else
            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
            Omph = Om_tE + Omth/iota
            dOmphdv = dOmthdv/iota
            dOmphdeta = dOmthdeta/iota
            ! Passing orbits have their own magnetic-drift switch, and Om_tB
            ! honours it. Gating this branch on magdrift instead would drop
            ! the passing drift whenever the two switches are set apart.
            if (magdrift_passing > 0) then
                call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
                Omph = Omph + OmtB
                dOmphdv = dOmphdv + dOmtBdv
                dOmphdeta = dOmphdeta + dOmtBdeta
            end if
        end if
    end subroutine Om_ph

    subroutine Om_ph_passing_from_omth(v, eta, Omth, dOmthdv, dOmthdeta, &
            Omph, dOmphdv, dOmphdeta)
        ! Complete the passing canonical toroidal frequency from one already
        ! evaluated poloidal frequency.  This keeps the optimized resonance
        ! residual consistent with Om_ph while preserving the independent
        ! passing magnetic-drift switch.
        real(dp), intent(in) :: v, eta, Omth, dOmthdv, dOmthdeta
        real(dp), intent(out) :: Omph, dOmphdv, dOmphdeta
        real(dp) :: OmtB, dOmtBdv, dOmtBdeta

        Omph = Om_tE + Omth / iota
        dOmphdv = dOmthdv / iota
        dOmphdeta = dOmthdeta / iota
        if (magdrift_passing > 0) then
            call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
            Omph = Omph + OmtB
            dOmphdv = dOmphdv + dOmtBdv
            dOmphdeta = dOmphdeta + dOmtBdeta
        end if
    end subroutine Om_ph_passing_from_omth

    subroutine Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        ! returns canonical poloidal frequency
        ! and derivatives w.r.t. v and eta
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: Omth, dOmthdv, dOmthdeta
        real(dp) :: splineval(3)
        real(dp) :: x, tau_model, dtau_dx, omega_model, domega_dx
        real(dp) :: drift_model, ddrift_dx
        logical :: fit_ok

        if (eta > etatp) then
            x = (eta - etatp)/etatp
            if (x <= separatrix_trapped%xscale) then
                call evaluate_separatrix_model(separatrix_trapped, x, tau_model, &
                    dtau_dx, omega_model, domega_dx, drift_model, ddrift_dx, fit_ok)
                if (.not. fit_ok) error stop "trapped separatrix frequency evaluation failed"
                Omth = sign_vpar*v*omega_model
                dOmthdv = sign_vpar*omega_model
                dOmthdeta = sign_vpar*v*domega_dx/etatp
                return
            else
                splineval = spline_val_0(Omth_spl_coeff, eta)
            end if
        else
            x = (etatp - eta)/etatp
            if (x <= separatrix_passing%xscale) then
                call evaluate_separatrix_model(separatrix_passing, x, tau_model, &
                    dtau_dx, omega_model, domega_dx, drift_model, ddrift_dx, fit_ok)
                if (.not. fit_ok) error stop "passing separatrix frequency evaluation failed"
                Omth = sign_vpar*v*omega_model
                dOmthdv = sign_vpar*omega_model
                dOmthdeta = -sign_vpar*v*domega_dx/etatp
                return
            else
                splineval = spline_val_0(Omth_pass_spl_coeff, eta)
            end if
        end if
        Omth = sign_vpar * splineval(1) * v
        dOmthdv = sign_vpar * splineval(1)
        dOmthdeta = sign_vpar * splineval(2) * v
    end subroutine Om_th

    subroutine d_Om_ds(v, eta, taub_estimate, dOmthds, dOmphds)
        real(dp), intent(in) :: v, eta, taub_estimate
        real(dp), intent(out) :: dOmthds, dOmphds
        real(dp) :: s0, ds, bounceavg(nvar)
        real(dp) :: taub, taub_est, Omth, Omph_noE

        call trace('d_Om_ds')

        ! store current flux surface values
        s0 = s

        ds = 2.0e-8_dp
        s = s0 - ds / 2.0_dp
        taub_est = bounce_time(v, eta, taub_estimate)
        taub = taub_est
        call bounce_fast(v, eta, taub, bounceavg, timestep)
        Omth = sign_vpar_htheta * 2.0_dp * pi / taub
        if (magdrift) then
            if (eta > etatp) then
                Omph_noE = bounceavg(3) * v**2
            else
                Omph_noE = bounceavg(3) * v**2 + Omth / iota
            end if
        else
            if (eta > etatp) then
                Omph_noE = 0.0_dp
            else
                Omph_noE = Omth / iota
            end if
        end if
        s = s0 + ds / 2.0_dp
        taub_est = bounce_time(v, eta, taub_estimate)
        taub = taub_est
        call bounce_fast(v, eta, taub, bounceavg, timestep)
        dOmthds = sign_vpar_htheta * (2.0_dp * pi / taub - sign_vpar_htheta * Omth) / ds
        if (magdrift) then
            if (eta > etatp) then
                dOmphds = dOm_tEds + (bounceavg(3) * v**2 - Omph_noE) / ds
            else
                dOmphds = dOm_tEds + (bounceavg(3) * v**2 + (2.0_dp * pi / taub) / iota - &
                                      Omph_noE) / ds
            end if
        else
            if (eta > etatp) then
                dOmphds = dOm_tEds
            else
                dOmphds = dOm_tEds + ((2.0_dp * pi / taub) / iota - Omph_noE) / ds
            end if
        end if

        ! re-set current flux surface values
        s = s0
        call trace('d_Om_ds complete')
    end subroutine d_Om_ds
end module neort_freq
