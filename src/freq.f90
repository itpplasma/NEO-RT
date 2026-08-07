module neort_freq
    use iso_fortran_env, only: dp => real64
    use logger, only: debug, trace, get_log_level, log_result, LOG_TRACE
    use util, only: pi, mi, qi
    use spline, only: spline_coeff, spline_val_0
    use neort_orbit, only: nvar, bounce_fast, bounce_time, timestep
    use neort_profiles, only: vth, Om_tE, dOm_tEds
    use neort_gc_wall_context, only: configured_wall_file, configured_wall_units
    use driftorbit, only: etamin, etamax, etatp, etadt, epsst_spl, epst_spl, epst, magdrift, &
        magdrift_passing, frequency_model, FREQUENCY_MODEL_GC_THIN, &
        FREQUENCY_MODEL_GC_FULL, &
        epssp_spl, epsp_spl, sign_vpar, sign_vpar_htheta, mph, nonlin, supban
    use shaing, only: omph_shaing
    use do_magfie_mod, only: iota, s, Bthcov, Bphcov, q, bfac
    use neort_gc_frequency_provider, only: GC_FREQUENCY_SUCCESS, &
        gc_frequency_context_t, gc_frequency_result_t, &
        initialize_gc_frequency_context, evaluate_gc_frequency
    use neort_gc_frequency_splines, only: GC_SPLINE_SUCCESS, &
        gc_spline_diagnostics_t, initialize_gc_spline_surface, &
        fit_gc_frequency_region, evaluate_gc_spline, gc_spline_q, &
        get_gc_spline_diagnostics
    use neort_gc_orbit_integrator, only: GC_ORBIT_TRAPPED, GC_ORBIT_PASSING
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

    real(dp) :: k_taub_p=0.0_dp, d_taub_p=0.0_dp, k_taub_t=0.0_dp, d_taub_t=0.0_dp ! extrapolation at tp bound
    real(dp) :: k_OmtB_p=0.0_dp, d_Omtb_p=0.0_dp, k_Omtb_t=0.0_dp, d_Omtb_t=0.0_dp ! extrapolation at tp bound

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
    !$omp threadprivate (k_taub_p, d_taub_p, k_taub_t, d_taub_t)
    !$omp threadprivate (k_OmtB_p, d_Omtb_p, k_Omtb_t, d_Omtb_t)

contains

    subroutine freq_thread_init()
        ! Initialize threadprivate variables for this thread
        ! Must be called once per thread before using frequency routines
        freq_trapped_initialized = .false.
        freq_passing_initialized = .false.
        k_taub_p = 0.0_dp
        d_taub_p = 0.0_dp
        k_taub_t = 0.0_dp
        d_taub_t = 0.0_dp
        k_OmtB_p = 0.0_dp
        d_Omtb_p = 0.0_dp
        k_Omtb_t = 0.0_dp
        d_Omtb_t = 0.0_dp
    end subroutine freq_thread_init

    subroutine init_canon_freq_trapped_spline
        ! Initialise splines for canonical frequencies of trapped orbits
        use neort_orbit, only: th0

        real(dp) :: etarange(netaspl), Om_tB_v(netaspl), Omth_v(netaspl)
        integer :: k, gc_status
        real(dp) :: aa, b
        real(dp) :: taub0, taub1, leta0, leta1, OmtB0, OmtB1
        real(dp) :: v, eta, taub, taub_est, bounceavg(nvar)

        ! Auto-initialize threadprivate state if not yet done for this thread
        if (freq_thread_init_state /= FREQ_INIT_SENTINEL) then
            call freq_thread_init()
            freq_thread_init_state = FREQ_INIT_SENTINEL
        end if

        call trace('init_canon_freq_trapped_spline')

        if (frequency_model == FREQUENCY_MODEL_GC_THIN) then
            call init_gc_frequency_trapped_spline()
            call trace('init_canon_freq_trapped_spline complete')
            return
        end if
        if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
            call initialize_gc_spline_surface(s, th0, bfac, Om_tE, mi, qi, &
                vth, gc_status, selected_frequency_model=FREQUENCY_MODEL_GC_FULL, &
                wall_path=configured_wall_file, wall_units=configured_wall_units)
            if (gc_status /= GC_SPLINE_SUCCESS) then
                write(0, *) 'GC full-orbit surface initialization failed:', &
                    gc_status, s
                error stop 'GC full-orbit surface initialization failed'
            end if
        end if

        taub0 = 0.0_dp
        taub1 = 0.0_dp
        leta0 = 0.0_dp
        leta1 = 0.0_dp
        OmtB0 = 0.0_dp
        OmtB1 = 0.0_dp

        v = vth
        ! Start the spline knots where the spline is actually evaluated
        ! (eta > etatp*(1+epst_spl)); knots closer to the trapped-passing
        ! boundary sit in the logarithmic taub singularity and make the global
        ! cubic spline ring. The region below is covered by the log extrapolation.
        etamin = (1.0_dp + epst_spl) * etatp
        etamax = etatp + (etadt - etatp) * (1.0_dp - epsst_spl)
     ! Allocate coefficient arrays for trapped region splines (safe for undefined allocation status)
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
            Omth_v(k + 1) = 2*pi/(v*taub)
            if (k == 0) then
                leta0 = log(eta - etatp)
                taub0 = v*taub
                if (magdrift) OmtB0 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
            if (k == 1) then
                leta1 = log(eta - etatp)
                taub1 = v*taub
                if (magdrift) OmtB1 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
        end do

        k_taub_t = (taub1 - taub0)/(leta1 - leta0)
        d_taub_t = taub0 - k_taub_t*leta0
        Omth_spl_coeff = spline_coeff(etarange, Omth_v)

        if (magdrift) then
            k_OmtB_t = (OmtB1 - OmtB0)/(leta1 - leta0)
            d_OmtB_t = OmtB0 - k_OmtB_t*leta0
            OmtB_spl_coeff = spline_coeff(etarange, Om_tB_v)
        end if

        call trace('init_canon_freq_trapped_spline complete')

    end subroutine init_canon_freq_trapped_spline

    subroutine init_canon_freq_passing_spline
        ! Initialise splines for canonical frequencies of passing orbits

        real(dp) :: etarange(netaspl_pass), Om_tB_v(netaspl_pass), Omth_v(netaspl_pass)
        real(dp) :: aa, b
        integer :: k
        real(dp) :: leta0, leta1, taub0, taub1, OmtB0, OmtB1
        real(dp) :: v, eta, taub, taub_est, bounceavg(nvar)

        ! Auto-initialize threadprivate state if not yet done for this thread
        if (freq_thread_init_state /= FREQ_INIT_SENTINEL) then
            call freq_thread_init()
            freq_thread_init_state = FREQ_INIT_SENTINEL
        end if

        call trace('init_canon_freq_passing_spline')

        if (frequency_model == FREQUENCY_MODEL_GC_THIN) then
            call init_gc_frequency_passing_spline()
            call trace('init_canon_freq_passing_spline complete')
            return
        end if

        taub0 = 0.0_dp
        taub1 = 0.0_dp
        leta0 = 0.0_dp
        leta1 = 0.0_dp
        OmtB0 = 0.0_dp
        OmtB1 = 0.0_dp

        v = vth
        etamin = etatp*epssp_spl
        etamax = etatp
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
            Omth_v(k + 1) = 2*pi/(v*taub)
            if (k == netaspl_pass - 2) then
                leta0 = log(etatp - eta)
                taub0 = v*taub
                if (magdrift_passing > 0) OmtB0 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
            if (k == netaspl_pass - 1) then
                leta1 = log(etatp - eta)
                taub1 = v*taub
                if (magdrift_passing > 0) OmtB1 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
        end do

        k_taub_p = (taub1 - taub0)/(leta1 - leta0)
        d_taub_p = taub0 - k_taub_p*leta0
        Omth_pass_spl_coeff = spline_coeff(etarange, Omth_v)

        if (magdrift_passing > 0) then
            k_OmtB_p = (OmtB1 - OmtB0)/(leta1 - leta0)
            d_OmtB_p = OmtB0 - k_OmtB_p*leta0
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
        real(dp) :: Omth, dOmthdv, dOmthdeta
        real(dp) :: OmE, dOmEdeta
        integer :: orbit_class

        if (frequency_model == FREQUENCY_MODEL_GC_THIN) then
            orbit_class = merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, eta > etatp)
            call evaluate_gc_components_or_stop(v, eta, orbit_class, Omth, &
                dOmthdv, dOmthdeta, OmtB, dOmtBdv, dOmtBdeta, OmE, &
                dOmEdeta)
            if ((orbit_class == GC_ORBIT_TRAPPED .and. .not. magdrift) &
                    .or. (orbit_class == GC_ORBIT_PASSING .and. &
                    (.not. magdrift .or. magdrift_passing <= 0))) then
                OmtB = 0.0_dp
                dOmtBdv = 0.0_dp
                dOmtBdeta = 0.0_dp
            end if
            return
        end if

        if (eta > etatp) then
            if (eta > etatp * (1 + epst_spl)) then
                splineval = spline_val_0(OmtB_spl_coeff, eta)
            else ! extrapolation
                call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
                splineval(1) = sign_vpar * (k_OmtB_t * log(eta - etatp) + d_OmtB_t) * Omth / v
                splineval(2) = sign_vpar * (Omth / v * k_OmtB_t / (eta - etatp) + &
                                            dOmthdeta / v * (k_OmtB_t * log(eta - etatp) &
                                                + d_OmtB_t))
            end if
        else if (magdrift_passing <= 0) then
            ! MARS-K has no passing magnetic-drift path; this reproduces that
            ! scope exactly rather than approximately.
            splineval = 0.0_dp
        else
            if (eta < etatp * (1 - epsp_spl)) then
                splineval = spline_val_0(OmtB_pass_spl_coeff, eta)
            else ! extrapolation
                call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
                splineval(1) = sign_vpar * (k_OmtB_p * log(etatp - eta) + d_OmtB_p) * Omth / v
                splineval(2) = sign_vpar * (Omth / v * k_OmtB_p / (eta - etatp) + &
                                            dOmthdeta / v * (k_OmtB_p * log(etatp - eta) &
                                                + d_OmtB_p))
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
        real(dp) :: OmtE_gc, dOmtEdeta, q_gc
        real(dp) :: deta, dv
        integer :: orbit_class

        if (frequency_model == FREQUENCY_MODEL_GC_THIN) then
            orbit_class = merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, eta > etatp)
            call evaluate_gc_components_or_stop(v, eta, orbit_class, Omth, &
                dOmthdv, dOmthdeta, OmtB, dOmtBdv, dOmtBdeta, OmtE_gc, &
                dOmtEdeta)
            Omph = OmtE_gc
            dOmphdv = 0.0_dp
            dOmphdeta = dOmtEdeta
            if (orbit_class == GC_ORBIT_PASSING) then
                q_gc = gc_spline_q()
                Omph = Omph + q_gc*Omth
                dOmphdv = dOmphdv + q_gc*dOmthdv
                dOmphdeta = dOmphdeta + q_gc*dOmthdeta
            end if
            if (magdrift .and. (orbit_class == GC_ORBIT_TRAPPED &
                    .or. magdrift_passing > 0)) then
                Omph = Omph + OmtB
                dOmphdv = dOmphdv + dOmtBdv
                dOmphdeta = dOmphdeta + dOmtBdeta
            end if
            return
        end if

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
            if (magdrift) then
                call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
                Omph = Omph + OmtB
                dOmphdv = dOmphdv + dOmtBdv
                dOmphdeta = dOmphdeta + dOmtBdeta
            end if
        end if
    end subroutine Om_ph

    subroutine Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        ! returns canonical poloidal frequency
        ! and derivatives w.r.t. v and eta
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: Omth, dOmthdv, dOmthdeta
        real(dp) :: splineval(3)
        real(dp) :: OmtB, dOmtBdv, dOmtBdeta, OmE, dOmEdeta
        integer :: orbit_class

        if (frequency_model == FREQUENCY_MODEL_GC_THIN) then
            orbit_class = merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, eta > etatp)
            call evaluate_gc_components_or_stop(v, eta, orbit_class, Omth, &
                dOmthdv, dOmthdeta, OmtB, dOmtBdv, dOmtBdeta, OmE, &
                dOmEdeta)
            return
        end if

        if (eta > etatp) then
            if (eta > etatp*(1 + epst_spl)) then
                splineval = spline_val_0(Omth_spl_coeff, eta)
            else  ! extrapolation
                splineval(1) = 2.0_dp * pi / (k_taub_t * log(eta - etatp) + d_taub_t)
                splineval(2) = -splineval(1)**2 / (2.0_dp * pi) * k_taub_t / (eta - etatp)
            end if
        else
            if (eta < etatp * (1 - epsp_spl)) then
                splineval = spline_val_0(Omth_pass_spl_coeff, eta)
            else  ! extrapolation
                splineval(1) = 2.0_dp * pi / (k_taub_p * log(etatp - eta) + d_taub_p)
                splineval(2) = -splineval(1)**2 / (2.0_dp * pi) * k_taub_p / (eta - etatp)
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

        if (frequency_model == FREQUENCY_MODEL_GC_THIN) then
            call d_gc_Om_ds(v, eta, taub_estimate, dOmthds, dOmphds)
            call trace('d_Om_ds complete')
            return
        end if

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

    subroutine init_gc_frequency_trapped_spline()
        use neort_orbit, only: th0

        real(dp) :: eta_grid(netaspl), period_estimate(netaspl)
        real(dp) :: aa, b, eta, taub_est
        integer :: k, status, failed_knot

        etamin = (1.0_dp + epst_spl)*etatp
        etamax = etatp + (etadt - etatp)*(1.0_dp - epsst_spl)
        call initialize_gc_spline_surface(s, th0, bfac, Om_tE, mi, qi, vth, &
            status, selected_frequency_model=FREQUENCY_MODEL_GC_THIN)
        if (status /= GC_SPLINE_SUCCESS) then
            write(0, *) 'GC thin-limit surface initialization failed:', status, s
            error stop 'GC thin-limit surface initialization failed'
        end if

        b = log(epst_spl)
        aa = (log(etamax/etamin - 1.0_dp) - b)/real(netaspl - 1, dp)
        taub_est = 0.0_dp
        do k = netaspl - 1, 0, -1
            eta = etamin*(1.0_dp + exp(aa*real(k, dp) + b))
            eta_grid(k + 1) = eta
            if (k == netaspl - 1) then
                taub_est = bounce_time(vth, eta)
            else
                taub_est = bounce_time(vth, eta, taub_estimate=taub_est)
            end if
            period_estimate(k + 1) = taub_est
        end do
        call fit_gc_frequency_region(eta_grid, period_estimate, etatp, &
            GC_ORBIT_TRAPPED, status, failed_knot)
        if (status /= GC_SPLINE_SUCCESS) then
            call get_gc_spline_diagnostics_for_failure()
            write(0, *) 'GC trapped spline failed:', status, failed_knot, &
                eta_grid(max(1, failed_knot))
            error stop 'GC trapped thin-limit spline failed'
        end if
        freq_trapped_initialized = .true.
        call log_gc_spline_diagnostics('trapped')
    end subroutine init_gc_frequency_trapped_spline

    subroutine init_gc_frequency_passing_spline()
        real(dp) :: eta_grid(netaspl_pass), period_estimate(netaspl_pass)
        real(dp) :: aa, b, eta, taub_est
        integer :: k, status, failed_knot

        etamin = etatp*epssp_spl
        etamax = etatp
        b = log((etamax - etamin)/etamax)
        aa = (log(epsp_spl) - b)/real(netaspl_pass - 1, dp)
        taub_est = 0.0_dp
        do k = netaspl_pass - 1, 0, -1
            eta = etamax*(1.0_dp - exp(aa*real(k, dp) + b))
            eta_grid(k + 1) = eta
            if (k == netaspl_pass - 1) then
                taub_est = bounce_time(vth, eta)
            else
                taub_est = bounce_time(vth, eta, taub_estimate=taub_est)
            end if
            period_estimate(k + 1) = taub_est
        end do
        call fit_gc_frequency_region(eta_grid, period_estimate, etatp, &
            GC_ORBIT_PASSING, status, failed_knot)
        if (status /= GC_SPLINE_SUCCESS) then
            call get_gc_spline_diagnostics_for_failure()
            write(0, *) 'GC passing spline failed:', status, failed_knot, &
                eta_grid(max(1, failed_knot))
            error stop 'GC passing thin-limit spline failed'
        end if
        freq_passing_initialized = .true.
        call log_gc_spline_diagnostics('trapped+passing')
    end subroutine init_gc_frequency_passing_spline

    subroutine evaluate_gc_components_or_stop(v, eta, orbit_class, Omth, &
            dOmthdv, dOmthdeta, OmtB, dOmtBdv, dOmtBdeta, OmtE, dOmtEdeta)
        real(dp), intent(in) :: v, eta
        integer, intent(in) :: orbit_class
        real(dp), intent(out) :: Omth, dOmthdv, dOmthdeta
        real(dp), intent(out) :: OmtB, dOmtBdv, dOmtBdeta, OmtE, dOmtEdeta
        integer :: status

        call evaluate_gc_spline(v, eta, int(sign_vpar), orbit_class, Omth, &
            dOmthdv, dOmthdeta, OmtB, dOmtBdv, dOmtBdeta, OmtE, &
            dOmtEdeta, status)
        if (status /= GC_SPLINE_SUCCESS) then
            write(0, *) 'GC frequency spline evaluation failed:', status, &
                s, v, eta, orbit_class
            error stop 'GC frequency spline evaluation failed'
        end if
    end subroutine evaluate_gc_components_or_stop

    subroutine d_gc_Om_ds(v, eta, period_estimate, dOmthds, dOmphds)
        use neort_orbit, only: th0

        real(dp), intent(in) :: v, eta, period_estimate
        real(dp), intent(out) :: dOmthds, dOmphds

        type(gc_frequency_context_t) :: context_minus, context_plus
        type(gc_frequency_result_t) :: result_minus, result_plus
        real(dp) :: radial_step, omega_phi_minus, omega_phi_plus
        integer :: attempt, orbit_class, direction, status_minus, status_plus

        orbit_class = merge(GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, eta > etatp)
        direction = int(sign_vpar)
        radial_step = min(1.0e-4_dp, 0.2_dp*min(s, 1.0_dp - s))
        do attempt = 0, 5
            call initialize_gc_frequency_context(s - radial_step, th0, bfac, &
                Om_tE - radial_step*dOm_tEds, mi, qi, v, context_minus, &
                status_minus, selected_frequency_model=FREQUENCY_MODEL_GC_THIN)
            call initialize_gc_frequency_context(s + radial_step, th0, bfac, &
                Om_tE + radial_step*dOm_tEds, mi, qi, v, context_plus, &
                status_plus, selected_frequency_model=FREQUENCY_MODEL_GC_THIN)
            if (status_minus == GC_FREQUENCY_SUCCESS .and. &
                    status_plus == GC_FREQUENCY_SUCCESS) then
                call evaluate_gc_frequency(context_minus, eta, direction, &
                    orbit_class, abs(period_estimate), result_minus, status_minus)
                call evaluate_gc_frequency(context_plus, eta, direction, &
                    orbit_class, abs(period_estimate), result_plus, status_plus)
            end if
            if (status_minus == GC_FREQUENCY_SUCCESS .and. &
                    status_plus == GC_FREQUENCY_SUCCESS) exit
            radial_step = 0.5_dp*radial_step
        end do
        if (status_minus /= GC_FREQUENCY_SUCCESS .or. &
                status_plus /= GC_FREQUENCY_SUCCESS) then
            write(0, *) 'GC radial frequency derivative failed:', &
                status_minus, status_plus, s, eta, radial_step
            error stop 'GC radial frequency derivative failed'
        end if

        omega_phi_minus = canonical_toroidal_frequency(result_minus, &
            orbit_class)
        omega_phi_plus = canonical_toroidal_frequency(result_plus, orbit_class)
        dOmthds = (result_plus%omega_b - result_minus%omega_b) &
            /(2.0_dp*radial_step)
        dOmphds = (omega_phi_plus - omega_phi_minus)/(2.0_dp*radial_step)
    end subroutine d_gc_Om_ds

    function canonical_toroidal_frequency(result, orbit_class) result(omega)
        type(gc_frequency_result_t), intent(in) :: result
        integer, intent(in) :: orbit_class
        real(dp) :: omega

        omega = result%omega_electric
        if (orbit_class == GC_ORBIT_PASSING) then
            omega = omega + result%q_fieldline*result%omega_b
        end if
        if (magdrift .and. (orbit_class == GC_ORBIT_TRAPPED &
                .or. magdrift_passing > 0)) then
            omega = omega + result%omega_magnetic
        end if
    end function canonical_toroidal_frequency

    subroutine log_gc_spline_diagnostics(label)
        character(*), intent(in) :: label
        type(gc_spline_diagnostics_t) :: value
        character(len=512) :: buffer

        call get_gc_spline_diagnostics(value)
        write(buffer, '(A,1X,A,1X,A,I0,1X,A,ES11.4,1X,A,I0,1X,A,ES11.4,1X,A,ES11.4,1X,A,I0)') &
            'GC thin-limit spline', trim(label), 'knots=', value%orbit_evaluations, &
            'cpu_s=', value%elapsed_seconds, 'max_refine=', &
            value%maximum_refinements, 'max_relerr=', &
            max(value%maximum_magnetic_relative_error, &
                value%maximum_electric_relative_error), 'min_lambda=', &
            value%minimum_lambda, 'boundary_extrap=', &
            value%boundary_extrapolations
        call log_result(buffer)
    end subroutine log_gc_spline_diagnostics

    subroutine get_gc_spline_diagnostics_for_failure()
        type(gc_spline_diagnostics_t) :: value

        call get_gc_spline_diagnostics(value)
        write(0, *) 'GC provider failure status/limits:', &
            value%failed_provider_status, value%failed_magnetic_status, &
            value%failed_total_status
        write(0, *) 'GC provider failure eta/errors/orders:', value%failed_eta, &
            value%failed_magnetic_error, value%failed_electric_error, &
            value%failed_magnetic_order, value%failed_total_order
    end subroutine get_gc_spline_diagnostics_for_failure
end module neort_freq
