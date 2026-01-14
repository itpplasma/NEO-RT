module neort_freq
    use iso_fortran_env, only: dp => real64
    use logger, only: debug, trace, get_log_level, LOG_TRACE
    use util, only: pi
    use spline, only: spline_coeff, spline_val_0
    use neort_orbit, only: nvar, bounce_fast, bounce_time, timestep
    use neort_profiles, only: vth, Om_tE, dOm_tEds
    use driftorbit, only: etamin, etamax, etatp, etadt, epsst_spl, epst_spl, epst, magdrift, &
        epssp_spl, epsp_spl, sign_vpar, sign_vpar_htheta, mph, nonlin
    use do_magfie_mod, only: iota, s, Bthcov, Bphcov, q
    use dvode_f90_m, only: vode_thread_init
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
        ! Also initializes VODE threadprivate state since freq uses bounce integration
        call vode_thread_init()
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

        real(dp) :: etarange(netaspl), Om_tB_v(netaspl), Omth_v(netaspl)
        integer :: k
        real(dp) :: aa, b
        real(dp) :: taub0, taub1, leta0, leta1, OmtB0, OmtB1
        real(dp) :: v, eta, taub, taub_est, bounceavg(nvar)

        ! Auto-initialize threadprivate state if not yet done for this thread
        if (freq_thread_init_state /= FREQ_INIT_SENTINEL) then
            call freq_thread_init()
            freq_thread_init_state = FREQ_INIT_SENTINEL
        end if

        call trace('init_canon_freq_trapped_spline')

        taub0 = 0.0_dp
        taub1 = 0.0_dp
        leta0 = 0.0_dp
        leta1 = 0.0_dp
        OmtB0 = 0.0_dp
        OmtB1 = 0.0_dp

        v = vth
        etamin = (1.0_dp + epst) * etatp
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
            if (magdrift) Om_tB_v(k + 1) = bounceavg(3)
            Omth_v(k + 1) = 2*pi/(v*taub)
            if (k == netaspl_pass - 2) then
                leta0 = log(etatp - eta)
                taub0 = v*taub
                if (magdrift) OmtB0 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
            if (k == netaspl_pass - 1) then
                leta1 = log(etatp - eta)
                taub1 = v*taub
                if (magdrift) OmtB1 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
        end do

        k_taub_p = (taub1 - taub0)/(leta1 - leta0)
        d_taub_p = taub0 - k_taub_p*leta0
        Omth_pass_spl_coeff = spline_coeff(etarange, Omth_v)

        if (magdrift) then
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
