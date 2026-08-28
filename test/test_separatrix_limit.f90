program test_separatrix_limit
    ! Behavioural check of the near-separatrix thin-orbit model against
    ! oracles that do not use the fit.
    !
    ! Conventions.  eta = v_perp**2/(v**2 B) is the pitch variable, etatp =
    ! 1/B_barrier the trapped-passing boundary of the selected orbit class,
    ! and x = |eta-etatp|/etatp the pitch distance from it.  Om_tB is the
    ! bounce-averaged toroidal magnetic drift frequency and scales as v**2;
    ! the model carries Om_tB/v**2.  sign_theta is the Boozer chart
    ! handedness and psi_pr the toroidal flux at the boundary, both signed.
    !
    ! Oracle 1 (endpoint).  A near-separatrix orbit spends almost all of its
    ! period at the barrier top, so the bounce average of any bounded
    ! integrand tends to that integrand evaluated at the barrier with
    ! eta*B = 1.  For the drift integrand of neort_orbit::timestep the second
    ! bracket carries a factor (1 - eta*B) and drops out, leaving
    !
    !   lim Om_tB/v**2 = -mi c q /(2 qi sign_theta psi_pr) * dlnB/ds |_barrier.
    !
    ! This is a direct do_magfie evaluation, independent of any orbit
    ! integration and of the asymptotic fit. The trapped and the passing
    ! branch approach the same barrier, so both must reproduce it.
    !
    ! Oracle 2 (structure).  The logarithmic coefficient of an orbit integral
    ! is (number of barrier passages per period)/lambda, with lambda the
    ! instability exponent of the saddle. On a field with a single maximum
    ! per period a passing orbit passes the barrier once and a trapped orbit
    ! turns near it twice, so tau_a(trapped) = 2 tau_a(passing).
    !
    ! Oracle 3 (invariance).  taub_estimate only sizes the search window of
    ! the turn-finding integration. A physical bounce time cannot depend on
    ! it.
    !
    ! Oracle 4 (direct orbit average).  Below the fitted window but above the
    ! distance where the parallel-velocity integration loses the pitch
    ! invariant, direct integration is still trustworthy, so the extrapolated
    ! model must reproduce it there.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: s, do_magfie, psi_pr, q, sign_theta
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_freq, only: separatrix_trapped, separatrix_passing, Om_th, Om_tB
    use neort_orbit, only: bounce_time, bounce_fast, timestep, nvar
    use neort_orbit_classes, only: find_periodic_extrema
    use neort_profiles, only: vth
    use driftorbit, only: etatp, sign_vpar
    use util, only: files_exist, pi, mi, qi, c
    implicit none

    ! The endpoint is reached only logarithmically, so the quantity compared
    ! here is the extrapolated limit num_a/tau_a, not a sampled value.  The
    ! measured deviation is a few 1e-3; the tolerance keeps an order of
    ! margin and still rejects a model with a truncated 1/log correction,
    ! which biases this limit by more than two percent.
    real(dp), parameter :: LIMIT_TOLERANCE = 1.0e-2_dp
    real(dp), parameter :: BARRIER_PASSAGE_TOLERANCE = 1.0e-2_dp
    real(dp), parameter :: WINDOW_TOLERANCE = 1.0e-4_dp
    ! Ten times below the innermost fitted sample the model is a genuine
    ! extrapolation, while the parallel-velocity integration still keeps the
    ! pitch invariant to far better than that distance, so the two must agree
    ! to the level of the sampling noise; the measured deviation is below 1e-4.
    real(dp), parameter :: DIRECT_TOLERANCE = 1.0e-3_dp

    real(dp) :: v, g_barrier, theta_barrier
    real(dp) :: trapped_limit, passing_limit
    integer :: maxima_count

    call neort_init('driftorbit.in', 'in_file', 'in_file_pert')
    if (files_exist('plasma.in', 'profile.in')) then
        call neort_prepare_splines('plasma.in', 'profile.in')
    end if
    call neort_setup_at_s(s)
    v = vth
    sign_vpar = 1.0_dp

    if (.not. separatrix_trapped%valid) error stop "trapped separatrix model is invalid"
    if (.not. separatrix_passing%valid) error stop "passing separatrix model is invalid"

    call barrier_drift(theta_barrier, g_barrier, maxima_count)
    trapped_limit = separatrix_trapped%numerator(1)/separatrix_trapped%tau(1)
    passing_limit = separatrix_passing%numerator(1)/separatrix_passing%tau(1)

    write (*, '(A,ES23.15)') 'barrier theta            = ', theta_barrier
    write (*, '(A,ES23.15)') 'analytic  v**2*Om_tB lim = ', v**2*g_barrier
    write (*, '(A,ES23.15)') 'trapped   v**2*Om_tB lim = ', v**2*trapped_limit
    write (*, '(A,ES23.15)') 'passing   v**2*Om_tB lim = ', v**2*passing_limit

    if (abs(g_barrier) <= tiny(1.0_dp)) then
        error stop "degenerate barrier drift makes this test vacuous"
    end if
    if (abs(trapped_limit/g_barrier - 1.0_dp) > LIMIT_TOLERANCE) then
        error stop "trapped separatrix limit disagrees with the barrier-top drift"
    end if
    if (abs(passing_limit/g_barrier - 1.0_dp) > LIMIT_TOLERANCE) then
        error stop "passing separatrix limit disagrees with the barrier-top drift"
    end if
    if (abs(trapped_limit/passing_limit - 1.0_dp) > LIMIT_TOLERANCE) then
        error stop "trapped and passing separatrix limits disagree"
    end if

    if (maxima_count == 1) then
        write (*, '(A,ES23.15)') 'tau_a(trapped)/tau_a(passing) = ', &
            separatrix_trapped%tau(1)/separatrix_passing%tau(1)
        if (abs(separatrix_trapped%tau(1)/(2.0_dp*separatrix_passing%tau(1)) &
                - 1.0_dp) > BARRIER_PASSAGE_TOLERANCE) then
            error stop "trapped log coefficient is not twice the passing one"
        end if
    end if

    call check_window_invariance(1.0e-5_dp)
    call check_window_invariance(1.0e-4_dp)
    call check_window_invariance(1.0e-2_dp)

    call check_against_direct_orbit(0.1_dp*separatrix_trapped%xscale, .true.)
    call check_against_direct_orbit(0.1_dp*separatrix_passing%xscale, .false.)

    write (*, '(A)') "PASS separatrix limit"

contains

    subroutine barrier_drift(theta_out, g_out, nmax_out)
        ! Locate the barrier top of the selected orbit class by bisecting
        ! dlnB/dtheta, then evaluate the drift integrand there at eta*B = 1.
        real(dp), intent(out) :: theta_out, g_out
        integer, intent(out) :: nmax_out
        integer, parameter :: nscan = 4096
        real(dp) :: theta(nscan), bmod_scan(nscan)
        logical :: minima(nscan), maxima(nscan)
        real(dp) :: bmod, sqrtg, xv(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: left, right, middle, dleft, dmiddle
        integer :: k, nmin, best

        do k = 1, nscan
            theta(k) = -pi + 2.0_dp*pi*real(k - 1, dp)/real(nscan, dp)
            xv = [s, 0.0_dp, theta(k)]
            call do_magfie(xv, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            bmod_scan(k) = bmod
        end do
        call find_periodic_extrema(bmod_scan, minima, maxima, nmin, nmax_out)
        best = maxloc(bmod_scan, 1)

        left = theta(modulo(best - 2, nscan) + 1)
        right = left + 4.0_dp*pi/real(nscan, dp)
        dleft = dlnb_dtheta(left)
        do k = 1, 200
            middle = 0.5_dp*(left + right)
            if (middle <= left .or. middle >= right) exit
            dmiddle = dlnb_dtheta(middle)
            if (dleft*dmiddle <= 0.0_dp) then
                right = middle
            else
                left = middle
                dleft = dmiddle
            end if
        end do
        theta_out = 0.5_dp*(left + right)

        xv = [s, 0.0_dp, theta_out]
        call do_magfie(xv, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        if (abs(bmod*etatp - 1.0_dp) > 1.0e-9_dp) then
            error stop "located barrier does not match etatp"
        end if
        g_out = -mi*c*q/(2.0_dp*qi*sign_theta*psi_pr)*hder(1)
    end subroutine barrier_drift

    real(dp) function dlnb_dtheta(theta_in) result(dlnb)
        real(dp), intent(in) :: theta_in
        real(dp) :: bmod, sqrtg, xv(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        xv = [s, 0.0_dp, theta_in]
        call do_magfie(xv, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        dlnb = hder(3)
    end function dlnb_dtheta

    subroutine check_window_invariance(x)
        real(dp), intent(in) :: x
        real(dp) :: reference, lo, hi, trial, factor
        integer :: k

        reference = bounce_time(v, etatp*(1.0_dp + x))
        lo = reference
        hi = reference
        do k = -6, 6
            factor = 1.15_dp**k
            trial = bounce_time(v, etatp*(1.0_dp + x), taub_estimate=reference*factor)
            lo = min(lo, trial)
            hi = max(hi, trial)
        end do
        write (*, '(A,ES12.5,A,ES12.5)') 'window spread at x = ', x, ' : ', &
            (hi - lo)/reference
        if ((hi - lo)/reference > WINDOW_TOLERANCE) then
            error stop "bounce time depends on the turn-search window size"
        end if
    end subroutine check_window_invariance

    subroutine check_against_direct_orbit(x, trapped)
        real(dp), intent(in) :: x
        logical, intent(in) :: trapped
        real(dp) :: eta, taub, bounceavg(nvar)
        real(dp) :: Omth, dOmthdv, dOmthdeta, OmtB, dOmtBdv, dOmtBdeta
        real(dp) :: direct_omega, direct_drift

        if (trapped) then
            eta = etatp*(1.0_dp + x)
        else
            eta = etatp*(1.0_dp - x)
        end if
        taub = bounce_time(v, eta)
        call bounce_fast(v, eta, taub, bounceavg, timestep)
        direct_omega = 2.0_dp*pi/taub
        direct_drift = v**2*bounceavg(3)

        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
        write (*, '(A,ES12.5,A,2ES12.5,A,2ES12.5)') 'direct vs model at x = ', x, &
            ' omega_b ', direct_omega, abs(Omth) , ' Om_tB ', direct_drift, OmtB
        if (abs(abs(Omth)/direct_omega - 1.0_dp) > DIRECT_TOLERANCE) then
            error stop "extrapolated bounce frequency disagrees with direct integration"
        end if
        if (abs(OmtB/direct_drift - 1.0_dp) > DIRECT_TOLERANCE) then
            error stop "extrapolated drift disagrees with direct integration"
        end if
    end subroutine check_against_direct_orbit

end program test_separatrix_limit
