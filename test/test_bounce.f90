program test_bounce_program
    use iso_fortran_env, only: dp => real64
    use util, only: qe, mu

    implicit none

    real(dp), parameter :: tol = 1.0e-13_dp
    real(dp) :: v, eta

    call setup
    !  The absolute oracle runs first: it is the only one of these that can
    !  see a defect the two integrations share.
    call test_bounce_against_quadrature
    call test_bounce
    call test_bounce_time
    call test_bounce_fast

contains

    subroutine setup
        use neort, only: init
        use driftorbit, only: do_magfie_init, etamin, etamax, &
            Om_tE, dOm_tEds, etatp, etadt, epst, sign_vpar, vth, M_t, dM_tds
        use do_magfie_mod, only: R0

        call setup_control
        call do_magfie_init("in_file")
        call init

        Om_tE = vth * M_t / R0
        dOm_tEds = vth * dM_tds / R0

        etamin = (1 + epst) * etatp
        etamax = (1 - epst) * etadt
        sign_vpar = 1

        v = vth
        eta = 0.5_dp * (etamin + etamax)

    end subroutine setup

    subroutine setup_control
        use driftorbit, only: s, M_t, qi, mi, vth, epsmn, m0, &
            mph, mth, magdrift, nopassing, pertfile, &
            nonlin, bfac, efac, inp_swi
        use neort_orbit, only: noshear
        real(dp) :: qs, ms

        s = 0.153_dp
        M_t = 0.1_dp
        qs = 1.0_dp
        ms = 2.014_dp
        vth = 37280978.0_dp
        epsmn = 1.0e-3_dp
        m0 = 0
        mph = 18
        mth = -1
        magdrift = .false.
        nopassing = .false.
        noshear = .true.
        pertfile = .false.
        nonlin = .false.
        bfac = 1.0_dp
        efac = 1.0_dp
        inp_swi = 8

        M_t = M_t*efac/bfac
        qi = qs*qe
        mi = ms*mu
    end subroutine setup_control

    subroutine test_bounce
        use neort_orbit, only: bounce, nvar

        real(dp) :: taub, bounceavg(nvar), bounceavg_tmp(nvar)

        bounceavg = 0.0_dp
        call bounce(v, eta, taub, bounceavg_tmp)

        bounceavg = 0.0_dp
        call bounce(v, eta, taub, bounceavg)

        if (maxval(abs(bounceavg - bounceavg_tmp)) > tol) then
            print *, 'test_bounce failed'
            error stop
        end if

        print *, 'test_bounce OK'
    end subroutine test_bounce


    subroutine test_bounce_time
        use neort_orbit, only: bounce, bounce_time, nvar

        real(dp) :: taub, taub_ref, bounceavg(nvar)

        call bounce(v, eta, taub_ref, bounceavg)
        taub = bounce_time(v, eta)

        if (abs(taub - taub_ref) / taub_ref > 1.0e-8_dp) then
            print *, 'test_bounce_time failed', taub, taub_ref
            error stop
        end if

        print *, 'test_bounce_time OK'
    end subroutine test_bounce_time


    !> Check the bounce time against a quadrature of the same orbit.
    !>
    !> `test_bounce` and `test_bounce_time` compare two integrations of the
    !> same equations against each other, and two integrations sharing a defect
    !> agree. They did: both were wrong by about 0.4% while differing by only
    !> 0.26%, because the ODE solver moved its clock onto a located event root
    !> without rebuilding the history it interpolates from, so everything past
    !> the half-bounce crossing was shifted. Only an oracle that does not
    !> integrate the orbit can see that.
    !>
    !> For a trapped orbit at fixed `v` and `eta`,
    !>
    !>     taub = 2 * int_{th-}^{th+} dtheta / (v sqrt(1 - eta B) h^theta)
    !>
    !> with the turning points where `1 - eta B` vanishes. The square root is
    !> integrable but not differentiable there, so `theta = thc + a sin(u)`
    !> maps the interval onto `u` in `[-pi/2, pi/2]` and cancels it: the
    !> `cos(u)` from `dtheta` divides the `cos(u)` the root contributes.
    !> Gauss-Legendre never lands on the endpoints, so the cancellation is
    !> never evaluated as 0/0.
    subroutine test_bounce_against_quadrature
        use neort_orbit, only: bounce_time

        real(dp) :: taub, taub_quad, coarse, error

        taub = bounce_time(v, eta)
        coarse = quadrature_bounce_time(96)
        taub_quad = quadrature_bounce_time(192)
        !  The quadrature must have converged before it is used as an oracle.
        if (abs(taub_quad - coarse) > 1.0e-9_dp*taub_quad) then
            print *, 'test_bounce_against_quadrature: oracle not converged', &
                coarse, taub_quad
            error stop
        end if
        error = abs(taub - taub_quad)/taub_quad
        print *, 'test_bounce_against_quadrature: taub, oracle, rel', &
            taub, taub_quad, error
        !  The defect this catches was 3.6e-3 here and up to 1.6e-2 elsewhere
        !  in eta; the integrator's own accuracy is several decades better.
        if (error > 1.0e-5_dp) then
            print *, 'test_bounce_against_quadrature failed'
            error stop
        end if
        print *, 'test_bounce_against_quadrature OK'
    end subroutine test_bounce_against_quadrature

    !> Gauss-Legendre bounce time on `n` nodes, independent of the ODE solver.
    function quadrature_bounce_time(n) result(taub)
        use fortnum_quadrature, only: gauss_legendre_ab
        use driftorbit, only: s
        use do_magfie_mod, only: do_magfie

        integer, intent(in) :: n
        real(dp) :: taub

        real(dp) :: u(n), w(n), thm, thp, thc, a, theta, vperp2
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: half_pi
        integer :: k

        half_pi = 2.0_dp*atan(1.0_dp)
        thp = turning_point(0.0_dp, 4.0_dp*half_pi)
        thm = turning_point(0.0_dp, -4.0_dp*half_pi)
        thc = 0.5_dp*(thp + thm)
        a = 0.5_dp*(thp - thm)
        call gauss_legendre_ab(n, -half_pi, half_pi, u, w)
        taub = 0.0_dp
        do k = 1, n
            theta = thc + a*sin(u(k))
            x(1) = s
            x(2) = 0.0_dp
            x(3) = theta
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            vperp2 = 1.0_dp - eta*bmod
            if (vperp2 <= 0.0_dp) then
                print *, 'quadrature_bounce_time: node outside the orbit', theta
                error stop
            end if
            taub = taub + w(k)*a*cos(u(k))/(v*sqrt(vperp2)*hctrvr(3))
        end do
        taub = 2.0_dp*taub
    end function quadrature_bounce_time

    !> Bisect for the turning point of `1 - eta B` between `th0` and `th1`.
    function turning_point(th0_in, th1_in) result(root)
        use driftorbit, only: s
        use do_magfie_mod, only: do_magfie

        real(dp), intent(in) :: th0_in, th1_in
        real(dp) :: root

        real(dp) :: lo, hi, mid, flo, fmid
        integer :: k

        lo = th0_in
        hi = th1_in
        flo = trapping_residual(lo)
        if (flo <= 0.0_dp) then
            print *, 'turning_point: the orbit does not start inside', flo
            error stop
        end if
        do k = 1, 200
            if (trapping_residual(hi) < 0.0_dp) exit
            hi = lo + 0.5_dp*(hi - lo)
            if (k == 200) then
                print *, 'turning_point: no sign change found'
                error stop
            end if
        end do
        do k = 1, 200
            mid = 0.5_dp*(lo + hi)
            fmid = trapping_residual(mid)
            if (fmid > 0.0_dp) then
                lo = mid
            else
                hi = mid
            end if
        end do
        root = 0.5_dp*(lo + hi)
    end function turning_point

    function trapping_residual(theta) result(f)
        use driftorbit, only: s
        use do_magfie_mod, only: do_magfie

        real(dp), intent(in) :: theta
        real(dp) :: f

        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0.0_dp
        x(3) = theta
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        f = 1.0_dp - eta*bmod
    end function trapping_residual

    subroutine test_bounce_fast
        use neort_orbit, only: bounce, bounce_fast, nvar, timestep
        real(dp) :: taub, bounceavg(nvar), bounceavg_tmp(nvar)

        !  `taub` was passed both as the intent(out) result and as the optional
        !  `taub_estimate`, which aliases an argument the callee writes with one
        !  it reads, and reads it before anything has written it.  The estimate
        !  is what sizes the integration step, so an undefined one made the
        !  search run to MXSTEP.  `bounce` finds its own bounce time without it.
        bounceavg_tmp = 0.0_dp
        call bounce(v, eta, taub, bounceavg_tmp)

        bounceavg = 0.0_dp
        call bounce_fast(v, eta, taub, bounceavg, timestep)

        !  Component 1 is the poloidal angle at the end of the orbit, not a
        !  bounce average.  Both codes return to th0, so it is near zero and a
        !  relative comparison against it measures nothing; it is checked
        !  against the orbit's own angular extent instead.
        !
        !  The bound cannot be tightened much further, and the reason is not
        !  accuracy.  `bounce` stops at the root, so it lands on th0 to the root
        !  resolution; `bounce_fast` integrates to a prescribed taub, so its
        !  terminal angle is off by thetadot times however well that taub was
        !  resolved.  The two therefore agree on where the orbit closes to about
        !  1e-5 of a full turn, and asking for more would be asking bounce_fast
        !  to be something it is not.
        print *, 'test_bounce_fast: terminal angles', &
            bounceavg(1), bounceavg_tmp(1)
        if (abs(bounceavg(1) - bounceavg_tmp(1)) > 1.0e-4_dp*8.0_dp*atan(1.0_dp)) then
            print *, 'test_bounce_fast failed: terminal angle'
            error stop
        end if
        if (maxval(abs((bounceavg(2:) - bounceavg_tmp(2:)) &
                       /(bounceavg(2:) + 1.0e-6_dp))) > 1.0e-3_dp) then
            print *, 'test_bounce_fast failed'
            error stop
        end if

        print *, 'test_bounce_fast OK'
    end subroutine test_bounce_fast
end program test_bounce_program
