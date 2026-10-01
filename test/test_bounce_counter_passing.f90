program test_bounce_counter_passing
    ! Transit time of co- and counter-passing orbits against a quadrature that
    ! does not integrate the orbit,
    !
    !     taub = int_0^{2 pi} dtheta / (v sqrt(1 - eta B) |h^theta|),
    !
    ! which is the same for both directions of motion. bounce() and
    ! bounce_time() used to watch only for theta reaching th0 + 2*pi, so an
    ! orbit with decreasing theta never completed its turn and returned the end
    ! of the search window (~1e3 times the transit time) instead.
    use iso_fortran_env, only: dp => real64
    use util, only: qe, mu, pi
    use driftorbit, only: do_magfie_init, s, M_t, qi, mi, vth, epsmn, m0, &
        mph, mth, magdrift, nopassing, pertfile, nonlin, bfac, efac, inp_swi, &
        etatp, sign_vpar
    use neort, only: init
    use neort_orbit, only: bounce, bounce_time, nvar, noshear

    implicit none

    ! Co-passing bounce() agrees with the quadrature to ~2e-7 here; the bound
    ! leaves room for that, not for a missed turn.
    real(dp), parameter :: tol = 1.0e-5_dp
    real(dp), parameter :: sgn(2) = [1.0_dp, -1.0_dp]
    real(dp) :: v, eta, taub_quad, coarse, taub, taub_t, bounceavg(nvar), err
    integer :: k, i, nfail

    s = 0.153_dp
    M_t = 0.1_dp
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
    qi = qe
    mi = 2.014_dp*mu
    call do_magfie_init("in_file")
    call init
    v = vth

    nfail = 0
    do k = 1, 3
        eta = 0.25_dp*k*etatp
        coarse = quadrature_transit_time(256)
        taub_quad = quadrature_transit_time(512)
        if (abs(taub_quad - coarse) > 1.0e-10_dp*taub_quad) then
            print *, 'oracle not converged', coarse, taub_quad
            error stop
        end if
        do i = 1, 2
            sign_vpar = sgn(i)
            call bounce(v, eta, taub, bounceavg)
            taub_t = bounce_time(v, eta)
            err = max(abs(taub - taub_quad), abs(taub_t - taub_quad))/taub_quad
            print '(A,F5.2,A,F5.1,3ES22.14,ES10.2)', 'eta/etatp=', 0.25_dp*k, &
                ' sign_vpar=', sgn(i), taub_quad, taub, taub_t, err
            if (err > tol) nfail = nfail + 1
        end do
    end do
    sign_vpar = 1.0_dp

    if (nfail > 0) then
        print *, 'test_bounce_counter_passing failed:', nfail
        error stop
    end if
    print *, 'test_bounce_counter_passing OK'

contains

    function quadrature_transit_time(n) result(taub)
        use fortnum_quadrature, only: gauss_legendre_ab
        use do_magfie_mod, only: do_magfie

        integer, intent(in) :: n
        real(dp) :: taub
        real(dp) :: u(n), w(n), x(3), bmod, sqrtg, hder(3), hcovar(3), &
            hctrvr(3), hcurl(3)
        integer :: j

        call gauss_legendre_ab(n, 0.0_dp, 2.0_dp*pi, u, w)
        taub = 0.0_dp
        do j = 1, n
            x(1) = s
            x(2) = 0.0_dp
            x(3) = u(j)
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            taub = taub + w(j)/(v*sqrt(1.0_dp - eta*bmod)*abs(hctrvr(3)))
        end do
    end function quadrature_transit_time

end program test_bounce_counter_passing
