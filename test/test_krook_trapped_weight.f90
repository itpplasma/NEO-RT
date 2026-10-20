program test_krook_trapped_weight
    ! Trapped Hamiltonian harmonics next to the trapped-passing boundary
    ! against a fixed-panel quadrature of t = int dtheta/|v_par*h^theta| with
    ! theta = theta_m + Delta*sin(phi) between the turning points and
    ! cumulative (spectral) bounce time, converged in the panel count. Time
    ! integration of v_parallel loses the energy invariant; the time spent
    ! near the X-point amplifies that error by 1/(eta*B_max - 1) and made the
    ! Krook pitch quadrature fail to converge.
    use iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_transport, only: transport_bounce, Omth
    use neort_orbit, only: nvar, th0
    use do_magfie_mod, only: do_magfie, s
    use driftorbit, only: etatp, vth, mth, mph, m0, epsmn, sign_vpar, &
        nopassing, comptorque, drift_phase, fow_radial
    implicit none
    integer, parameter :: n = 16
    real(dp), parameter :: near_turning = 0.05_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    ! eta/etatp - 1: the time integration of v_parallel failed the first
    ! three by 3e-4, 3e-5 and 3e-6; the last two cover deeper orbits.
    real(dp), parameter :: gaps(5) = [1.0e-8_dp, 1.0e-7_dp, 1.0e-6_dp, 1.0e-3_dp, 1.0e-1_dp]
    real(dp) :: x(n), w(n), cumulative(n, n)
    real(dp) :: eta, period, coarse_period, average(nvar), relative_error, worst
    complex(dp) :: coarse, fine
    integer :: k, kstate

    call gauss_legendre(n, x, w)
    call integration_matrix()
    call neort_init('driftorbit.in', 'in_file')
    nopassing = .false.
    comptorque = .false.
    drift_phase = .false.
    fow_radial = .false.
    call neort_prepare_splines('plasma.in', 'profile.in')
    call neort_setup_at_s(0.5_dp)
    ! H = (2 - eta*B)*exp(-i*mth*Omth*t): the bounce harmonic is sensitive
    ! to the time spent near the X-point. Odd mth vanish by up-down symmetry.
    mth = 2
    mph = 0
    m0 = 0
    epsmn = 1.0_dp
    sign_vpar = 1.0_dp
    worst = 0.0_dp
    do k = 1, size(gaps)
        eta = etatp*(1.0_dp + gaps(k))
        call orbit_reference(eta, 256, coarse_period, coarse)
        call orbit_reference(eta, 512, period, fine)
        if (abs(period - coarse_period) > 1.0e-9_dp*period .or. &
            abs(fine - coarse) > 1.0e-9_dp*abs(fine)) &
            error stop 'trapped pitch-angle reference has not converged'
        Omth = 2.0_dp*pi/period
        call transport_bounce(vth, eta, period, average, kstate)
        if (kstate /= 2) error stop 'trapped transport bounce failed'
        relative_error = abs(cmplx(average(3), average(4), dp) - fine)/abs(fine)
        worst = max(worst, relative_error)
        print *, 'eta/etatp - 1, harmonic relative error:', &
            eta/etatp - 1.0_dp, relative_error
        if (relative_error > 1.0e-7_dp) &
            error stop 'trapped Hamiltonian harmonic disagrees with pitch-angle oracle'
    end do
    print *, 'trapped weight maximum relative error:', worst
    print *, 'test_krook_trapped_weight PASSED'
contains
    subroutine integration_matrix()
        ! cumulative(i, j) = integral of the j-th Lagrange basis polynomial
        ! of the nodes x from -1 to x(i).
        integer :: i, j, l
        real(dp) :: xi
        cumulative = 0.0_dp
        do i = 1, n
            do l = 1, n
                xi = -1.0_dp + (x(i) + 1.0_dp)*(x(l) + 1.0_dp)/2.0_dp
                do j = 1, n
                    cumulative(i, j) = cumulative(i, j) + &
                        (x(i) + 1.0_dp)/2.0_dp*w(l)*lagrange(j, xi)
                end do
            end do
        end do
    end subroutine integration_matrix

    real(dp) function lagrange(j, xi)
        integer, intent(in) :: j
        real(dp), intent(in) :: xi
        integer :: l
        lagrange = 1.0_dp
        do l = 1, n
            if (l /= j) lagrange = lagrange*(xi - x(l))/(x(j) - x(l))
        end do
    end function lagrange

    subroutine orbit_reference(eta_, panels, period_, harmonic)
        ! One bounce from theta = th0 with increasing theta, as in bounce_fast.
        ! Panels end at the turning points phi = pi/2 and 3*pi/2.
        real(dp), intent(in) :: eta_
        integer, intent(in) :: panels
        real(dp), intent(out) :: period_
        complex(dp), intent(out) :: harmonic
        real(dp) :: lower, upper, middle, half, phi0, edges(4), h, a, t
        real(dp) :: dtdphi(n), hamiltonian(n), times(n)
        integer :: segment, panel, count, pass
        lower = turning_point(eta_, -1.0_dp)
        upper = turning_point(eta_, 1.0_dp)
        middle = (lower + upper)/2.0_dp
        half = (upper - lower)/2.0_dp
        phi0 = asin((th0 - middle)/half)
        edges = [phi0, pi/2.0_dp, 3.0_dp*pi/2.0_dp, phi0 + 2.0_dp*pi]
        do pass = 1, 2
            t = 0.0_dp
            harmonic = (0.0_dp, 0.0_dp)
            do segment = 1, 3
                count = max(1, nint(panels*(edges(segment + 1) - edges(segment))/(2.0_dp*pi)))
                h = (edges(segment + 1) - edges(segment))/count
                do panel = 1, count
                    a = edges(segment) + (panel - 1)*h
                    call panel_values(eta_, lower, upper, a, h, dtdphi, hamiltonian)
                    times = t + h/2.0_dp*matmul(cumulative, dtdphi)
                    if (pass == 2) harmonic = harmonic + h/2.0_dp*sum(w*dtdphi* &
                        hamiltonian*exp(-cmplx(0.0_dp, mth*2.0_dp*pi*times/period_, dp)))
                    t = t + h/2.0_dp*sum(w*dtdphi)
                end do
            end do
            period_ = t
        end do
        harmonic = harmonic/period_
    end subroutine orbit_reference

    subroutine panel_values(eta_, lower, upper, a, h, dtdphi, hamiltonian)
        real(dp), intent(in) :: eta_, lower, upper, a, h
        real(dp), intent(out) :: dtdphi(n), hamiltonian(n)
        real(dp) :: phi, half, gap, position(3), bmod, sqrtg
        real(dp) :: hder(3), hcovar(3), hctrvr(3), hcurl(3)
        integer :: k
        half = (upper - lower)/2.0_dp
        position(1:2) = [s, 0.0_dp]
        do k = 1, n
            phi = a + h*(x(k) + 1.0_dp)/2.0_dp
            position(3) = (lower + upper)/2.0_dp + half*sin(phi)
            call do_magfie(position, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            gap = 1.0_dp - eta_*bmod
            ! Next to a turning point theta_t, 1 - eta*B is a small difference
            ! of O(1) numbers. Integrate dB/dtheta from theta_t instead, which
            ! keeps its relative precision and vanishes exactly at theta_t.
            if (upper - position(3) < near_turning) then
                gap = turning_gap(eta_, position(3), upper)
            else if (position(3) - lower < near_turning) then
                gap = turning_gap(eta_, position(3), lower)
            end if
            dtdphi(k) = half*abs(cos(phi))/(vth*sqrt(gap)*abs(hctrvr(3)))
            hamiltonian(k) = 2.0_dp - eta_*bmod
        end do
    end subroutine panel_values

    real(dp) function turning_gap(eta_, theta, turning)
        ! 1 - eta*B(theta) = eta*integral of dB/dtheta from theta to turning.
        real(dp), intent(in) :: eta_, theta, turning
        real(dp) :: position(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        integer :: k
        turning_gap = 0.0_dp
        position(1:2) = [s, 0.0_dp]
        do k = 1, n
            position(3) = theta + (turning - theta)*(x(k) + 1.0_dp)/2.0_dp
            call do_magfie(position, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            turning_gap = turning_gap + w(k)*bmod*hder(3)
        end do
        turning_gap = eta_*turning_gap*(turning - theta)/2.0_dp
    end function turning_gap

    real(dp) function turning_point(eta_, direction)
        ! First root of 1 - eta*B from th0 in the given direction.
        real(dp), intent(in) :: eta_, direction
        real(dp) :: inside, outside, middle_
        integer :: k
        inside = th0
        do k = 1, 4096
            outside = th0 + direction*pi*k/4096.0_dp
            if (passing_gap(eta_, outside) < 0.0_dp) exit
            inside = outside
        end do
        if (k > 4096) error stop 'trapped orbit has no turning point'
        do
            middle_ = (inside + outside)/2.0_dp
            if (middle_ == inside .or. middle_ == outside) exit
            if (passing_gap(eta_, middle_) < 0.0_dp) then
                outside = middle_
            else
                inside = middle_
            end if
        end do
        turning_point = inside
    end function turning_point

    real(dp) function passing_gap(eta_, theta)
        real(dp), intent(in) :: eta_, theta
        real(dp) :: position(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        position = [s, 0.0_dp, theta]
        call do_magfie(position, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        passing_gap = 1.0_dp - eta_*bmod
    end function passing_gap
end program test_krook_trapped_weight
