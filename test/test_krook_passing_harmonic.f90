program test_krook_passing_harmonic
    ! Passing Hamiltonian harmonics against an angular reference on the closed
    ! orbit. The caller's transit time comes from the interpolated frequency
    ! Om_th, which differs from the orbit period by the spline error (about
    ! 5e-7 relative for ITER C12 at m = -11). The harmonic
    ! exp(-i*(mth + q*mph)*Omth*t) is a strongly cancelling average; integrated
    ! over that span with that frequency it is wrong by about its own
    ! size, and its noise in eta stopped the Krook pitch quadrature.
    ! The orbit integration controls the error relative to the integrand
    ! scale max|H| (transport_bounce), so that is the scale of the gate: a
    ! supplied period bias of 5e-7 gave errors of 4e-9*max|H| before.
    use iso_fortran_env, only: dp => real64
    use fortnum_quadrature, only: gauss_legendre
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_transport, only: transport_bounce, Omth
    use neort_orbit, only: nvar, th0
    use do_magfie_mod, only: do_magfie, s
    use driftorbit, only: etatp, vth, mth, mph, m0, epsmn, sign_vpar, &
        nopassing, comptorque, drift_phase, fow_radial, Bmin
    implicit none
    integer, parameter :: n = 16
    real(dp), parameter :: pi = acos(-1.0_dp)
    ! Relative errors of the supplied transit time; zero is the exact period.
    real(dp), parameter :: biases(3) = [0.0_dp, 5.0e-7_dp, -1.0e-6_dp]
    real(dp), parameter :: depths(3) = [0.5_dp, 0.97_dp, 0.999_dp]
    real(dp), parameter :: tolerance = 1.0e-10_dp
    real(dp) :: x(n), w(n), cumulative(n, n)
    real(dp) :: eta, period, coarse_period, average(nvar), error, worst
    complex(dp) :: coarse, fine
    integer :: k, l, direction, kstate

    call gauss_legendre(n, x, w)
    call integration_matrix()
    call neort_init('driftorbit.in', 'in_file')
    nopassing = .false.
    comptorque = .false.
    drift_phase = .false.
    fow_radial = .false.
    call neort_prepare_splines('plasma.in', 'profile.in')
    call neort_setup_at_s(0.5_dp)
    ! H = (2 - eta*B)*exp(-i*mth*Omth*t), a cancelling transit harmonic.
    mth = 3
    mph = 0
    m0 = 0
    epsmn = 1.0_dp
    worst = 0.0_dp
    do k = 1, size(depths)
        eta = etatp*depths(k)
        do direction = -1, 1, 2
            sign_vpar = real(direction, dp)
            call orbit_reference(eta, 256, coarse_period, coarse)
            call orbit_reference(eta, 512, period, fine)
            if (abs(period - coarse_period) > 1.0e-10_dp*period .or. &
                abs(fine - coarse) > 1.0e-9_dp*abs(fine)) &
                error stop 'passing angular reference has not converged'
            do l = 1, size(biases)
                Omth = sign_vpar*2.0_dp*pi/(period*(1.0_dp + biases(l)))
                call transport_bounce(vth, eta, 2.0_dp*pi/abs(Omth), average, kstate)
                if (kstate /= 2) error stop 'passing transport bounce failed'
                ! max|H| = 2 - eta*Bmin.
                error = abs(cmplx(average(3), average(4), dp) - fine)/ &
                    (2.0_dp - eta*Bmin)
                worst = max(worst, error)
                print '(a,es10.2,i3,es10.2,a,es12.4,a,es10.2,a,es10.2)', &
                    ' eta/etatp, direction, bias:', depths(k), direction, biases(l), &
                    ' |harmonic|', abs(fine), ' error/max|H|', error, &
                    ' relative', error*(2.0_dp - eta*Bmin)/abs(fine)
                if (error > tolerance) &
                    error stop 'passing harmonic disagrees with closed-orbit reference'
            end do
        end do
    end do
    print *, 'passing harmonic maximum error/max|H|:', worst
    print *, 'test_krook_passing_harmonic PASSED'
contains
    subroutine integration_matrix()
        ! cumulative(i, j) = integral of the j-th Lagrange basis polynomial
        ! of the nodes x from -1 to x(i).
        integer :: i, j, p
        real(dp) :: xi
        cumulative = 0.0_dp
        do i = 1, n
            do p = 1, n
                xi = -1.0_dp + (x(i) + 1.0_dp)*(x(p) + 1.0_dp)/2.0_dp
                do j = 1, n
                    cumulative(i, j) = cumulative(i, j) + &
                        (x(i) + 1.0_dp)/2.0_dp*w(p)*lagrange(j, xi)
                end do
            end do
        end do
    end subroutine integration_matrix

    real(dp) function lagrange(j, xi)
        integer, intent(in) :: j
        real(dp), intent(in) :: xi
        integer :: p
        lagrange = 1.0_dp
        do p = 1, n
            if (p /= j) lagrange = lagrange*(xi - x(p))/(x(j) - x(p))
        end do
    end function lagrange

    subroutine orbit_reference(eta_, panels, period_, harmonic)
        ! One transit from theta = th0 in the direction sign_vpar, with
        ! t(theta) = integral of dtheta/(v*sqrt(1 - eta*B)*|h^theta|) and the
        ! harmonic frequency 2*pi/period of this closed orbit.
        real(dp), intent(in) :: eta_
        integer, intent(in) :: panels
        real(dp), intent(out) :: period_
        complex(dp), intent(out) :: harmonic
        real(dp) :: h, a, t, omega, dtdphi(n), hamiltonian(n), times(n)
        integer :: panel, pass
        h = 2.0_dp*pi/real(panels, dp)
        omega = 0.0_dp
        do pass = 1, 2
            t = 0.0_dp
            harmonic = (0.0_dp, 0.0_dp)
            do panel = 1, panels
                a = h*real(panel - 1, dp)
                call panel_values(eta_, a, h, dtdphi, hamiltonian)
                times = t + h/2.0_dp*matmul(cumulative, dtdphi)
                harmonic = harmonic + h/2.0_dp*sum(w*dtdphi*hamiltonian* &
                    exp(-cmplx(0.0_dp, 1.0_dp, dp)*mth*omega*times))
                t = t + h/2.0_dp*sum(w*dtdphi)
            end do
            period_ = t
            omega = sign_vpar*2.0_dp*pi/period_
        end do
        harmonic = harmonic/period_
    end subroutine orbit_reference

    subroutine panel_values(eta_, a, h, dtdphi, hamiltonian)
        real(dp), intent(in) :: eta_, a, h
        real(dp), intent(out) :: dtdphi(n), hamiltonian(n)
        real(dp) :: position(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        integer :: k
        position(1:2) = [s, 0.0_dp]
        do k = 1, n
            position(3) = th0 + sign_vpar*(a + h*(x(k) + 1.0_dp)/2.0_dp)
            call do_magfie(position, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            dtdphi(k) = 1.0_dp/(vth*sqrt(1.0_dp - eta_*bmod)*abs(hctrvr(3)))
            hamiltonian(k) = 2.0_dp - eta_*bmod
        end do
    end subroutine panel_values
end program test_krook_passing_harmonic
