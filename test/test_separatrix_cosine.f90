program test_separatrix_cosine
    ! Independent single-well oracle: B(theta)=1-e*cos(theta), v=h^theta=1.
    ! theta/2 reduces the transit integral to the complete elliptic K;
    ! sin(theta/2)=sqrt(m)*sin(u) does the same for a trapped bounce.
    ! K is evaluated by the arithmetic-geometric mean, not the orbit ODE or fit.
    ! A bounded weight B(theta) also checks the numerator/period ratio, using
    ! a separate Gauss-Legendre integral for elliptic E.
    use iso_fortran_env, only: dp => real64
    use neort_separatrix, only: separatrix_model_t, fit_separatrix_model, &
        evaluate_separatrix_model
    use fortnum_quadrature, only: gauss_legendre_ab
    use util, only: pi
    implicit none
    real(dp), parameter :: modulation(3) = [0.05_dp, 0.3_dp, 0.7_dp]
    real(dp) :: x(16), periods(16), drifts(16), nodes(256), weights(256)
    real(dp) :: tau, dtau, omega, domega, drift, ddrift, tau_ref, drift_ref, distance
    real(dp) :: a_pass, a_trap
    type(separatrix_model_t) :: model
    integer :: ecase, branch, k
    logical :: ok

    call gauss_legendre_ab(size(nodes), 0.0_dp, pi/2, nodes, weights)
    do ecase = 1, size(modulation)
        do branch = -1, 1, 2
            do k = 1, size(x)
                x(k) = 1.0e-6_dp*exp(real(k - 1, dp)*log(100.0_dp)/(size(x)-1))
                call exact_integrals(x(k), modulation(ecase), branch, periods(k), drifts(k))
            end do
            call fit_separatrix_model(x, periods, drifts, model, ok)
            if (.not. ok) error stop 'cosine-well fit failed'
            if (branch < 0) then
                a_pass = model%tau(1)
            else
                a_trap = model%tau(1)
            end if
            if (abs(model%numerator(1)/model%tau(1)/(1+modulation(ecase))-1) &
                    > 2.0e-6_dp) error stop 'cosine barrier-average limit'
            ! Extrapolate four to nine decades below the orbit samples.
            do k = 10, 15
                distance = 10.0_dp**(-k)
                call exact_integrals(distance, modulation(ecase), branch, tau_ref, drift_ref)
                call evaluate_separatrix_model(model, distance, tau, dtau, omega, &
                    domega, drift, ddrift, ok)
                if (.not. ok) error stop 'cosine extrapolation failed'
                if (abs(tau/tau_ref-1) > 2.0e-6_dp) error stop 'cosine elliptic period'
                if (abs(omega*tau_ref/(2*pi)-1) > 2.0e-6_dp) &
                    error stop 'cosine elliptic frequency'
                if (abs(drift/drift_ref-1) > 2.0e-6_dp) &
                    error stop 'cosine bounded drift average'
                if (dtau >= 0 .or. domega <= 0) error stop 'cosine derivative signs'
            end do
        end do
        if (abs(a_trap/(2*a_pass)-1) > 2.0e-6_dp) &
            error stop 'cosine barrier passage multiplicity'
    end do
    print *, 'test_separatrix_cosine PASSED'
contains
    subroutine exact_integrals(distance, e, branch, tau, drift)
        real(dp), intent(in) :: distance, e
        integer, intent(in) :: branch
        real(dp), intent(out) :: tau, drift
        real(dp) :: eta, a, b, complementary, elliptic_k, elliptic_e, gap_average
        real(dp) :: arithmetic, geometric, next_a
        integer :: j
        eta = (1 + branch*distance)/(1+e)
        a = 1 - eta*(1-e)
        b = 2*eta*e
        ! Compute 1-m directly from distance, avoiding cancellation near m=1.
        if (branch < 0) then
            complementary = distance/a
        else
            complementary = distance/b
        end if
        arithmetic = 1
        geometric = sqrt(complementary)
        do j = 1, 30
            next_a = (arithmetic + geometric)/2
            geometric = sqrt(arithmetic*geometric)
            arithmetic = next_a
            if (abs(arithmetic-geometric) < epsilon(arithmetic)*arithmetic) exit
        end do
        elliptic_k = pi/(2*arithmetic)
        elliptic_e = sum(weights*sqrt(cos(nodes)**2 + complementary*sin(nodes)**2))
        if (branch < 0) then
            tau = 4*elliptic_k/sqrt(a)
            gap_average = a*elliptic_e/elliptic_k
        else
            tau = 8*elliptic_k/sqrt(b)
            gap_average = b*(elliptic_e/elliptic_k-complementary)
        end if
        drift = (1-gap_average)/eta
    end subroutine exact_integrals
end program test_separatrix_cosine
