program test_thin_orbit_limit
    !! Closed-form oracle for the strict thin-orbit limiting expression
    !!
    !!   F(lambda) = (Delta phi(lambda) - Delta phi(0))/tau(lambda)
    !!   Omega_t^(1) = F'(0).
    !!
    !! The synthetic return map contains even and odd powers and a varying
    !! period.  It therefore detects a one-sided quotient, omission of the
    !! period denominator, or an incorrect Richardson stencil.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_thin_orbit_limit, only: THIN_LIMIT_SUCCESS, &
        THIN_LIMIT_BASELINE_ERROR, THIN_LIMIT_TOPOLOGY_CHANGE, &
        orbit_return_t, thin_limit_result_t, estimate_thin_limit
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: omega_exact = -3.75_dp
    real(dp), parameter :: topological_delta = 7.5_dp
    real(dp), parameter :: lambda(3) = [0.4_dp, 0.2_dp, 0.1_dp]
    type(orbit_return_t) :: base, plus(3), minus(3)
    type(thin_limit_result_t) :: result
    integer :: k

    base = make_return(0.0_dp)
    do k = 1, size(lambda)
        plus(k) = make_return(lambda(k))
        minus(k) = make_return(-lambda(k))
    end do

    call estimate_thin_limit(topological_delta, base, lambda, plus, minus, result)
    if (result%status /= THIN_LIMIT_SUCCESS) then
        write(*,*) 'thin-limit status:', result%status
        error stop 'closed polynomial thin limit was rejected'
    end if
    if (abs(result%omega - omega_exact) > 2.0e-13_dp) then
        write(*,*) 'thin-limit value/reference:', result%omega, omega_exact
        error stop 'Richardson thin-limit oracle failed'
    end if
    if (abs(result%observed_order - 2.0_dp) > 2.0e-11_dp) then
        write(*,*) 'thin-limit observed order:', result%observed_order
        error stop 'centered thin-limit order is not quadratic'
    end if

    base%delta_phi = base%delta_phi + 1.0e-5_dp
    call estimate_thin_limit(topological_delta, base, lambda, plus, minus, result)
    if (result%status /= THIN_LIMIT_BASELINE_ERROR) then
        error stop 'nonzero zero-width residual was not rejected'
    end if
    base = make_return(0.0_dp)

    plus(2)%winding = plus(2)%winding + 1
    call estimate_thin_limit(topological_delta, base, lambda, plus, minus, result)
    if (result%status /= THIN_LIMIT_TOPOLOGY_CHANGE) then
        error stop 'lambda-ladder topology change was not rejected'
    end if

    call pass_test

contains

    function make_return(scale) result(value)
        real(dp), intent(in) :: scale
        type(orbit_return_t) :: value
        real(dp) :: period, residual_frequency

        period = 2.3_dp + 0.4_dp*scale - 0.2_dp*scale**2
        residual_frequency = omega_exact*scale + 1.1_dp*scale**2 &
            - 0.8_dp*scale**3 + 0.6_dp*scale**4
        value%period = period
        value%delta_phi = topological_delta + period*residual_frequency
        value%orbit_class = 2
        value%winding = 1
        value%status = 0
    end function make_return

end program test_thin_orbit_limit
