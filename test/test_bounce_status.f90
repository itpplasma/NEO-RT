program test_bounce_status
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf
    use iso_fortran_env, only: dp => real64
    use neort_orbit, only: bounce_fast_toleranced, nvar

    implicit none

    real(dp) :: bounceavg(nvar)
    integer :: status

    call bounce_fast_toleranced(1.0_dp, 0.5_dp, 0.0_dp, bounceavg, &
        constant_rhs, status, 1.0e-9_dp, 1.0e-10_dp)
    if (status == 2) error stop "invalid bounce interval accepted"
    if (any(bounceavg /= 0.0_dp)) error stop "invalid bounce output not cleared"

    call bounce_fast_toleranced(ieee_value(0.0_dp, ieee_quiet_nan), 0.5_dp, &
        1.0_dp, bounceavg, constant_rhs, status, 1.0e-9_dp, 1.0e-10_dp)
    if (status == 2) error stop "nonfinite orbit input accepted"
    if (any(bounceavg /= 0.0_dp)) error stop "nonfinite orbit output not cleared"

    call bounce_fast_toleranced(1.0_dp, &
        ieee_value(0.0_dp, ieee_positive_inf), 1.0_dp, bounceavg, &
        constant_rhs, status, 1.0e-9_dp, 1.0e-10_dp)
    if (status == 2) error stop "nonfinite pitch input accepted"
    if (any(bounceavg /= 0.0_dp)) error stop "nonfinite pitch output not cleared"

    call bounce_fast_toleranced(1.0_dp, 0.5_dp, 1.0_dp, bounceavg, &
        constant_rhs, status, ieee_value(0.0_dp, ieee_quiet_nan), 1.0e-10_dp)
    if (status == 2) error stop "nonfinite tolerance accepted"
    if (any(bounceavg /= 0.0_dp)) error stop "nonfinite tolerance output not cleared"
    print *, "test_bounce_status PASSED"

contains

    subroutine constant_rhs(v, eta, neq, t, y, ydot)
        real(dp), intent(in) :: v, eta, t, y(neq)
        integer, intent(in) :: neq
        real(dp), intent(out) :: ydot(neq)

        ydot = 0.0_dp
    end subroutine constant_rhs

end program test_bounce_status
