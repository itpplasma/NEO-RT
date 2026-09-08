program test_bounce_status
    use iso_fortran_env, only: dp => real64
    use neort_orbit, only: bounce_fast_toleranced, nvar

    implicit none

    real(dp) :: bounceavg(nvar)
    integer :: status

    call bounce_fast_toleranced(1.0_dp, 0.5_dp, 0.0_dp, bounceavg, &
        constant_rhs, status, 1.0e-9_dp, 1.0e-10_dp)
    if (status == 2) error stop "invalid bounce interval accepted"
    if (any(bounceavg /= 0.0_dp)) error stop "invalid bounce output not cleared"
    print *, "test_bounce_status PASSED"

contains

    subroutine constant_rhs(v, eta, neq, t, y, ydot)
        real(dp), intent(in) :: v, eta, t, y(neq)
        integer, intent(in) :: neq
        real(dp), intent(out) :: ydot(neq)

        ydot = 0.0_dp
    end subroutine constant_rhs

end program test_bounce_status
