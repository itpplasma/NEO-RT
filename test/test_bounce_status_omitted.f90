program test_bounce_status_omitted
    use iso_fortran_env, only: dp => real64
    use neort_orbit, only: bounce_fast, nvar

    implicit none

    real(dp) :: bounceavg(nvar)

    call bounce_fast(1.0_dp, 0.5_dp, 0.0_dp, bounceavg, constant_rhs)
    error stop "non-success bounce status was not rejected without an output"

contains

    subroutine constant_rhs(v, eta, neq, t, y, ydot)
        real(dp), intent(in) :: v, eta, t, y(neq)
        integer, intent(in) :: neq
        real(dp), intent(out) :: ydot(neq)

        ydot = 0.0_dp
    end subroutine constant_rhs

end program test_bounce_status_omitted
