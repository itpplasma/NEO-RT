module acceptance_fixture
    use iso_fortran_env, only: dp => real64
    implicit none
    integer :: mode = 1
contains
    subroutine nonroot(v, eta, res, derivative)
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: res, derivative
        if (v <= 0) error stop "invalid test velocity"
        if (mode == 1) then
            ! A discontinuity between adjacent floats has a sign change but no zero.
            res = merge(-1.0_dp, 1.0_dp, eta <= 1.0_dp)
            derivative = 1.0_dp
        else if (mode == 3) then
            ! A continuous root half way between floats has a nonzero rounding
            ! residual even though its pitch is known to one ulp.
            derivative = 1.0e16_dp
            res = derivative*((eta - 1.0_dp) - 1.5_dp*epsilon(eta))
        else
            res = 0.0_dp
            derivative = 0.0_dp
        end if
    end subroutine nonroot
end module acceptance_fixture

program test_resonance_acceptance
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_next_after
    use neort_resonance, only: driftorbit_root, resonance_override
    use acceptance_fixture, only: mode, nonroot
    use driftorbit, only: vth
    use logger, only: set_log_level
    implicit none
    real(dp) :: root(2), upper
    call set_log_level(-1)
    vth = 1
    resonance_override => nonroot
    upper = ieee_next_after(1.0_dp, 2.0_dp)
    root = driftorbit_root(1.0_dp, 1.0e-12_dp, 1.0_dp, upper)
    if (root(1) >= 0 .or. root(2) /= 0) &
        error stop "machine-width nonroot accepted"
    mode = 2
    root = driftorbit_root(1.0_dp, 1.0e-12_dp, 1.0_dp, upper)
    if (root(1) >= 0 .or. root(2) /= 0) &
        error stop "zero-Jacobian endpoint accepted"
    mode = 3
    upper = ieee_next_after(upper, 2.0_dp)
    root = driftorbit_root(1.0_dp, 1.0e-12_dp, 1.0_dp, upper)
    if (root(1) < 1.0_dp .or. root(1) > upper) &
        error stop "continuous machine-width root rejected"
    if (abs(root(1) - (1.0_dp + epsilon(upper))) > epsilon(upper)) &
        error stop "continuous root exceeds pitch rounding error"
    if (root(2) /= 1.0e16_dp) error stop "continuous root Jacobian lost"
    nullify(resonance_override)
    print *, "test_resonance_acceptance PASSED"
end program test_resonance_acceptance
