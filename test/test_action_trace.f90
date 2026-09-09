program test_action_trace
    use iso_fortran_env, only: dp => real64
    use neort_action_trace_contract, only: action_trace_weight, &
        complex_orbit_action, action_modulus_squared
    implicit none

    real(dp), parameter :: expected_signed = -0.5_dp
    real(dp), parameter :: expected_absolute_jacobian = 0.5_dp
    real(dp) :: actual
    complex(dp) :: action, rotated

    action = complex_orbit_action(3.0_dp, -4.0_dp, 2.0_dp)
    if (abs(action - cmplx(6.0_dp, -8.0_dp, dp)) > 1.0e-14_dp) then
        error stop "complex action does not preserve amplitude and signed phase"
    end if
    if (abs(action_modulus_squared(action) - 100.0_dp) > 1.0e-14_dp) then
        error stop "complex action norm does not reproduce the independent oracle"
    end if
    rotated = cmplx(0.0_dp, 1.0_dp, dp)*action
    if (abs(rotated - cmplx(8.0_dp, 6.0_dp, dp)) > 1.0e-14_dp) then
        error stop "complex action does not retain a unit phase rotation"
    end if
    if (abs(action_modulus_squared(rotated) - &
            action_modulus_squared(action)) > 1.0e-14_dp) then
        error stop "action norm changed under a unit phase rotation"
    end if

    actual = action_trace_weight(0.25_dp, -8.0_dp, 2.0_dp, 0.5_dp)
    if (abs(actual - expected_signed) > 1.0e-14_dp) then
        error stop "action trace weight does not preserve native sign"
    end if

    actual = action_trace_weight(0.25_dp, 8.0_dp, -2.0_dp, 0.5_dp)
    if (abs(actual - expected_absolute_jacobian) > 1.0e-14_dp) then
        error stop "action trace weight does not use absolute Jacobian"
    end if

    print *, "test_action_trace ... OK"
end program test_action_trace
