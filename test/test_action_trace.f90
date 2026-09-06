program test_action_trace
    use iso_fortran_env, only: dp => real64
    use neort_action_trace_contract, only: action_trace_weight
    implicit none

    real(dp), parameter :: expected_signed = -0.5_dp
    real(dp), parameter :: expected_absolute_jacobian = 0.5_dp
    real(dp) :: actual

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
