program test_resonance_jacobian
    use iso_fortran_env, only: dp => real64, int64
    use ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan, ieee_positive_inf
    use neort_resonance, only: valid_resonance_jacobian

    implicit none

    real(dp) :: nan_value, inf_value, reciprocal_overflow
    integer(int64) :: subnormal_bits

    nan_value = ieee_value(0.0_dp, ieee_quiet_nan)
    inf_value = ieee_value(0.0_dp, ieee_positive_inf)
    ! Construct the least positive subnormal by its IEEE bit pattern.  An
    ! arithmetic tiny()/4 expression would itself trip the underflow trap in
    ! Debug builds before the predicate can reject it.
    subnormal_bits = 1_int64
    reciprocal_overflow = transfer(subnormal_bits, reciprocal_overflow)

    if (.not. valid_resonance_jacobian(1.0_dp)) error stop "positive Jacobian rejected"
    if (.not. valid_resonance_jacobian(-1.0_dp)) error stop "negative Jacobian rejected"
    if (valid_resonance_jacobian(0.0_dp)) error stop "zero Jacobian accepted"
    if (valid_resonance_jacobian(-0.0_dp)) error stop "negative zero Jacobian accepted"
    if (valid_resonance_jacobian(nan_value)) error stop "NaN Jacobian accepted"
    if (valid_resonance_jacobian(inf_value)) error stop "infinite Jacobian accepted"
    if (valid_resonance_jacobian(reciprocal_overflow)) error stop "overflowing Jacobian reciprocal accepted"
    if (.not. ieee_is_finite(nan_value)) then
        print *, "test_resonance_jacobian PASSED"
    else
        error stop "NaN fixture is not nonfinite"
    end if
end program test_resonance_jacobian
