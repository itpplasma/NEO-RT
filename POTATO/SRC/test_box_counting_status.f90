program test_box_counting_status
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use box_counting_status_mod, only: box_count_ok, &
        box_count_nonpositive_bounce, box_count_vode_mxstep, &
        box_count_vode_failed, box_count_nonfinite, &
        box_count_negative_residence, box_count_residence_mismatch, &
        box_count_vode_status, box_count_residence_status
    implicit none

    real(8) :: tau(3), nan_value

    ! Independent arithmetic oracle: 1/4 + 1/2 + 1/8 inside plus 1/8 outside
    ! is exactly one complete bounce.
    tau = [0.25d0, 0.5d0, 0.125d0]
    if (box_count_residence_status(1d0, tau, 0.125d0) /= box_count_ok) &
        error stop "complete residence rejected"

    if (box_count_residence_status(-1d0, tau, 0d0) /= &
        box_count_nonpositive_bounce) &
        error stop "nonpositive bounce accepted"

    tau = [0.25d0, -0.5d0, 1.25d0]
    if (box_count_residence_status(1d0, tau, 0d0) /= &
        box_count_negative_residence) &
        error stop "negative residence accepted"

    tau = [0.2d0, 0.2d0, 0.2d0]
    if (box_count_residence_status(1d0, tau, 0d0) /= &
        box_count_residence_mismatch) &
        error stop "incomplete residence accepted"

    nan_value = ieee_value(0d0, ieee_quiet_nan)
    tau = [0.25d0, nan_value, 0.75d0]
    if (box_count_residence_status(1d0, tau, 0d0) /= box_count_nonfinite) &
        error stop "nonfinite residence accepted"

    if (box_count_vode_status(2) /= box_count_ok) &
        error stop "successful VODE state rejected"
    if (box_count_vode_status(-1) /= box_count_vode_mxstep) &
        error stop "VODE mxstep state misclassified"
    if (box_count_vode_status(-3) /= box_count_vode_failed) &
        error stop "generic VODE failure state accepted"
end program test_box_counting_status
