program test_gc_variational_event
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_variational_event_diagnostic_symbolic, only: &
        evaluate_neort_gc_variational_event_diagnostic
    use neort_gc_variational_event_symbolic, only: &
        evaluate_neort_gc_variational_event
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: parameter_value = 0.35_dp
    real(dp), parameter :: finite_difference_step = 1.0e-4_dp
    real(dp), parameter :: transversality_floor = 1.0e-6_dp
    real(dp), parameter :: tolerance = 2.0e-11_dp
    real(dp) :: time_event, time_derivative, observable_derivative
    real(dp) :: f_x, f_y, c_x, c_y, c_t, c_lambda
    real(dp) :: s_x, s_y, y_x, y_y, y_t, y_lambda_explicit
    real(dp) :: d, event_numerator, tau_lambda, y_explicit, y_transport
    real(dp) :: y_lambda, d_reverse, numerator_reverse, tau_reverse
    real(dp) :: y_lambda_reverse, tau_residual, y_residual, abs_d, d_squared
    real(dp) :: diagnostic_d, diagnostic_numerator
    real(dp) :: diagnostic_y_explicit, diagnostic_y_transport
    real(dp) :: diagnostic_d_reverse, diagnostic_numerator_reverse
    real(dp) :: diagnostic_abs_d, diagnostic_d_squared
    real(dp) :: diagnostic_margin, diagnostic_margin_reverse
    real(dp) :: exact_d, exact_numerator
    real(dp) :: near_d
    real(dp) :: c_x_zero, c_y_zero, c_t_zero, c_lambda_zero

    time_event = exact_event_time(parameter_value)
    f_x = 1.0_dp + parameter_value
    f_y = 2.0_dp*time_event*(1.0_dp+parameter_value) + parameter_value
    c_x = 1.0_dp
    c_y = 0.5_dp
    c_t = 0.2_dp
    c_lambda = 0.1_dp
    s_x = time_event
    s_y = time_event**2 + time_event
    y_x = 0.0_dp
    y_y = 1.0_dp
    y_t = 3.0_dp
    y_lambda_explicit = 2.0_dp

    call evaluate_neort_gc_variational_event_diagnostic(f_x, f_y, c_x, c_y, &
        c_t, c_lambda, s_x, s_y, y_x, y_y, y_t, y_lambda_explicit, &
        transversality_floor, diagnostic_d, diagnostic_numerator, &
        diagnostic_y_explicit, diagnostic_y_transport, diagnostic_d_reverse, &
        diagnostic_numerator_reverse, diagnostic_abs_d, diagnostic_d_squared, &
        diagnostic_margin, diagnostic_margin_reverse)
    if (diagnostic_margin <= 0.0_dp) then
        error stop 'regular event was rejected by generated margin'
    end if
    call require_close('regular diagnostic margin reversal', &
        diagnostic_margin_reverse, diagnostic_margin)

    call evaluate_neort_gc_variational_event(f_x, f_y, c_x, c_y, c_t, &
        c_lambda, s_x, s_y, y_x, y_y, y_t, y_lambda_explicit, d, &
        event_numerator, tau_lambda, y_explicit, y_transport, y_lambda, &
        d_reverse, numerator_reverse, tau_reverse, y_lambda_reverse, &
        tau_residual, y_residual, abs_d, d_squared)
    exact_d = f_x+0.5_dp*f_y+0.2_dp
    exact_numerator = s_x+0.5_dp*s_y+0.1_dp
    time_derivative = exact_event_time_derivative(parameter_value)
    observable_derivative = exact_observable_derivative(parameter_value)

    call require_close('event transversality', d, exact_d)
    call require_close('event numerator', event_numerator, exact_numerator)
    call require_close('event time derivative', tau_lambda, time_derivative)
    call require_close('explicit Y derivative', y_explicit, 2.0_dp)
    call require_close('terminal transport derivative', y_transport, f_y+3.0_dp)
    call require_close('event-corrected Y derivative', y_lambda, &
        observable_derivative)
    call require_close('D reversal', d_reverse, -d)
    call require_close('numerator reversal', numerator_reverse, -event_numerator)
    call require_close('tau reversal', tau_reverse, tau_lambda)
    call require_close('Y reversal', y_lambda_reverse, y_lambda)
    call require_close('tau identity residual', tau_residual, 0.0_dp)
    call require_close('Y identity residual', y_residual, 0.0_dp)
    call require_close('absolute D', abs_d, abs(d))
    call require_close('squared D', d_squared, d*d)

    call require_close('tau high-quality finite difference', tau_lambda, &
        five_point_time_derivative(parameter_value, finite_difference_step))
    call require_close('Y high-quality finite difference', y_lambda, &
        five_point_observable_derivative(parameter_value, &
        finite_difference_step))

    c_x_zero = 1.4_dp
    c_y_zero = -1.0_dp
    c_t_zero = -1.16_dp
    c_lambda_zero = 0.3_dp
    call evaluate_neort_gc_variational_event_diagnostic(1.4_dp, 0.8_dp, &
        c_x_zero, c_y_zero, c_t_zero, c_lambda_zero, 0.7_dp, -0.2_dp, &
        0.4_dp, 1.0_dp, 2.0_dp, 0.5_dp, transversality_floor, &
        diagnostic_d, diagnostic_numerator, diagnostic_y_explicit, &
        diagnostic_y_transport, diagnostic_d_reverse, &
        diagnostic_numerator_reverse, diagnostic_abs_d, diagnostic_d_squared, &
        diagnostic_margin, diagnostic_margin_reverse)
    call require_absolute_close('exact grazing D exposed', diagnostic_d, &
        0.0_dp, 1.0e-14_dp)
    call require_close('exact grazing numerator exposed', &
        diagnostic_numerator, 1.48_dp)
    call require_absolute_close('exact grazing D magnitude', diagnostic_abs_d, &
        0.0_dp, 1.0e-14_dp)
    call require_absolute_close('exact grazing D square', diagnostic_d_squared, &
        0.0_dp, 1.0e-28_dp)
    if (diagnostic_margin >= 0.0_dp) then
        error stop 'exact grazing event passed generated margin'
    end if
    call require_close('exact grazing margin reversal', &
        diagnostic_margin_reverse, diagnostic_margin)
    call require_absolute_close('exact grazing reversed D', &
        diagnostic_d_reverse, -diagnostic_d, 1.0e-14_dp)
    call require_close('exact grazing reversed numerator', &
        diagnostic_numerator_reverse, -diagnostic_numerator)
    call require_finite('exact grazing D', diagnostic_d)
    call require_finite('exact grazing numerator', diagnostic_numerator)
    call require_finite('exact grazing explicit Y', diagnostic_y_explicit)
    call require_finite('exact grazing transport Y', diagnostic_y_transport)
    call require_finite('exact grazing reversed D', diagnostic_d_reverse)
    call require_finite('exact grazing reversed numerator', &
        diagnostic_numerator_reverse)
    call require_finite('exact grazing absolute D', diagnostic_abs_d)
    call require_finite('exact grazing squared D', diagnostic_d_squared)
    call require_finite('exact grazing margin', diagnostic_margin)
    call require_finite('exact grazing reversed margin', &
        diagnostic_margin_reverse)

    call evaluate_neort_gc_variational_event_diagnostic(0.0_dp, 0.0_dp, &
        1.0_dp, 0.0_dp, transversality_floor, 0.3_dp, 0.7_dp, -0.2_dp, &
        0.4_dp, 1.0_dp, 2.0_dp, 0.5_dp, transversality_floor, &
        diagnostic_d, diagnostic_numerator, diagnostic_y_explicit, &
        diagnostic_y_transport, diagnostic_d_reverse, &
        diagnostic_numerator_reverse, diagnostic_abs_d, diagnostic_d_squared, &
        diagnostic_margin, diagnostic_margin_reverse)
    call require_close('margin boundary D', diagnostic_d, &
        transversality_floor)
    call require_absolute_close('margin equality at absolute D floor', &
        diagnostic_margin, 0.0_dp, 1.0e-30_dp)
    call require_absolute_close('reversed margin equality', &
        diagnostic_margin_reverse, 0.0_dp, 1.0e-30_dp)

    near_d = 1.0e-10_dp
    call evaluate_neort_gc_variational_event_diagnostic(1.4_dp, 0.8_dp, &
        c_x_zero, c_y_zero, &
        -(c_x_zero*1.4_dp+c_y_zero*0.8_dp)+near_d, c_lambda_zero, &
        0.7_dp, -0.2_dp, 0.4_dp, 1.0_dp, 2.0_dp, 0.5_dp, &
        transversality_floor, diagnostic_d, diagnostic_numerator, &
        diagnostic_y_explicit, diagnostic_y_transport, diagnostic_d_reverse, &
        diagnostic_numerator_reverse, diagnostic_abs_d, diagnostic_d_squared, &
        diagnostic_margin, diagnostic_margin_reverse)
    call require_absolute_close('near-grazing D exposed', diagnostic_d, &
        near_d, 1.0e-14_dp)
    call require_close('near-grazing numerator exposed', &
        diagnostic_numerator, 1.48_dp)
    call require_absolute_close('near-grazing magnitude exposed', &
        diagnostic_abs_d, abs(near_d), 1.0e-14_dp)
    call require_absolute_close('near-grazing margin value', &
        diagnostic_margin, near_d**2-transversality_floor**2, 1.0e-18_dp)
    if (diagnostic_margin >= 0.0_dp) then
        error stop 'near-grazing event passed generated margin'
    end if
    call require_close('near-grazing margin reversal', &
        diagnostic_margin_reverse, diagnostic_margin)
    call require_absolute_close('near-grazing reversed D', &
        diagnostic_d_reverse, -diagnostic_d, 1.0e-14_dp)
    call require_close('near-grazing reversed numerator', &
        diagnostic_numerator_reverse, -diagnostic_numerator)

    write (*, '(A)') 'test_gc_variational_event OK'
    call pass_test

contains

    pure real(dp) function exact_event_time(lambda_value)
        real(dp), intent(in) :: lambda_value
        real(dp) :: coefficient_a, coefficient_b, coefficient_c
        real(dp) :: discriminant

        coefficient_a = 0.5_dp*(1.0_dp+lambda_value)
        coefficient_b = 1.2_dp+1.5_dp*lambda_value
        coefficient_c = 0.1_dp*lambda_value-1.0_dp
        discriminant = coefficient_b**2-4.0_dp*coefficient_a*coefficient_c
        exact_event_time = (-coefficient_b+sqrt(discriminant))/ &
            (2.0_dp*coefficient_a)
    end function exact_event_time

    pure real(dp) function exact_event_time_derivative(lambda_value)
        real(dp), intent(in) :: lambda_value
        real(dp) :: t_value, f_x_value, f_y_value

        t_value = exact_event_time(lambda_value)
        f_x_value = 1.0_dp+lambda_value
        f_y_value = 2.0_dp*t_value*f_x_value+lambda_value
        exact_event_time_derivative = -(t_value+0.5_dp*(t_value**2+t_value)+ &
            0.1_dp)/(f_x_value+0.5_dp*f_y_value+0.2_dp)
    end function exact_event_time_derivative

    pure real(dp) function exact_observable(lambda_value)
        real(dp), intent(in) :: lambda_value
        real(dp) :: t_value

        t_value = exact_event_time(lambda_value)
        exact_observable = (1.0_dp+lambda_value)*t_value**2+ &
            lambda_value*t_value+3.0_dp*t_value+2.0_dp*lambda_value
    end function exact_observable

    pure real(dp) function exact_observable_derivative(lambda_value)
        real(dp), intent(in) :: lambda_value
        real(dp) :: t_value, t_derivative

        t_value = exact_event_time(lambda_value)
        t_derivative = exact_event_time_derivative(lambda_value)
        exact_observable_derivative = t_value**2+t_value+2.0_dp+ &
            (2.0_dp*t_value*(1.0_dp+lambda_value)+lambda_value+3.0_dp)* &
            t_derivative
    end function exact_observable_derivative

    pure real(dp) function five_point_time_derivative(lambda_value, step)
        real(dp), intent(in) :: lambda_value, step

        five_point_time_derivative = (exact_event_time(lambda_value-2.0_dp*step) &
            -8.0_dp*exact_event_time(lambda_value-step) + &
            8.0_dp*exact_event_time(lambda_value+step) - &
            exact_event_time(lambda_value+2.0_dp*step))/(12.0_dp*step)
    end function five_point_time_derivative

    pure real(dp) function five_point_observable_derivative(lambda_value, step)
        real(dp), intent(in) :: lambda_value, step

        five_point_observable_derivative = &
            (exact_observable(lambda_value-2.0_dp*step) - &
            8.0_dp*exact_observable(lambda_value-step) + &
            8.0_dp*exact_observable(lambda_value+step) - &
            exact_observable(lambda_value+2.0_dp*step))/(12.0_dp*step)
    end function five_point_observable_derivative

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual-expected) > tolerance*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'variational event oracle failed'
        end if
    end subroutine require_close

    subroutine require_absolute_close(label, actual, expected, absolute_tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected, absolute_tolerance

        if (abs(actual-expected) > absolute_tolerance) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'variational event absolute oracle failed'
        end if
    end subroutine require_absolute_close

    subroutine require_finite(label, value)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value

        if (.not. ieee_is_finite(value)) then
            write (*, '(A,1X,ES24.16)') trim(label), value
            error stop 'variational event diagnostic is not finite'
        end if
    end subroutine require_finite

end program test_gc_variational_event
