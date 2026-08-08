program test_gc_homoclinic_limit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_homoclinic_coefficients_symbolic, only: &
        evaluate_neort_gc_homoclinic_coefficients
    use neort_gc_homoclinic_diagnostic_symbolic, only: &
        evaluate_neort_gc_homoclinic_diagnostic
    use neort_gc_homoclinic_finite_part_symbolic, only: &
        evaluate_neort_gc_homoclinic_finite_part
    use util_for_test, only: pass_test
    implicit none

    integer, parameter :: NUMBER_OF_EPSILONS = 4
    integer, parameter :: TIME_INTEGRAND = 1
    integer, parameter :: PHI_INTEGRAND = 2
    integer, parameter :: REVERSED_PHI_INTEGRAND = 3
    real(dp), parameter :: SADDLE_RATE = 1.7_dp
    real(dp), parameter :: TOROIDAL_RATE = -0.7_dp
    real(dp), parameter :: LOCAL_DOMAIN = 1.3_dp
    real(dp), parameter :: SADDLE_RATE_FLOOR = 1.0e-4_dp
    real(dp), parameter :: EPSILON_FLOOR = 1.0e-8_dp
    real(dp), parameter :: QUADRATURE_TOLERANCE = 2.0e-12_dp
    real(dp), parameter :: EPSILONS(NUMBER_OF_EPSILONS) = &
        [3.0e-2_dp, 1.0e-2_dp, 3.0e-3_dp, 1.0e-3_dp]
    real(dp) :: hessian_determinant, saddle_discriminant
    real(dp) :: saddle_rate_value, c_tau, c_phi
    real(dp) :: one_leg_tau_coefficient, two_leg_tau_coefficient
    real(dp) :: one_leg_phi_coefficient, two_leg_phi_coefficient
    real(dp) :: c_tau_section_reverse, c_phi_section_reverse
    real(dp) :: c_tau_toroidal_reverse, c_phi_toroidal_reverse
    real(dp) :: absolute_phi_reversal_defect
    real(dp) :: diagnostic_hessian_determinant
    real(dp) :: diagnostic_saddle_discriminant
    real(dp) :: nondegeneracy_margin, section_reversed_margin
    real(dp) :: epsilon_margin, section_reversed_epsilon_margin
    real(dp) :: saddle_rate_floor_margin, epsilon_floor_margin
    real(dp) :: tau_one(NUMBER_OF_EPSILONS)
    real(dp) :: tau_two(NUMBER_OF_EPSILONS)
    real(dp) :: phi_one(NUMBER_OF_EPSILONS)
    real(dp) :: phi_two(NUMBER_OF_EPSILONS)
    real(dp) :: phi_one_reversed(NUMBER_OF_EPSILONS)
    real(dp) :: tau_finite_one(NUMBER_OF_EPSILONS)
    real(dp) :: tau_finite_two(NUMBER_OF_EPSILONS)
    real(dp) :: phi_finite_one(NUMBER_OF_EPSILONS)
    real(dp) :: phi_finite_two(NUMBER_OF_EPSILONS)
    real(dp) :: phi_finite_one_reverse(NUMBER_OF_EPSILONS)
    real(dp) :: phi_finite_two_reverse(NUMBER_OF_EPSILONS)
    real(dp) :: c_tau_one_estimate(NUMBER_OF_EPSILONS-1)
    real(dp) :: c_tau_two_estimate(NUMBER_OF_EPSILONS-1)
    real(dp) :: c_phi_one_estimate(NUMBER_OF_EPSILONS-1)
    real(dp) :: c_phi_two_estimate(NUMBER_OF_EPSILONS-1)
    real(dp) :: c_phi_reverse_estimate(NUMBER_OF_EPSILONS-1)
    real(dp) :: logarithmic_interval
    integer :: i

    call evaluate_neort_gc_homoclinic_diagnostic(-SADDLE_RATE**2, 0.0_dp, &
        1.0_dp, SADDLE_RATE_FLOOR, EPSILONS(1), EPSILON_FLOOR, &
        diagnostic_hessian_determinant, diagnostic_saddle_discriminant, &
        nondegeneracy_margin, section_reversed_margin, epsilon_margin, &
        section_reversed_epsilon_margin, saddle_rate_floor_margin, &
        epsilon_floor_margin)
    if (nondegeneracy_margin <= 0.0_dp) then
        error stop 'manufactured saddle failed generated nondegeneracy gate'
    end if
    if (epsilon_margin <= 0.0_dp) then
        error stop 'manufactured epsilon failed generated domain gate'
    end if
    if (saddle_rate_floor_margin <= 0.0_dp) then
        error stop 'saddle-rate floor is not positive'
    end if
    if (epsilon_floor_margin <= 0.0_dp) then
        error stop 'epsilon floor is not positive'
    end if
    call require_close('section keeps nondegeneracy margin', &
        section_reversed_margin, nondegeneracy_margin)
    call require_close('section keeps epsilon margin', &
        section_reversed_epsilon_margin, epsilon_margin)

    call evaluate_neort_gc_homoclinic_coefficients(-SADDLE_RATE**2, 0.0_dp, &
        1.0_dp, TOROIDAL_RATE, hessian_determinant, saddle_discriminant, &
        saddle_rate_value, c_tau, c_phi, one_leg_tau_coefficient, &
        two_leg_tau_coefficient, one_leg_phi_coefficient, &
        two_leg_phi_coefficient, c_tau_section_reverse, &
        c_phi_section_reverse, c_tau_toroidal_reverse, &
        c_phi_toroidal_reverse, absolute_phi_reversal_defect)
    call require_close('quadratic Hessian determinant', hessian_determinant, &
        -SADDLE_RATE**2)
    call require_close('quadratic saddle discriminant', saddle_discriminant, &
        SADDLE_RATE**2)
    call require_close('quadratic saddle rate', saddle_rate_value, SADDLE_RATE)
    call require_close('one-leg C tau', c_tau, 1.0_dp/SADDLE_RATE)
    call require_close('one-leg C phi', c_phi, TOROIDAL_RATE/SADDLE_RATE)
    call require_close('one-leg time coefficient', one_leg_tau_coefficient, &
        c_tau)
    call require_close('two-leg time coefficient', two_leg_tau_coefficient, &
        2.0_dp*c_tau)
    call require_close('one-leg phi coefficient', one_leg_phi_coefficient, &
        c_phi)
    call require_close('two-leg phi coefficient', two_leg_phi_coefficient, &
        2.0_dp*c_phi)
    call require_close('section reversal keeps C tau', c_tau_section_reverse, &
        c_tau)
    call require_close('section reversal keeps C phi', c_phi_section_reverse, &
        c_phi)
    call require_close('toroidal reversal keeps C tau', &
        c_tau_toroidal_reverse, c_tau)
    call require_close('toroidal reversal flips C phi', &
        c_phi_toroidal_reverse, -c_phi)
    call require_close('absolute C phi reversal defect', &
        absolute_phi_reversal_defect, 2.0_dp*abs(c_phi))
    if (absolute_phi_reversal_defect <= 0.0_dp) then
        error stop 'absolute C phi replacement did not expose sign loss'
    end if

    do i = 1, NUMBER_OF_EPSILONS
        tau_one(i) = integrate_manufactured_saddle(0.0_dp, LOCAL_DOMAIN, &
            EPSILONS(i), TIME_INTEGRAND)
        tau_two(i) = integrate_manufactured_saddle(-LOCAL_DOMAIN, &
            LOCAL_DOMAIN, EPSILONS(i), TIME_INTEGRAND)
        phi_one(i) = integrate_manufactured_saddle(0.0_dp, LOCAL_DOMAIN, &
            EPSILONS(i), PHI_INTEGRAND)
        phi_two(i) = integrate_manufactured_saddle(-LOCAL_DOMAIN, &
            LOCAL_DOMAIN, EPSILONS(i), PHI_INTEGRAND)
        phi_one_reversed(i) = integrate_manufactured_saddle(0.0_dp, &
            LOCAL_DOMAIN, EPSILONS(i), REVERSED_PHI_INTEGRAND)
        call evaluate_neort_gc_homoclinic_finite_part(tau_one(i), tau_two(i), &
            phi_one(i), phi_two(i), EPSILONS(i), c_tau, c_phi, &
            tau_finite_one(i), tau_finite_two(i), phi_finite_one(i), &
            phi_finite_two(i), phi_finite_one_reverse(i), &
            phi_finite_two_reverse(i))
        call require_close('numerical toroidal-coordinate reversal', &
            phi_one_reversed(i), -phi_one(i))
        call require_close('generated one-leg finite-part reversal', &
            phi_finite_one_reverse(i), -phi_finite_one(i))
        call require_close('generated two-leg finite-part reversal', &
            phi_finite_two_reverse(i), -phi_finite_two(i))
    end do

    do i = 1, NUMBER_OF_EPSILONS-1
        logarithmic_interval = log(EPSILONS(i)/EPSILONS(i+1))
        c_tau_one_estimate(i) = &
            (tau_one(i+1)-tau_one(i))/logarithmic_interval
        c_tau_two_estimate(i) = &
            (tau_two(i+1)-tau_two(i))/logarithmic_interval
        c_phi_one_estimate(i) = &
            (phi_one(i+1)-phi_one(i))/logarithmic_interval
        c_phi_two_estimate(i) = &
            (phi_two(i+1)-phi_two(i))/logarithmic_interval
        c_phi_reverse_estimate(i) = &
            (phi_one_reversed(i+1)-phi_one_reversed(i))/ &
            logarithmic_interval
        call require_close('numerical two-to-one time coefficient', &
            c_tau_two_estimate(i), 2.0_dp*c_tau_one_estimate(i))
        call require_close('numerical two-to-one phi coefficient', &
            c_phi_two_estimate(i), 2.0_dp*c_phi_one_estimate(i))
    end do

    if (abs(c_tau_one_estimate(3)-c_tau) >= &
            abs(c_tau_one_estimate(1)-c_tau)) then
        error stop 'numerical C tau did not converge'
    end if
    if (abs(c_phi_one_estimate(3)-c_phi) >= &
            abs(c_phi_one_estimate(1)-c_phi)) then
        error stop 'numerical C phi did not converge'
    end if
    call require_absolute_close('converged one-leg C tau', &
        c_tau_one_estimate(3), c_tau, 2.0e-6_dp)
    call require_absolute_close('converged two-leg C tau', &
        c_tau_two_estimate(3), 2.0_dp*c_tau, 4.0e-6_dp)
    call require_absolute_close('converged signed C phi', &
        c_phi_one_estimate(3), c_phi, 2.0e-6_dp)
    call require_absolute_close('converged reversed C phi', &
        c_phi_reverse_estimate(3), c_phi_toroidal_reverse, 2.0e-6_dp)
    call require_absolute_close('one-leg tau finite-part convergence', &
        tau_finite_one(4), tau_finite_one(3), 2.0e-6_dp)
    call require_absolute_close('two-leg tau finite-part convergence', &
        tau_finite_two(4), tau_finite_two(3), 4.0e-6_dp)
    call require_absolute_close('one-leg phi finite-part convergence', &
        phi_finite_one(4), phi_finite_one(3), 2.0e-6_dp)
    call require_absolute_close('two-leg phi finite-part convergence', &
        phi_finite_two(4), phi_finite_two(3), 4.0e-6_dp)

    call check_degenerate_diagnostics

    write (*, '(A)') 'test_gc_homoclinic_limit OK'
    call pass_test

contains

    subroutine check_degenerate_diagnostics
        real(dp) :: near_h_qq

        call evaluate_neort_gc_homoclinic_diagnostic(0.0_dp, 0.0_dp, 1.0_dp, &
            SADDLE_RATE_FLOOR, EPSILONS(4), EPSILON_FLOOR, &
            diagnostic_hessian_determinant, diagnostic_saddle_discriminant, &
            nondegeneracy_margin, section_reversed_margin, epsilon_margin, &
            section_reversed_epsilon_margin, saddle_rate_floor_margin, &
            epsilon_floor_margin)
        call require_close('degenerate discriminant', &
            diagnostic_saddle_discriminant, 0.0_dp)
        if (nondegeneracy_margin >= 0.0_dp) then
            error stop 'exactly degenerate saddle passed generated margin'
        end if
        call require_close('degenerate section margin', &
            section_reversed_margin, nondegeneracy_margin)
        call require_finite('degenerate determinant', &
            diagnostic_hessian_determinant)
        call require_finite('degenerate discriminant', &
            diagnostic_saddle_discriminant)
        call require_finite('degenerate margin', nondegeneracy_margin)

        near_h_qq = -(0.5_dp*SADDLE_RATE_FLOOR)**2
        call evaluate_neort_gc_homoclinic_diagnostic(near_h_qq, 0.0_dp, &
            1.0_dp, SADDLE_RATE_FLOOR, EPSILONS(4), EPSILON_FLOOR, &
            diagnostic_hessian_determinant, diagnostic_saddle_discriminant, &
            nondegeneracy_margin, section_reversed_margin, epsilon_margin, &
            section_reversed_epsilon_margin, saddle_rate_floor_margin, &
            epsilon_floor_margin)
        call require_absolute_close('near-degenerate discriminant', &
            diagnostic_saddle_discriminant, &
            (0.5_dp*SADDLE_RATE_FLOOR)**2, 1.0e-20_dp)
        if (nondegeneracy_margin >= 0.0_dp) then
            error stop 'near-degenerate saddle passed generated margin'
        end if
        call require_close('near-degenerate section margin', &
            section_reversed_margin, nondegeneracy_margin)

        call evaluate_neort_gc_homoclinic_diagnostic(-SADDLE_RATE**2, &
            0.0_dp, 1.0_dp, SADDLE_RATE_FLOOR, EPSILON_FLOOR, &
            EPSILON_FLOOR, diagnostic_hessian_determinant, &
            diagnostic_saddle_discriminant, nondegeneracy_margin, &
            section_reversed_margin, epsilon_margin, &
            section_reversed_epsilon_margin, saddle_rate_floor_margin, &
            epsilon_floor_margin)
        call require_absolute_close('epsilon margin at floor', epsilon_margin, &
            0.0_dp, 0.0_dp)
        call require_close('section keeps epsilon floor boundary', &
            section_reversed_epsilon_margin, epsilon_margin)
    end subroutine check_degenerate_diagnostics

    real(dp) function integrate_manufactured_saddle(lower, upper, &
            epsilon_value, integrand_kind)
        real(dp), intent(in) :: lower, upper, epsilon_value
        integer, intent(in) :: integrand_kind
        real(dp) :: midpoint, f_lower, f_midpoint, f_upper, whole

        midpoint = 0.5_dp*(lower+upper)
        f_lower = manufactured_integrand(lower, epsilon_value, integrand_kind)
        f_midpoint = manufactured_integrand(midpoint, epsilon_value, &
            integrand_kind)
        f_upper = manufactured_integrand(upper, epsilon_value, integrand_kind)
        whole = (upper-lower)*(f_lower+4.0_dp*f_midpoint+f_upper)/6.0_dp
        integrate_manufactured_saddle = adaptive_simpson(lower, upper, &
            f_lower, f_midpoint, f_upper, whole, epsilon_value, &
            integrand_kind, QUADRATURE_TOLERANCE, 40)
    end function integrate_manufactured_saddle

    recursive real(dp) function adaptive_simpson(lower, upper, f_lower, &
            f_midpoint, f_upper, whole, epsilon_value, integrand_kind, &
            local_tolerance, depth) result(integral)
        real(dp), intent(in) :: lower, upper, f_lower, f_midpoint, f_upper
        real(dp), intent(in) :: whole, epsilon_value, local_tolerance
        integer, intent(in) :: integrand_kind, depth
        real(dp) :: midpoint, left_midpoint, right_midpoint
        real(dp) :: f_left_midpoint, f_right_midpoint
        real(dp) :: left, right, correction

        midpoint = 0.5_dp*(lower+upper)
        left_midpoint = 0.5_dp*(lower+midpoint)
        right_midpoint = 0.5_dp*(midpoint+upper)
        f_left_midpoint = manufactured_integrand(left_midpoint, epsilon_value, &
            integrand_kind)
        f_right_midpoint = manufactured_integrand(right_midpoint, &
            epsilon_value, integrand_kind)
        left = (midpoint-lower)*(f_lower+4.0_dp*f_left_midpoint+f_midpoint)/ &
            6.0_dp
        right = (upper-midpoint)*(f_midpoint+4.0_dp*f_right_midpoint+f_upper)/ &
            6.0_dp
        correction = left+right-whole
        if (abs(correction) <= 15.0_dp*local_tolerance) then
            integral = left+right+correction/15.0_dp
            return
        end if
        if (depth <= 0) then
            error stop 'manufactured saddle quadrature did not converge'
        end if
        integral = adaptive_simpson(lower, midpoint, f_lower, &
            f_left_midpoint, f_midpoint, left, epsilon_value, integrand_kind, &
            0.5_dp*local_tolerance, depth-1) + &
            adaptive_simpson(midpoint, upper, f_midpoint, f_right_midpoint, &
            f_upper, right, epsilon_value, integrand_kind, &
            0.5_dp*local_tolerance, depth-1)
    end function adaptive_simpson

    real(dp) function manufactured_integrand(coordinate, epsilon_value, &
            integrand_kind)
        real(dp), intent(in) :: coordinate, epsilon_value
        integer, intent(in) :: integrand_kind
        real(dp) :: inverse_speed

        ! Independent energy-shell oracle for
        ! H=(p**2-SADDLE_RATE**2*q**2)/2=epsilon**2/2.
        inverse_speed = 1.0_dp/sqrt((SADDLE_RATE*coordinate)**2+ &
            epsilon_value**2)
        select case (integrand_kind)
        case (TIME_INTEGRAND)
            manufactured_integrand = inverse_speed
        case (PHI_INTEGRAND)
            manufactured_integrand = TOROIDAL_RATE*inverse_speed
        case (REVERSED_PHI_INTEGRAND)
            manufactured_integrand = -TOROIDAL_RATE*inverse_speed
        case default
            error stop 'unknown manufactured saddle integrand'
        end select
    end function manufactured_integrand

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual-expected) > 5.0e-9_dp*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'homoclinic limit oracle failed'
        end if
    end subroutine require_close

    subroutine require_absolute_close(label, actual, expected, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected, tolerance

        if (abs(actual-expected) > tolerance) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'homoclinic absolute oracle failed'
        end if
    end subroutine require_absolute_close

    subroutine require_finite(label, value)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value

        if (.not. ieee_is_finite(value)) then
            write (*, '(A,1X,ES24.16)') trim(label), value
            error stop 'homoclinic diagnostic is not finite'
        end if
    end subroutine require_finite

end program test_gc_homoclinic_limit
