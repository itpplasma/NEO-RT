program test_gc_full_fow_simple_root_force
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_full_fow_simple_root_symbolic, only: &
        evaluate_neort_full_fow_simple_root_force
    use neort_gc_nonlocal_resonance_integral, only: &
        evaluate_gc_nonlocal_root_contribution
    use neort_gc_nonlocal_resonance_types, only: GC_NONLOCAL_SAMPLE_VALID, &
        GC_NONLOCAL_SINGULAR_RESONANCE, GC_NONLOCAL_SUCCESS, &
        gc_nonlocal_orbit_sample_t
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: dpsi_star_dx = -1.75_dp
    real(dp), parameter :: hm_real = 1.25_dp, hm_imag = -0.5_dp
    real(dp), parameter :: tau_b = 0.9_dp
    real(dp), parameter :: positive_force = 2.75_dp
    real(dp), parameter :: positive_derivative = 2.4_dp
    ! This is only a forbidden alternative for the behavioral oracle: the
    ! kernel API has no mode input, so no mode-squared factor is available.
    real(dp), parameter :: mode_for_no_n_squared = 3.0_dp
    real(dp) :: positive_weight, negative_derivative_weight
    real(dp) :: positive_force_result, negative_force_result
    real(dp) :: positive_force_oracle, negative_force_oracle
    real(dp) :: negative_derivative_force_result
    real(dp) :: runtime_contribution(2), runtime_reversed(2)
    type(gc_nonlocal_orbit_sample_t) :: runtime_sample
    integer :: runtime_status

    call check_case('positive force', positive_derivative, positive_force, &
        positive_force_result, positive_force_oracle, positive_weight)
    call check_case('negative force', positive_derivative, -positive_force, &
        negative_force_result, negative_force_oracle)
    call check_case('negative residual derivative', -positive_derivative, &
        positive_force, negative_derivative_force_result, &
        positive_force_oracle, negative_derivative_weight)

    call require_close('force sign reversal', negative_force_result, &
        -positive_force_result)
    call require_close('negative force independent oracle', negative_force_result, &
        negative_force_oracle)
    call require_close('residual derivative sign invariance', &
        negative_derivative_weight, positive_weight)
    call require_close('force derivative sign invariance', &
        negative_derivative_force_result, positive_force_result)
    call require_close('force sign remains signed', positive_force_result, &
        abs(positive_force_result))
    call require_not_close('no hidden n squared factor', positive_force_result, &
        mode_for_no_n_squared**2*positive_force_oracle)

    runtime_sample%status = GC_NONLOCAL_SAMPLE_VALID
    runtime_sample%dpsi_star_dx = dpsi_star_dx
    runtime_sample%tau_b = tau_b
    runtime_sample%h_m = cmplx(hm_real, hm_imag, kind=dp)
    runtime_sample%nforce = 2
    runtime_sample%thermodynamic_force(1:2) = &
        [positive_force, -positive_force]
    call evaluate_gc_nonlocal_root_contribution(runtime_sample, &
        positive_derivative, 2, 1.0e-13_dp, runtime_contribution, &
        runtime_status)
    call require_status('runtime generated-root status', runtime_status, &
        GC_NONLOCAL_SUCCESS)
    call require_close('runtime positive signed force', &
        runtime_contribution(1), positive_force_oracle)
    call require_close('runtime negative signed force', &
        runtime_contribution(2), negative_force_oracle)
    call evaluate_gc_nonlocal_root_contribution(runtime_sample, &
        -positive_derivative, 2, 1.0e-13_dp, runtime_reversed, &
        runtime_status)
    call require_status('runtime derivative-reversal status', runtime_status, &
        GC_NONLOCAL_SUCCESS)
    call require_close('runtime derivative-sign invariance', &
        runtime_reversed(1), runtime_contribution(1))
    call evaluate_gc_nonlocal_root_contribution(runtime_sample, &
        1.0e-14_dp, 2, 1.0e-13_dp, runtime_reversed, runtime_status)
    call require_status('runtime tangent root is rejected', runtime_status, &
        GC_NONLOCAL_SINGULAR_RESONANCE)

    write (*, '(A)') 'test_gc_full_fow_simple_root_force OK'
    call pass_test

contains

    subroutine check_case(label, frequency_residual_derivative, force_value, &
            force_result, expected_force_contribution, weight_result)
        character(*), intent(in) :: label
        real(dp), intent(in) :: frequency_residual_derivative, force_value
        real(dp), intent(out) :: force_result, expected_force_contribution
        real(dp), intent(out), optional :: weight_result
        real(dp) :: abs_Hm_squared, simple_root_weight, force_contribution
        real(dp) :: expected_abs_Hm_squared, expected_simple_root_weight

        call evaluate_neort_full_fow_simple_root_force(dpsi_star_dx, hm_real, &
            hm_imag, tau_b, frequency_residual_derivative, force_value, &
            abs_Hm_squared, simple_root_weight, force_contribution)

        expected_abs_Hm_squared = hm_real**2 + hm_imag**2
        ! Evaluate the complete positive delta-root weight independently;
        ! force_value supplies its sign and is not folded into the measure.
        expected_simple_root_weight = abs(dpsi_star_dx)*expected_abs_Hm_squared*&
            tau_b/abs(frequency_residual_derivative)
        expected_force_contribution = force_value*expected_simple_root_weight

        call require_close(trim(label)//' |Hm| squared', abs_Hm_squared, &
            expected_abs_Hm_squared)
        call require_close(trim(label)//' simple-root weight', &
            simple_root_weight, expected_simple_root_weight)
        call require_close(trim(label)//' signed force contribution', &
            force_contribution, expected_force_contribution)
        force_result = force_contribution
        if (present(weight_result)) weight_result = simple_root_weight
    end subroutine check_case

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual - expected) > 2.0e-12_dp*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'simple-root force oracle failed'
        end if
    end subroutine require_close

    subroutine require_not_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(actual), abs(expected)))
        if (abs(actual - expected) <= 2.0e-12_dp*scale) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'simple-root force distinction failed'
        end if
    end subroutine require_not_close

    subroutine require_status(label, actual, expected)
        character(*), intent(in) :: label
        integer, intent(in) :: actual, expected

        if (actual /= expected) then
            write (*, '(A,2(1X,I0))') trim(label), actual, expected
            error stop 'simple-root status oracle failed'
        end if
    end subroutine require_status

end program test_gc_full_fow_simple_root_force
