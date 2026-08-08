program test_gc_profile_potential_map
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_t
    use neort_profile_potential_map_interval_symbolic, only: &
        evaluate_neort_profile_potential_map_interval
    use neort_profile_potential_map_symbolic, only: &
        evaluate_neort_profile_potential_map
    implicit none

    real(dp) :: phi, derivative, second_derivative
    real(dp) :: reversed_phi, reversed_derivative, reversed_second
    type(gc_outward_interval_t) :: phi_box, derivative_box, second_box

    call evaluate_neort_profile_potential_map(4.0_dp, 2.0_dp, 5.0_dp, &
        7.0_dp, 3.0_dp, 9.0_dp, 2.0_dp, phi, derivative, &
        second_derivative)
    call require_close(phi, 12.0_dp, 'quadratic potential value')
    call require_close(derivative, 3.5_dp, 'local electric-field slope')
    call require_close(second_derivative, 1.0_dp, &
        'potential segment curvature')

    ! The independently integrated right endpoint is Phi(5)=16.  Reversing
    ! the segment data must describe the same physical quadratic at psi=4.
    call evaluate_neort_profile_potential_map(4.0_dp, 5.0_dp, 2.0_dp, &
        16.0_dp, 9.0_dp, 3.0_dp, 2.0_dp, reversed_phi, &
        reversed_derivative, reversed_second)
    call require_close(reversed_phi, phi, 'segment reversal potential')
    call require_close(reversed_derivative, derivative, &
        'segment reversal derivative')
    call require_close(reversed_second, second_derivative, &
        'segment reversal curvature')

    call evaluate_neort_profile_potential_map_interval( &
        gc_outward_interval(3.0_dp, 4.0_dp), &
        gc_outward_interval(2.0_dp, 2.0_dp), &
        gc_outward_interval(5.0_dp, 5.0_dp), &
        gc_outward_interval(7.0_dp, 7.0_dp), &
        gc_outward_interval(3.0_dp, 3.0_dp), &
        gc_outward_interval(9.0_dp, 9.0_dp), &
        gc_outward_interval(2.0_dp, 2.0_dp), phi_box, derivative_box, &
        second_box)
    call require_contains(phi_box, 9.0_dp, 12.0_dp, &
        'interval potential image')
    call require_contains(derivative_box, 2.5_dp, 3.5_dp, &
        'interval potential derivative')
    call require_contains(second_box, 1.0_dp, 1.0_dp, &
        'interval potential curvature')

    write (*, '(A)') 'test_gc_profile_potential_map: PASS'

contains

    subroutine require_close(actual, expected, message)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: message

        if (abs(actual-expected) > 2.0e-13_dp*max(1.0_dp, abs(expected))) then
            error stop message
        end if
    end subroutine require_close

    subroutine require_contains(interval, lower, upper, message)
        type(gc_outward_interval_t), intent(in) :: interval
        real(dp), intent(in) :: lower, upper
        character(*), intent(in) :: message

        if (interval%lo > lower .or. interval%hi < upper) error stop message
    end subroutine require_contains

end program test_gc_profile_potential_map
