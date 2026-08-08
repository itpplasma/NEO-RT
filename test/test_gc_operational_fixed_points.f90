program test_gc_operational_fixed_points
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_operational_fixed_points, only: &
        GC_FIXED_POINT_DEGENERATE, GC_FIXED_POINT_O, &
        GC_FIXED_POINT_SUCCESS, GC_FIXED_POINT_X, &
        find_gc_operational_fixed_points, gc_operational_fixed_point_result_t
    use neort_gc_operational_scalar_roots, only: &
        gc_operational_root_options_t
    use util_for_test, only: pass_test
    implicit none

    type(gc_operational_root_options_t) :: options
    type(gc_operational_fixed_point_result_t) :: result

    options%initial_intervals = 5
    call find_gc_operational_fixed_points(cubic_canonical, &
        manufactured_classifier, -2.0_dp, 2.0_dp, options, result)
    call require(result%status == GC_FIXED_POINT_SUCCESS, &
        'interior fixed-point search failed')
    call require(result%complete, 'fixed-point set is incomplete')
    call require(result%npoints == 2, 'interior stationary points were lost')
    call require(result%n_o_points == 1 .and. result%n_x_points == 1, &
        'O/X counts are wrong')
    call require(result%points(1)%kind == GC_FIXED_POINT_O .and. &
        result%points(2)%kind == GC_FIXED_POINT_X, &
        'flow discriminant signs were not retained')
    call require(maxval(abs([result%points%x]-[-1.0_dp, 1.0_dp])) < &
        1.0e-10_dp, 'fixed-point locations are wrong')
    call require(maxval(abs([result%points%canonical_momentum]- &
        [2.0_dp, -2.0_dp])) < 1.0e-10_dp, &
        'fixed-point canonical momenta are wrong')

    call find_gc_operational_fixed_points(endpoint_stationary, &
        always_o_classifier, 0.0_dp, 1.0_dp, options, result)
    call require(result%status == GC_FIXED_POINT_SUCCESS, &
        'endpoint-stationary search failed')
    call require(result%npoints == 0, &
        'allowed-region endpoint was misclassified as a fixed point')

    call find_gc_operational_fixed_points(cubic_canonical, &
        degenerate_classifier, -2.0_dp, 2.0_dp, options, result)
    call require(result%status == GC_FIXED_POINT_DEGENERATE, &
        'zero flow discriminant was assigned an O/X label')
    call require(.not. result%complete, &
        'degenerate fixed-point set was marked complete')

    call pass_test

contains

    subroutine cubic_canonical(x, p_star, first, second, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: p_star, first, second
        integer, intent(out) :: status

        p_star = x**3-3.0_dp*x
        first = 3.0_dp*x**2-3.0_dp
        second = 6.0_dp*x
        status = 0
    end subroutine cubic_canonical

    subroutine endpoint_stationary(x, p_star, first, second, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: p_star, first, second
        integer, intent(out) :: status

        p_star = x**2
        first = 2.0_dp*x
        second = 2.0_dp
        status = 0
    end subroutine endpoint_stationary

    subroutine manufactured_classifier(x, discriminant, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: discriminant
        integer, intent(out) :: status

        discriminant = 4.0_dp*x
        status = 0
    end subroutine manufactured_classifier

    subroutine always_o_classifier(x, discriminant, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: discriminant
        integer, intent(out) :: status

        discriminant = -4.0_dp
        status = 0
        associate (unused_x => x)
        end associate
    end subroutine always_o_classifier

    subroutine degenerate_classifier(x, discriminant, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: discriminant
        integer, intent(out) :: status

        discriminant = 0.0_dp
        status = 0
        associate (unused_x => x)
        end associate
    end subroutine degenerate_classifier

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_operational_fixed_points
