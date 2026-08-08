program test_gc_operational_scalar_roots
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_operational_scalar_roots, only: &
        GC_OPERATIONAL_ROOT_CALLBACK_FAILURE, &
        GC_OPERATIONAL_ROOT_DEGENERATE, GC_OPERATIONAL_ROOT_SUCCESS, &
        find_gc_operational_scalar_roots, gc_operational_root_options_t, &
        gc_operational_root_result_t
    use util_for_test, only: pass_test
    implicit none

    type(gc_operational_root_options_t) :: options
    type(gc_operational_root_result_t) :: result

    options%initial_intervals = 7
    call find_gc_operational_scalar_roots(four_root_polynomial, 0.0_dp, &
        5.0_dp, options, result)
    call require(result%status == GC_OPERATIONAL_ROOT_SUCCESS, &
        'four-root search failed')
    call require(result%complete .and. result%coarse_fine_agree, &
        'four-root coarse/fine gate did not close')
    call require(result%nroots == 4, 'four-root search lost a component')
    call require(maxval(abs([result%roots%x]-[1.0_dp, 2.0_dp, 3.0_dp, &
        4.0_dp])) < 1.0e-10_dp, 'four-root locations are wrong')
    call require(all(result%roots%simple), 'simple roots were not retained')

    options%initial_intervals = 1
    call find_gc_operational_scalar_roots(narrow_pair, 0.0_dp, 1.0_dp, &
        options, result)
    call require(result%status == GC_OPERATIONAL_ROOT_SUCCESS, &
        'same-cell root-pair search failed')
    call require(result%nroots == 2, &
        'derivative extremum did not split a same-cell root pair')
    call require(maxval(abs([result%roots%x]-[0.49_dp, 0.51_dp])) < &
        1.0e-10_dp, 'same-cell root-pair locations are wrong')

    options%initial_intervals = 5
    call find_gc_operational_scalar_roots(endpoint_polynomial, 0.0_dp, &
        1.0_dp, options, result)
    call require(result%status == GC_OPERATIONAL_ROOT_SUCCESS, &
        'endpoint-root search failed')
    call require(result%nroots == 3, 'endpoint roots were omitted')
    call require(result%roots(1)%endpoint_root .and. &
        result%roots(3)%endpoint_root, 'endpoint roles were not retained')

    call find_gc_operational_scalar_roots(tangent_polynomial, 0.0_dp, &
        1.0_dp, options, result)
    call require(result%status == GC_OPERATIONAL_ROOT_DEGENERATE, &
        'tangent root was promoted to a simple operational root')
    call require(.not. result%complete, &
        'tangent-root topology was marked complete')

    options%initial_intervals = 10
    call find_gc_operational_scalar_roots(interior_hole, 0.0_dp, 1.0_dp, &
        options, result)
    call require(result%status == GC_OPERATIONAL_ROOT_CALLBACK_FAILURE, &
        'interior invalid-domain hole was accepted')
    call require(result%nroots == 0 .and. .not. result%complete, &
        'partial roots survived a callback failure')

    call pass_test

contains

    subroutine four_root_polynomial(x, value, derivative, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value, derivative
        integer, intent(out) :: status

        value = (x-1.0_dp)*(x-2.0_dp)*(x-3.0_dp)*(x-4.0_dp)
        derivative = 4.0_dp*x**3-30.0_dp*x**2+70.0_dp*x-50.0_dp
        status = 0
    end subroutine four_root_polynomial

    subroutine narrow_pair(x, value, derivative, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value, derivative
        integer, intent(out) :: status

        value = (x-0.49_dp)*(x-0.51_dp)
        derivative = 2.0_dp*x-1.0_dp
        status = 0
    end subroutine narrow_pair

    subroutine endpoint_polynomial(x, value, derivative, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value, derivative
        integer, intent(out) :: status

        value = x*(x-0.5_dp)*(x-1.0_dp)
        derivative = 3.0_dp*x**2-3.0_dp*x+0.5_dp
        status = 0
    end subroutine endpoint_polynomial

    subroutine tangent_polynomial(x, value, derivative, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value, derivative
        integer, intent(out) :: status

        value = (x-0.5_dp)**2
        derivative = 2.0_dp*(x-0.5_dp)
        status = 0
    end subroutine tangent_polynomial

    subroutine interior_hole(x, value, derivative, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value, derivative
        integer, intent(out) :: status

        value = x-0.5_dp
        derivative = 1.0_dp
        status = 0
        if (x >= 0.4_dp .and. x <= 0.6_dp) status = 1
    end subroutine interior_hole

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_operational_scalar_roots
