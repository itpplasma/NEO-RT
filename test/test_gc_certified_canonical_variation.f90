program test_gc_certified_canonical_variation
    use, intrinsic :: ieee_arithmetic, only: ieee_next_after
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_certified_interval_roots, only: &
        gc_interval_callback_result_t, &
        gc_interval_root_box_t, gc_interval_root_options_t, &
        gc_interval_root_result_t, gc_interval_t
    use neort_gc_certified_canonical_variation, only: &
        GC_CANONICAL_VARIATION_INVALID_ROOT_EVIDENCE, &
        GC_CANONICAL_VARIATION_ROOT_FAILURE, GC_CANONICAL_VARIATION_SUCCESS, &
        gc_canonical_variation_options_t, gc_canonical_variation_result_t, &
        assemble_gc_canonical_total_variation, certify_gc_canonical_total_variation
    implicit none

    integer, parameter :: cut_id = 77
    integer, parameter :: value_certificate_id = 700
    integer, parameter :: root_certificate_id = 701
    real(dp), parameter :: domain_lo = 0.0_dp
    real(dp), parameter :: domain_hi = 4.0_dp

    type(gc_canonical_variation_options_t) :: options, unresolved_options
    type(gc_canonical_variation_result_t) :: result, missing_result, malformed_result
    type(gc_interval_root_result_t) :: missing_roots, malformed_roots

    options%root_options%initial_partition = 1
    options%root_options%max_depth = 80
    options%root_options%max_boxes = 20000
    options%root_options%max_roots = 16
    options%root_options%max_stationary_iterations = 160
    options%root_options%x_tolerance = 1.0e-8_dp
    options%root_options%expected_enclosure_certificate_id = root_certificate_id
    options%root_options%expected_stationary_certificate_id = 0
    options%expected_cut_id = cut_id
    options%expected_value_certificate_id = value_certificate_id

    call certify_gc_canonical_total_variation(evaluate_psi, verify_psi, &
        evaluate_derivative, verify_derivative, verify_stationary, domain_lo, domain_hi, &
        options, result)
    call require(result%status == GC_CANONICAL_VARIATION_SUCCESS .and. result%certified, &
        'polynomial variation was not certified')
    call require(result%nroots == 2 .and. result%nspans == 3, &
        'the two stationary roots did not partition the cut into three spans')
    call require(result%total_variation_enclosure%lo <= 12.0_dp .and. &
        result%total_variation_enclosure%hi >= 12.0_dp, &
        'certified variation enclosure does not contain the oracle value 12')
    call require(result%endpoint_difference_enclosure%lo <= 4.0_dp .and. &
        result%endpoint_difference_enclosure%hi >= 4.0_dp, &
        'endpoint-difference diagnostic does not contain the oracle value 4')
    call require(result%total_variation_enclosure%lo > &
        result%endpoint_difference_enclosure%hi, &
        'variation collapsed to the endpoint difference')
    call require(abs(result%total_variation_estimate - 12.0_dp) < 1.0e-5_dp, &
        'variation estimate is not the monotone-decomposition value')
    call require(abs(result%endpoint_difference_estimate - 4.0_dp) < 1.0e-5_dp, &
        'endpoint-difference estimate is not 4')

    !! A missing derivative root leaves a zero of d psi*/dR inside a span;
    !! the monotonicity certificate must reject that evidence.
    missing_roots = result%derivative_roots
    deallocate(missing_roots%roots)
    allocate(missing_roots%roots(1))
    missing_roots%roots(1) = result%derivative_roots%roots(1)
    missing_roots%nroots = 1
    call assemble_gc_canonical_total_variation(evaluate_psi, verify_psi, domain_lo, &
        domain_hi, options, missing_roots, missing_result)
    call require(missing_result%status /= GC_CANONICAL_VARIATION_SUCCESS .and. &
        .not. missing_result%certified, 'missing derivative root was accepted')

    !! A malformed certificate must fail before any numerical assembly is used.
    malformed_roots = result%derivative_roots
    malformed_roots%roots(1)%classification_certified = .false.
    call assemble_gc_canonical_total_variation(evaluate_psi, verify_psi, domain_lo, &
        domain_hi, options, malformed_roots, malformed_result)
    call require(malformed_result%status == GC_CANONICAL_VARIATION_INVALID_ROOT_EVIDENCE .and. &
        .not. malformed_result%certified, 'malformed root evidence was accepted')

    !! Force the interval-root isolator to stop before it can close the roots.
    !! The top-level variation certificate must preserve that unresolved status.
    unresolved_options = options
    unresolved_options%root_options%initial_partition = 1
    unresolved_options%root_options%max_depth = 0
    call certify_gc_canonical_total_variation(evaluate_psi, verify_psi, &
        evaluate_derivative, verify_derivative, verify_stationary, domain_lo, domain_hi, &
        unresolved_options, malformed_result)
    call require(malformed_result%status == GC_CANONICAL_VARIATION_ROOT_FAILURE .and. &
        .not. malformed_result%certified, 'unresolved roots were accepted')

    write (*, '(A)') 'test_gc_certified_canonical_variation OK'

contains

    subroutine evaluate_psi(lo, hi, value)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(out) :: value
        type(gc_interval_t) :: r, r2, r3

        call initialize_value(lo, hi, value, value_certificate_id)
        if (lo > hi) return
        r = gc_interval_t(lo, hi)
        r2 = multiply_interval(r, r)
        r3 = multiply_interval(r2, r)
        value%f = add_interval(add_interval(r3, scale_interval(r2, -6.0_dp)), &
            scale_interval(r, 9.0_dp))
        value%df = derivative_interval(r)
        value%d2f = linear_interval(r, 6.0_dp, -12.0_dp)
    end subroutine evaluate_psi

    subroutine evaluate_derivative(lo, hi, value)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(out) :: value
        type(gc_interval_t) :: r

        call initialize_value(lo, hi, value, root_certificate_id)
        if (lo > hi) return
        r = gc_interval_t(lo, hi)
        value%f = derivative_interval(r)
        value%df = linear_interval(r, 6.0_dp, -12.0_dp)
        value%d2f = gc_interval_t(6.0_dp, 6.0_dp)
    end subroutine evaluate_derivative

    function derivative_interval(r) result(value)
        type(gc_interval_t), intent(in) :: r
        type(gc_interval_t) :: value
        type(gc_interval_t) :: left_factor, right_factor

        left_factor = subtract_interval(r, gc_interval_t(1.0_dp, 1.0_dp))
        right_factor = subtract_interval(r, gc_interval_t(3.0_dp, 3.0_dp))
        value = scale_interval(multiply_interval(left_factor, right_factor), 3.0_dp)
    end function derivative_interval

    subroutine initialize_value(lo, hi, value, certificate_id)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(out) :: value
        integer, intent(in) :: certificate_id

        value = gc_interval_callback_result_t()
        value%query_lo = lo
        value%query_hi = hi
        value%status = 0
        value%cut_id = cut_id
        value%enclosure_certificate_id = certificate_id
        value%stationary_certificate_id = 0
    end subroutine initialize_value

    subroutine verify_psi(lo, hi, value, expected_id, status)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(in) :: value
        integer, intent(in) :: expected_id
        integer, intent(out) :: status
        type(gc_interval_t) :: expected_f, expected_df, expected_d2f

        status = 1
        if (expected_id /= value_certificate_id) return
        if (value%cut_id /= cut_id) return
        if (value%query_lo /= lo .or. value%query_hi /= hi) return
        expected_f = exact_psi_range(lo, hi)
        expected_df = exact_derivative_range(lo, hi)
        expected_d2f = exact_linear_range(lo, hi, 6.0_dp, -12.0_dp)
        if (.not. encloses(value%f, expected_f)) return
        if (.not. encloses(value%df, expected_df)) return
        if (.not. encloses(value%d2f, expected_d2f)) return
        status = 0
    end subroutine verify_psi

    subroutine verify_derivative(lo, hi, value, expected_id, status)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(in) :: value
        integer, intent(in) :: expected_id
        integer, intent(out) :: status
        type(gc_interval_t) :: expected_f, expected_df, expected_d2f

        status = 1
        if (expected_id /= root_certificate_id) return
        if (value%cut_id /= cut_id) return
        if (value%query_lo /= lo .or. value%query_hi /= hi) return
        expected_f = exact_derivative_range(lo, hi)
        expected_df = exact_linear_range(lo, hi, 6.0_dp, -12.0_dp)
        expected_d2f = gc_interval_t(6.0_dp, 6.0_dp)
        if (.not. encloses(value%f, expected_f)) return
        if (.not. encloses(value%df, expected_df)) return
        if (.not. encloses(value%d2f, expected_d2f)) return
        status = 0
    end subroutine verify_derivative

    subroutine verify_stationary(lo, hi, point, value, expected_id, expected_stationary_id, &
            status)
        real(dp), intent(in) :: lo, hi, point
        type(gc_interval_callback_result_t), intent(out) :: value
        integer, intent(in) :: expected_id, expected_stationary_id
        integer, intent(out) :: status

        value = gc_interval_callback_result_t()
        value%query_lo = point
        value%query_hi = point
        value%cut_id = cut_id
        value%enclosure_certificate_id = expected_id
        value%stationary_certificate_id = expected_stationary_id
        value%stationary_point = point
        status = 1
        if (lo > hi) return
        if (point < lo .or. point > hi) return
        !! This fixture has no tangent derivative roots, so this callback is
        !! deliberately fail-closed if the isolator ever asks for one.
    end subroutine verify_stationary

    function exact_psi_range(lo, hi) result(value)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_t) :: value

        value = range_from_values([direct_psi(lo), direct_psi(hi)])
        if (lo <= 1.0_dp .and. hi >= 1.0_dp) call include_value(value, 4.0_dp)
        if (lo <= 3.0_dp .and. hi >= 3.0_dp) call include_value(value, 0.0_dp)
    end function exact_psi_range

    function exact_derivative_range(lo, hi) result(value)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_t) :: value

        value = range_from_values([direct_derivative(lo), direct_derivative(hi)])
        if (lo <= 2.0_dp .and. hi >= 2.0_dp) call include_value(value, -3.0_dp)
    end function exact_derivative_range

    function exact_linear_range(lo, hi, slope, intercept) result(value)
        real(dp), intent(in) :: lo, hi, slope, intercept
        type(gc_interval_t) :: value

        value = range_from_values([slope*lo + intercept, slope*hi + intercept])
    end function exact_linear_range

    function direct_psi(x) result(value)
        real(dp), intent(in) :: x
        real(dp) :: value

        value = x*x*x - 6.0_dp*x*x + 9.0_dp*x
    end function direct_psi

    function direct_derivative(x) result(value)
        real(dp), intent(in) :: x
        real(dp) :: value

        value = 3.0_dp*(x - 1.0_dp)*(x - 3.0_dp)
    end function direct_derivative

    function range_from_values(values) result(value)
        real(dp), intent(in) :: values(:)
        type(gc_interval_t) :: value

        value%lo = down(minval(values))
        value%hi = up(maxval(values))
        if (minval(values) == 0.0_dp .and. maxval(values) == 0.0_dp) then
            value = gc_interval_t(0.0_dp, 0.0_dp)
        end if
    end function range_from_values

    subroutine include_value(interval, number)
        type(gc_interval_t), intent(inout) :: interval
        real(dp), intent(in) :: number

        interval%lo = min(interval%lo, down(number))
        interval%hi = max(interval%hi, up(number))
    end subroutine include_value

    function linear_interval(interval, slope, intercept) result(value)
        type(gc_interval_t), intent(in) :: interval
        real(dp), intent(in) :: slope, intercept
        type(gc_interval_t) :: value

        value = exact_linear_range(interval%lo, interval%hi, slope, intercept)
    end function linear_interval

    function scale_interval(interval, factor) result(value)
        type(gc_interval_t), intent(in) :: interval
        real(dp), intent(in) :: factor
        type(gc_interval_t) :: value

        value = multiply_interval(interval, gc_interval_t(factor, factor))
    end function scale_interval

    function add_interval(left, right) result(value)
        type(gc_interval_t), intent(in) :: left, right
        type(gc_interval_t) :: value
        real(dp) :: lo, hi

        lo = left%lo + right%lo
        hi = left%hi + right%hi
        if (lo == 0.0_dp .and. hi == 0.0_dp) then
            value = gc_interval_t(0.0_dp, 0.0_dp)
        else
            value = gc_interval_t(down(lo), up(hi))
        end if
    end function add_interval

    function subtract_interval(left, right) result(value)
        type(gc_interval_t), intent(in) :: left, right
        type(gc_interval_t) :: value
        real(dp) :: lo, hi

        lo = left%lo - right%hi
        hi = left%hi - right%lo
        if (lo == 0.0_dp .and. hi == 0.0_dp) then
            value = gc_interval_t(0.0_dp, 0.0_dp)
        else
            value = gc_interval_t(down(lo), up(hi))
        end if
    end function subtract_interval

    function multiply_interval(left, right) result(value)
        type(gc_interval_t), intent(in) :: left, right
        type(gc_interval_t) :: value
        real(dp) :: products(4), lo, hi

        products = [left%lo*right%lo, left%lo*right%hi, &
            left%hi*right%lo, left%hi*right%hi]
        lo = minval(products)
        hi = maxval(products)
        if (lo == 0.0_dp .and. hi == 0.0_dp) then
            value = gc_interval_t(0.0_dp, 0.0_dp)
        else
            value = gc_interval_t(down(lo), up(hi))
        end if
    end function multiply_interval

    logical function encloses(actual, expected)
        type(gc_interval_t), intent(in) :: actual, expected

        encloses = actual%lo <= expected%lo .and. actual%hi >= expected%hi
    end function encloses

    real(dp) function down(number)
        real(dp), intent(in) :: number

        if (number == 0.0_dp) then
            down = 0.0_dp
        else
            down = ieee_next_after(number, -huge(number))
        end if
    end function down

    real(dp) function up(number)
        real(dp), intent(in) :: number

        if (number == 0.0_dp) then
            up = 0.0_dp
        else
            up = ieee_next_after(number, huge(number))
        end if
    end function up

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*, '(A)') trim(message)
            error stop 1
        end if
    end subroutine require

end program test_gc_certified_canonical_variation
