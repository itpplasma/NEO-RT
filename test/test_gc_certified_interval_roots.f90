program test_gc_certified_interval_roots
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after
    use neort_gc_certified_interval_roots
    implicit none

    type(gc_interval_root_options_t) :: options
    type(gc_interval_root_result_t) :: result, coarse, fine, repeat
    integer :: callback_mode, callback_calls

    options%max_depth = 80
    options%max_boxes = 20000
    options%max_roots = 16
    options%max_stationary_iterations = 160
    options%expected_enclosure_certificate_id = 901
    options%x_tolerance = 1.0e-8_dp

    call run_case(1, -3.0_dp, 3.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_SUCCESS, &
        'x^2+1 was not certified root-free')
    call require(result%nroots == 0 .and. result%coverage_certified, &
        'root-free coverage certificate missing')

    call run_case(12, -3.0_dp, 3.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_SUCCESS .and. result%coverage_certified, &
        'acceptance incorrectly depended on callback certified boolean')

    call run_case(2, -3.0_dp, 3.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_SUCCESS .and. result%nroots == 2, &
        'simple roots were unresolved or missed')
    call require(all_simple_evidence(result), 'simple-root evidence missing')
    call require(all_widths_within(result, options%x_tolerance), &
        'simple-root enclosure exceeds x_tolerance')
    call require(result%roots(1)%hi < -1.9_dp .and. result%roots(1)%lo > -2.1_dp, &
        'root near -2 is not enclosed')
    call require(result%roots(2)%lo < 1.1_dp .and. result%roots(2)%hi > 0.9_dp, &
        'root near 1 is not enclosed')

    options%expected_stationary_certificate_id = 303
    call run_case(3, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_SUCCESS .and. result%nroots == 1, &
        'tangent root was unresolved')
    call require(result%roots(1)%kind == GC_INTERVAL_ROOT_TANGENT, &
        'double root was not classified tangent')
    call require(result%roots(1)%multiplicity_lower == 2 .and. &
        result%roots(1)%multiplicity_upper == 2 .and. &
        result%roots(1)%transversality_kind == GC_INTERVAL_ROOT_EXTREMAL .and. &
        result%roots(1)%stationary_certified .and. &
        (result%roots(1)%bracket_certified .or. &
            result%roots(1)%interval_newton_certified) .and. &
        .not. result%roots(1)%transversality_certified .and. &
        result%roots(1)%classification_certified .and. &
        excludes_zero_for_test(result%roots(1)%second_derivative_enclosure), &
        'tangent evidence missing')
    call require(result%roots(1)%lo == 0.3_dp .and. result%roots(1)%hi == 0.3_dp, &
        'exact stationary certificate did not return its point enclosure')

    call run_case(14, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified, 'forged stationary certificate was accepted')

    options%expected_stationary_certificate_id = 0
    call run_case(13, -3.0_dp, 3.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified, 'wrong enclosure with valid ID was accepted')

    call run_case(10, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_UNRESOLVED .and. &
        .not. result%coverage_certified, 'tangent without certificate was accepted')

    call run_case(15, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified, 'callback failure did not dominate unresolved')

    options%initial_partition = 2
    call run_case(16, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_SUCCESS .and. result%nroots == 1, &
        'shared endpoint root was duplicated or lost')
    call require(result%roots(1)%lo == 0.0_dp .and. result%roots(1)%hi == 0.0_dp .and. &
        result%roots(1)%left_endpoint_root .and. result%roots(1)%right_endpoint_root, &
        'shared endpoint ownership was not canonicalized')
    options%initial_partition = 1

    options%max_stationary_iterations = 1
    options%max_depth = 0
    call run_case(17, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_UNRESOLVED .and. &
        .not. result%coverage_certified, 'wide root escaped the x_tolerance gate')
    options%max_stationary_iterations = 160
    options%max_depth = 80

    call run_case(17, -1.0_dp, 1.0_dp, coarse)
    options%initial_partition = 8
    call run_case(17, -1.0_dp, 1.0_dp, fine)
    call require(coarse%status == GC_INTERVAL_ROOT_SUCCESS .and. &
        fine%status == GC_INTERVAL_ROOT_SUCCESS, 'partition-invariance setup failed')
    call require(partition_invariant_enclosures(coarse, fine, options%x_tolerance), &
        'simple-root enclosure depends on the initial partition')
    options%initial_partition = 1

    call run_case(4, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_SUCCESS .and. result%nroots == 2, &
        'nearby simple roots were merged or missed')

    options%max_depth = 8
    call run_case(5, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_UNRESOLVED .and. &
        .not. result%coverage_certified, 'identically-zero interval was accepted')
    options%max_depth = 80

    call run_case(6, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified, 'non-certified callback was accepted')

    call run_case(7, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified, 'callback cut identity drift was accepted')

    call run_case(8, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified, 'missing enclosure ID was accepted')

    call run_case(11, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified, 'unregistered enclosure ID was accepted')

    call run_case(9, -1.0_dp, 1.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified, 'callback interval identity was accepted')

    options%max_boxes = 1
    call run_case(2, -3.0_dp, 3.0_dp, result)
    call require(result%status == GC_INTERVAL_ROOT_CAPACITY, 'capacity failure was not explicit')
    options%max_boxes = 20000

    call run_case(2, -3.0_dp, 3.0_dp, coarse)
    options%initial_partition = 8
    call run_case(2, -3.0_dp, 3.0_dp, fine)
    call run_case(2, -3.0_dp, 3.0_dp, repeat)
    call require(coarse%status == GC_INTERVAL_ROOT_SUCCESS .and. &
        fine%status == GC_INTERVAL_ROOT_SUCCESS .and. &
        repeat%status == GC_INTERVAL_ROOT_SUCCESS, 'determinism setup failed')
    call require(coarse%nroots == fine%nroots, 'root count depends on initial partition')
    call require(all_kinds_equal(coarse, fine), 'root kind depends on initial partition')
    call require(same_canonical_evidence(fine, repeat), &
        'full canonical evidence is not repeatable')

    write (*, '(A)') 'test_gc_certified_interval_roots OK'

contains

    subroutine run_case(mode, lo, hi, output)
        integer, intent(in) :: mode
        real(dp), intent(in) :: lo, hi
        type(gc_interval_root_result_t), intent(out) :: output
        callback_mode = mode
        callback_calls = 0
        call isolate_gc_interval_roots(analytic_callback, analytic_enclosure_verifier, &
            analytic_stationary_verifier, lo, hi, options, output)
    end subroutine run_case

    subroutine analytic_callback(lo, hi, value)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(out) :: value
        real(dp) :: a, b, x

        callback_calls = callback_calls + 1
        value = gc_interval_callback_result_t()
        value%query_lo = lo
        value%query_hi = hi
        value%status = 0
        value%certified = .true.
        value%cut_id = 17
        value%enclosure_certificate_id = 901

        if (callback_mode == 15) then
            if (callback_calls == 1) then
                value%f = gc_interval_t(0.0_dp, 0.0_dp)
                value%df = gc_interval_t(0.0_dp, 0.0_dp)
                value%d2f = gc_interval_t(0.0_dp, 0.0_dp)
            else
                value%status = 99
            end if
            return
        end if

        select case (callback_mode)
        case (1, 12, 13)
            value%f = square_plus_constant(lo, hi, 1.0_dp)
            value%df = gc_interval_t(down(2.0_dp*lo), up(2.0_dp*hi))
            value%d2f = gc_interval_t(2.0_dp, 2.0_dp)
        case (2)
            value%f = quadratic_interval(lo, hi, 1.0_dp, 1.0_dp, -2.0_dp)
            value%df = gc_interval_t(down(2.0_dp*lo + 1.0_dp), &
                up(2.0_dp*hi + 1.0_dp))
            value%d2f = gc_interval_t(2.0_dp, 2.0_dp)
        case (3, 10, 14)
            value%f = square_plus_constant(lo - 0.3_dp, hi - 0.3_dp, 0.0_dp)
            value%df = gc_interval_t(down(2.0_dp*(lo - 0.3_dp)), &
                up(2.0_dp*(hi - 0.3_dp)))
            value%d2f = gc_interval_t(2.0_dp, 2.0_dp)
            if (callback_mode == 3) then
                value%stationary_certificate_id = 303
                value%stationary_point = 0.3_dp
            else if (callback_mode == 14) then
                value%stationary_certificate_id = 303
                value%stationary_point = ieee_next_after(0.3_dp, huge(0.3_dp))
            end if
        case (4)
            a = 0.2_dp
            b = 0.21_dp
            value%f = quadratic_interval(lo, hi, 1.0_dp, -(a + b), a*b)
            value%df = gc_interval_t(down(2.0_dp*lo - (a + b)), &
                up(2.0_dp*hi - (a + b)))
            value%d2f = gc_interval_t(2.0_dp, 2.0_dp)
        case (5)
            value%f = gc_interval_t(0.0_dp, 0.0_dp)
            value%df = gc_interval_t(0.0_dp, 0.0_dp)
            value%d2f = gc_interval_t(0.0_dp, 0.0_dp)
        case (6)
            value%certified = .false.
            value%status = 99
        case (7)
            if (hi - lo < 1.5_dp) value%cut_id = 18
            value%f = gc_interval_t(0.0_dp, 0.0_dp)
            value%df = gc_interval_t(0.0_dp, 0.0_dp)
            value%d2f = gc_interval_t(0.0_dp, 0.0_dp)
        case (8)
            value%enclosure_certificate_id = 0
            value%f = square_plus_constant(lo, hi, 1.0_dp)
            value%df = gc_interval_t(down(2.0_dp*lo), up(2.0_dp*hi))
            value%d2f = gc_interval_t(2.0_dp, 2.0_dp)
        case (9)
            value%query_lo = lo + 1.0_dp
            value%f = square_plus_constant(lo, hi, 1.0_dp)
            value%df = gc_interval_t(down(2.0_dp*lo), up(2.0_dp*hi))
            value%d2f = gc_interval_t(2.0_dp, 2.0_dp)
        case (11)
            value%enclosure_certificate_id = 902
            value%f = square_plus_constant(lo, hi, 1.0_dp)
            value%df = gc_interval_t(down(2.0_dp*lo), up(2.0_dp*hi))
            value%d2f = gc_interval_t(2.0_dp, 2.0_dp)
        case (16)
            value%f = gc_interval_t(down(lo), up(hi))
            value%df = gc_interval_t(1.0_dp, 1.0_dp)
            value%d2f = gc_interval_t(0.0_dp, 0.0_dp)
        case (17)
            value%f = gc_interval_t(down(lo - 0.3_dp), up(hi - 0.3_dp))
            value%df = gc_interval_t(1.0_dp, 1.0_dp)
            value%d2f = gc_interval_t(0.0_dp, 0.0_dp)
        case default
            value%status = 99
        end select

        if (callback_mode == 12) value%certified = .false.

        if (lo == hi) then
            x = lo
            select case (callback_mode)
            case (1, 12, 13)
                value%f = point_interval(x*x + 1.0_dp)
            case (2)
                value%f = point_hull((x - 1.0_dp)*(x + 2.0_dp), &
                    x*x + x - 2.0_dp)
            case (3, 10, 14)
                value%f = point_interval((x - 0.3_dp)**2)
            case (4)
                value%f = point_hull((x - 0.2_dp)*(x - 0.21_dp), &
                    x*x - 0.41_dp*x + 0.2_dp*0.21_dp)
            case (16)
                value%f = point_interval(x)
            case (17)
                value%f = point_interval(x - 0.3_dp)
            end select
        end if
        if (callback_mode == 13) value%f%lo = value%f%lo + 1.0_dp
    end subroutine analytic_callback

    subroutine analytic_enclosure_verifier(lo, hi, value, expected_id, status)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(in) :: value
        integer, intent(in) :: expected_id
        integer, intent(out) :: status
        type(gc_interval_t) :: expected_f, expected_df, expected_d2f

        status = 1
        if (expected_id /= 901 .or. value%enclosure_certificate_id /= expected_id) return
        if (value%query_lo /= lo .or. value%query_hi /= hi) return
        if (.not. valid_test_interval(value%f) .or. &
                .not. valid_test_interval(value%df) .or. &
                .not. valid_test_interval(value%d2f)) return
        expected_f = reference_f(lo, hi, callback_mode)
        expected_df = reference_df(lo, hi, callback_mode)
        expected_d2f = reference_d2f(callback_mode)
        if (.not. contains_oracle_interval(value%f, expected_f) .or. &
                .not. contains_oracle_interval(value%df, expected_df) .or. &
                .not. contains_oracle_interval(value%d2f, expected_d2f)) return
        status = 0
    end subroutine analytic_enclosure_verifier

    subroutine analytic_stationary_verifier(lo, hi, point, value, expected_enclosure_id, &
            expected_stationary_id, status)
        real(dp), intent(in) :: lo, hi, point
        type(gc_interval_callback_result_t), intent(out) :: value
        integer, intent(in) :: expected_enclosure_id, expected_stationary_id
        integer, intent(out) :: status

        value = gc_interval_callback_result_t()
        value%query_lo = point
        value%query_hi = point
        value%cut_id = 17
        value%enclosure_certificate_id = expected_enclosure_id
        value%stationary_certificate_id = expected_stationary_id
        value%stationary_point = point
        value%status = 1
        status = 1
        if (expected_enclosure_id /= 901 .or. expected_stationary_id /= 303) return
        if (point < lo .or. point > hi) return
        if (point /= 0.3_dp) return
        value%f = point_interval((point - 0.3_dp)**2)
        value%df = point_interval(2.0_dp*(point - 0.3_dp))
        value%d2f = point_interval(2.0_dp)
        value%status = 0
        status = 0
    end subroutine analytic_stationary_verifier

    function reference_f(lo, hi, mode) result(value)
        real(dp), intent(in) :: lo, hi
        integer, intent(in) :: mode
        type(gc_interval_t) :: value
        select case (mode)
        case (1, 12, 13)
            value = reference_square(lo, hi, 1.0_dp)
        case (2)
            value = reference_quadratic(lo, hi, 1.0_dp, 1.0_dp, -2.0_dp)
        case (3, 10, 14)
            value = reference_square(lo - 0.3_dp, hi - 0.3_dp, 0.0_dp)
        case (4)
            value = reference_quadratic(lo, hi, 1.0_dp, -0.41_dp, 0.2_dp*0.21_dp)
        case (5, 15)
            value = gc_interval_t(0.0_dp, 0.0_dp)
        case (7)
            value = gc_interval_t(0.0_dp, 0.0_dp)
        case (8, 9, 11)
            value = reference_square(lo, hi, 1.0_dp)
        case (16)
            if (lo == hi) then
                value = reference_point_interval(lo)
            else
                value = gc_interval_t(down(lo), up(hi))
            end if
        case (17)
            if (lo == hi) then
                value = reference_point_interval(lo - 0.3_dp)
            else
                value = gc_interval_t(down(lo - 0.3_dp), up(hi - 0.3_dp))
            end if
        case default
            value = gc_interval_t(0.0_dp, 0.0_dp)
        end select
    end function reference_f

    function reference_df(lo, hi, mode) result(value)
        real(dp), intent(in) :: lo, hi
        integer, intent(in) :: mode
        type(gc_interval_t) :: value
        select case (mode)
        case (1, 8, 9, 11, 12, 13)
            value = gc_interval_t(down(2.0_dp*lo), up(2.0_dp*hi))
        case (2)
            value = gc_interval_t(down(2.0_dp*lo + 1.0_dp), up(2.0_dp*hi + 1.0_dp))
        case (3, 10, 14)
            value = gc_interval_t(down(2.0_dp*(lo - 0.3_dp)), &
                up(2.0_dp*(hi - 0.3_dp)))
        case (4)
            value = gc_interval_t(down(2.0_dp*lo - 0.41_dp), &
                up(2.0_dp*hi - 0.41_dp))
        case (5, 7, 15)
            value = gc_interval_t(0.0_dp, 0.0_dp)
        case (16, 17)
            value = gc_interval_t(1.0_dp, 1.0_dp)
        case default
            value = gc_interval_t(0.0_dp, 0.0_dp)
        end select
    end function reference_df

    function reference_d2f(mode) result(value)
        integer, intent(in) :: mode
        type(gc_interval_t) :: value
        select case (mode)
        case (1, 2, 3, 4, 8, 9, 10, 11, 12, 13, 14)
            value = gc_interval_t(2.0_dp, 2.0_dp)
        case default
            value = gc_interval_t(0.0_dp, 0.0_dp)
        end select
    end function reference_d2f

    function reference_square(lo, hi, constant) result(value)
        real(dp), intent(in) :: lo, hi, constant
        type(gc_interval_t) :: value
        real(dp) :: smallest, largest
        if (lo == hi) then
            value = reference_point_interval(lo*lo + constant)
            return
        end if
        smallest = min(abs(lo), abs(hi))
        largest = max(abs(lo), abs(hi))
        if (lo <= 0.0_dp .and. hi >= 0.0_dp) smallest = 0.0_dp
        value%lo = smallest*smallest + constant
        if (value%lo /= 0.0_dp) value%lo = down(value%lo)
        value%hi = up(largest*largest + constant)
    end function reference_square

    function reference_quadratic(lo, hi, qa, qb, qc) result(value)
        real(dp), intent(in) :: lo, hi, qa, qb, qc
        type(gc_interval_t) :: value
        real(dp) :: vertex, v1, v2, vv
        if (lo == hi) then
            value = reference_point_interval(qa*lo*lo + qb*lo + qc)
            return
        end if
        vertex = -qb/(2.0_dp*qa)
        v1 = qa*lo*lo + qb*lo + qc
        v2 = qa*hi*hi + qb*hi + qc
        value%lo = down(min(v1, v2))
        value%hi = up(max(v1, v2))
        if (lo <= vertex .and. vertex <= hi) then
            vv = qa*vertex*vertex + qb*vertex + qc
            value%lo = min(value%lo, down(vv))
            value%hi = max(value%hi, up(vv))
        end if
    end function reference_quadratic

    function reference_point_interval(number) result(value)
        real(dp), intent(in) :: number
        type(gc_interval_t) :: value
        if (number == 0.0_dp) then
            value = gc_interval_t(0.0_dp, 0.0_dp)
        else
            value = gc_interval_t(down(number), up(number))
        end if
    end function reference_point_interval

    function square_plus_constant(lo, hi, constant) result(value)
        real(dp), intent(in) :: lo, hi, constant
        type(gc_interval_t) :: value
        real(dp) :: smallest, largest
        smallest = min(abs(lo), abs(hi))
        largest = max(abs(lo), abs(hi))
        if (lo <= 0.0_dp .and. hi >= 0.0_dp) smallest = 0.0_dp
        value%lo = smallest*smallest + constant
        if (value%lo /= 0.0_dp) value%lo = down(value%lo)
        value%hi = up(largest*largest + constant)
    end function square_plus_constant

    function quadratic_interval(lo, hi, qa, qb, qc) result(value)
        real(dp), intent(in) :: lo, hi, qa, qb, qc
        type(gc_interval_t) :: value
        real(dp) :: vertex, v(2), vv
        vertex = -qb/(2.0_dp*qa)
        v = [qa*lo*lo + qb*lo + qc, qa*hi*hi + qb*hi + qc]
        value%lo = down(minval(v))
        value%hi = up(maxval(v))
        if (lo <= vertex .and. vertex <= hi) then
            vv = qa*vertex*vertex + qb*vertex + qc
            value%lo = min(value%lo, down(vv))
            value%hi = max(value%hi, up(vv))
        end if
    end function quadratic_interval

    function point_interval(number) result(value)
        real(dp), intent(in) :: number
        type(gc_interval_t) :: value
        if (number == 0.0_dp) then
            value = gc_interval_t(0.0_dp, 0.0_dp)
        else
            value = gc_interval_t(down(number), up(number))
        end if
    end function point_interval

    function point_hull(first, second) result(value)
        !! Enclose both independently evaluated forms of the manufactured
        !! polynomial. Floating-point distributivity is not assumed.
        real(dp), intent(in) :: first, second
        type(gc_interval_t) :: value
        if (first == 0.0_dp .and. second == 0.0_dp) then
            value = gc_interval_t(0.0_dp, 0.0_dp)
        else
            value = gc_interval_t(down(min(first, second)), &
                up(max(first, second)))
        end if
    end function point_hull

    real(dp) function down(number)
        real(dp), intent(in) :: number
        down = ieee_next_after(number, -huge(number))
    end function down

    real(dp) function up(number)
        real(dp), intent(in) :: number
        up = ieee_next_after(number, huge(number))
    end function up

    logical function valid_test_interval(interval)
        type(gc_interval_t), intent(in) :: interval
        valid_test_interval = ieee_is_finite(interval%lo) .and. &
            ieee_is_finite(interval%hi) .and. interval%lo <= interval%hi
    end function valid_test_interval

    logical function contains_oracle_interval(actual, expected)
        type(gc_interval_t), intent(in) :: actual, expected
        contains_oracle_interval = actual%lo <= expected%lo .and. &
            actual%hi >= expected%hi
    end function contains_oracle_interval

    logical function excludes_zero_for_test(interval)
        type(gc_interval_t), intent(in) :: interval
        excludes_zero_for_test = interval%hi < 0.0_dp .or. interval%lo > 0.0_dp
    end function excludes_zero_for_test

    logical function all_widths_within(value, tolerance)
        type(gc_interval_root_result_t), intent(in) :: value
        real(dp), intent(in) :: tolerance
        integer :: k
        all_widths_within = .true.
        do k = 1, value%nroots
            if (value%roots(k)%hi - value%roots(k)%lo > tolerance) then
                all_widths_within = .false.
            end if
        end do
    end function all_widths_within

    logical function partition_invariant_enclosures(left, right, tolerance)
        type(gc_interval_root_result_t), intent(in) :: left, right
        real(dp), intent(in) :: tolerance
        integer :: k

        partition_invariant_enclosures = left%nroots == right%nroots
        if (.not. partition_invariant_enclosures) return
        do k = 1, left%nroots
            if (left%roots(k)%kind /= right%roots(k)%kind) then
                partition_invariant_enclosures = .false.
            end if
            if (left%roots(k)%hi - left%roots(k)%lo > tolerance) then
                partition_invariant_enclosures = .false.
            end if
            if (right%roots(k)%hi - right%roots(k)%lo > tolerance) then
                partition_invariant_enclosures = .false.
            end if
            if (left%roots(k)%hi < right%roots(k)%lo) then
                partition_invariant_enclosures = .false.
            end if
            if (right%roots(k)%hi < left%roots(k)%lo) then
                partition_invariant_enclosures = .false.
            end if
        end do
    end function partition_invariant_enclosures

    logical function all_simple_evidence(value)
        type(gc_interval_root_result_t), intent(in) :: value
        integer :: k
        all_simple_evidence = .true.
        do k = 1, value%nroots
            all_simple_evidence = all_simple_evidence .and. &
                value%roots(k)%kind == GC_INTERVAL_ROOT_SIMPLE .and. &
                value%roots(k)%multiplicity_lower == 1 .and. &
                value%roots(k)%multiplicity_upper == 1 .and. &
                value%roots(k)%derivative_excludes_zero .and. &
                value%roots(k)%transversality_certified .and. &
                value%roots(k)%classification_certified .and. &
                (value%roots(k)%bracket_certified .or. &
                    value%roots(k)%interval_newton_certified) .and. &
                value%roots(k)%enclosure_certificate_id == 901
        end do
    end function all_simple_evidence

    logical function all_kinds_equal(left, right)
        type(gc_interval_root_result_t), intent(in) :: left, right
        integer :: k
        all_kinds_equal = left%nroots == right%nroots
        if (.not. all_kinds_equal) return
        do k = 1, left%nroots
            all_kinds_equal = all_kinds_equal .and. &
                left%roots(k)%kind == right%roots(k)%kind
        end do
    end function all_kinds_equal

    logical function same_canonical_evidence(left, right)
        type(gc_interval_root_result_t), intent(in) :: left, right
        integer :: k
        same_canonical_evidence = left%status == right%status .and. &
            left%nroots == right%nroots .and. &
            (left%coverage_certified .eqv. right%coverage_certified)
        if (.not. same_canonical_evidence) return
        do k = 1, left%nroots
            same_canonical_evidence = same_canonical_evidence .and. &
                same_root_evidence(left%roots(k), right%roots(k))
        end do
    end function same_canonical_evidence

    logical function same_root_evidence(left, right)
        type(gc_interval_root_box_t), intent(in) :: left, right
        same_root_evidence = left%lo == right%lo .and. left%hi == right%hi .and. &
            left%kind == right%kind .and. left%cut_id == right%cut_id .and. &
            left%boundary_role == right%boundary_role .and. &
            left%multiplicity_lower == right%multiplicity_lower .and. &
            left%multiplicity_upper == right%multiplicity_upper .and. &
            (left%derivative_excludes_zero .eqv. right%derivative_excludes_zero) .and. &
            (left%stationary_certified .eqv. right%stationary_certified) .and. &
            (left%bracket_certified .eqv. right%bracket_certified) .and. &
            (left%classification_certified .eqv. right%classification_certified) .and. &
            (left%transversality_certified .eqv. right%transversality_certified) .and. &
            left%transversality_kind == right%transversality_kind .and. &
            (left%interval_newton_certified .eqv. right%interval_newton_certified) .and. &
            (left%left_endpoint_root .eqv. right%left_endpoint_root) .and. &
            (left%right_endpoint_root .eqv. right%right_endpoint_root) .and. &
            left%enclosure_certificate_id == right%enclosure_certificate_id .and. &
            left%stationary_certificate_id == right%stationary_certificate_id .and. &
            left%derivative_enclosure%lo == right%derivative_enclosure%lo .and. &
            left%derivative_enclosure%hi == right%derivative_enclosure%hi .and. &
            left%second_derivative_enclosure%lo == right%second_derivative_enclosure%lo .and. &
            left%second_derivative_enclosure%hi == right%second_derivative_enclosure%hi
    end function same_root_evidence

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) error stop message
    end subroutine require

end program test_gc_certified_interval_roots
