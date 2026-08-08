module neort_gc_certified_canonical_variation
    !! Certified total variation of a canonical coordinate on one allowed cut.
    !!
    !! The derivative of psi_star is treated as the scalar root function.  Its
    !! certified root boxes partition the allowed component into monotone
    !! spans.  Variation on a span is enclosed from endpoint values; variation
    !! inside a root box is bounded by width times the supremum of the
    !! certified derivative enclosure.  All additions are rounded outwards.
    !!
    !! Root isolation remains the responsibility of
    !! neort_gc_certified_interval_roots.  This module consumes its evidence
    !! without weakening it and also exposes the assembly step for callers
    !! that have already isolated the derivative roots.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_certified_interval_roots, only: &
        GC_INTERVAL_ROOT_EXTREMAL, GC_INTERVAL_ROOT_SIMPLE, GC_INTERVAL_ROOT_TANGENT, &
        GC_INTERVAL_ROOT_SUCCESS, GC_INTERVAL_ROOT_TRANSVERSE, gc_interval_callback_result_t, &
        gc_interval_enclosure_verifier_i, gc_interval_root_box_t, &
        gc_interval_root_callback_i, gc_interval_root_options_t, &
        gc_interval_root_result_t, gc_interval_stationary_verifier_i, &
        gc_interval_t, isolate_gc_interval_roots

    implicit none
    private

    integer, parameter, public :: GC_CANONICAL_VARIATION_SUCCESS = 0
    integer, parameter, public :: GC_CANONICAL_VARIATION_INVALID_INPUT = 1
    integer, parameter, public :: GC_CANONICAL_VARIATION_ROOT_FAILURE = 2
    integer, parameter, public :: GC_CANONICAL_VARIATION_CALLBACK_FAILURE = 3
    integer, parameter, public :: GC_CANONICAL_VARIATION_INVALID_ROOT_EVIDENCE = 4
    integer, parameter, public :: GC_CANONICAL_VARIATION_NONMONOTONE_SPAN = 5
    integer, parameter, public :: GC_CANONICAL_VARIATION_NONFINITE = 6

    type, public :: gc_canonical_variation_options_t
        type(gc_interval_root_options_t) :: root_options
        integer :: expected_cut_id = 0
        integer :: expected_value_certificate_id = 0
    end type gc_canonical_variation_options_t

    type, public :: gc_canonical_variation_root_t
        type(gc_interval_root_box_t) :: evidence
        !! The returned root evidence enlarged by one representable number on
        !! each available side.  The enlargement closes the small pieces
        !! removed when a root endpoint is excluded from a monotone span.
        type(gc_interval_t) :: bound_coordinate
        type(gc_interval_t) :: derivative_enclosure
        type(gc_interval_t) :: variation_bound
        logical :: variation_certified = .false.
    end type gc_canonical_variation_root_t

    type, public :: gc_canonical_variation_span_t
        type(gc_interval_t) :: coordinate
        type(gc_interval_t) :: derivative_enclosure
        type(gc_interval_t) :: psi_left
        type(gc_interval_t) :: psi_right
        type(gc_interval_t) :: variation
        integer :: derivative_sign = 0
        logical :: monotonicity_certified = .false.
    end type gc_canonical_variation_span_t

    type, public :: gc_canonical_variation_result_t
        integer :: status = GC_CANONICAL_VARIATION_INVALID_INPUT
        logical :: certified = .false.
        integer :: nroots = 0
        integer :: nspans = 0
        type(gc_interval_root_result_t) :: derivative_roots
        type(gc_interval_t) :: total_variation_enclosure
        type(gc_interval_t) :: endpoint_difference_enclosure
        real(dp) :: total_variation_estimate = 0.0_dp
        real(dp) :: endpoint_difference_estimate = 0.0_dp
        type(gc_canonical_variation_root_t), allocatable :: roots(:)
        type(gc_canonical_variation_span_t), allocatable :: spans(:)
    end type gc_canonical_variation_result_t

    abstract interface
        subroutine gc_canonical_variation_callback_i(lo, hi, value)
            import dp, gc_interval_callback_result_t
            real(dp), intent(in) :: lo, hi
            type(gc_interval_callback_result_t), intent(out) :: value
        end subroutine gc_canonical_variation_callback_i
    end interface

    public :: gc_canonical_variation_callback_i
    public :: assemble_gc_canonical_total_variation
    public :: certify_gc_canonical_total_variation

contains

    subroutine certify_gc_canonical_total_variation(evaluate_canonical, verify_canonical, &
            evaluate_derivative, verify_derivative, verify_derivative_stationary, &
            domain_lo, domain_hi, options, result)
        procedure(gc_canonical_variation_callback_i) :: evaluate_canonical
        procedure(gc_interval_enclosure_verifier_i) :: verify_canonical
        procedure(gc_interval_root_callback_i) :: evaluate_derivative
        procedure(gc_interval_enclosure_verifier_i) :: verify_derivative
        procedure(gc_interval_stationary_verifier_i) :: verify_derivative_stationary
        real(dp), intent(in) :: domain_lo, domain_hi
        type(gc_canonical_variation_options_t), intent(in) :: options
        type(gc_canonical_variation_result_t), intent(out) :: result

        type(gc_interval_root_result_t) :: roots

        call initialize_result(result)
        if (.not. valid_variation_input(domain_lo, domain_hi, options)) then
            return
        end if

        call isolate_gc_interval_roots(evaluate_derivative, verify_derivative, &
            verify_derivative_stationary, domain_lo, domain_hi, options%root_options, roots)
        result%derivative_roots = roots
        if (roots%status /= GC_INTERVAL_ROOT_SUCCESS .or. &
                .not. roots%coverage_certified) then
            result%status = GC_CANONICAL_VARIATION_ROOT_FAILURE
            return
        end if

        call assemble_gc_canonical_total_variation(evaluate_canonical, verify_canonical, &
            domain_lo, domain_hi, options, roots, result)
    end subroutine certify_gc_canonical_total_variation

    subroutine assemble_gc_canonical_total_variation(evaluate_canonical, verify_canonical, &
            domain_lo, domain_hi, options, roots, result)
        procedure(gc_canonical_variation_callback_i) :: evaluate_canonical
        procedure(gc_interval_enclosure_verifier_i) :: verify_canonical
        real(dp), intent(in) :: domain_lo, domain_hi
        type(gc_canonical_variation_options_t), intent(in) :: options
        type(gc_interval_root_result_t), intent(in) :: roots
        type(gc_canonical_variation_result_t), intent(out) :: result

        type(gc_canonical_variation_root_t), allocatable :: root_data(:)
        type(gc_canonical_variation_span_t), allocatable :: span_data(:)
        type(gc_interval_callback_result_t) :: left_value, right_value
        type(gc_interval_t) :: total, endpoint_difference, updated_total
        type(gc_interval_t), allocatable :: bound_coordinates(:)
        real(dp) :: left, right
        integer :: i, nspans, local_status

        call initialize_result(result)
        result%derivative_roots = roots
        if (.not. valid_variation_input(domain_lo, domain_hi, options)) return
        if (.not. valid_root_result(roots, domain_lo, domain_hi, options)) then
            result%status = GC_CANONICAL_VARIATION_INVALID_ROOT_EVIDENCE
            return
        end if

        allocate(root_data(roots%nroots), bound_coordinates(roots%nroots))
        do i = 1, roots%nroots
            bound_coordinates(i) = expanded_root_box(roots%roots(i), domain_lo, domain_hi)
            if (.not. valid_interval(bound_coordinates(i))) then
                result%status = GC_CANONICAL_VARIATION_INVALID_ROOT_EVIDENCE
                return
            end if
            if (i > 1) then
                if (bound_coordinates(i)%lo <= bound_coordinates(i-1)%hi) then
                    result%status = GC_CANONICAL_VARIATION_INVALID_ROOT_EVIDENCE
                    return
                end if
            end if

            call evaluate_canonical_checked(evaluate_canonical, verify_canonical, &
                bound_coordinates(i)%lo, bound_coordinates(i)%hi, options, left_value, &
                local_status)
            if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
                result%status = local_status
                return
            end if
            root_data(i)%evidence = roots%roots(i)
            root_data(i)%bound_coordinate = bound_coordinates(i)
            root_data(i)%derivative_enclosure = left_value%df
            call root_box_variation_bound(bound_coordinates(i), left_value%df, &
                root_data(i)%variation_bound, local_status)
            if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
                result%status = local_status
                return
            end if
            root_data(i)%variation_certified = .true.
        end do

        allocate(span_data(roots%nroots+1))
        left = domain_lo
        nspans = 0
        do i = 1, roots%nroots
            right = bound_coordinates(i)%lo
            if (right > left) then
                call certify_monotone_span(evaluate_canonical, verify_canonical, left, &
                    right, options, span_data(nspans+1), local_status)
                if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
                    result%status = local_status
                    return
                end if
                nspans = nspans+1
            end if
            left = bound_coordinates(i)%hi
        end do
        right = domain_hi
        if (right > left) then
            call certify_monotone_span(evaluate_canonical, verify_canonical, left, right, &
                options, span_data(nspans+1), local_status)
            if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
                result%status = local_status
                return
            end if
            nspans = nspans+1
        end if
        if (left > domain_hi) then
            result%status = GC_CANONICAL_VARIATION_INVALID_ROOT_EVIDENCE
            return
        end if

        call evaluate_canonical_checked(evaluate_canonical, verify_canonical, domain_lo, &
            domain_lo, options, left_value, local_status)
        if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
            result%status = local_status
            return
        end if
        call evaluate_canonical_checked(evaluate_canonical, verify_canonical, domain_hi, &
            domain_hi, options, right_value, local_status)
        if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
            result%status = local_status
            return
        end if
        endpoint_difference = interval_abs_difference(right_value%f, left_value%f, &
            local_status)
        if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
            result%status = local_status
            return
        end if

        total = gc_interval_t(0.0_dp, 0.0_dp)
        do i = 1, nspans
            call add_nonnegative_intervals(total, span_data(i)%variation, updated_total, &
                local_status)
            if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
                result%status = local_status
                return
            end if
            total = updated_total
        end do
        do i = 1, roots%nroots
            call add_nonnegative_intervals(total, root_data(i)%variation_bound, updated_total, &
                local_status)
            if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
                result%status = local_status
                return
            end if
            total = updated_total
        end do

        if (.not. valid_interval(total) .or. .not. valid_interval(endpoint_difference) .or. &
                .not. ieee_is_finite(interval_midpoint(total)) .or. &
                .not. ieee_is_finite(interval_midpoint(endpoint_difference))) then
            result%status = GC_CANONICAL_VARIATION_NONFINITE
            return
        end if
        result%roots = root_data
        result%spans = span_data(:nspans)
        result%nroots = roots%nroots
        result%nspans = nspans
        result%total_variation_enclosure = total
        result%endpoint_difference_enclosure = endpoint_difference
        result%total_variation_estimate = interval_midpoint(total)
        result%endpoint_difference_estimate = interval_midpoint(endpoint_difference)
        result%status = GC_CANONICAL_VARIATION_SUCCESS
        result%certified = .true.
    end subroutine assemble_gc_canonical_total_variation

    subroutine certify_monotone_span(evaluate, verify, left, right, options, span, status)
        procedure(gc_canonical_variation_callback_i) :: evaluate
        procedure(gc_interval_enclosure_verifier_i) :: verify
        real(dp), intent(in) :: left, right
        type(gc_canonical_variation_options_t), intent(in) :: options
        type(gc_canonical_variation_span_t), intent(out) :: span
        integer, intent(out) :: status

        type(gc_interval_callback_result_t) :: value, left_value, right_value
        integer :: local_status

        span = gc_canonical_variation_span_t()
        status = GC_CANONICAL_VARIATION_INVALID_INPUT
        if (.not. ieee_is_finite(left) .or. .not. ieee_is_finite(right)) return
        if (right <= left) return

        call evaluate_canonical_checked(evaluate, verify, left, right, options, value, &
            local_status)
        if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
            status = local_status
            return
        end if
        if (.not. excludes_zero(value%df)) then
            status = GC_CANONICAL_VARIATION_NONMONOTONE_SPAN
            return
        end if

        call evaluate_canonical_checked(evaluate, verify, left, left, options, left_value, &
            local_status)
        if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
            status = local_status
            return
        end if
        call evaluate_canonical_checked(evaluate, verify, right, right, options, right_value, &
            local_status)
        if (local_status /= GC_CANONICAL_VARIATION_SUCCESS) then
            status = local_status
            return
        end if

        span%coordinate = gc_interval_t(left, right)
        span%derivative_enclosure = value%df
        span%psi_left = left_value%f
        span%psi_right = right_value%f
        span%variation = interval_abs_difference(right_value%f, left_value%f, status)
        if (status /= GC_CANONICAL_VARIATION_SUCCESS) return
        span%derivative_sign = merge(1, -1, value%df%lo > 0.0_dp)
        span%monotonicity_certified = .true.
        status = GC_CANONICAL_VARIATION_SUCCESS
    end subroutine certify_monotone_span

    subroutine evaluate_canonical_checked(evaluate, verify, lo, hi, options, value, status)
        procedure(gc_canonical_variation_callback_i) :: evaluate
        procedure(gc_interval_enclosure_verifier_i) :: verify
        real(dp), intent(in) :: lo, hi
        type(gc_canonical_variation_options_t), intent(in) :: options
        type(gc_interval_callback_result_t), intent(out) :: value
        integer, intent(out) :: status

        integer :: verifier_status

        value = gc_interval_callback_result_t()
        status = GC_CANONICAL_VARIATION_CALLBACK_FAILURE
        call evaluate(lo, hi, value)
        if (value%status /= 0) return
        if (value%query_lo /= lo .or. value%query_hi /= hi) return
        if (value%cut_id /= options%expected_cut_id) return
        if (value%enclosure_certificate_id /= options%expected_value_certificate_id) return
        if (value%stationary_certificate_id < 0) return
        if (value%boundary_role < 0 .or. value%boundary_role > 3) return
        if (value%stationary_certificate_id > 0 .and. &
                .not. ieee_is_finite(value%stationary_point)) return
        if (.not. valid_interval(value%f)) return
        if (.not. valid_interval(value%df)) return
        if (.not. valid_interval(value%d2f)) return
        call verify(lo, hi, value, options%expected_value_certificate_id, verifier_status)
        if (verifier_status /= 0) return
        status = GC_CANONICAL_VARIATION_SUCCESS
    end subroutine evaluate_canonical_checked

    subroutine root_box_variation_bound(coordinate, derivative, bound, status)
        type(gc_interval_t), intent(in) :: coordinate, derivative
        type(gc_interval_t), intent(out) :: bound
        integer, intent(out) :: status

        real(dp) :: width, supremum, upper

        bound = gc_interval_t()
        status = GC_CANONICAL_VARIATION_NONFINITE
        if (.not. valid_interval(coordinate)) return
        if (.not. valid_interval(derivative)) return
        width = coordinate%hi - coordinate%lo
        if (.not. ieee_is_finite(width)) return
        supremum = max(abs(derivative%lo), abs(derivative%hi))
        if (.not. ieee_is_finite(supremum)) return
        upper = width*supremum
        if (.not. ieee_is_finite(upper)) return
        if (upper < 0.0_dp) return
        if (width == 0.0_dp .or. supremum == 0.0_dp) then
            bound = gc_interval_t(0.0_dp, 0.0_dp)
        else
            bound = gc_interval_t(0.0_dp, round_up(upper))
        end if
        status = GC_CANONICAL_VARIATION_SUCCESS
    end subroutine root_box_variation_bound

    function expanded_root_box(root, domain_lo, domain_hi) result(coordinate)
        type(gc_interval_root_box_t), intent(in) :: root
        real(dp), intent(in) :: domain_lo, domain_hi
        type(gc_interval_t) :: coordinate

        coordinate%lo = root%lo
        coordinate%hi = root%hi
        if (coordinate%lo > domain_lo) then
            coordinate%lo = ieee_next_after(coordinate%lo, -huge(coordinate%lo))
        end if
        if (coordinate%hi < domain_hi) then
            coordinate%hi = ieee_next_after(coordinate%hi, huge(coordinate%hi))
        end if
        coordinate%lo = max(coordinate%lo, domain_lo)
        coordinate%hi = min(coordinate%hi, domain_hi)
    end function expanded_root_box

    logical function valid_variation_input(domain_lo, domain_hi, options)
        real(dp), intent(in) :: domain_lo, domain_hi
        type(gc_canonical_variation_options_t), intent(in) :: options

        valid_variation_input = .false.
        if (.not. ieee_is_finite(domain_lo) .or. .not. ieee_is_finite(domain_hi)) return
        if (domain_hi <= domain_lo) return
        if (options%expected_cut_id <= 0) return
        if (options%expected_value_certificate_id <= 0) return
        if (options%root_options%expected_enclosure_certificate_id <= 0) return
        if (options%root_options%expected_stationary_certificate_id < 0) return
        valid_variation_input = .true.
    end function valid_variation_input

    logical function valid_root_result(roots, domain_lo, domain_hi, options)
        type(gc_interval_root_result_t), intent(in) :: roots
        real(dp), intent(in) :: domain_lo, domain_hi
        type(gc_canonical_variation_options_t), intent(in) :: options

        integer :: i
        type(gc_interval_root_box_t) :: root

        valid_root_result = .false.
        if (roots%status /= GC_INTERVAL_ROOT_SUCCESS) return
        if (.not. roots%coverage_certified) return
        if (roots%nroots < 0) return
        if (roots%nroots > options%root_options%max_roots) return
        if (.not. allocated(roots%roots)) return
        if (size(roots%roots) /= roots%nroots) return
        do i = 1, roots%nroots
            root = roots%roots(i)
            if (.not. ieee_is_finite(root%lo) .or. .not. ieee_is_finite(root%hi)) return
            if (root%lo > root%hi) return
            if (root%lo < domain_lo .or. root%hi > domain_hi) return
            if (root%cut_id /= options%expected_cut_id) return
            if (root%enclosure_certificate_id /= &
                    options%root_options%expected_enclosure_certificate_id) return
            if (root%stationary_certificate_id /= &
                    options%root_options%expected_stationary_certificate_id) return
            if (root%boundary_role < 0 .or. root%boundary_role > 3) return
            if (.not. root%classification_certified) return
            if (.not. valid_interval(root%derivative_enclosure)) return
            if (.not. valid_interval(root%second_derivative_enclosure)) return
            select case (root%kind)
            case (GC_INTERVAL_ROOT_SIMPLE)
                if (.not. root%derivative_excludes_zero) return
                if (.not. root%transversality_certified) return
                if (root%transversality_kind /= GC_INTERVAL_ROOT_TRANSVERSE) return
                if (root%multiplicity_lower /= 1 .or. root%multiplicity_upper /= 1) return
                if (.not. (root%bracket_certified .or. root%interval_newton_certified)) return
                if (.not. excludes_zero(root%derivative_enclosure)) return
            case (GC_INTERVAL_ROOT_TANGENT)
                if (.not. root%stationary_certified) return
                if (root%transversality_kind /= GC_INTERVAL_ROOT_EXTREMAL) return
                if (root%multiplicity_lower /= 2 .or. root%multiplicity_upper /= 2) return
                if (.not. (root%bracket_certified .or. root%interval_newton_certified)) return
                if (.not. contains_zero(root%derivative_enclosure)) return
                if (.not. excludes_zero(root%second_derivative_enclosure)) return
            case default
                return
            end select
            if (i > 1) then
                if (root%lo <= roots%roots(i-1)%hi) return
            end if
        end do
        valid_root_result = .true.
    end function valid_root_result

    function interval_abs_difference(right, left, status) result(value)
        type(gc_interval_t), intent(in) :: right, left
        integer, intent(out) :: status
        type(gc_interval_t) :: difference, value
        real(dp) :: raw_lo, raw_hi

        value = gc_interval_t()
        status = GC_CANONICAL_VARIATION_NONFINITE
        if (.not. valid_interval(right)) return
        if (.not. valid_interval(left)) return
        raw_lo = right%lo-left%hi
        raw_hi = right%hi-left%lo
        if (.not. ieee_is_finite(raw_lo) .or. .not. ieee_is_finite(raw_hi)) return
        difference%lo = round_down(raw_lo)
        difference%hi = round_up(raw_hi)
        if (.not. valid_interval(difference)) return
        if (difference%lo >= 0.0_dp) then
            value%lo = max(0.0_dp, difference%lo)
            value%hi = round_up(difference%hi)
        else if (difference%hi <= 0.0_dp) then
            value%lo = max(0.0_dp, round_down(-difference%hi))
            value%hi = round_up(-difference%lo)
        else
            value%lo = 0.0_dp
            value%hi = round_up(max(abs(difference%lo), abs(difference%hi)))
        end if
        if (.not. valid_interval(value)) return
        status = GC_CANONICAL_VARIATION_SUCCESS
    end function interval_abs_difference

    subroutine add_nonnegative_intervals(first, second, sum, status)
        type(gc_interval_t), intent(in) :: first, second
        type(gc_interval_t), intent(out) :: sum
        integer, intent(out) :: status
        real(dp) :: raw_lo, raw_hi

        sum = gc_interval_t()
        status = GC_CANONICAL_VARIATION_NONFINITE
        if (.not. valid_interval(first)) return
        if (.not. valid_interval(second)) return
        if (first%lo < 0.0_dp .or. second%lo < 0.0_dp) return
        raw_lo = first%lo+second%lo
        raw_hi = first%hi+second%hi
        if (.not. ieee_is_finite(raw_lo) .or. .not. ieee_is_finite(raw_hi)) return
        sum%lo = max(0.0_dp, round_down(raw_lo))
        sum%hi = round_up(raw_hi)
        if (.not. valid_interval(sum)) return
        status = GC_CANONICAL_VARIATION_SUCCESS
    end subroutine add_nonnegative_intervals

    subroutine initialize_result(result)
        type(gc_canonical_variation_result_t), intent(out) :: result

        result = gc_canonical_variation_result_t()
        allocate(result%roots(0), result%spans(0))
        result%total_variation_enclosure = gc_interval_t(0.0_dp, 0.0_dp)
        result%endpoint_difference_enclosure = gc_interval_t(0.0_dp, 0.0_dp)
    end subroutine initialize_result

    logical function valid_interval(interval)
        type(gc_interval_t), intent(in) :: interval

        valid_interval = .false.
        if (.not. ieee_is_finite(interval%lo)) return
        if (.not. ieee_is_finite(interval%hi)) return
        if (interval%lo > interval%hi) return
        valid_interval = .true.
    end function valid_interval

    logical function excludes_zero(interval)
        type(gc_interval_t), intent(in) :: interval

        excludes_zero = interval%hi < 0.0_dp .or. interval%lo > 0.0_dp
    end function excludes_zero

    logical function contains_zero(interval)
        type(gc_interval_t), intent(in) :: interval

        contains_zero = interval%lo <= 0.0_dp .and. interval%hi >= 0.0_dp
    end function contains_zero

    real(dp) function round_down(number)
        real(dp), intent(in) :: number

        round_down = ieee_next_after(number, -huge(number))
    end function round_down

    real(dp) function round_up(number)
        real(dp), intent(in) :: number

        round_up = ieee_next_after(number, huge(number))
    end function round_up

    real(dp) function interval_midpoint(interval)
        type(gc_interval_t), intent(in) :: interval

        interval_midpoint = 0.5_dp*interval%lo+0.5_dp*interval%hi
    end function interval_midpoint

end module neort_gc_certified_canonical_variation
