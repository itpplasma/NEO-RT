module neort_gc_certified_interval_roots
    !! Fail-closed isolation of scalar roots from rigorous interval callbacks.
    !!
    !! This module deliberately contains no model-specific algebra.  The
    !! callback is the boundary between this algorithm and a generated
    !! (for example Fortsym-generated) field/kernel evaluator.  Every
    !! enclosure used to discard a box or certify a root must be accepted by
    !! the pinned outward-enclosure verifier.
    !!
    !! A sampled endpoint value is not an interval certificate.  The
    !! callback must return outward enclosures on every requested interval,
    !! including point intervals used by interval Newton steps.  Acceptance
    !! is delegated to a required verifier procedure pinned to a registry
    !! certificate ID; the callback's certified flag is diagnostic only.
    !! Production verifiers must be generated directed-rounding/interval
    !! kernels or their independently reviewed registry adapters.
    !!
    !! A stationary point plus 0 in an enclosure for f is only a candidate:
    !! it does not prove f(stationary point)=0.  Tangent roots therefore also
    !! require a nonzero stationary certificate ID, an exact stationary point,
    !! and exact generated-kernel evaluations f(point)=f'(point)=0.  A loose
    !! boolean assertion is deliberately not part of this contract.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_next_after

    implicit none
    private

    integer, parameter, public :: GC_INTERVAL_ROOT_SUCCESS = 0
    integer, parameter, public :: GC_INTERVAL_ROOT_INVALID_INPUT = 1
    integer, parameter, public :: GC_INTERVAL_ROOT_CALLBACK_FAILURE = 2
    integer, parameter, public :: GC_INTERVAL_ROOT_UNRESOLVED = 3
    integer, parameter, public :: GC_INTERVAL_ROOT_CAPACITY = 4
    integer, parameter, public :: GC_INTERVAL_ROOT_NONFINITE = 5

    integer, parameter, public :: GC_INTERVAL_ROOT_SIMPLE = 1
    integer, parameter, public :: GC_INTERVAL_ROOT_TANGENT = 2
    integer, parameter, public :: GC_INTERVAL_ROOT_INTERIOR = 0
    integer, parameter, public :: GC_INTERVAL_ROOT_LEFT_BOUNDARY = 1
    integer, parameter, public :: GC_INTERVAL_ROOT_RIGHT_BOUNDARY = 2
    integer, parameter, public :: GC_INTERVAL_ROOT_BOTH_BOUNDARIES = 3
    integer, parameter, public :: GC_INTERVAL_ROOT_TRANSVERSE = 1
    integer, parameter, public :: GC_INTERVAL_ROOT_EXTREMAL = 2

    type, public :: gc_interval_t
        real(dp) :: lo = 0.0_dp
        real(dp) :: hi = 0.0_dp
    end type gc_interval_t

    type, public :: gc_interval_callback_result_t
        type(gc_interval_t) :: f
        type(gc_interval_t) :: df
        type(gc_interval_t) :: d2f
        !! Echoes of the requested interval.  They make accidental callback
        !! state/cache mix-ups detectable instead of silently accepted.
        real(dp) :: query_lo = 0.0_dp
        real(dp) :: query_hi = 0.0_dp
        integer :: status = 1
        !! Diagnostic only.  Production acceptance never consults this bit.
        logical :: certified = .false.
        !! A cut/branch label is provenance, not an invariant label.  Roots
        !! with different labels are never coalesced, even when their
        !! coordinate intervals or mapped invariants coincide.
        integer :: cut_id = 0
        integer :: boundary_role = GC_INTERVAL_ROOT_INTERIOR
        !! Nonzero ID of the generated outward-enclosure certificate.  The
        !! ID must remain unchanged for every child and point query.
        integer :: enclosure_certificate_id = 0
        !! Optional analytic stationary-root/factor certificate.  It is
        !! checked only through the exact point data below, never as a bare
        !! logical assertion.
        integer :: stationary_certificate_id = 0
        real(dp) :: stationary_point = 0.0_dp
        !! Required only to certify an even/tangent root.  It must be an
        !! analytic proof attached to the stationary-value evaluation, not
        !! an assertion based on a floating-point tolerance.
    end type gc_interval_callback_result_t

    abstract interface
        subroutine gc_interval_root_callback_i(lo, hi, result)
            import dp, gc_interval_callback_result_t
            real(dp), intent(in) :: lo, hi
            type(gc_interval_callback_result_t), intent(out) :: result
        end subroutine gc_interval_root_callback_i

        subroutine gc_interval_enclosure_verifier_i(lo, hi, value, expected_id, status)
            import dp, gc_interval_callback_result_t
            real(dp), intent(in) :: lo, hi
            type(gc_interval_callback_result_t), intent(in) :: value
            integer, intent(in) :: expected_id
            integer, intent(out) :: status
        end subroutine gc_interval_enclosure_verifier_i

        subroutine gc_interval_stationary_verifier_i(lo, hi, point, value, &
                expected_enclosure_id, expected_stationary_id, status)
            import dp, gc_interval_callback_result_t
            real(dp), intent(in) :: lo, hi, point
            type(gc_interval_callback_result_t), intent(out) :: value
            integer, intent(in) :: expected_enclosure_id
            integer, intent(in) :: expected_stationary_id
            integer, intent(out) :: status
        end subroutine gc_interval_stationary_verifier_i
    end interface

    type, public :: gc_interval_root_options_t
        integer :: initial_partition = 1
        integer :: max_depth = 64
        integer :: max_boxes = 100000
        integer :: max_roots = 4096
        integer :: max_stationary_iterations = 128
        !! Registry IDs supplied by the generated Fortsym kernel package.
        !! The first callback is checked against these values; they are not
        !! learned from callback output.
        integer :: expected_enclosure_certificate_id = 0
        integer :: expected_stationary_certificate_id = 0
        real(dp) :: x_tolerance = 1.0e-10_dp
    end type gc_interval_root_options_t

    type, public :: gc_interval_root_box_t
        real(dp) :: lo = 0.0_dp
        real(dp) :: hi = 0.0_dp
        integer :: kind = 0
        integer :: cut_id = 0
        integer :: boundary_role = GC_INTERVAL_ROOT_INTERIOR
        integer :: multiplicity_lower = 0
        integer :: multiplicity_upper = 0
        logical :: derivative_excludes_zero = .false.
        logical :: stationary_certified = .false.
        logical :: bracket_certified = .false.
        logical :: classification_certified = .false.
        logical :: transversality_certified = .false.
        integer :: transversality_kind = 0
        logical :: interval_newton_certified = .false.
        logical :: left_endpoint_root = .false.
        logical :: right_endpoint_root = .false.
        integer :: enclosure_certificate_id = 0
        integer :: stationary_certificate_id = 0
        type(gc_interval_t) :: derivative_enclosure
        type(gc_interval_t) :: second_derivative_enclosure
    end type gc_interval_root_box_t

    type, public :: gc_interval_root_result_t
        integer :: status = GC_INTERVAL_ROOT_INVALID_INPUT
        integer :: nroots = 0
        integer :: boxes_visited = 0
        integer :: boxes_discarded = 0
        integer :: unresolved_boxes = 0
        logical :: coverage_certified = .false.
        type(gc_interval_root_box_t), allocatable :: roots(:)
    end type gc_interval_root_result_t

    public :: isolate_gc_interval_roots, gc_interval_root_callback_i, &
        gc_interval_enclosure_verifier_i, gc_interval_stationary_verifier_i

    type :: work_box_t
        real(dp) :: lo = 0.0_dp
        real(dp) :: hi = 0.0_dp
        integer :: depth = 0
        integer :: cut_id = 0
        integer :: enclosure_certificate_id = 0
        integer :: stationary_certificate_id = 0
        logical :: identity_initialized = .false.
    end type work_box_t

contains

    subroutine isolate_gc_interval_roots(evaluate, verify, verify_stationary, domain_lo, &
            domain_hi, options, result)
        procedure(gc_interval_root_callback_i) :: evaluate
        procedure(gc_interval_enclosure_verifier_i) :: verify
        procedure(gc_interval_stationary_verifier_i) :: verify_stationary
        real(dp), intent(in) :: domain_lo, domain_hi
        type(gc_interval_root_options_t), intent(in) :: options
        type(gc_interval_root_result_t), intent(out) :: result

        type(work_box_t), allocatable :: queue(:)
        type(gc_interval_root_box_t), allocatable :: candidates(:)
        type(gc_interval_callback_result_t) :: value
        integer :: head, tail, k, ninitial, istat
        real(dp) :: width, lo, hi
        logical :: valid, resolved, found, needs_split
        type(gc_interval_root_box_t) :: root

        result = gc_interval_root_result_t()
        allocate(result%roots(0))

        if (.not. valid_options(domain_lo, domain_hi, options)) then
            result%status = GC_INTERVAL_ROOT_INVALID_INPUT
            return
        end if
        if (options%initial_partition > options%max_boxes) then
            result%status = GC_INTERVAL_ROOT_CAPACITY
            return
        end if

        allocate(queue(options%max_boxes), stat=istat)
        if (istat /= 0) then
            result%status = GC_INTERVAL_ROOT_CAPACITY
            return
        end if
        allocate(candidates(options%max_roots), stat=istat)
        if (istat /= 0) then
            result%status = GC_INTERVAL_ROOT_CAPACITY
            return
        end if

        ninitial = options%initial_partition
        width = (domain_hi - domain_lo)/real(ninitial, dp)
        head = 1
        tail = ninitial
        do k = 1, ninitial
            queue(k)%lo = domain_lo + real(k - 1, dp)*width
            if (k == ninitial) then
                queue(k)%hi = domain_hi
            else
                queue(k)%hi = domain_lo + real(k, dp)*width
            end if
            queue(k)%depth = 0
            queue(k)%enclosure_certificate_id = options%expected_enclosure_certificate_id
            queue(k)%stationary_certificate_id = options%expected_stationary_certificate_id
            !! Cut/certificate identity comes only from the callback; the
            !! initial partition index is never used as provenance.
            queue(k)%identity_initialized = .false.
        end do
        result%status = GC_INTERVAL_ROOT_SUCCESS

        do while (head <= tail)
            lo = queue(head)%lo
            hi = queue(head)%hi
            result%boxes_visited = result%boxes_visited + 1

            call evaluate_checked(evaluate, verify, lo, hi, value, valid, &
                queue(head)%cut_id, queue(head)%enclosure_certificate_id, &
                queue(head)%stationary_certificate_id, queue(head)%identity_initialized, &
                options%expected_enclosure_certificate_id, options%expected_stationary_certificate_id)
            if (.not. valid) then
                result%unresolved_boxes = result%unresolved_boxes + 1
                call promote_status(result%status, callback_failure_status(value))
                head = head + 1
                cycle
            end if
            if (.not. queue(head)%identity_initialized) then
                queue(head)%cut_id = value%cut_id
                queue(head)%enclosure_certificate_id = options%expected_enclosure_certificate_id
                queue(head)%stationary_certificate_id = options%expected_stationary_certificate_id
                queue(head)%identity_initialized = .true.
            end if

            !! This is the only unconditional no-root discard in the module.
            if (.not. contains_zero(value%f)) then
                result%boxes_discarded = result%boxes_discarded + 1
                head = head + 1
                cycle
            end if

            found = .false.
            needs_split = .true.
            if (.not. contains_zero(value%df)) then
                call certify_simple_root(evaluate, verify, lo, hi, value, options, queue(head)%cut_id, &
                    queue(head)%enclosure_certificate_id, queue(head)%stationary_certificate_id, root, &
                    found, valid)
                needs_split = .not. found
            else if (.not. contains_zero(value%d2f)) then
                call certify_tangent_root(evaluate, verify, verify_stationary, lo, hi, value, options, &
                    queue(head)%cut_id, &
                    queue(head)%enclosure_certificate_id, queue(head)%stationary_certificate_id, root, &
                    found, valid)
                needs_split = .not. found
            end if

            if (found) then
                call append_candidate(candidates, result%nroots, root, options, resolved)
                if (.not. resolved) then
                    call promote_status(result%status, GC_INTERVAL_ROOT_CAPACITY)
                    result%unresolved_boxes = result%unresolved_boxes + 1
                end if
                head = head + 1
                cycle
            end if

            if (.not. valid) then
                call promote_status(result%status, GC_INTERVAL_ROOT_CALLBACK_FAILURE)
                result%unresolved_boxes = result%unresolved_boxes + 1
                head = head + 1
                cycle
            end if
            if (.not. needs_split .and. valid) then
                !! The current box was not certified as a root, but it also
                !! must not be discarded while its f enclosure contains zero.
                needs_split = .true.
            end if
            if (queue(head)%depth >= options%max_depth .or. &
                    hi - lo <= options%x_tolerance) then
                result%unresolved_boxes = result%unresolved_boxes + 1
                call promote_status(result%status, GC_INTERVAL_ROOT_UNRESOLVED)
                head = head + 1
                cycle
            end if
            if (tail + 2 > size(queue)) then
                call promote_status(result%status, GC_INTERVAL_ROOT_CAPACITY)
                result%unresolved_boxes = result%unresolved_boxes + 1
                head = head + 1
                cycle
            end if
            call split_box(queue, head, tail)
        end do

        if (result%unresolved_boxes > 0) then
            call promote_status(result%status, GC_INTERVAL_ROOT_UNRESOLVED)
        end if

        call order_and_validate_candidates(candidates, result%nroots, result, options)
        result%coverage_certified = result%status == GC_INTERVAL_ROOT_SUCCESS
        if (.not. result%coverage_certified) then
            !! A partial list is diagnostic only and must not be consumed as
            !! a certified root set.
            deallocate(result%roots)
            allocate(result%roots(0))
            result%nroots = 0
        end if
    end subroutine isolate_gc_interval_roots

    logical function valid_options(domain_lo, domain_hi, options)
        real(dp), intent(in) :: domain_lo, domain_hi
        type(gc_interval_root_options_t), intent(in) :: options
        valid_options = .true.
        valid_options = valid_options .and. ieee_is_finite(domain_lo)
        valid_options = valid_options .and. ieee_is_finite(domain_hi)
        valid_options = valid_options .and. domain_hi > domain_lo
        valid_options = valid_options .and. options%initial_partition >= 1
        valid_options = valid_options .and. options%max_depth >= 0
        valid_options = valid_options .and. options%max_boxes >= options%initial_partition
        valid_options = valid_options .and. options%max_roots >= 1
        valid_options = valid_options .and. options%max_stationary_iterations >= 1
        valid_options = valid_options .and. options%expected_enclosure_certificate_id > 0
        valid_options = valid_options .and. options%expected_stationary_certificate_id >= 0
        valid_options = valid_options .and. ieee_is_finite(options%x_tolerance)
        valid_options = valid_options .and. options%x_tolerance > 0.0_dp
    end function valid_options

    subroutine evaluate_checked(evaluate, verify, lo, hi, result, valid, expected_cut_id, &
            expected_enclosure_certificate_id, expected_stationary_certificate_id, enforce_identity, &
            registry_enclosure_certificate_id, registry_stationary_certificate_id)
        procedure(gc_interval_root_callback_i) :: evaluate
        procedure(gc_interval_enclosure_verifier_i) :: verify
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(out) :: result
        logical, intent(out) :: valid
        integer, intent(in), optional :: expected_cut_id
        integer, intent(in), optional :: expected_enclosure_certificate_id
        integer, intent(in), optional :: expected_stationary_certificate_id
        logical, intent(in), optional :: enforce_identity
        integer, intent(in), optional :: registry_enclosure_certificate_id
        integer, intent(in), optional :: registry_stationary_certificate_id
        integer :: verifier_status, verifier_expected_id

        call evaluate(lo, hi, result)
        valid = result%status == 0
        if (.not. valid) return
        if (result%cut_id <= 0 .or. result%enclosure_certificate_id <= 0 .or. &
                result%stationary_certificate_id < 0) then
            valid = .false.
            return
        end if
        if (result%stationary_certificate_id > 0 .and. &
                .not. ieee_is_finite(result%stationary_point)) then
            valid = .false.
            return
        end if
        if (present(registry_enclosure_certificate_id)) then
            if (result%enclosure_certificate_id /= registry_enclosure_certificate_id) then
                valid = .false.
                return
            end if
        end if
        if (present(registry_stationary_certificate_id)) then
            if (result%stationary_certificate_id /= registry_stationary_certificate_id) then
                valid = .false.
                return
            end if
        end if
        if (.not. ieee_is_finite(result%query_lo) .or. .not. ieee_is_finite(result%query_hi)) then
            valid = .false.
            return
        end if
        if (result%query_lo /= lo .or. result%query_hi /= hi) then
            valid = .false.
            return
        end if
        verifier_expected_id = 0
        if (present(registry_enclosure_certificate_id)) then
            verifier_expected_id = registry_enclosure_certificate_id
        else if (present(expected_enclosure_certificate_id)) then
            verifier_expected_id = expected_enclosure_certificate_id
        end if
        if (verifier_expected_id <= 0) then
            valid = .false.
            return
        end if
        call verify(lo, hi, result, verifier_expected_id, verifier_status)
        if (verifier_status /= 0) then
            valid = .false.
            return
        end if
        if (present(enforce_identity)) then
            if (enforce_identity) then
                if (.not. present(expected_cut_id)) then
                    valid = .false.
                    return
                end if
                if (.not. present(expected_enclosure_certificate_id)) then
                    valid = .false.
                    return
                end if
                if (.not. present(expected_stationary_certificate_id)) then
                    valid = .false.
                    return
                end if
                if (result%cut_id /= expected_cut_id .or. &
                        result%enclosure_certificate_id /= expected_enclosure_certificate_id .or. &
                        result%stationary_certificate_id /= expected_stationary_certificate_id) then
                    valid = .false.
                    return
                end if
            end if
        end if
        if (result%boundary_role < GC_INTERVAL_ROOT_INTERIOR .or. &
                result%boundary_role > GC_INTERVAL_ROOT_BOTH_BOUNDARIES) then
            valid = .false.
            return
        end if
        valid = valid_interval(result%f) .and. valid_interval(result%df) .and. &
            valid_interval(result%d2f)
    end subroutine evaluate_checked

    integer function callback_failure_status(value)
        type(gc_interval_callback_result_t), intent(in) :: value
        callback_failure_status = GC_INTERVAL_ROOT_CALLBACK_FAILURE
        if ((.not. ieee_is_finite(value%query_lo) .or. &
                .not. ieee_is_finite(value%query_hi) .or. &
                .not. valid_interval(value%f) .or. .not. valid_interval(value%df) .or. &
                .not. valid_interval(value%d2f))) &
            callback_failure_status = GC_INTERVAL_ROOT_NONFINITE
    end function callback_failure_status

    integer function status_priority(status)
        integer, intent(in) :: status
        select case (status)
        case (GC_INTERVAL_ROOT_SUCCESS)
            status_priority = 0
        case (GC_INTERVAL_ROOT_UNRESOLVED)
            status_priority = 1
        case (GC_INTERVAL_ROOT_CAPACITY)
            status_priority = 2
        case (GC_INTERVAL_ROOT_CALLBACK_FAILURE)
            status_priority = 3
        case (GC_INTERVAL_ROOT_NONFINITE)
            status_priority = 4
        case (GC_INTERVAL_ROOT_INVALID_INPUT)
            status_priority = 5
        case default
            status_priority = 6
        end select
    end function status_priority

    subroutine promote_status(current, candidate)
        integer, intent(inout) :: current
        integer, intent(in) :: candidate
        if (status_priority(candidate) > status_priority(current)) current = candidate
    end subroutine promote_status

    logical function valid_interval(interval)
        type(gc_interval_t), intent(in) :: interval
        valid_interval = ieee_is_finite(interval%lo) .and. ieee_is_finite(interval%hi) .and. &
            interval%lo <= interval%hi
    end function valid_interval

    logical function contains_zero(interval)
        type(gc_interval_t), intent(in) :: interval
        contains_zero = interval%lo <= 0.0_dp .and. interval%hi >= 0.0_dp
    end function contains_zero

    logical function excludes_zero(interval)
        type(gc_interval_t), intent(in) :: interval
        excludes_zero = interval%hi < 0.0_dp .or. interval%lo > 0.0_dp
    end function excludes_zero

    logical function exact_zero(interval)
        type(gc_interval_t), intent(in) :: interval
        exact_zero = interval%lo == 0.0_dp .and. interval%hi == 0.0_dp
    end function exact_zero

    subroutine certify_simple_root(evaluate, verify, lo, hi, value, options, source_cut_id, &
            source_enclosure_certificate_id, source_stationary_certificate_id, root, found, valid)
        procedure(gc_interval_root_callback_i) :: evaluate
        procedure(gc_interval_enclosure_verifier_i) :: verify
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(in) :: value
        type(gc_interval_root_options_t), intent(in) :: options
        integer, intent(in) :: source_cut_id
        integer, intent(in) :: source_enclosure_certificate_id
        integer, intent(in) :: source_stationary_certificate_id
        type(gc_interval_root_box_t), intent(out) :: root
        logical, intent(out) :: found, valid

        type(gc_interval_callback_result_t) :: left_value, right_value, mid_value, current_value
        type(gc_interval_callback_result_t) :: final_left_value, final_right_value
        type(gc_interval_t) :: current, newton, intersection
        real(dp) :: mid
        logical :: left_ok, right_ok, mid_ok, current_ok, bracketed, shrunk
        logical :: bracket_proved, newton_proved, strict_inclusion, progressed
        integer :: iteration

        root = gc_interval_root_box_t()
        root%cut_id = source_cut_id
        root%enclosure_certificate_id = source_enclosure_certificate_id
        root%stationary_certificate_id = source_stationary_certificate_id
        root%derivative_enclosure = value%df
        root%second_derivative_enclosure = value%d2f
        root%derivative_excludes_zero = excludes_zero(value%df)
        root%multiplicity_lower = 1
        root%multiplicity_upper = 1
        root%transversality_certified = .false.
        root%classification_certified = .false.
        root%transversality_kind = GC_INTERVAL_ROOT_TRANSVERSE
        found = .false.
        valid = .true.
        current = gc_interval_t(lo, hi)

        call evaluate_checked(evaluate, verify, lo, lo, left_value, left_ok, source_cut_id, &
            source_enclosure_certificate_id, source_stationary_certificate_id, .true.)
        call evaluate_checked(evaluate, verify, hi, hi, right_value, right_ok, source_cut_id, &
            source_enclosure_certificate_id, source_stationary_certificate_id, .true.)
        valid = left_ok .and. right_ok
        if (.not. valid) return
        bracketed = certified_sign_bracket(left_value%f, right_value%f)
        bracket_proved = bracketed
        newton_proved = .false.

        do iteration = 1, options%max_stationary_iterations
            if (current%hi - current%lo <= options%x_tolerance) then
                found = bracket_proved .or. newton_proved
                if (found) exit
            end if
            call evaluate_checked(evaluate, verify, current%lo, current%hi, current_value, current_ok, &
                source_cut_id, source_enclosure_certificate_id, source_stationary_certificate_id, .true.)
            if (.not. current_ok) then
                valid = .false.
                return
            end if
            if (.not. excludes_zero(current_value%df)) then
                valid = .false.
                return
            end if
            if (.not. bracketed) then
                call evaluate_checked(evaluate, verify, current%lo, current%lo, left_value, left_ok, &
                    source_cut_id, source_enclosure_certificate_id, &
                    source_stationary_certificate_id, .true.)
                if (.not. left_ok) then
                    valid = .false.
                    return
                end if
                call evaluate_checked(evaluate, verify, current%hi, current%hi, right_value, right_ok, &
                    source_cut_id, source_enclosure_certificate_id, &
                    source_stationary_certificate_id, .true.)
                if (.not. right_ok) then
                    valid = .false.
                    return
                end if
                bracketed = certified_sign_bracket(left_value%f, right_value%f)
                if (bracketed) bracket_proved = .true.
            end if
            mid = midpoint(current%lo, current%hi)
            call evaluate_checked(evaluate, verify, mid, mid, mid_value, mid_ok, source_cut_id, &
                source_enclosure_certificate_id, source_stationary_certificate_id, .true.)
            if (.not. mid_ok) then
                valid = .false.
                return
            end if
            progressed = .false.
            if (bracketed) then
                if (exact_zero(mid_value%f)) then
                    current = gc_interval_t(mid, mid)
                    left_value = mid_value
                    right_value = mid_value
                    bracket_proved = .true.
                    progressed = .true.
                else if (.not. contains_zero(mid_value%f)) then
                    if (same_sign_as_left(mid_value%f, left_value%f)) then
                        current%lo = mid
                        left_value = mid_value
                        progressed = .true.
                    else if (same_sign_as_right(mid_value%f, right_value%f)) then
                        current%hi = mid
                        right_value = mid_value
                        progressed = .true.
                    end if
                end if
            end if
            if (.not. progressed) then
                call interval_newton(mid, mid_value%f, current_value%df, newton, shrunk)
                strict_inclusion = shrunk .and. newton%lo > current%lo .and. &
                    newton%hi < current%hi
                if (strict_inclusion) newton_proved = .true.
                if (shrunk) then
                    intersection%lo = max(current%lo, newton%lo)
                    intersection%hi = min(current%hi, newton%hi)
                    if (valid_interval(intersection)) then
                        if (intersection%hi - intersection%lo < current%hi - current%lo) then
                            current = intersection
                            bracketed = .false.
                            progressed = .true.
                        end if
                    end if
                end if
            end if
            if (.not. progressed) exit
        end do

        if (.not. found) then
            if (current%hi - current%lo <= options%x_tolerance) then
                found = bracket_proved .or. newton_proved
            end if
        end if
        if (found) then
            call evaluate_checked(evaluate, verify, current%lo, current%hi, current_value, current_ok, &
                source_cut_id, source_enclosure_certificate_id, source_stationary_certificate_id, .true.)
            if (.not. current_ok) then
                valid = .false.
                found = .false.
                return
            end if
            if (.not. excludes_zero(current_value%df)) then
                valid = .false.
                found = .false.
                return
            end if
            call evaluate_checked(evaluate, verify, current%lo, current%lo, final_left_value, left_ok, &
                source_cut_id, source_enclosure_certificate_id, &
                source_stationary_certificate_id, .true.)
            if (.not. left_ok) then
                valid = .false.
                found = .false.
                return
            end if
            call evaluate_checked(evaluate, verify, current%hi, current%hi, final_right_value, right_ok, &
                source_cut_id, source_enclosure_certificate_id, &
                source_stationary_certificate_id, .true.)
            if (.not. right_ok) then
                valid = .false.
                found = .false.
                return
            end if
            root%lo = current%lo
            root%hi = current%hi
            root%kind = GC_INTERVAL_ROOT_SIMPLE
            root%bracket_certified = bracket_proved
            root%interval_newton_certified = newton_proved
            root%left_endpoint_root = .false.
            if (root%lo == lo) root%left_endpoint_root = exact_zero(final_left_value%f)
            root%right_endpoint_root = .false.
            if (root%hi == hi) root%right_endpoint_root = exact_zero(final_right_value%f)
            root%transversality_certified = .true.
            root%classification_certified = .true.
            root%boundary_role = merge(value%boundary_role, &
                boundary_role(lo, hi, root%lo, root%hi), value%boundary_role /= 0)
            root%derivative_enclosure = current_value%df
        end if
    end subroutine certify_simple_root

    subroutine certify_tangent_root(evaluate, verify, verify_stationary, lo, hi, value, options, &
            source_cut_id, &
            source_enclosure_certificate_id, source_stationary_certificate_id, root, found, valid)
        procedure(gc_interval_root_callback_i) :: evaluate
        procedure(gc_interval_enclosure_verifier_i) :: verify
        procedure(gc_interval_stationary_verifier_i) :: verify_stationary
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(in) :: value
        type(gc_interval_root_options_t), intent(in) :: options
        integer, intent(in) :: source_cut_id
        integer, intent(in) :: source_enclosure_certificate_id
        integer, intent(in) :: source_stationary_certificate_id
        type(gc_interval_root_box_t), intent(out) :: root
        logical, intent(out) :: found, valid

        type(gc_interval_callback_result_t) :: left_value, right_value
        type(gc_interval_callback_result_t) :: mid_value, stationary_value, current_value, certificate_value
        type(gc_interval_t) :: current, newton
        real(dp) :: mid
        logical :: left_ok, right_ok, mid_ok, current_ok, stationary_ok, certificate_ok
        logical :: bracketed, shrunk, strict_inclusion, stationary_closed
        integer :: iteration, stationary_status

        root = gc_interval_root_box_t()
        root%cut_id = source_cut_id
        root%enclosure_certificate_id = source_enclosure_certificate_id
        root%stationary_certificate_id = source_stationary_certificate_id
        root%derivative_enclosure = value%df
        root%second_derivative_enclosure = value%d2f
        root%multiplicity_lower = 2
        root%multiplicity_upper = 2
        root%transversality_certified = .false.
        root%classification_certified = .false.
        root%transversality_kind = GC_INTERVAL_ROOT_EXTREMAL
        found = .false.
        valid = .true.
        current = gc_interval_t(lo, hi)
        call evaluate_checked(evaluate, verify, lo, lo, left_value, left_ok, source_cut_id, &
            source_enclosure_certificate_id, source_stationary_certificate_id, .true.)
        call evaluate_checked(evaluate, verify, hi, hi, right_value, right_ok, source_cut_id, &
            source_enclosure_certificate_id, source_stationary_certificate_id, .true.)
        valid = left_ok .and. right_ok
        if (.not. valid) return
        bracketed = certified_sign_bracket(left_value%df, right_value%df)
        stationary_closed = .false.

        do iteration = 1, options%max_stationary_iterations
            call evaluate_checked(evaluate, verify, current%lo, current%hi, current_value, current_ok, &
                source_cut_id, source_enclosure_certificate_id, source_stationary_certificate_id, .true.)
            if (.not. current_ok) then
                valid = .false.
                return
            end if
            if (.not. excludes_zero(current_value%d2f)) then
                valid = .false.
                return
            end if
            mid = midpoint(current%lo, current%hi)
            call evaluate_checked(evaluate, verify, mid, mid, mid_value, mid_ok, source_cut_id, &
                source_enclosure_certificate_id, source_stationary_certificate_id, .true.)
            if (.not. mid_ok) then
                valid = .false.
                return
            end if
            if (bracketed .and. .not. contains_zero(mid_value%df)) then
                if (same_sign_as_left(mid_value%df, left_value%df)) then
                    current%lo = mid
                    left_value = mid_value
                else
                    current%hi = mid
                    right_value = mid_value
                end if
            else
                call interval_newton(mid, mid_value%df, current_value%d2f, newton, shrunk)
                strict_inclusion = shrunk .and. newton%lo > current%lo .and. &
                    newton%hi < current%hi
                if (strict_inclusion) then
                    current = newton
                    stationary_closed = .true.
                end if
            end if
            if (.not. bracketed) bracketed = certified_sign_bracket(&
                left_value%df, right_value%df)
            if (current%hi - current%lo <= options%x_tolerance .or. stationary_closed) exit
        end do

        !! A small interval is not a stationary-point existence proof.  The
        !! derivative bracket or strict interval-Newton inclusion must close
        !! the stationary isolation before the exact certificate is examined.
        if (.not. bracketed .and. .not. stationary_closed) return
        !! Evaluate the whole isolated stationary box.  The exact generated
        !! stationary certificate below, not a sampled point value, closes
        !! the tangent test.
        call evaluate_checked(evaluate, verify, current%lo, current%hi, stationary_value, stationary_ok, &
            source_cut_id, source_enclosure_certificate_id, source_stationary_certificate_id, .true.)
        if (.not. stationary_ok) then
            valid = .false.
            return
        end if
        if (.not. contains_zero(stationary_value%f)) return
        if (stationary_value%stationary_certificate_id <= 0) return
        if (stationary_value%stationary_certificate_id /= source_stationary_certificate_id) return
        if (stationary_value%stationary_point < current%lo .or. &
                stationary_value%stationary_point > current%hi) return
        call verify_stationary(current%lo, current%hi, stationary_value%stationary_point, &
            certificate_value, source_enclosure_certificate_id, source_stationary_certificate_id, &
            stationary_status)
        certificate_ok = stationary_status == 0
        if (.not. certificate_ok) then
            valid = .false.
            return
        end if
        if (.not. valid_interval(certificate_value%f) .or. &
                .not. valid_interval(certificate_value%df)) then
            valid = .false.
            return
        end if
        if (certificate_value%query_lo /= stationary_value%stationary_point .or. &
                certificate_value%query_hi /= stationary_value%stationary_point) then
            valid = .false.
            return
        end if
        if (certificate_value%stationary_point /= stationary_value%stationary_point .or. &
                certificate_value%stationary_certificate_id /= &
                    stationary_value%stationary_certificate_id) return
        if (.not. exact_zero(certificate_value%f) .or. &
                .not. exact_zero(certificate_value%df)) return
        if (.not. excludes_zero(stationary_value%d2f)) return

        root%lo = certificate_value%stationary_point
        root%hi = certificate_value%stationary_point
        root%kind = GC_INTERVAL_ROOT_TANGENT
        root%stationary_certified = .true.
        root%bracket_certified = bracketed
        root%derivative_enclosure = stationary_value%df
        root%second_derivative_enclosure = stationary_value%d2f
        root%stationary_certificate_id = stationary_value%stationary_certificate_id
        root%interval_newton_certified = stationary_closed
        root%left_endpoint_root = root%lo == lo
        root%right_endpoint_root = root%hi == hi
        root%classification_certified = .true.
        root%boundary_role = merge(stationary_value%boundary_role, &
            boundary_role(lo, hi, root%lo, root%hi), stationary_value%boundary_role /= 0)
        found = .true.
    end subroutine certify_tangent_root

    subroutine interval_newton(mid, function_interval, derivative_interval, newton, useful)
        real(dp), intent(in) :: mid
        type(gc_interval_t), intent(in) :: function_interval, derivative_interval
        type(gc_interval_t), intent(out) :: newton
        logical, intent(out) :: useful

        type(gc_interval_t) :: quotient
        real(dp) :: values(4)

        newton = gc_interval_t()
        useful = .false.
        if (.not. excludes_zero(derivative_interval)) return
        values = [function_interval%lo/derivative_interval%lo, &
            function_interval%lo/derivative_interval%hi, &
            function_interval%hi/derivative_interval%lo, &
            function_interval%hi/derivative_interval%hi]
        if (.not. all(ieee_is_finite(values))) return
        quotient%lo = outward_lower(minval(values))
        quotient%hi = outward_upper(maxval(values))
        newton%lo = outward_lower(mid - quotient%hi)
        newton%hi = outward_upper(mid - quotient%lo)
        useful = valid_interval(newton)
    end subroutine interval_newton

    logical function certified_sign_bracket(left, right)
        type(gc_interval_t), intent(in) :: left, right
        certified_sign_bracket = (left%hi < 0.0_dp .and. right%lo > 0.0_dp) .or. &
            (left%lo > 0.0_dp .and. right%hi < 0.0_dp) .or. &
            (left%lo == 0.0_dp .and. left%hi == 0.0_dp) .or. &
            (right%lo == 0.0_dp .and. right%hi == 0.0_dp)
    end function certified_sign_bracket

    logical function same_sign_as_left(value, left)
        type(gc_interval_t), intent(in) :: value, left
        same_sign_as_left = (left%hi < 0.0_dp .and. value%hi < 0.0_dp) .or. &
            (left%lo > 0.0_dp .and. value%lo > 0.0_dp)
    end function same_sign_as_left

    logical function same_sign_as_right(value, right)
        type(gc_interval_t), intent(in) :: value, right
        same_sign_as_right = (right%hi < 0.0_dp .and. value%hi < 0.0_dp) .or. &
            (right%lo > 0.0_dp .and. value%lo > 0.0_dp)
    end function same_sign_as_right

    real(dp) function midpoint(lo, hi)
        real(dp), intent(in) :: lo, hi
        midpoint = lo + 0.5_dp*(hi - lo)
    end function midpoint

    real(dp) function outward_lower(value)
        real(dp), intent(in) :: value
        outward_lower = ieee_next_after(value, -huge(value))
    end function outward_lower

    real(dp) function outward_upper(value)
        real(dp), intent(in) :: value
        outward_upper = ieee_next_after(value, huge(value))
    end function outward_upper

    subroutine split_box(queue, head, tail)
        type(work_box_t), intent(inout) :: queue(:)
        integer, intent(inout) :: head
        integer, intent(inout) :: tail
        real(dp) :: mid

        mid = midpoint(queue(head)%lo, queue(head)%hi)
        queue(tail + 1) = queue(head)
        queue(tail + 2) = queue(head)
        queue(tail + 1)%lo = queue(head)%lo
        queue(tail + 1)%hi = mid
        queue(tail + 2)%lo = mid
        queue(tail + 2)%hi = queue(head)%hi
        queue(tail + 1)%depth = queue(head)%depth + 1
        queue(tail + 2)%depth = queue(head)%depth + 1
        tail = tail + 2
        head = head + 1
    end subroutine split_box

    subroutine append_candidate(candidates, n, root, options, ok)
        type(gc_interval_root_box_t), intent(inout) :: candidates(:)
        integer, intent(inout) :: n
        type(gc_interval_root_box_t), intent(in) :: root
        type(gc_interval_root_options_t), intent(in) :: options
        logical, intent(out) :: ok

        ok = .false.
        if (root%hi < root%lo .or. root%kind == 0) return
        if (n >= min(size(candidates), options%max_roots)) return
        n = n + 1
        candidates(n) = root
        ok = .true.
    end subroutine append_candidate

    subroutine order_and_validate_candidates(candidates, n, result, options)
        type(gc_interval_root_box_t), intent(in) :: candidates(:)
        integer, intent(in) :: n
        type(gc_interval_root_result_t), intent(inout) :: result
        type(gc_interval_root_options_t), intent(in) :: options

        type(gc_interval_root_box_t), allocatable :: ordered(:)
        type(gc_interval_root_box_t) :: temporary
        integer :: i, j, count

        if (allocated(result%roots)) deallocate(result%roots)
        if (n == 0) then
            allocate(result%roots(0))
            return
        end if
        allocate(ordered(n))
        ordered = candidates(1:n)
        do i = 2, n
            temporary = ordered(i)
            j = i - 1
            do while (j >= 1)
                if (.not. root_precedes(temporary, ordered(j))) exit
                ordered(j + 1) = ordered(j)
                j = j - 1
            end do
            ordered(j + 1) = temporary
        end do

        count = 0
        do i = 1, n
            if (count > 0) then
                if (shared_endpoint_duplicate(ordered(count), ordered(i))) then
                    call merge_shared_endpoint_root(ordered(count), ordered(i))
                    cycle
                end if
                if (ordered(i)%lo < ordered(count)%hi) then
                    call promote_status(result%status, GC_INTERVAL_ROOT_UNRESOLVED)
                    result%unresolved_boxes = result%unresolved_boxes + 1
                    cycle
                end if
            end if
            if (count < size(ordered)) then
                count = count + 1
                ordered(count) = ordered(i)
            else
                call promote_status(result%status, GC_INTERVAL_ROOT_CAPACITY)
                result%unresolved_boxes = result%unresolved_boxes + 1
            end if
        end do
        if (count > options%max_roots) then
            call promote_status(result%status, GC_INTERVAL_ROOT_CAPACITY)
            allocate(result%roots(0))
            return
        end if
        allocate(result%roots(count))
        if (count > 0) result%roots = ordered(1:count)
        result%nroots = count
    end subroutine order_and_validate_candidates

    logical function shared_endpoint_duplicate(left, right)
        type(gc_interval_root_box_t), intent(in) :: left, right
        shared_endpoint_duplicate = left%hi == right%lo .and. &
            left%right_endpoint_root .and. right%left_endpoint_root .and. &
            left%kind == right%kind .and. left%cut_id == right%cut_id .and. &
            left%enclosure_certificate_id == right%enclosure_certificate_id .and. &
            left%stationary_certificate_id == right%stationary_certificate_id
    end function shared_endpoint_duplicate

    subroutine merge_shared_endpoint_root(left, right)
        type(gc_interval_root_box_t), intent(inout) :: left
        type(gc_interval_root_box_t), intent(in) :: right
        real(dp) :: point

        point = left%hi
        left%lo = point
        left%hi = point
        left%boundary_role = ior(left%boundary_role, right%boundary_role)
        left%left_endpoint_root = .true.
        left%right_endpoint_root = .true.
        left%bracket_certified = left%bracket_certified .or. right%bracket_certified
        left%interval_newton_certified = left%interval_newton_certified .or. &
            right%interval_newton_certified
        left%stationary_certified = left%stationary_certified .or. right%stationary_certified
        left%classification_certified = left%classification_certified .and. &
            right%classification_certified
        left%transversality_certified = left%transversality_certified .and. &
            right%transversality_certified
    end subroutine merge_shared_endpoint_root

    logical function root_precedes(left, right)
        type(gc_interval_root_box_t), intent(in) :: left, right
        root_precedes = left%lo < right%lo .or. &
            (left%lo == right%lo .and. left%hi < right%hi) .or. &
            (left%lo == right%lo .and. left%hi == right%hi .and. &
                left%cut_id < right%cut_id)
    end function root_precedes

    integer function boundary_role(source_lo, source_hi, root_lo, root_hi)
        real(dp), intent(in) :: source_lo, source_hi, root_lo, root_hi
        boundary_role = GC_INTERVAL_ROOT_INTERIOR
        if (root_lo == source_lo) boundary_role = GC_INTERVAL_ROOT_LEFT_BOUNDARY
        if (root_hi == source_hi) boundary_role = ior(boundary_role, &
            GC_INTERVAL_ROOT_RIGHT_BOUNDARY)
    end function boundary_role

end module neort_gc_certified_interval_roots
