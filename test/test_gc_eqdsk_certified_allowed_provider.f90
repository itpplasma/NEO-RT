program test_gc_eqdsk_certified_allowed_provider
    !! Contract-level coverage for an Eq. 13-style certified allowed provider.
    !!
    !! Keep the fixture independent of a graph or field representation and
    !! exercise both the provider's fixed-invariant certificate seam and the
    !! interval/root-isolation contract it must satisfy.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_next_after
    use neort_gc_eqdsk_certified_allowed_provider, only: &
        GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_FAILURE, &
        verify_gc_eqdsk_fixed_invariant_stationary_certificate
    use neort_gc_certified_interval_roots, only: &
        GC_INTERVAL_ROOT_CALLBACK_FAILURE, GC_INTERVAL_ROOT_SIMPLE, &
        GC_INTERVAL_ROOT_SUCCESS, GC_INTERVAL_ROOT_UNRESOLVED, &
        gc_interval_callback_result_t, gc_interval_enclosure_verifier_i, &
        gc_interval_root_callback_i, gc_interval_root_options_t, &
        gc_interval_root_result_t, gc_interval_stationary_verifier_i, &
        isolate_gc_interval_roots, gc_interval_t
    implicit none

    integer, parameter :: MODE_SIMPLE = 1
    integer, parameter :: MODE_TANGENT = 2
    integer, parameter :: MODE_UNRESOLVED = 3
    integer, parameter :: MODE_FORGED_QUERY = 4
    integer, parameter :: MODE_FORGED_CERTIFICATE = 5
    integer, parameter :: MODE_FORGED_IDENTITY = 6
    integer, parameter :: CUT_ID = 2718
    integer, parameter :: ENCLOSURE_ID = 9142
    integer, parameter :: STATIONARY_ID = 7731
    integer, parameter :: FORGED_CUT_ID = 9991
    integer, parameter :: FORGED_ENCLOSURE_ID = 9992
    real(dp), parameter :: SIMPLE_ROOT = 0.25_dp
    real(dp), parameter :: TANGENT_ROOT = 0.5_dp

    type(gc_interval_root_options_t) :: options
    type(gc_interval_root_result_t) :: result
    integer :: mode
    integer :: verifier_calls, point_verifier_calls
    integer :: stationary_verifier_calls

    options%initial_partition = 1
    options%max_depth = 16
    options%max_boxes = 512
    options%max_roots = 16
    options%max_stationary_iterations = 8
    options%expected_enclosure_certificate_id = ENCLOSURE_ID
    options%expected_stationary_certificate_id = STATIONARY_ID
    options%x_tolerance = 1.0e-6_dp

    mode = MODE_SIMPLE
    verifier_calls = 0
    point_verifier_calls = 0
    stationary_verifier_calls = 0
    call isolate_gc_interval_roots(manufactured_callback, independent_verifier, &
        rejecting_stationary_verifier, 0.0_dp, 1.0_dp, options, result)
    call require(result%status == GC_INTERVAL_ROOT_SUCCESS .and. &
        result%coverage_certified, 'certified simple-root coverage did not close')
    call require(result%nroots == 1, 'certified simple root was not isolated')
    call require(result%roots(1)%kind == GC_INTERVAL_ROOT_SIMPLE, &
        'simple root was not classified as simple')
    call require(result%roots(1)%cut_id == CUT_ID .and. &
        result%roots(1)%enclosure_certificate_id == ENCLOSURE_ID .and. &
        result%roots(1)%stationary_certificate_id == STATIONARY_ID, &
        'simple-root certificate identity was not retained')
    call require(result%roots(1)%lo <= SIMPLE_ROOT .and. &
        result%roots(1)%hi >= SIMPLE_ROOT, &
        'simple-root enclosure does not contain the manufactured root')
    call require(verifier_calls > 0 .and. point_verifier_calls > 0, &
        'exact interval and point queries were not independently verified')

    mode = MODE_TANGENT
    stationary_verifier_calls = 0
    call isolate_gc_interval_roots(manufactured_callback, independent_verifier, &
        rejecting_stationary_verifier, 0.25_dp, 0.75_dp, options, result)
    call require(result%status == GC_INTERVAL_ROOT_UNRESOLVED .and. &
        .not. result%coverage_certified .and. result%nroots == 0, &
        'tangent candidate was not failed closed')
    call require(stationary_verifier_calls > 0, &
        'tangent candidate did not reach the exact stationary verifier')

    mode = MODE_UNRESOLVED
    options%max_depth = 0
    call isolate_gc_interval_roots(manufactured_callback, independent_verifier, &
        rejecting_stationary_verifier, 0.0_dp, 1.0_dp, options, result)
    call require(result%status == GC_INTERVAL_ROOT_UNRESOLVED .and. &
        result%unresolved_boxes > 0 .and. .not. result%coverage_certified .and. &
        result%nroots == 0, 'unresolved interval coverage was consumed')

    options%max_depth = 16
    mode = MODE_FORGED_QUERY
    call isolate_gc_interval_roots(manufactured_callback, independent_verifier, &
        rejecting_stationary_verifier, 0.0_dp, 1.0_dp, options, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified .and. result%nroots == 0, &
        'forged callback query identity was accepted')

    mode = MODE_FORGED_CERTIFICATE
    call isolate_gc_interval_roots(manufactured_callback, independent_verifier, &
        rejecting_stationary_verifier, 0.0_dp, 1.0_dp, options, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified .and. result%nroots == 0, &
        'forged enclosure certificate identity was accepted')

    mode = MODE_FORGED_IDENTITY
    call isolate_gc_interval_roots(manufactured_callback, independent_verifier, &
        rejecting_stationary_verifier, 0.0_dp, 1.0_dp, options, result)
    call require(result%status == GC_INTERVAL_ROOT_CALLBACK_FAILURE .and. &
        .not. result%coverage_certified .and. result%nroots == 0, &
        'forged cut identity was accepted')

    call check_fixed_invariant_stationary_certificate_seam()

    write (*, '(A)') 'test_gc_eqdsk_certified_allowed_provider OK'

contains

    subroutine check_fixed_invariant_stationary_certificate_seam()
        type(gc_interval_callback_result_t) :: candidate, point_value
        integer :: status

        candidate = gc_interval_callback_result_t()
        candidate%query_lo = 0.25_dp
        candidate%query_hi = 0.75_dp
        candidate%cut_id = CUT_ID
        candidate%status = 0
        candidate%enclosure_certificate_id = ENCLOSURE_ID
        candidate%stationary_certificate_id = STATIONARY_ID
        candidate%stationary_point = TANGENT_ROOT
        candidate%f = gc_interval_t(-1.0_dp, 1.0_dp)
        candidate%df = gc_interval_t(-1.0_dp, 1.0_dp)
        candidate%d2f = gc_interval_t(1.0_dp, 3.0_dp)

        point_value = candidate
        point_value%query_lo = TANGENT_ROOT
        point_value%query_hi = TANGENT_ROOT
        point_value%f = gc_interval_t(0.0_dp, 0.0_dp)
        point_value%df = gc_interval_t(0.0_dp, 0.0_dp)
        point_value%d2f = gc_interval_t(2.0_dp, 2.0_dp)
        call verify_gc_eqdsk_fixed_invariant_stationary_certificate( &
            candidate%query_lo, candidate%query_hi, candidate, point_value, &
            ENCLOSURE_ID, STATIONARY_ID, status)
        call require(status == GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_FAILURE, &
            'unregistered stationary identity was accepted')

        point_value%f = gc_interval_t(-epsilon(1.0_dp), epsilon(1.0_dp))
        call verify_gc_eqdsk_fixed_invariant_stationary_certificate( &
            candidate%query_lo, candidate%query_hi, candidate, point_value, &
            ENCLOSURE_ID, STATIONARY_ID, status)
        call require(status == GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_FAILURE, &
            'non-exact stationary value was promoted to a certificate')

        point_value%f = gc_interval_t(0.0_dp, 0.0_dp)
        point_value%stationary_certificate_id = 0
        call verify_gc_eqdsk_fixed_invariant_stationary_certificate( &
            candidate%query_lo, candidate%query_hi, candidate, point_value, &
            ENCLOSURE_ID, STATIONARY_ID, status)
        call require(status == GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_FAILURE, &
            'missing stationary identity was accepted')
    end subroutine check_fixed_invariant_stationary_certificate_seam

    subroutine manufactured_callback(lo, hi, value)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(out) :: value
        real(dp) :: a, b

        value = gc_interval_callback_result_t()
        value%query_lo = lo
        value%query_hi = hi
        value%status = 0
        value%certified = .true.
        value%cut_id = CUT_ID
        value%enclosure_certificate_id = ENCLOSURE_ID
        value%stationary_certificate_id = STATIONARY_ID
        value%stationary_point = TANGENT_ROOT

        select case (mode)
        case (MODE_SIMPLE, MODE_FORGED_QUERY, MODE_FORGED_CERTIFICATE, &
                MODE_FORGED_IDENTITY)
            value%f = callback_linear_interval(lo, hi)
            if (lo == hi) then
                value%df = point_interval(1.0_dp)
            else
                value%df = gc_interval_t(1.0_dp, 1.0_dp)
            end if
            value%d2f = gc_interval_t(0.0_dp, 0.0_dp)
        case (MODE_TANGENT)
            value%stationary_point = TANGENT_ROOT
            value%f = callback_square_interval(lo, hi, TANGENT_ROOT)
            if (lo == hi) then
                value%df = point_interval(2.0_dp*(lo - TANGENT_ROOT))
            else
                value%df = gc_interval_t( &
                    lower_outward(2.0_dp*(lo - TANGENT_ROOT)), &
                    upper_outward(2.0_dp*(hi - TANGENT_ROOT)))
            end if
            value%d2f = gc_interval_t(2.0_dp, 2.0_dp)
        case (MODE_UNRESOLVED)
            value%f = gc_interval_t(0.0_dp, 1.0_dp)
            value%df = gc_interval_t(-1.0_dp, 1.0_dp)
            value%d2f = gc_interval_t(-1.0_dp, 1.0_dp)
        case default
            value%status = 1
        end select

        if (mode == MODE_FORGED_QUERY .and. lo == hi) then
            value%query_hi = hi + 1.0e-3_dp
        else if (mode == MODE_FORGED_CERTIFICATE .and. lo == hi) then
            value%enclosure_certificate_id = FORGED_ENCLOSURE_ID
        else if (mode == MODE_FORGED_IDENTITY .and. lo == hi) then
            value%cut_id = FORGED_CUT_ID
        end if

        ! Keep the point enclosure exact; the interval cases above are
        ! outward enclosures of the same manufactured functions.
        if (mode == MODE_SIMPLE .or. mode == MODE_FORGED_QUERY .or. &
                mode == MODE_FORGED_CERTIFICATE .or. mode == MODE_FORGED_IDENTITY) then
            if (lo == hi) value%f = point_interval(lo - SIMPLE_ROOT)
        else if (mode == MODE_TANGENT) then
            if (lo == hi) then
                a = lo - TANGENT_ROOT
                value%f = point_interval(a*a)
            end if
        end if
    end subroutine manufactured_callback

    subroutine independent_verifier(lo, hi, value, expected_id, status)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_callback_result_t), intent(in) :: value
        integer, intent(in) :: expected_id
        integer, intent(out) :: status
        type(gc_interval_t) :: expected_f, expected_df, expected_d2f

        verifier_calls = verifier_calls + 1
        if (lo == hi) point_verifier_calls = point_verifier_calls + 1
        status = 1
        if (expected_id /= ENCLOSURE_ID) return
        if (value%query_lo /= lo .or. value%query_hi /= hi) return
        if (value%cut_id /= CUT_ID .or. &
                value%enclosure_certificate_id /= ENCLOSURE_ID .or. &
                value%stationary_certificate_id /= STATIONARY_ID) return

        select case (mode)
        case (MODE_SIMPLE, MODE_FORGED_QUERY, MODE_FORGED_CERTIFICATE, &
                MODE_FORGED_IDENTITY)
            if (lo == hi) then
                expected_f = point_interval(lo - SIMPLE_ROOT)
            else
                expected_f = verifier_linear_interval(lo, hi)
            end if
            expected_df = point_interval(1.0_dp)
            expected_d2f = point_interval(0.0_dp)
        case (MODE_TANGENT)
            if (lo == hi) then
                expected_f = point_interval((lo - TANGENT_ROOT)**2)
            else
                expected_f = verifier_square_interval(lo, hi, TANGENT_ROOT)
            end if
            if (lo == hi) then
                expected_df = point_interval(2.0_dp*(lo - TANGENT_ROOT))
            else
                expected_df = gc_interval_t( &
                    lower_outward(2.0_dp*(lo - TANGENT_ROOT)), &
                    upper_outward(2.0_dp*(hi - TANGENT_ROOT)))
            end if
            expected_d2f = point_interval(2.0_dp)
        case (MODE_UNRESOLVED)
            expected_f = gc_interval_t(0.0_dp, 1.0_dp)
            expected_df = gc_interval_t(-1.0_dp, 1.0_dp)
            expected_d2f = gc_interval_t(-1.0_dp, 1.0_dp)
        case default
            return
        end select

        if (.not. encloses(value%f, expected_f) .or. &
                .not. encloses(value%df, expected_df) .or. &
                .not. encloses(value%d2f, expected_d2f)) return
        status = 0
    end subroutine independent_verifier

    subroutine rejecting_stationary_verifier(lo, hi, point, value, &
            expected_enclosure_id, expected_stationary_id, status)
        real(dp), intent(in) :: lo, hi, point
        type(gc_interval_callback_result_t), intent(out) :: value
        integer, intent(in) :: expected_enclosure_id, expected_stationary_id
        integer, intent(out) :: status

        stationary_verifier_calls = stationary_verifier_calls + 1
        value = gc_interval_callback_result_t()
        value%query_lo = point
        value%query_hi = point
        value%cut_id = CUT_ID
        value%enclosure_certificate_id = expected_enclosure_id
        value%stationary_certificate_id = expected_stationary_id
        value%stationary_point = point
        value%status = 0
        value%f = gc_interval_t(0.0_dp, 1.0_dp)
        value%df = point_interval(0.0_dp)
        value%d2f = point_interval(2.0_dp)
        status = 1
        if (expected_enclosure_id /= ENCLOSURE_ID .or. &
                expected_stationary_id /= STATIONARY_ID) return
        if (point < lo .or. point > hi) return
        ! Deliberately non-exact f(point): a stationary candidate is not a
        ! tangent-root certificate unless the exact generated value is zero.
        status = 0
    end subroutine rejecting_stationary_verifier

    function callback_linear_interval(lo, hi) result(value)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_t) :: value

        value = gc_interval_t(lower_outward(lo - SIMPLE_ROOT), &
            upper_outward(hi - SIMPLE_ROOT))
    end function callback_linear_interval

    function verifier_linear_interval(lo, hi) result(value)
        real(dp), intent(in) :: lo, hi
        type(gc_interval_t) :: value

        value%lo = ieee_next_after(lo - SIMPLE_ROOT, -huge(0.0_dp))
        value%hi = ieee_next_after(hi - SIMPLE_ROOT, huge(0.0_dp))
    end function verifier_linear_interval

    function callback_square_interval(lo, hi, root) result(value)
        real(dp), intent(in) :: lo, hi, root
        type(gc_interval_t) :: value
        real(dp) :: left, right, smallest, largest

        left = lo - root
        right = hi - root
        smallest = min(abs(left), abs(right))
        largest = max(abs(left), abs(right))
        if (left <= 0.0_dp .and. right >= 0.0_dp) smallest = 0.0_dp
        value = gc_interval_t(lower_outward(smallest*smallest), &
            upper_outward(largest*largest))
    end function callback_square_interval

    function verifier_square_interval(lo, hi, root) result(value)
        real(dp), intent(in) :: lo, hi, root
        type(gc_interval_t) :: value
        real(dp) :: left, right, smallest, largest

        left = lo - root
        right = hi - root
        smallest = min(abs(left), abs(right))
        largest = max(abs(left), abs(right))
        if (left <= 0.0_dp .and. right >= 0.0_dp) smallest = 0.0_dp
        value%lo = ieee_next_after(smallest*smallest, -huge(0.0_dp))
        value%hi = ieee_next_after(largest*largest, huge(0.0_dp))
    end function verifier_square_interval

    pure function point_interval(x) result(value)
        real(dp), intent(in) :: x
        type(gc_interval_t) :: value

        value = gc_interval_t(x, x)
    end function point_interval

    pure real(dp) function lower_outward(x)
        real(dp), intent(in) :: x

        lower_outward = ieee_next_after(x, -huge(x))
    end function lower_outward

    pure real(dp) function upper_outward(x)
        real(dp), intent(in) :: x

        upper_outward = ieee_next_after(x, huge(x))
    end function upper_outward

    logical function encloses(actual, expected)
        type(gc_interval_t), intent(in) :: actual, expected

        encloses = actual%lo <= expected%lo .and. actual%hi >= expected%hi
    end function encloses

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*, '(A)') trim(message)
            error stop 1
        end if
    end subroutine require

end program test_gc_eqdsk_certified_allowed_provider
