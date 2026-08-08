module neort_gc_operational_scalar_roots
    !! Data-oriented scalar root search for the operational POTATO class flow.
    !!
    !! Function values and derivatives belong to generated physics callbacks.
    !! This module only samples, brackets, bisects, compares two resolutions,
    !! and fails closed on invalid domains or degenerate roots.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: GC_OPERATIONAL_ROOT_SUCCESS = 0
    integer, parameter, public :: GC_OPERATIONAL_ROOT_INVALID_INPUT = 1
    integer, parameter, public :: GC_OPERATIONAL_ROOT_CALLBACK_FAILURE = 2
    integer, parameter, public :: GC_OPERATIONAL_ROOT_NONFINITE = 3
    integer, parameter, public :: GC_OPERATIONAL_ROOT_NONCONVERGED = 4
    integer, parameter, public :: GC_OPERATIONAL_ROOT_CAPACITY = 5
    integer, parameter, public :: GC_OPERATIONAL_ROOT_DEGENERATE = 6
    integer, parameter, public :: GC_OPERATIONAL_ROOT_REFINEMENT_MISMATCH = 7
    integer, parameter, public :: GC_OPERATIONAL_ROOT_UNRESOLVED_SEPARATION = 8

    type, public :: gc_operational_root_options_t
        integer :: initial_intervals = 100
        integer :: maximum_iterations = 100
        integer :: maximum_roots = 4096
        real(dp) :: relative_tolerance = 1.0e-12_dp
    end type gc_operational_root_options_t

    type, public :: gc_operational_root_t
        real(dp) :: x = 0.0_dp
        real(dp) :: bracket_lo = 0.0_dp
        real(dp) :: bracket_hi = 0.0_dp
        real(dp) :: residual = 0.0_dp
        real(dp) :: derivative = 0.0_dp
        logical :: bracketed = .false.
        logical :: endpoint_root = .false.
        logical :: simple = .false.
    end type gc_operational_root_t

    type, public :: gc_operational_root_result_t
        integer :: status = GC_OPERATIONAL_ROOT_INVALID_INPUT
        integer :: nroots = 0
        integer :: coarse_nroots = 0
        integer :: evaluations = 0
        logical :: coarse_fine_agree = .false.
        logical :: complete = .false.
        type(gc_operational_root_t), allocatable :: roots(:)
    end type gc_operational_root_result_t

    abstract interface
        subroutine gc_operational_scalar_i(x, value, derivative, status)
            import :: dp
            real(dp), intent(in) :: x
            real(dp), intent(out) :: value, derivative
            integer, intent(out) :: status
        end subroutine gc_operational_scalar_i
    end interface

    public :: find_gc_operational_scalar_roots
    public :: gc_operational_scalar_i

contains

    subroutine find_gc_operational_scalar_roots(evaluate, domain_lo, domain_hi, &
            options, result)
        procedure(gc_operational_scalar_i) :: evaluate
        real(dp), intent(in) :: domain_lo, domain_hi
        type(gc_operational_root_options_t), intent(in) :: options
        type(gc_operational_root_result_t), intent(out) :: result

        type(gc_operational_root_result_t) :: coarse, fine
        real(dp) :: coordinate_tolerance
        integer :: i

        result = gc_operational_root_result_t()
        allocate(result%roots(0))
        if (.not. valid_options(domain_lo, domain_hi, options)) return

        call scan_roots(evaluate, domain_lo, domain_hi, &
            options%initial_intervals, options, coarse)
        if (coarse%status /= GC_OPERATIONAL_ROOT_SUCCESS) then
            result = coarse
            return
        end if
        call scan_roots(evaluate, domain_lo, domain_hi, &
            2*options%initial_intervals, options, fine)
        if (fine%status /= GC_OPERATIONAL_ROOT_SUCCESS) then
            result = fine
            result%evaluations = result%evaluations+coarse%evaluations
            return
        end if

        result = fine
        result%coarse_nroots = coarse%nroots
        result%evaluations = fine%evaluations+coarse%evaluations
        if (coarse%nroots /= fine%nroots) then
            result%status = GC_OPERATIONAL_ROOT_REFINEMENT_MISMATCH
            result%complete = .false.
            return
        end if
        coordinate_tolerance = root_coordinate_tolerance(domain_lo, &
            domain_hi, options)
        do i = 1, fine%nroots
            if (abs(coarse%roots(i)%x-fine%roots(i)%x) > &
                    4.0_dp*coordinate_tolerance) then
                result%status = GC_OPERATIONAL_ROOT_REFINEMENT_MISMATCH
                result%complete = .false.
                return
            end if
        end do
        result%coarse_fine_agree = .true.
        result%complete = .true.
        result%status = GC_OPERATIONAL_ROOT_SUCCESS
    end subroutine find_gc_operational_scalar_roots

    subroutine scan_roots(evaluate, domain_lo, domain_hi, interval_count, &
            options, result)
        procedure(gc_operational_scalar_i) :: evaluate
        real(dp), intent(in) :: domain_lo, domain_hi
        integer, intent(in) :: interval_count
        type(gc_operational_root_options_t), intent(in) :: options
        type(gc_operational_root_result_t), intent(out) :: result

        real(dp), allocatable :: x(:), value(:), derivative(:)
        real(dp) :: extremum_x, extremum_value, extremum_derivative
        real(dp) :: value_scale, degeneracy_tolerance
        integer :: i, local_status

        result = gc_operational_root_result_t()
        allocate(result%roots(0))
        result%status = GC_OPERATIONAL_ROOT_SUCCESS
        allocate(x(0:interval_count), value(0:interval_count), &
            derivative(0:interval_count))
        do i = 0, interval_count
            x(i) = domain_lo+(domain_hi-domain_lo)*real(i, dp) &
                /real(interval_count, dp)
            if (i == 0) x(i) = domain_lo
            if (i == interval_count) x(i) = domain_hi
            call evaluate_checked(evaluate, x(i), value(i), derivative(i), &
                result, local_status)
            if (local_status /= GC_OPERATIONAL_ROOT_SUCCESS) then
                call fail_result(result, local_status)
                return
            end if
        end do

        value_scale = max(1.0_dp, maxval(abs(value)))
        degeneracy_tolerance = 256.0_dp*epsilon(1.0_dp)*value_scale
        do i = 1, interval_count
            if (value(i-1) == 0.0_dp) then
                call append_exact_root(x(i-1), value(i-1), derivative(i-1), &
                    domain_lo, domain_hi, options, result, local_status)
                if (local_status /= GC_OPERATIONAL_ROOT_SUCCESS) then
                    call fail_result(result, local_status)
                    return
                end if
            end if

            if (derivative(i-1) == 0.0_dp .and. derivative(i) == 0.0_dp) then
                if (abs(value(i-1)) <= degeneracy_tolerance .or. &
                        abs(value(i)) <= degeneracy_tolerance) then
                    call fail_result(result, GC_OPERATIONAL_ROOT_DEGENERATE)
                    return
                end if
            else if (derivative(i-1)*derivative(i) <= 0.0_dp) then
                call refine_derivative_root(evaluate, x(i-1), x(i), &
                    derivative(i-1), derivative(i), options, result, &
                    extremum_x, extremum_value, extremum_derivative, &
                    local_status)
                if (local_status /= GC_OPERATIONAL_ROOT_SUCCESS) then
                    call fail_result(result, local_status)
                    return
                end if
                if (abs(extremum_value) <= degeneracy_tolerance) then
                    call fail_result(result, GC_OPERATIONAL_ROOT_DEGENERATE)
                    return
                end if
                call refine_function_bracket(evaluate, x(i-1), extremum_x, &
                    value(i-1), extremum_value, domain_lo, domain_hi, &
                    options, result, local_status)
                if (local_status /= GC_OPERATIONAL_ROOT_SUCCESS) then
                    call fail_result(result, local_status)
                    return
                end if
                call refine_function_bracket(evaluate, extremum_x, x(i), &
                    extremum_value, value(i), domain_lo, domain_hi, options, &
                    result, local_status)
                if (local_status /= GC_OPERATIONAL_ROOT_SUCCESS) then
                    call fail_result(result, local_status)
                    return
                end if
                associate (unused_extremum_derivative => extremum_derivative)
                end associate
            else
                call refine_function_bracket(evaluate, x(i-1), x(i), &
                    value(i-1), value(i), domain_lo, domain_hi, options, &
                    result, local_status)
                if (local_status /= GC_OPERATIONAL_ROOT_SUCCESS) then
                    call fail_result(result, local_status)
                    return
                end if
            end if
        end do
        if (value(interval_count) == 0.0_dp) then
            call append_exact_root(x(interval_count), value(interval_count), &
                derivative(interval_count), domain_lo, domain_hi, options, &
                result, local_status)
            if (local_status /= GC_OPERATIONAL_ROOT_SUCCESS) then
                call fail_result(result, local_status)
                return
            end if
        end if
        result%complete = .true.
        result%status = GC_OPERATIONAL_ROOT_SUCCESS
    end subroutine scan_roots

    subroutine refine_function_bracket(evaluate, x_lo, x_hi, f_lo, f_hi, &
            domain_lo, domain_hi, options, result, status)
        procedure(gc_operational_scalar_i) :: evaluate
        real(dp), intent(in) :: x_lo, x_hi, f_lo, f_hi, domain_lo, domain_hi
        type(gc_operational_root_options_t), intent(in) :: options
        type(gc_operational_root_result_t), intent(inout) :: result
        integer, intent(out) :: status

        type(gc_operational_root_t) :: root
        real(dp) :: a, b, fa, fb, midpoint, fm, dfm, tolerance
        integer :: iteration, local_status

        status = GC_OPERATIONAL_ROOT_SUCCESS
        if (x_hi < x_lo) then
            status = GC_OPERATIONAL_ROOT_INVALID_INPUT
            return
        end if
        if (f_lo == 0.0_dp) then
            call append_exact_root(x_lo, f_lo, 0.0_dp, domain_lo, domain_hi, &
                options, result, status, evaluate)
            return
        end if
        if (f_hi == 0.0_dp) then
            call append_exact_root(x_hi, f_hi, 0.0_dp, domain_lo, domain_hi, &
                options, result, status, evaluate)
            return
        end if
        if (f_lo*f_hi > 0.0_dp .or. x_hi == x_lo) return

        a = x_lo
        b = x_hi
        fa = f_lo
        fb = f_hi
        tolerance = root_coordinate_tolerance(domain_lo, domain_hi, options)
        do iteration = 1, options%maximum_iterations
            midpoint = 0.5_dp*(a+b)
            call evaluate_checked(evaluate, midpoint, fm, dfm, result, &
                local_status)
            if (local_status /= GC_OPERATIONAL_ROOT_SUCCESS) then
                status = local_status
                return
            end if
            if (fm == 0.0_dp .or. b-a <= tolerance) exit
            if (opposite_sign(fa, fm)) then
                b = midpoint
                fb = fm
            else
                a = midpoint
                fa = fm
            end if
        end do
        if (iteration > options%maximum_iterations) then
            status = GC_OPERATIONAL_ROOT_NONCONVERGED
            return
        end if
        root%x = midpoint
        root%bracket_lo = a
        root%bracket_hi = b
        root%residual = fm
        root%derivative = dfm
        root%bracketed = opposite_sign(fa, fb) .or. fm == 0.0_dp
        root%endpoint_root = midpoint == domain_lo .or. midpoint == domain_hi
        root%simple = derivative_is_resolved(dfm, fm)
        if (.not. root%simple) then
            status = GC_OPERATIONAL_ROOT_DEGENERATE
            return
        end if
        call append_root(root, options, result, status)
    end subroutine refine_function_bracket

    subroutine refine_derivative_root(evaluate, x_lo, x_hi, df_lo, df_hi, &
            options, result, root_x, root_value, root_derivative, status)
        procedure(gc_operational_scalar_i) :: evaluate
        real(dp), intent(in) :: x_lo, x_hi, df_lo, df_hi
        type(gc_operational_root_options_t), intent(in) :: options
        type(gc_operational_root_result_t), intent(inout) :: result
        real(dp), intent(out) :: root_x, root_value, root_derivative
        integer, intent(out) :: status

        real(dp) :: a, b, da, db, midpoint, fm, dm, tolerance
        integer :: iteration, local_status

        root_x = 0.0_dp
        root_value = 0.0_dp
        root_derivative = 0.0_dp
        status = GC_OPERATIONAL_ROOT_SUCCESS
        if (df_lo == 0.0_dp) then
            root_x = x_lo
            call evaluate_checked(evaluate, root_x, root_value, &
                root_derivative, result, status)
            return
        end if
        if (df_hi == 0.0_dp) then
            root_x = x_hi
            call evaluate_checked(evaluate, root_x, root_value, &
                root_derivative, result, status)
            return
        end if
        if (df_lo*df_hi > 0.0_dp) then
            status = GC_OPERATIONAL_ROOT_INVALID_INPUT
            return
        end if
        a = x_lo
        b = x_hi
        da = df_lo
        db = df_hi
        tolerance = options%relative_tolerance*max(1.0_dp, abs(x_hi-x_lo))
        do iteration = 1, options%maximum_iterations
            midpoint = 0.5_dp*(a+b)
            call evaluate_checked(evaluate, midpoint, fm, dm, result, &
                local_status)
            if (local_status /= GC_OPERATIONAL_ROOT_SUCCESS) then
                status = local_status
                return
            end if
            if (dm == 0.0_dp .or. b-a <= tolerance) then
                root_x = midpoint
                root_value = fm
                root_derivative = dm
                return
            end if
            if (opposite_sign(da, dm)) then
                b = midpoint
                db = dm
            else
                a = midpoint
                da = dm
            end if
        end do
        associate (unused_db => db)
        end associate
        status = GC_OPERATIONAL_ROOT_NONCONVERGED
    end subroutine refine_derivative_root

    subroutine append_exact_root(x, value, derivative, domain_lo, domain_hi, &
            options, result, status, evaluate)
        real(dp), intent(in) :: x, value, derivative, domain_lo, domain_hi
        type(gc_operational_root_options_t), intent(in) :: options
        type(gc_operational_root_result_t), intent(inout) :: result
        integer, intent(out) :: status
        procedure(gc_operational_scalar_i), optional :: evaluate

        type(gc_operational_root_t) :: root
        real(dp) :: checked_value, checked_derivative
        integer :: local_status

        checked_value = value
        checked_derivative = derivative
        if (present(evaluate)) then
            call evaluate_checked(evaluate, x, checked_value, &
                checked_derivative, result, local_status)
            if (local_status /= GC_OPERATIONAL_ROOT_SUCCESS) then
                status = local_status
                return
            end if
        end if
        root%x = x
        root%bracket_lo = x
        root%bracket_hi = x
        root%residual = checked_value
        root%derivative = checked_derivative
        root%bracketed = .true.
        root%endpoint_root = x == domain_lo .or. x == domain_hi
        root%simple = derivative_is_resolved(checked_derivative, checked_value)
        if (.not. root%simple) then
            status = GC_OPERATIONAL_ROOT_DEGENERATE
            return
        end if
        call append_root(root, options, result, status)
    end subroutine append_exact_root

    subroutine append_root(root, options, result, status)
        type(gc_operational_root_t), intent(in) :: root
        type(gc_operational_root_options_t), intent(in) :: options
        type(gc_operational_root_result_t), intent(inout) :: result
        integer, intent(out) :: status

        type(gc_operational_root_t), allocatable :: enlarged(:)
        real(dp) :: separation
        integer :: i

        status = GC_OPERATIONAL_ROOT_SUCCESS
        do i = 1, result%nroots
            if (root%x == result%roots(i)%x) return
            separation = 128.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
                abs(root%x), abs(result%roots(i)%x))
            if (abs(root%x-result%roots(i)%x) <= separation) then
                status = GC_OPERATIONAL_ROOT_UNRESOLVED_SEPARATION
                return
            end if
        end do
        if (result%nroots >= options%maximum_roots) then
            status = GC_OPERATIONAL_ROOT_CAPACITY
            return
        end if
        allocate(enlarged(result%nroots+1))
        if (result%nroots > 0) enlarged(1:result%nroots) = result%roots
        enlarged(result%nroots+1) = root
        call move_alloc(enlarged, result%roots)
        result%nroots = result%nroots+1
    end subroutine append_root

    subroutine evaluate_checked(evaluate, x, value, derivative, result, status)
        procedure(gc_operational_scalar_i) :: evaluate
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value, derivative
        type(gc_operational_root_result_t), intent(inout) :: result
        integer, intent(out) :: status

        call evaluate(x, value, derivative, status)
        result%evaluations = result%evaluations+1
        if (status /= 0) then
            status = GC_OPERATIONAL_ROOT_CALLBACK_FAILURE
            return
        end if
        if (.not. all(ieee_is_finite([value, derivative]))) then
            status = GC_OPERATIONAL_ROOT_NONFINITE
            return
        end if
        status = GC_OPERATIONAL_ROOT_SUCCESS
    end subroutine evaluate_checked

    subroutine fail_result(result, status)
        type(gc_operational_root_result_t), intent(inout) :: result
        integer, intent(in) :: status

        result%status = status
        result%nroots = 0
        result%complete = .false.
        result%coarse_fine_agree = .false.
        if (allocated(result%roots)) deallocate(result%roots)
        allocate(result%roots(0))
    end subroutine fail_result

    pure logical function valid_options(domain_lo, domain_hi, options)
        real(dp), intent(in) :: domain_lo, domain_hi
        type(gc_operational_root_options_t), intent(in) :: options

        valid_options = all(ieee_is_finite([domain_lo, domain_hi, &
            options%relative_tolerance])) .and. domain_hi > domain_lo .and. &
            options%initial_intervals >= 1 .and. &
            options%initial_intervals <= huge(options%initial_intervals)/2 .and. &
            options%maximum_iterations >= 1 .and. &
            options%maximum_roots >= 1 .and. options%relative_tolerance > 0.0_dp
    end function valid_options

    pure real(dp) function root_coordinate_tolerance(domain_lo, domain_hi, &
            options)
        real(dp), intent(in) :: domain_lo, domain_hi
        type(gc_operational_root_options_t), intent(in) :: options

        root_coordinate_tolerance = max(options%relative_tolerance* &
            abs(domain_hi-domain_lo), 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, abs(domain_lo), abs(domain_hi)))
    end function root_coordinate_tolerance

    pure logical function derivative_is_resolved(derivative, value)
        real(dp), intent(in) :: derivative, value

        derivative_is_resolved = abs(derivative) > 256.0_dp*epsilon(1.0_dp) &
            *max(1.0_dp, abs(value))
    end function derivative_is_resolved

    pure logical function opposite_sign(left, right)
        real(dp), intent(in) :: left, right

        opposite_sign = (left < 0.0_dp .and. right > 0.0_dp) .or. &
            (left > 0.0_dp .and. right < 0.0_dp)
    end function opposite_sign

end module neort_gc_operational_scalar_roots
