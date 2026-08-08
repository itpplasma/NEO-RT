module neort_gc_operational_fixed_points
    !! Interior canonical-momentum fixed points for one allowed cut interval.
    !!
    !! The canonical jet and flow discriminant are supplied by generated
    !! physics callbacks.  This module owns root orchestration, endpoint
    !! exclusion, classification, and explicit failure on degenerate points.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_operational_scalar_roots, only: &
        GC_OPERATIONAL_ROOT_SUCCESS, find_gc_operational_scalar_roots, &
        gc_operational_root_options_t, gc_operational_root_t, &
        gc_operational_root_result_t
    implicit none
    private

    integer, parameter, public :: GC_FIXED_POINT_SUCCESS = 0
    integer, parameter, public :: GC_FIXED_POINT_INVALID_INPUT = 1
    integer, parameter, public :: GC_FIXED_POINT_CANONICAL_FAILURE = 2
    integer, parameter, public :: GC_FIXED_POINT_ROOT_FAILURE = 3
    integer, parameter, public :: GC_FIXED_POINT_CLASSIFIER_FAILURE = 4
    integer, parameter, public :: GC_FIXED_POINT_DEGENERATE = 5
    integer, parameter, public :: GC_FIXED_POINT_NONFINITE = 6

    integer, parameter, public :: GC_FIXED_POINT_O = 1
    integer, parameter, public :: GC_FIXED_POINT_X = 2

    type, public :: gc_operational_fixed_point_t
        integer :: point_id = 0
        integer :: kind = 0
        real(dp) :: x = 0.0_dp
        real(dp) :: canonical_momentum = 0.0_dp
        real(dp) :: first_derivative = 0.0_dp
        real(dp) :: second_derivative = 0.0_dp
        real(dp) :: flow_discriminant = 0.0_dp
        type(gc_operational_root_t) :: stationary_root
    end type gc_operational_fixed_point_t

    type, public :: gc_operational_fixed_point_result_t
        integer :: status = GC_FIXED_POINT_INVALID_INPUT
        integer :: npoints = 0
        integer :: n_o_points = 0
        integer :: n_x_points = 0
        logical :: complete = .false.
        type(gc_operational_fixed_point_t), allocatable :: points(:)
    end type gc_operational_fixed_point_result_t

    abstract interface
        subroutine gc_operational_canonical_jet_i(x, canonical_momentum, &
                first_derivative, second_derivative, status)
            import :: dp
            real(dp), intent(in) :: x
            real(dp), intent(out) :: canonical_momentum
            real(dp), intent(out) :: first_derivative, second_derivative
            integer, intent(out) :: status
        end subroutine gc_operational_canonical_jet_i

        subroutine gc_operational_fixed_point_classifier_i(x, discriminant, &
                status)
            import :: dp
            real(dp), intent(in) :: x
            real(dp), intent(out) :: discriminant
            integer, intent(out) :: status
        end subroutine gc_operational_fixed_point_classifier_i
    end interface

    public :: find_gc_operational_fixed_points
    public :: gc_operational_canonical_jet_i
    public :: gc_operational_fixed_point_classifier_i

contains

    subroutine find_gc_operational_fixed_points(evaluate_canonical, classify, &
            x_lo, x_hi, root_options, result)
        procedure(gc_operational_canonical_jet_i) :: evaluate_canonical
        procedure(gc_operational_fixed_point_classifier_i) :: classify
        real(dp), intent(in) :: x_lo, x_hi
        type(gc_operational_root_options_t), intent(in) :: root_options
        type(gc_operational_fixed_point_result_t), intent(out) :: result

        type(gc_operational_root_result_t) :: roots
        type(gc_operational_fixed_point_t), allocatable :: points(:)
        real(dp) :: p_star, first_derivative, second_derivative
        real(dp) :: discriminant, endpoint_tolerance, discriminant_tolerance
        integer :: i, count, local_status

        result = gc_operational_fixed_point_result_t()
        allocate(result%points(0))
        if (.not. all(ieee_is_finite([x_lo, x_hi]))) return
        if (x_hi <= x_lo) return

        call find_gc_operational_scalar_roots(stationary_value, x_lo, x_hi, &
            root_options, roots)
        if (roots%status /= GC_OPERATIONAL_ROOT_SUCCESS .or. &
                .not. roots%complete) then
            result%status = GC_FIXED_POINT_ROOT_FAILURE
            return
        end if

        endpoint_tolerance = max(root_options%relative_tolerance* &
            abs(x_hi-x_lo), 64.0_dp*epsilon(1.0_dp)* &
            max(1.0_dp, abs(x_lo), abs(x_hi)))
        allocate(points(roots%nroots))
        count = 0
        do i = 1, roots%nroots
            if (roots%roots(i)%x <= x_lo+endpoint_tolerance .or. &
                    roots%roots(i)%x >= x_hi-endpoint_tolerance) cycle
            call evaluate_canonical(roots%roots(i)%x, p_star, &
                first_derivative, second_derivative, local_status)
            if (local_status /= 0) then
                result%status = GC_FIXED_POINT_CANONICAL_FAILURE
                return
            end if
            if (.not. all(ieee_is_finite([p_star, first_derivative, &
                    second_derivative]))) then
                result%status = GC_FIXED_POINT_NONFINITE
                return
            end if
            call classify(roots%roots(i)%x, discriminant, local_status)
            if (local_status /= 0) then
                result%status = GC_FIXED_POINT_CLASSIFIER_FAILURE
                return
            end if
            if (.not. ieee_is_finite(discriminant)) then
                result%status = GC_FIXED_POINT_NONFINITE
                return
            end if
            discriminant_tolerance = 256.0_dp*epsilon(1.0_dp)* &
                max(1.0_dp, abs(discriminant))
            if (abs(discriminant) <= discriminant_tolerance) then
                result%status = GC_FIXED_POINT_DEGENERATE
                return
            end if
            count = count+1
            points(count)%point_id = count
            points(count)%x = roots%roots(i)%x
            points(count)%canonical_momentum = p_star
            points(count)%first_derivative = first_derivative
            points(count)%second_derivative = second_derivative
            points(count)%flow_discriminant = discriminant
            points(count)%stationary_root = roots%roots(i)
            if (discriminant < 0.0_dp) then
                points(count)%kind = GC_FIXED_POINT_O
                result%n_o_points = result%n_o_points+1
            else
                points(count)%kind = GC_FIXED_POINT_X
                result%n_x_points = result%n_x_points+1
            end if
        end do

        if (allocated(result%points)) deallocate(result%points)
        allocate(result%points(count))
        if (count > 0) result%points = points(1:count)
        result%npoints = count
        result%complete = .true.
        result%status = GC_FIXED_POINT_SUCCESS

    contains

        subroutine stationary_value(x, value, derivative, status)
            real(dp), intent(in) :: x
            real(dp), intent(out) :: value, derivative
            integer, intent(out) :: status
            real(dp) :: unused_canonical

            call evaluate_canonical(x, unused_canonical, value, derivative, &
                status)
        end subroutine stationary_value

    end subroutine find_gc_operational_fixed_points

end module neort_gc_operational_fixed_points
