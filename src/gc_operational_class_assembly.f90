module neort_gc_operational_class_assembly
    !! Assemble the sorted, POTATO-compatible classes for one cut interval.
    !!
    !! This module is deliberately data-oriented.  The allowed-interval
    !! endpoints, fixed-point values, and separatrix partner values are
    !! supplied by their respective operational stages.  Assembly only
    !! validates, orders, and combines those values; it does not evaluate a
    !! field, differentiate a quantity, or invent a topology.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_operational_fixed_points, only: &
        GC_FIXED_POINT_O, GC_FIXED_POINT_SUCCESS, GC_FIXED_POINT_X, &
        gc_operational_fixed_point_t, gc_operational_fixed_point_result_t
    use neort_gc_operational_partner_crossings, only: &
        GC_PARTNER_BOUNDARY_USUAL, GC_PARTNER_BOUNDARY_X, GC_PARTNER_SUCCESS, &
        gc_operational_partner_pair_t, gc_operational_partner_result_t, &
        gc_operational_separatrix_boundary_t
    implicit none
    private

    integer, parameter, public :: GC_CLASS_ASSEMBLY_SUCCESS = 0
    integer, parameter, public :: GC_CLASS_ASSEMBLY_INVALID_INPUT = 1
    integer, parameter, public :: GC_CLASS_ASSEMBLY_INCOMPLETE_FIXED_POINTS = 2
    integer, parameter, public :: GC_CLASS_ASSEMBLY_INCOMPLETE_PARTNERS = 3
    integer, parameter, public :: GC_CLASS_ASSEMBLY_INVALID_BOUNDARY = 4
    integer, parameter, public :: GC_CLASS_ASSEMBLY_DUPLICATE_BOUNDARY = 5
    integer, parameter, public :: GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR = 6
    integer, parameter, public :: GC_CLASS_ASSEMBLY_ENDPOINT_MISMATCH = 7
    integer, parameter, public :: GC_CLASS_ASSEMBLY_NONFINITE = 8

    integer, parameter, public :: GC_CLASS_BOUNDARY_REGULAR = 1
    integer, parameter, public :: GC_CLASS_BOUNDARY_TURNING = 2
    integer, parameter, public :: GC_CLASS_BOUNDARY_SEPARATRIX = 3
    integer, parameter, public :: GC_CLASS_BOUNDARY_X = 4

    type, public :: gc_operational_allowed_interval_t
        !! x is physical R, canonicalized to increase from HFS to LFS.
        !! sigma is the parallel-velocity branch and is not an orientation.
        integer :: component_id = 0
        integer :: sigma = 0
        real(dp) :: x_lo = 0.0_dp
        real(dp) :: x_hi = 0.0_dp
        integer :: left_kind = 0
        integer :: right_kind = 0
        real(dp) :: canonical_lo = 0.0_dp
        real(dp) :: canonical_hi = 0.0_dp
    end type gc_operational_allowed_interval_t

    type, public :: gc_operational_class_interval_t
        integer :: class_id = 0
        integer :: component_id = 0
        integer :: sigma = 0
        real(dp) :: x_lo = 0.0_dp
        real(dp) :: x_hi = 0.0_dp
        integer :: left_kind = 0
        integer :: right_kind = 0
        integer :: ifuntype = 0
        real(dp) :: canonical_lo = 0.0_dp
        real(dp) :: canonical_hi = 0.0_dp
        real(dp) :: canonical_total_variation = 0.0_dp
        integer :: left_boundary_id = 0
        integer :: right_boundary_id = 0
    end type gc_operational_class_interval_t

    type, public :: gc_operational_class_assembly_result_t
        integer :: status = GC_CLASS_ASSEMBLY_INVALID_INPUT
        integer :: nclasses = 0
        logical :: complete = .false.
        type(gc_operational_class_interval_t), allocatable :: classes(:)
    end type gc_operational_class_assembly_result_t

    public :: assemble_gc_operational_classes

contains

    subroutine assemble_gc_operational_classes(allowed, fixed_points, partners, &
            result)
        type(gc_operational_allowed_interval_t), intent(in) :: allowed
        type(gc_operational_fixed_point_result_t), intent(in) :: fixed_points
        type(gc_operational_partner_result_t), intent(in) :: partners
        type(gc_operational_class_assembly_result_t), intent(out) :: result

        type(gc_operational_fixed_point_t), allocatable :: points(:)
        type(gc_operational_separatrix_boundary_t), allocatable :: &
            boundaries(:)
        type(gc_operational_class_interval_t), allocatable :: classes(:)
        real(dp) :: coordinate_tolerance
        integer :: i, j, status, npoints

        result = gc_operational_class_assembly_result_t()
        allocate(result%classes(0))

        status = validate_allowed(allowed)
        if (status /= GC_CLASS_ASSEMBLY_SUCCESS) then
            result%status = status
            return
        end if
        coordinate_tolerance = coordinate_scale(allowed)*256.0_dp* &
            epsilon(1.0_dp)
        call copy_and_validate_points(fixed_points, allowed, points, npoints, &
            status)
        if (status /= GC_CLASS_ASSEMBLY_SUCCESS) then
            result%status = status
            return
        end if
        call copy_and_validate_boundaries(partners, allowed, points, npoints, &
            coordinate_tolerance, boundaries, status)
        if (status /= GC_CLASS_ASSEMBLY_SUCCESS) then
            result%status = status
            return
        end if

        do i = 1, npoints
            if (near(points(i)%x, allowed%x_lo, coordinate_tolerance)) then
                result%status = GC_CLASS_ASSEMBLY_INVALID_BOUNDARY
                return
            end if
            if (near(points(i)%x, allowed%x_hi, coordinate_tolerance)) then
                result%status = GC_CLASS_ASSEMBLY_INVALID_BOUNDARY
                return
            end if
            do j = 1, size(boundaries)
                if (near(points(i)%x, boundaries(j)%x, &
                    coordinate_tolerance)) then
                    if (.not. canonical_near(points(i)%canonical_momentum, &
                        boundaries(j)%canonical_momentum)) then
                        result%status = GC_CLASS_ASSEMBLY_ENDPOINT_MISMATCH
                        return
                    end if
                end if
            end do
        end do

        allocate(classes(size(boundaries)+1))
        do i = 1, size(classes)
            classes(i)%class_id = i
            classes(i)%component_id = allowed%component_id
            classes(i)%sigma = allowed%sigma
            if (i == 1) then
                classes(i)%x_lo = allowed%x_lo
                classes(i)%canonical_lo = allowed%canonical_lo
                classes(i)%left_kind = allowed%left_kind
                classes(i)%left_boundary_id = 0
            else
                classes(i)%x_lo = boundaries(i-1)%x
                classes(i)%canonical_lo = boundaries(i-1)%canonical_momentum
                classes(i)%left_kind = boundaries(i-1)%kind
                classes(i)%left_boundary_id = boundaries(i-1)%boundary_id
            end if
            if (i == size(classes)) then
                classes(i)%x_hi = allowed%x_hi
                classes(i)%canonical_hi = allowed%canonical_hi
                classes(i)%right_kind = allowed%right_kind
                classes(i)%right_boundary_id = 0
            else
                classes(i)%x_hi = boundaries(i)%x
                classes(i)%canonical_hi = boundaries(i)%canonical_momentum
                classes(i)%right_kind = boundaries(i)%kind
                classes(i)%right_boundary_id = boundaries(i)%boundary_id
            end if
            classes(i)%ifuntype = 10*classes(i)%left_kind + &
                classes(i)%right_kind
            call class_total_variation(classes(i), points, npoints, &
                coordinate_tolerance, classes(i)%canonical_total_variation, &
                status)
            if (status /= GC_CLASS_ASSEMBLY_SUCCESS) then
                result%status = status
                return
            end if
        end do

        result%classes = classes
        result%nclasses = size(classes)
        result%complete = .true.
        result%status = GC_CLASS_ASSEMBLY_SUCCESS
    end subroutine assemble_gc_operational_classes

    integer function validate_allowed(allowed)
        type(gc_operational_allowed_interval_t), intent(in) :: allowed

        validate_allowed = GC_CLASS_ASSEMBLY_INVALID_INPUT
        if (allowed%component_id <= 0) return
        if (allowed%sigma /= 1 .and. allowed%sigma /= -1) return
        if (.not. all(ieee_is_finite([allowed%x_lo, allowed%x_hi, &
            allowed%canonical_lo, allowed%canonical_hi]))) then
            validate_allowed = GC_CLASS_ASSEMBLY_NONFINITE
            return
        end if
        if (allowed%x_hi <= allowed%x_lo) return
        if (allowed%left_kind < GC_CLASS_BOUNDARY_REGULAR .or. &
            allowed%left_kind > GC_CLASS_BOUNDARY_TURNING) return
        if (allowed%right_kind < GC_CLASS_BOUNDARY_REGULAR .or. &
            allowed%right_kind > GC_CLASS_BOUNDARY_TURNING) return
        validate_allowed = GC_CLASS_ASSEMBLY_SUCCESS
    end function validate_allowed

    subroutine copy_and_validate_points(input, allowed, points, npoints, status)
        type(gc_operational_fixed_point_result_t), intent(in) :: input
        type(gc_operational_allowed_interval_t), intent(in) :: allowed
        type(gc_operational_fixed_point_t), allocatable, intent(out) :: points(:)
        integer, intent(out) :: npoints, status

        integer :: actual_o_count, actual_x_count, i, j
        type(gc_operational_fixed_point_t) :: held
        real(dp) :: tolerance

        npoints = 0
        status = GC_CLASS_ASSEMBLY_INCOMPLETE_FIXED_POINTS
        allocate(points(0))
        if (input%status /= GC_FIXED_POINT_SUCCESS) return
        if (.not. input%complete) return
        if (input%npoints < 0 .or. input%n_o_points < 0 .or. &
            input%n_x_points < 0) then
            status = GC_CLASS_ASSEMBLY_INVALID_INPUT
            return
        end if
        if (allocated(input%points)) then
            if (size(input%points) /= input%npoints) then
                status = GC_CLASS_ASSEMBLY_INVALID_INPUT
                return
            end if
        else
            if (input%npoints /= 0) then
                status = GC_CLASS_ASSEMBLY_INVALID_INPUT
                return
            end if
        end if
        if (input%n_o_points+input%n_x_points /= input%npoints) then
            status = GC_CLASS_ASSEMBLY_INVALID_INPUT
            return
        end if
        deallocate(points)
        allocate(points(input%npoints))
        if (input%npoints > 0) points = input%points
        tolerance = coordinate_scale(allowed)*256.0_dp*epsilon(1.0_dp)
        actual_o_count = 0
        actual_x_count = 0
        do i = 1, size(points)
            if (points(i)%point_id <= 0) then
                status = GC_CLASS_ASSEMBLY_INVALID_INPUT
                return
            end if
            if (points(i)%kind /= GC_FIXED_POINT_O .and. &
                points(i)%kind /= GC_FIXED_POINT_X) then
                status = GC_CLASS_ASSEMBLY_INVALID_INPUT
                return
            end if
            if (points(i)%kind == GC_FIXED_POINT_O) &
                actual_o_count = actual_o_count+1
            if (points(i)%kind == GC_FIXED_POINT_X) &
                actual_x_count = actual_x_count+1
            if (.not. all(ieee_is_finite([points(i)%x, &
                points(i)%canonical_momentum, points(i)%first_derivative, &
                points(i)%second_derivative, &
                points(i)%flow_discriminant]))) then
                status = GC_CLASS_ASSEMBLY_NONFINITE
                return
            end if
            if (points(i)%x < allowed%x_lo-tolerance .or. &
                points(i)%x > allowed%x_hi+tolerance) then
                status = GC_CLASS_ASSEMBLY_INVALID_INPUT
                return
            end if
        end do
        if (actual_o_count /= input%n_o_points .or. &
            actual_x_count /= input%n_x_points) then
            status = GC_CLASS_ASSEMBLY_INVALID_INPUT
            return
        end if
        do i = 1, size(points)-1
            do j = i+1, size(points)
                if (points(i)%point_id == points(j)%point_id) then
                    status = GC_CLASS_ASSEMBLY_INVALID_INPUT
                    return
                end if
            end do
        end do
        do i = 2, size(points)
            held = points(i)
            j = i-1
            do while (j >= 1)
                if (points(j)%x <= held%x) exit
                points(j+1) = points(j)
                j = j-1
            end do
            points(j+1) = held
        end do
        do i = 2, size(points)
            if (near(points(i-1)%x, points(i)%x, tolerance)) then
                status = GC_CLASS_ASSEMBLY_DUPLICATE_BOUNDARY
                return
            end if
        end do
        npoints = size(points)
        status = GC_CLASS_ASSEMBLY_SUCCESS
    end subroutine copy_and_validate_points

    subroutine copy_and_validate_boundaries(input, allowed, points, npoints, &
            coordinate_tolerance, boundaries, status)
        type(gc_operational_partner_result_t), intent(in) :: input
        type(gc_operational_allowed_interval_t), intent(in) :: allowed
        type(gc_operational_fixed_point_t), intent(in) :: points(:)
        integer, intent(in) :: npoints
        real(dp), intent(in) :: coordinate_tolerance
        type(gc_operational_separatrix_boundary_t), allocatable, &
            intent(out) :: boundaries(:)
        integer, intent(out) :: status

        integer :: i, j, count, x_index, x_boundary_count
        integer :: expected_boundary_count
        real(dp) :: residual_limit
        type(gc_operational_separatrix_boundary_t) :: held
        type(gc_operational_partner_pair_t), allocatable :: pairs(:)

        status = GC_CLASS_ASSEMBLY_INCOMPLETE_PARTNERS
        allocate(boundaries(0))
        if (input%status /= GC_PARTNER_SUCCESS) return
        if (.not. input%complete) return
        if (input%nboundaries < 0 .or. input%npairs < 0) then
            status = GC_CLASS_ASSEMBLY_INVALID_INPUT
            return
        end if
        if (.not. all(ieee_is_finite([input%residual_absolute_tolerance, &
            input%residual_relative_tolerance]))) then
            status = GC_CLASS_ASSEMBLY_NONFINITE
            return
        end if
        if (input%residual_absolute_tolerance <= 0.0_dp .or. &
            input%residual_relative_tolerance < 0.0_dp) then
            status = GC_CLASS_ASSEMBLY_INVALID_INPUT
            return
        end if
        if (input%npairs /= count_x_points(points)) then
            status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
            return
        end if
        if (allocated(input%boundaries)) then
            if (size(input%boundaries) /= input%nboundaries) then
                status = GC_CLASS_ASSEMBLY_INVALID_INPUT
                return
            end if
        else
            if (input%nboundaries /= 0) then
                status = GC_CLASS_ASSEMBLY_INVALID_INPUT
                return
            end if
        end if
        if (allocated(input%pairs)) then
            if (size(input%pairs) /= input%npairs) then
                status = GC_CLASS_ASSEMBLY_INVALID_INPUT
                return
            end if
        else
            if (input%npairs /= 0) then
                status = GC_CLASS_ASSEMBLY_INVALID_INPUT
                return
            end if
        end if
        if (input%nboundaries == 0) then
            if (input%npairs /= 0) then
                status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                return
            end if
            status = GC_CLASS_ASSEMBLY_SUCCESS
            return
        end if
        if (input%npairs == 0) then
            status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
            return
        end if
        pairs = input%pairs
        deallocate(boundaries)
        allocate(boundaries(input%nboundaries))
        boundaries = input%boundaries
        do i = 1, size(pairs)
            if (pairs(i)%pair_id /= i .or. pairs(i)%x_point_id <= 0 .or. &
                pairs(i)%regular_crossing_count < 1 .or. &
                .not. ieee_is_finite(pairs(i)%canonical_momentum)) then
                status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                return
            end if
            x_index = find_point(points, pairs(i)%x_point_id)
            if (x_index == 0) then
                status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                return
            end if
            if (points(x_index)%kind /= GC_FIXED_POINT_X) then
                status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                return
            end if
            do j = 1, i-1
                if (pairs(j)%x_point_id == pairs(i)%x_point_id) then
                    status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                    return
                end if
            end do
        end do
        expected_boundary_count = size(pairs)
        do i = 1, size(pairs)
            expected_boundary_count = expected_boundary_count + &
                pairs(i)%regular_crossing_count
        end do
        if (size(boundaries) /= expected_boundary_count) then
            status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
            return
        end if
        do i = 1, size(boundaries)
            if (boundaries(i)%boundary_id < 1 .or. &
                boundaries(i)%boundary_id > size(boundaries)) then
                status = GC_CLASS_ASSEMBLY_INVALID_BOUNDARY
                return
            end if
            do j = 1, i-1
                if (boundaries(j)%boundary_id == boundaries(i)%boundary_id) then
                    status = GC_CLASS_ASSEMBLY_INVALID_BOUNDARY
                    return
                end if
            end do
            if (boundaries(i)%kind /= GC_PARTNER_BOUNDARY_USUAL .and. &
                boundaries(i)%kind /= GC_PARTNER_BOUNDARY_X) then
                status = GC_CLASS_ASSEMBLY_INVALID_BOUNDARY
                return
            end if
            if (.not. all(ieee_is_finite([boundaries(i)%x, &
                boundaries(i)%canonical_momentum, &
                boundaries(i)%canonical_residual]))) then
                status = GC_CLASS_ASSEMBLY_NONFINITE
                return
            end if
            if (boundaries(i)%x <= allowed%x_lo+coordinate_tolerance .or. &
                boundaries(i)%x >= allowed%x_hi-coordinate_tolerance) then
                status = GC_CLASS_ASSEMBLY_INVALID_BOUNDARY
                return
            end if
            if (boundaries(i)%pair_id < 1 .or. &
                boundaries(i)%pair_id > size(pairs)) then
                status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                return
            end if
            if (boundaries(i)%x_point_id /= &
                pairs(boundaries(i)%pair_id)%x_point_id) then
                status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                return
            end if
            if (.not. canonical_near(boundaries(i)%canonical_momentum, &
                pairs(boundaries(i)%pair_id)%canonical_momentum)) then
                status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                return
            end if
            residual_limit = input%residual_absolute_tolerance + &
                input%residual_relative_tolerance*max(1.0_dp, &
                abs(pairs(boundaries(i)%pair_id)%canonical_momentum))
            if (abs(boundaries(i)%canonical_residual) > residual_limit) then
                status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                return
            end if
            if (boundaries(i)%kind == GC_PARTNER_BOUNDARY_X) then
                x_index = 0
                do j = 1, npoints
                    if (points(j)%point_id == boundaries(i)%x_point_id) then
                        x_index = j
                        exit
                    end if
                end do
                if (x_index == 0) then
                    status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                    return
                end if
                if (points(x_index)%kind /= GC_FIXED_POINT_X .or. &
                    .not. near(points(x_index)%x, boundaries(i)%x, &
                    coordinate_tolerance)) then
                    status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                    return
                end if
            end if
        end do
        do i = 1, size(pairs)
            count = 0
            x_boundary_count = 0
            do j = 1, size(boundaries)
                if (boundaries(j)%pair_id == pairs(i)%pair_id) then
                    if (boundaries(j)%kind == GC_PARTNER_BOUNDARY_USUAL) &
                        count = count+1
                    if (boundaries(j)%kind == GC_PARTNER_BOUNDARY_X) &
                        x_boundary_count = x_boundary_count+1
                end if
            end do
            if (count /= pairs(i)%regular_crossing_count .or. &
                x_boundary_count /= 1) then
                status = GC_CLASS_ASSEMBLY_INCONSISTENT_PAIR
                return
            end if
        end do
        do i = 2, size(boundaries)
            held = boundaries(i)
            j = i-1
            do while (j >= 1)
                if (boundaries(j)%x <= held%x) exit
                boundaries(j+1) = boundaries(j)
                j = j-1
            end do
            boundaries(j+1) = held
        end do
        do i = 2, size(boundaries)
            if (near(boundaries(i-1)%x, boundaries(i)%x, &
                coordinate_tolerance)) then
                status = GC_CLASS_ASSEMBLY_DUPLICATE_BOUNDARY
                return
            end if
        end do
        status = GC_CLASS_ASSEMBLY_SUCCESS
    end subroutine copy_and_validate_boundaries

    subroutine class_total_variation(class, points, npoints, tolerance, variation, &
            status)
        type(gc_operational_class_interval_t), intent(in) :: class
        type(gc_operational_fixed_point_t), intent(in) :: points(:)
        integer, intent(in) :: npoints
        real(dp), intent(in) :: tolerance
        real(dp), intent(out) :: variation
        integer, intent(out) :: status

        real(dp) :: previous
        integer :: i

        variation = 0.0_dp
        status = GC_CLASS_ASSEMBLY_SUCCESS
        previous = class%canonical_lo
        do i = 1, npoints
            if (points(i)%x <= class%x_lo+tolerance) cycle
            if (points(i)%x >= class%x_hi-tolerance) exit
            variation = variation+abs(points(i)%canonical_momentum-previous)
            previous = points(i)%canonical_momentum
        end do
        variation = variation+abs(class%canonical_hi-previous)
        if (.not. ieee_is_finite(variation)) then
            status = GC_CLASS_ASSEMBLY_NONFINITE
        end if
    end subroutine class_total_variation

    pure real(dp) function coordinate_scale(allowed)
        type(gc_operational_allowed_interval_t), intent(in) :: allowed

        coordinate_scale = max(1.0_dp, abs(allowed%x_lo), abs(allowed%x_hi), &
            abs(allowed%x_hi-allowed%x_lo))
    end function coordinate_scale

    pure logical function near(first, second, tolerance)
        real(dp), intent(in) :: first, second, tolerance

        near = abs(first-second) <= tolerance
    end function near

    pure logical function canonical_near(first, second)
        real(dp), intent(in) :: first, second

        canonical_near = abs(first-second) <= &
            1.0e-10_dp*max(1.0_dp, abs(first), abs(second))
    end function canonical_near

    pure integer function count_x_points(points)
        type(gc_operational_fixed_point_t), intent(in) :: points(:)
        integer :: i

        count_x_points = 0
        do i = 1, size(points)
            if (points(i)%kind == GC_FIXED_POINT_X) &
                count_x_points = count_x_points+1
        end do
    end function count_x_points

    pure integer function find_point(points, point_id)
        type(gc_operational_fixed_point_t), intent(in) :: points(:)
        integer, intent(in) :: point_id
        integer :: i

        find_point = 0
        do i = 1, size(points)
            if (points(i)%point_id == point_id) then
                find_point = i
                return
            end if
        end do
    end function find_point

end module neort_gc_operational_class_assembly
