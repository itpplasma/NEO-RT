module neort_gc_operational_partner_crossings
    !! Same-canonical-momentum partners for operational X points.
    !!
    !! Each X point is inserted explicitly and excluded from the simple-root
    !! search.  Every other crossing of its canonical-momentum level is found
    !! by the converged scalar root engine.  Missing partners, repeated X
    !! levels, endpoint coincidences, and unresolved roots fail closed.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_operational_fixed_points, only: &
        GC_FIXED_POINT_X, gc_operational_canonical_jet_i, &
        gc_operational_fixed_point_result_t
    use neort_gc_operational_scalar_roots, only: &
        GC_OPERATIONAL_ROOT_SUCCESS, find_gc_operational_scalar_roots, &
        gc_operational_root_options_t, gc_operational_root_result_t
    implicit none
    private

    integer, parameter, public :: GC_PARTNER_SUCCESS = 0
    integer, parameter, public :: GC_PARTNER_INVALID_INPUT = 1
    integer, parameter, public :: GC_PARTNER_CANONICAL_FAILURE = 2
    integer, parameter, public :: GC_PARTNER_ROOT_FAILURE = 3
    integer, parameter, public :: GC_PARTNER_MISSING = 4
    integer, parameter, public :: GC_PARTNER_REPEATED_X_LEVEL = 5
    integer, parameter, public :: GC_PARTNER_ENDPOINT_COINCIDENCE = 6
    integer, parameter, public :: GC_PARTNER_UNRESOLVED_SEPARATION = 7
    integer, parameter, public :: GC_PARTNER_NONFINITE = 8

    integer, parameter, public :: GC_PARTNER_BOUNDARY_USUAL = 3
    integer, parameter, public :: GC_PARTNER_BOUNDARY_X = 4

    type, public :: gc_operational_partner_options_t
        type(gc_operational_root_options_t) :: root
        real(dp) :: x_exclusion_fraction = 1.0e-8_dp
    end type gc_operational_partner_options_t

    type, public :: gc_operational_partner_pair_t
        integer :: pair_id = 0
        integer :: x_point_id = 0
        integer :: regular_crossing_count = 0
        real(dp) :: canonical_momentum = 0.0_dp
    end type gc_operational_partner_pair_t

    type, public :: gc_operational_separatrix_boundary_t
        integer :: boundary_id = 0
        integer :: kind = 0
        integer :: pair_id = 0
        integer :: x_point_id = 0
        real(dp) :: x = 0.0_dp
        real(dp) :: canonical_momentum = 0.0_dp
        real(dp) :: canonical_residual = 0.0_dp
    end type gc_operational_separatrix_boundary_t

    type, public :: gc_operational_partner_result_t
        integer :: status = GC_PARTNER_INVALID_INPUT
        integer :: npairs = 0
        integer :: nboundaries = 0
        logical :: complete = .false.
        type(gc_operational_partner_pair_t), allocatable :: pairs(:)
        type(gc_operational_separatrix_boundary_t), allocatable :: boundaries(:)
    end type gc_operational_partner_result_t

    public :: find_gc_operational_partner_crossings

contains

    subroutine find_gc_operational_partner_crossings(evaluate_canonical, &
            fixed_points, x_lo, x_hi, options, result)
        procedure(gc_operational_canonical_jet_i) :: evaluate_canonical
        type(gc_operational_fixed_point_result_t), intent(in) :: fixed_points
        real(dp), intent(in) :: x_lo, x_hi
        type(gc_operational_partner_options_t), intent(in) :: options
        type(gc_operational_partner_result_t), intent(out) :: result

        type(gc_operational_partner_pair_t), allocatable :: pairs(:)
        type(gc_operational_separatrix_boundary_t), allocatable :: boundaries(:)
        real(dp) :: exclusion, level_tolerance, separation_tolerance
        real(dp) :: target_x, target_level, active_target_level
        integer :: i, j, pair_count, boundary_count, regular_count
        integer :: actual_x_count
        integer :: local_status

        result = gc_operational_partner_result_t()
        allocate(result%pairs(0), result%boundaries(0))
        if (.not. all(ieee_is_finite([x_lo, x_hi, &
                options%x_exclusion_fraction]))) return
        if (x_hi <= x_lo) return
        if (options%x_exclusion_fraction <= 0.0_dp .or. &
                options%x_exclusion_fraction >= 0.5_dp) return
        if (.not. fixed_points%complete) return
        if (fixed_points%npoints < 0 .or. fixed_points%n_x_points < 0) return
        if (fixed_points%npoints > 0) then
            if (.not. allocated(fixed_points%points)) return
            if (size(fixed_points%points) /= fixed_points%npoints) return
        end if

        actual_x_count = 0
        do i = 1, fixed_points%npoints
            if (fixed_points%points(i)%kind == GC_FIXED_POINT_X) then
                actual_x_count = actual_x_count+1
            end if
        end do
        if (actual_x_count /= fixed_points%n_x_points) return
        pair_count = actual_x_count
        allocate(pairs(pair_count))
        allocate(boundaries(max(1, fixed_points%n_x_points)))
        boundary_count = 0
        pair_count = 0
        exclusion = max(options%x_exclusion_fraction*abs(x_hi-x_lo), &
            64.0_dp*epsilon(1.0_dp)*max(1.0_dp, abs(x_lo), abs(x_hi)))
        separation_tolerance = 8.0_dp*exclusion

        do i = 1, fixed_points%npoints
            if (fixed_points%points(i)%kind /= GC_FIXED_POINT_X) cycle
            target_x = fixed_points%points(i)%x
            target_level = fixed_points%points(i)%canonical_momentum
            if (.not. all(ieee_is_finite([target_x, target_level]))) then
                result%status = GC_PARTNER_NONFINITE
                return
            end if
            if (target_x <= x_lo+exclusion .or. &
                    target_x >= x_hi-exclusion) then
                result%status = GC_PARTNER_ENDPOINT_COINCIDENCE
                return
            end if
            do j = 1, i-1
                if (fixed_points%points(j)%kind /= GC_FIXED_POINT_X) cycle
                level_tolerance = 256.0_dp*epsilon(1.0_dp)*max(1.0_dp, &
                    abs(target_level), &
                    abs(fixed_points%points(j)%canonical_momentum))
                if (abs(target_level- &
                        fixed_points%points(j)%canonical_momentum) <= &
                        level_tolerance) then
                    result%status = GC_PARTNER_REPEATED_X_LEVEL
                    return
                end if
            end do

            pair_count = pair_count+1
            pairs(pair_count)%pair_id = pair_count
            pairs(pair_count)%x_point_id = fixed_points%points(i)%point_id
            pairs(pair_count)%canonical_momentum = target_level
            call append_boundary(GC_PARTNER_BOUNDARY_X, pair_count, &
                fixed_points%points(i)%point_id, target_x, target_level, &
                0.0_dp, boundaries, boundary_count, separation_tolerance, &
                local_status)
            if (local_status /= GC_PARTNER_SUCCESS) then
                result%status = local_status
                return
            end if

            regular_count = 0
            call search_segment(x_lo, target_x-exclusion, target_x, &
                target_level, pair_count, fixed_points%points(i)%point_id, &
                regular_count, local_status)
            if (local_status /= GC_PARTNER_SUCCESS) then
                result%status = local_status
                return
            end if
            call search_segment(target_x+exclusion, x_hi, target_x, &
                target_level, pair_count, fixed_points%points(i)%point_id, &
                regular_count, local_status)
            if (local_status /= GC_PARTNER_SUCCESS) then
                result%status = local_status
                return
            end if
            if (regular_count < 1) then
                result%status = GC_PARTNER_MISSING
                return
            end if
            pairs(pair_count)%regular_crossing_count = regular_count
        end do

        call sort_boundaries(boundaries, boundary_count)
        if (allocated(result%pairs)) deallocate(result%pairs)
        if (allocated(result%boundaries)) deallocate(result%boundaries)
        allocate(result%pairs(pair_count), result%boundaries(boundary_count))
        if (pair_count > 0) result%pairs = pairs(1:pair_count)
        if (boundary_count > 0) then
            result%boundaries = boundaries(1:boundary_count)
            do i = 1, boundary_count
                result%boundaries(i)%boundary_id = i
            end do
        end if
        result%npairs = pair_count
        result%nboundaries = boundary_count
        result%complete = .true.
        result%status = GC_PARTNER_SUCCESS

    contains

        subroutine search_segment(segment_lo, segment_hi, excluded_x, &
                target, pair_id, x_point_id, count, status)
            real(dp), intent(in) :: segment_lo, segment_hi, excluded_x, target
            integer, intent(in) :: pair_id, x_point_id
            integer, intent(inout) :: count
            integer, intent(out) :: status

            type(gc_operational_root_result_t) :: roots
            real(dp) :: p_star, first, second
            integer :: k, callback_status

            status = GC_PARTNER_SUCCESS
            if (segment_hi <= segment_lo) return
            active_target_level = target
            call find_gc_operational_scalar_roots(partner_residual, segment_lo, &
                segment_hi, options%root, roots)
            if (roots%status /= GC_OPERATIONAL_ROOT_SUCCESS .or. &
                    .not. roots%complete) then
                status = GC_PARTNER_ROOT_FAILURE
                return
            end if
            do k = 1, roots%nroots
                if (abs(roots%roots(k)%x-excluded_x) <= &
                        separation_tolerance) then
                    status = GC_PARTNER_UNRESOLVED_SEPARATION
                    return
                end if
                if (roots%roots(k)%x <= x_lo+separation_tolerance .or. &
                        roots%roots(k)%x >= x_hi-separation_tolerance) then
                    status = GC_PARTNER_ENDPOINT_COINCIDENCE
                    return
                end if
                call evaluate_canonical(roots%roots(k)%x, p_star, first, &
                    second, callback_status)
                if (callback_status /= 0) then
                    status = GC_PARTNER_CANONICAL_FAILURE
                    return
                end if
                if (.not. all(ieee_is_finite([p_star, first, second]))) then
                    status = GC_PARTNER_NONFINITE
                    return
                end if
                call append_boundary(GC_PARTNER_BOUNDARY_USUAL, pair_id, &
                    x_point_id, roots%roots(k)%x, target, p_star-target, &
                    boundaries, &
                    boundary_count, separation_tolerance, status)
                if (status /= GC_PARTNER_SUCCESS) return
                count = count+1
            end do
        end subroutine search_segment

        subroutine partner_residual(x, value, derivative, callback_status)
            real(dp), intent(in) :: x
            real(dp), intent(out) :: value, derivative
            integer, intent(out) :: callback_status
            real(dp) :: canonical_momentum, second_derivative

            call evaluate_canonical(x, canonical_momentum, derivative, &
                second_derivative, callback_status)
            if (callback_status /= 0) return
            value = canonical_momentum-active_target_level
        end subroutine partner_residual

    end subroutine find_gc_operational_partner_crossings

    subroutine append_boundary(kind, pair_id, x_point_id, x, p_star, residual, &
            boundaries, count, separation_tolerance, status)
        integer, intent(in) :: kind, pair_id, x_point_id
        real(dp), intent(in) :: x, p_star, residual, separation_tolerance
        type(gc_operational_separatrix_boundary_t), allocatable, &
            intent(inout) :: boundaries(:)
        integer, intent(inout) :: count
        integer, intent(out) :: status

        type(gc_operational_separatrix_boundary_t), allocatable :: enlarged(:)
        integer :: i, capacity

        status = GC_PARTNER_SUCCESS
        if (.not. all(ieee_is_finite([x, p_star, residual, &
                separation_tolerance]))) then
            status = GC_PARTNER_NONFINITE
            return
        end if
        do i = 1, count
            if (abs(boundaries(i)%x-x) <= separation_tolerance) then
                status = GC_PARTNER_UNRESOLVED_SEPARATION
                return
            end if
        end do
        if (count >= size(boundaries)) then
            capacity = max(2*size(boundaries), count+1)
            allocate(enlarged(capacity))
            if (count > 0) enlarged(1:count) = boundaries(1:count)
            call move_alloc(enlarged, boundaries)
        end if
        count = count+1
        boundaries(count)%kind = kind
        boundaries(count)%pair_id = pair_id
        boundaries(count)%x_point_id = x_point_id
        boundaries(count)%x = x
        boundaries(count)%canonical_momentum = p_star
        boundaries(count)%canonical_residual = residual
    end subroutine append_boundary

    subroutine sort_boundaries(boundaries, count)
        type(gc_operational_separatrix_boundary_t), intent(inout) :: &
            boundaries(:)
        integer, intent(in) :: count

        type(gc_operational_separatrix_boundary_t) :: temporary
        integer :: i, j, first

        do i = 1, count-1
            first = i
            do j = i+1, count
                if (boundaries(j)%x < boundaries(first)%x) first = j
            end do
            if (first /= i) then
                temporary = boundaries(i)
                boundaries(i) = boundaries(first)
                boundaries(first) = temporary
            end if
        end do
    end subroutine sort_boundaries

end module neort_gc_operational_partner_crossings
