module neort_gc_eqdsk_composite_r_ownership
    !! Deterministic ownership of a physical-R interval on a composite cut.
    !!
    !! The composite cut has three computational source atlases.  Their
    !! physical ownership is fixed to
    !!
    !!     [R_left,R_axis_lo), [R_axis_lo,R_axis_hi), [R_axis_hi,R_right].
    !!
    !! Thus an exact first seam belongs to the axis branch and an exact second
    !! seam belongs to the outboard branch.  The source graph certificate and
    !! source interval are copied into every emitted piece; they are not
    !! replaced by the ownership partition's certificate.
    !!
    !! This module performs no field evaluation and makes no topology claim.
    !! It only validates the geometric partition and decomposes closed query
    !! intervals.  A query ending on a seam emits a zero-width point piece for
    !! the branch that owns that seam, so closed-interval coverage is explicit.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_eqdsk_composite_cut_atlas, only: &
        EQDSK_COMPOSITE_ATLAS_CERTIFICATE_ID, &
        EQDSK_COMPOSITE_ATLAS_SUCCESS, &
        eqdsk_composite_cut_atlas_t, &
        get_eqdsk_composite_cut_radius_bounds, &
        validate_eqdsk_composite_cut_atlas
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID
    implicit none
    private

    integer, parameter, public :: EQDSK_R_OWNERSHIP_SUCCESS = 0
    integer, parameter, public :: EQDSK_R_OWNERSHIP_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_R_OWNERSHIP_INVALID_ATLAS = 2
    integer, parameter, public :: EQDSK_R_OWNERSHIP_INVALID_PARTITION = 3
    integer, parameter, public :: EQDSK_R_OWNERSHIP_OUT_OF_RANGE = 4
    integer, parameter, public :: EQDSK_R_OWNERSHIP_INVALID_QUERY = 5
    integer, parameter, public :: EQDSK_R_OWNERSHIP_INCONSISTENT = 6

    integer, parameter, public :: EQDSK_R_BRANCH_INBOARD = 1
    integer, parameter, public :: EQDSK_R_BRANCH_AXIS = 2
    integer, parameter, public :: EQDSK_R_BRANCH_OUTBOARD = 3
    integer, parameter, public :: EQDSK_R_BRANCH_COUNT = 3

    ! These IDs are runtime theorem IDs, not generated registry indices.
    integer, parameter, public :: EQDSK_R_OWNERSHIP_CERTIFICATE_ID = 130018
    integer, parameter, public :: EQDSK_R_SEAM_INBOARD_AXIS_ID = 130019
    integer, parameter, public :: EQDSK_R_SEAM_AXIS_OUTBOARD_ID = 130020
    integer, parameter, public :: EQDSK_R_LEFT_ENDPOINT_ID = 130021
    integer, parameter, public :: EQDSK_R_RIGHT_ENDPOINT_ID = 130022

    type, public :: eqdsk_composite_r_seam_t
        real(dp) :: radius = 0.0_dp
        integer :: seam_id = 0
        integer :: left_branch_id = 0
        integer :: right_branch_id = 0
    end type eqdsk_composite_r_seam_t

    type, public :: eqdsk_composite_r_piece_t
        integer :: ordinal = 0
        integer :: source_branch_id = 0
        integer :: source_graph_certificate_id = 0
        real(dp) :: source_r_lo = 0.0_dp
        real(dp) :: source_r_hi = 0.0_dp
        real(dp) :: r_lo = 0.0_dp
        real(dp) :: r_hi = 0.0_dp
        logical :: left_closed = .false.
        logical :: right_closed = .false.
        logical :: point_only = .false.
        integer :: left_boundary_id = 0
        integer :: right_boundary_id = 0
    end type eqdsk_composite_r_piece_t

    type, public :: eqdsk_composite_r_partition_t
        integer :: certificate_id = 0
        integer :: branch_count = 0
        integer :: seam_count = 0
        real(dp) :: physical_r_lo = 0.0_dp
        real(dp) :: physical_r_hi = 0.0_dp
        type(eqdsk_composite_r_seam_t) :: seams(2)
        type(eqdsk_composite_r_piece_t) :: branches(3)
    end type eqdsk_composite_r_partition_t

    type, public :: eqdsk_composite_r_interval_t
        integer :: certificate_id = 0
        integer :: piece_count = 0
        real(dp) :: query_r_lo = 0.0_dp
        real(dp) :: query_r_hi = 0.0_dp
        type(eqdsk_composite_r_piece_t) :: pieces(3)
    end type eqdsk_composite_r_interval_t

    public :: build_eqdsk_composite_r_partition
    public :: build_eqdsk_composite_r_partition_from_bounds
    public :: decompose_eqdsk_composite_r_interval
    public :: locate_eqdsk_composite_r_point
    public :: validate_eqdsk_composite_r_interval
    public :: validate_eqdsk_composite_r_partition

contains

    subroutine build_eqdsk_composite_r_partition(atlas, partition, status)
        type(eqdsk_composite_cut_atlas_t), intent(inout) :: atlas
        type(eqdsk_composite_r_partition_t), intent(out) :: partition
        integer, intent(out) :: status

        real(dp) :: physical_bounds(4), source_bounds(2, 3)
        integer :: source_certificate_ids(3), local_status

        partition = eqdsk_composite_r_partition_t()
        status = EQDSK_R_OWNERSHIP_INVALID_ATLAS
        call validate_eqdsk_composite_cut_atlas(atlas, local_status)
        if (local_status /= EQDSK_COMPOSITE_ATLAS_SUCCESS) return
        call get_eqdsk_composite_cut_radius_bounds(atlas, physical_bounds(1), &
            physical_bounds(4), local_status)
        if (local_status /= EQDSK_COMPOSITE_ATLAS_SUCCESS) return
        physical_bounds(2) = atlas%axis%r_lo
        physical_bounds(3) = atlas%axis%r_hi

        source_bounds(:, 1) = [atlas%inboard_graph%requested_r_lo, &
            atlas%inboard_graph%requested_r_hi]
        source_bounds(:, 2) = [atlas%axis_graph%requested_r_lo, &
            atlas%axis_graph%requested_r_hi]
        source_bounds(:, 3) = [atlas%outboard_graph%requested_r_lo, &
            atlas%outboard_graph%requested_r_hi]
        source_certificate_ids = [atlas%inboard_graph%certificate_id, &
            atlas%axis_graph%certificate_id, &
            atlas%outboard_graph%certificate_id]

        call build_eqdsk_composite_r_partition_from_bounds(physical_bounds, &
            source_bounds, source_certificate_ids, partition, status)
    end subroutine build_eqdsk_composite_r_partition

    subroutine build_eqdsk_composite_r_partition_from_bounds(physical_bounds, &
            source_bounds, source_certificate_ids, partition, status)
        !! Manufactured/adapter-facing constructor for the ownership layer.
        !! `physical_bounds` are [R_left,R_seam1,R_seam2,R_right].
        !! `source_bounds(:,branch)` retain the source graph domain.
        real(dp), intent(in) :: physical_bounds(4), source_bounds(2, 3)
        integer, intent(in) :: source_certificate_ids(3)
        type(eqdsk_composite_r_partition_t), intent(out) :: partition
        integer, intent(out) :: status

        integer :: branch

        partition = eqdsk_composite_r_partition_t()
        status = EQDSK_R_OWNERSHIP_INVALID_INPUT
        if (.not. all(ieee_is_finite(physical_bounds))) return
        if (.not. all(ieee_is_finite(source_bounds))) return
        if (any(source_certificate_ids /= EQDSK_CUT_GRAPH_CERTIFICATE_ID)) return
        if (physical_bounds(1) <= 0.0_dp .or. &
                .not. strictly_increasing(physical_bounds)) return

        partition%certificate_id = EQDSK_R_OWNERSHIP_CERTIFICATE_ID
        partition%branch_count = EQDSK_R_BRANCH_COUNT
        partition%seam_count = 2
        partition%physical_r_lo = physical_bounds(1)
        partition%physical_r_hi = physical_bounds(4)

        partition%seams(1)%radius = physical_bounds(2)
        partition%seams(1)%seam_id = EQDSK_R_SEAM_INBOARD_AXIS_ID
        partition%seams(1)%left_branch_id = EQDSK_R_BRANCH_INBOARD
        partition%seams(1)%right_branch_id = EQDSK_R_BRANCH_AXIS
        partition%seams(2)%radius = physical_bounds(3)
        partition%seams(2)%seam_id = EQDSK_R_SEAM_AXIS_OUTBOARD_ID
        partition%seams(2)%left_branch_id = EQDSK_R_BRANCH_AXIS
        partition%seams(2)%right_branch_id = EQDSK_R_BRANCH_OUTBOARD

        do branch = 1, EQDSK_R_BRANCH_COUNT
            partition%branches(branch)%ordinal = branch
            partition%branches(branch)%source_branch_id = branch
            partition%branches(branch)%source_graph_certificate_id = &
                source_certificate_ids(branch)
            partition%branches(branch)%source_r_lo = source_bounds(1, branch)
            partition%branches(branch)%source_r_hi = source_bounds(2, branch)
        end do

        partition%branches(1)%r_lo = physical_bounds(1)
        partition%branches(1)%r_hi = physical_bounds(2)
        partition%branches(1)%left_closed = .true.
        partition%branches(1)%right_closed = .false.
        partition%branches(1)%left_boundary_id = EQDSK_R_LEFT_ENDPOINT_ID
        partition%branches(1)%right_boundary_id = &
            EQDSK_R_SEAM_INBOARD_AXIS_ID

        partition%branches(2)%r_lo = physical_bounds(2)
        partition%branches(2)%r_hi = physical_bounds(3)
        partition%branches(2)%left_closed = .true.
        partition%branches(2)%right_closed = .false.
        partition%branches(2)%left_boundary_id = &
            EQDSK_R_SEAM_INBOARD_AXIS_ID
        partition%branches(2)%right_boundary_id = &
            EQDSK_R_SEAM_AXIS_OUTBOARD_ID

        partition%branches(3)%r_lo = physical_bounds(3)
        partition%branches(3)%r_hi = physical_bounds(4)
        partition%branches(3)%left_closed = .true.
        partition%branches(3)%right_closed = .true.
        partition%branches(3)%left_boundary_id = &
            EQDSK_R_SEAM_AXIS_OUTBOARD_ID
        partition%branches(3)%right_boundary_id = &
            EQDSK_R_RIGHT_ENDPOINT_ID

        call validate_eqdsk_composite_r_partition(partition, status)
    end subroutine build_eqdsk_composite_r_partition_from_bounds

    subroutine validate_eqdsk_composite_r_partition(partition, status)
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        integer, intent(out) :: status

        integer :: branch

        status = EQDSK_R_OWNERSHIP_INVALID_PARTITION
        if (partition%certificate_id /= &
                EQDSK_R_OWNERSHIP_CERTIFICATE_ID) return
        if (partition%branch_count /= EQDSK_R_BRANCH_COUNT .or. &
                partition%seam_count /= 2) return
        if (.not. all(ieee_is_finite([partition%physical_r_lo, &
                partition%physical_r_hi, partition%seams%radius]))) return
        if (partition%physical_r_lo <= 0.0_dp .or. &
                partition%physical_r_lo >= partition%seams(1)%radius .or. &
                partition%seams(1)%radius >= partition%seams(2)%radius .or. &
                partition%seams(2)%radius >= partition%physical_r_hi) return

        if (partition%seams(1)%seam_id /= EQDSK_R_SEAM_INBOARD_AXIS_ID .or. &
                partition%seams(2)%seam_id /= &
                EQDSK_R_SEAM_AXIS_OUTBOARD_ID) return
        if (partition%seams(1)%left_branch_id /= &
                EQDSK_R_BRANCH_INBOARD .or. &
                partition%seams(1)%right_branch_id /= &
                EQDSK_R_BRANCH_AXIS) return
        if (partition%seams(2)%left_branch_id /= &
                EQDSK_R_BRANCH_AXIS .or. &
                partition%seams(2)%right_branch_id /= &
                EQDSK_R_BRANCH_OUTBOARD) return

        do branch = 1, EQDSK_R_BRANCH_COUNT
            if (.not. valid_partition_piece(partition, partition%branches(branch), &
                    branch)) return
        end do
        if (partition%branches(1)%r_lo /= partition%physical_r_lo .or. &
                partition%branches(1)%r_hi /= &
                partition%branches(2)%r_lo .or. &
                partition%branches(2)%r_hi /= &
                partition%branches(3)%r_lo .or. &
                partition%branches(3)%r_hi /= partition%physical_r_hi) return
        if (partition%branches(1)%right_boundary_id /= &
                partition%branches(2)%left_boundary_id .or. &
                partition%branches(2)%right_boundary_id /= &
                partition%branches(3)%left_boundary_id) return
        if (partition%branches(1)%source_r_hi /= &
                partition%branches(2)%source_r_lo .or. &
                partition%branches(2)%source_r_hi /= &
                partition%branches(3)%source_r_lo) return
        if (partition%branches(1)%source_r_lo > &
                partition%branches(1)%r_lo .or. &
                partition%branches(1)%source_r_hi < &
                partition%branches(1)%r_hi) return
        if (partition%branches(2)%source_r_lo > &
                partition%branches(2)%r_lo .or. &
                partition%branches(2)%source_r_hi < &
                partition%branches(2)%r_hi) return
        if (partition%branches(3)%source_r_lo > &
                partition%branches(3)%r_lo .or. &
                partition%branches(3)%source_r_hi < &
                partition%branches(3)%r_hi) return
        status = EQDSK_R_OWNERSHIP_SUCCESS
    end subroutine validate_eqdsk_composite_r_partition

    subroutine decompose_eqdsk_composite_r_interval(partition, r_lo, r_hi, &
            decomposition, status)
        !! Decompose the closed query [r_lo,r_hi] into owned source pieces.
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        real(dp), intent(in) :: r_lo, r_hi
        type(eqdsk_composite_r_interval_t), intent(out) :: decomposition
        integer, intent(out) :: status

        integer :: branch, count, local_status, seam_id, source_id
        real(dp) :: piece_lo, piece_hi

        decomposition = eqdsk_composite_r_interval_t()
        status = EQDSK_R_OWNERSHIP_INVALID_PARTITION
        call validate_eqdsk_composite_r_partition(partition, local_status)
        if (local_status /= EQDSK_R_OWNERSHIP_SUCCESS) return
        status = EQDSK_R_OWNERSHIP_INVALID_QUERY
        if (.not. all(ieee_is_finite([r_lo, r_hi]))) return
        if (r_lo > r_hi) return
        if (r_lo < partition%physical_r_lo .or. &
                r_hi > partition%physical_r_hi) then
            status = EQDSK_R_OWNERSHIP_OUT_OF_RANGE
            return
        end if

        decomposition%certificate_id = partition%certificate_id
        decomposition%query_r_lo = r_lo
        decomposition%query_r_hi = r_hi

        if (r_lo == r_hi) then
            call locate_eqdsk_composite_r_point(partition, r_lo, branch, &
                seam_id, source_id, status)
            if (status /= EQDSK_R_OWNERSHIP_SUCCESS) return
            decomposition%piece_count = 1
            call make_point_piece(partition, branch, r_lo, &
                decomposition%pieces(1))
            call validate_eqdsk_composite_r_interval(partition, decomposition, &
                status)
            return
        end if

        count = 0
        do branch = 1, EQDSK_R_BRANCH_COUNT
            piece_lo = max(r_lo, partition%branches(branch)%r_lo)
            piece_hi = min(r_hi, partition%branches(branch)%r_hi)
            if (piece_hi > piece_lo) then
                count = count+1
                if (count > EQDSK_R_BRANCH_COUNT) then
                    decomposition = eqdsk_composite_r_interval_t()
                    status = EQDSK_R_OWNERSHIP_INCONSISTENT
                    return
                end if
                decomposition%pieces(count) = &
                    partition%branches(branch)
                decomposition%pieces(count)%r_lo = piece_lo
                decomposition%pieces(count)%r_hi = piece_hi
                decomposition%pieces(count)%point_only = .false.
                decomposition%pieces(count)%left_closed = &
                    (piece_lo > partition%branches(branch)%r_lo) .or. &
                    partition%branches(branch)%left_closed
                decomposition%pieces(count)%right_closed = &
                    (piece_hi < partition%branches(branch)%r_hi) .or. &
                    partition%branches(branch)%right_closed
                if (piece_lo == partition%branches(branch)%r_lo) then
                    decomposition%pieces(count)%left_boundary_id = &
                        partition%branches(branch)%left_boundary_id
                else
                    decomposition%pieces(count)%left_boundary_id = 0
                end if
                if (piece_hi == partition%branches(branch)%r_hi) then
                    decomposition%pieces(count)%right_boundary_id = &
                        partition%branches(branch)%right_boundary_id
                else
                    decomposition%pieces(count)%right_boundary_id = 0
                end if
            else if (piece_hi == piece_lo) then
                if (piece_hi == r_hi) then
                    if (partition%branches(branch)%left_closed) then
                        if (piece_hi == partition%branches(branch)%r_lo) then
                            count = count+1
                            if (count > EQDSK_R_BRANCH_COUNT) then
                                decomposition = eqdsk_composite_r_interval_t()
                                status = EQDSK_R_OWNERSHIP_INCONSISTENT
                                return
                            end if
                            call make_point_piece(partition, branch, piece_hi, &
                                decomposition%pieces(count))
                        end if
                    end if
                end if
            end if
        end do
        decomposition%piece_count = count
        if (count == 0) then
            decomposition = eqdsk_composite_r_interval_t()
            status = EQDSK_R_OWNERSHIP_INCONSISTENT
            return
        end if
        call validate_eqdsk_composite_r_interval(partition, decomposition, status)
    end subroutine decompose_eqdsk_composite_r_interval

    subroutine validate_eqdsk_composite_r_interval(partition, decomposition, &
            status)
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        type(eqdsk_composite_r_interval_t), intent(in) :: decomposition
        integer, intent(out) :: status

        integer :: i, local_status

        status = EQDSK_R_OWNERSHIP_INVALID_PARTITION
        call validate_eqdsk_composite_r_partition(partition, local_status)
        if (local_status /= EQDSK_R_OWNERSHIP_SUCCESS) return
        if (decomposition%certificate_id /= &
                EQDSK_R_OWNERSHIP_CERTIFICATE_ID) return
        if (decomposition%piece_count < 1 .or. &
                decomposition%piece_count > EQDSK_R_BRANCH_COUNT) return
        if (.not. all(ieee_is_finite([decomposition%query_r_lo, &
                decomposition%query_r_hi]))) return
        if (decomposition%query_r_lo > decomposition%query_r_hi) return
        if (decomposition%query_r_lo < partition%physical_r_lo .or. &
                decomposition%query_r_hi > partition%physical_r_hi) then
            status = EQDSK_R_OWNERSHIP_OUT_OF_RANGE
            return
        end if

        if (decomposition%query_r_lo == decomposition%query_r_hi) then
            if (decomposition%piece_count /= 1) return
            if (.not. decomposition%pieces(1)%point_only) return
            if (decomposition%pieces(1)%r_lo /= &
                    decomposition%query_r_lo .or. &
                    decomposition%pieces(1)%r_hi /= &
                    decomposition%query_r_hi) return
            if (.not. valid_interval_piece(partition, &
                    decomposition%pieces(1))) return
            call expected_point_branch(partition, decomposition%query_r_lo, &
                i, local_status)
            if (local_status /= EQDSK_R_OWNERSHIP_SUCCESS) return
            if (decomposition%pieces(1)%source_branch_id /= i) return
            status = EQDSK_R_OWNERSHIP_SUCCESS
            return
        end if

        if (decomposition%pieces(1)%r_lo /= decomposition%query_r_lo .or. &
                .not. decomposition%pieces(1)%left_closed) return
        if (decomposition%pieces(decomposition%piece_count)%r_hi /= &
                decomposition%query_r_hi .or. &
                .not. decomposition%pieces(decomposition%piece_count)%right_closed) return
        do i = 1, decomposition%piece_count
            if (.not. valid_interval_piece(partition, &
                    decomposition%pieces(i))) return
        end do
        do i = 2, decomposition%piece_count
            if (decomposition%pieces(i)%source_branch_id <= &
                    decomposition%pieces(i-1)%source_branch_id) return
            if (decomposition%pieces(i)%r_lo /= &
                    decomposition%pieces(i-1)%r_hi) return
            if (decomposition%pieces(i-1)%point_only) return
            if (decomposition%pieces(i)%point_only) then
                if (.not. decomposition%pieces(i)%left_closed) return
                if (decomposition%pieces(i-1)%right_closed) return
            else
                if (decomposition%pieces(i-1)%right_closed .eqv. &
                        decomposition%pieces(i)%left_closed) return
            end if
        end do
        status = EQDSK_R_OWNERSHIP_SUCCESS
    end subroutine validate_eqdsk_composite_r_interval

    subroutine locate_eqdsk_composite_r_point(partition, radius, branch_id, &
            seam_id, source_graph_certificate_id, status)
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        real(dp), intent(in) :: radius
        integer, intent(out) :: branch_id, seam_id
        integer, intent(out) :: source_graph_certificate_id
        integer, intent(out) :: status

        integer :: local_status

        branch_id = 0
        seam_id = 0
        source_graph_certificate_id = 0
        call validate_eqdsk_composite_r_partition(partition, local_status)
        if (local_status /= EQDSK_R_OWNERSHIP_SUCCESS) then
            status = local_status
            return
        end if
        if (.not. ieee_is_finite(radius)) then
            status = EQDSK_R_OWNERSHIP_INVALID_INPUT
            return
        end if
        if (radius < partition%physical_r_lo .or. &
                radius > partition%physical_r_hi) then
            status = EQDSK_R_OWNERSHIP_OUT_OF_RANGE
            return
        end if

        if (radius < partition%seams(1)%radius) then
            branch_id = EQDSK_R_BRANCH_INBOARD
        else if (radius < partition%seams(2)%radius) then
            branch_id = EQDSK_R_BRANCH_AXIS
        else
            branch_id = EQDSK_R_BRANCH_OUTBOARD
        end if
        source_graph_certificate_id = &
            partition%branches(branch_id)%source_graph_certificate_id
        if (radius == partition%physical_r_lo) then
            seam_id = EQDSK_R_LEFT_ENDPOINT_ID
        else if (radius == partition%seams(1)%radius) then
            seam_id = EQDSK_R_SEAM_INBOARD_AXIS_ID
        else if (radius == partition%seams(2)%radius) then
            seam_id = EQDSK_R_SEAM_AXIS_OUTBOARD_ID
        else if (radius == partition%physical_r_hi) then
            seam_id = EQDSK_R_RIGHT_ENDPOINT_ID
        end if
        status = EQDSK_R_OWNERSHIP_SUCCESS
    end subroutine locate_eqdsk_composite_r_point

    pure logical function strictly_increasing(values)
        real(dp), intent(in) :: values(:)

        integer :: i

        strictly_increasing = .true.
        do i = 2, size(values)
            if (values(i) <= values(i-1)) then
                strictly_increasing = .false.
                return
            end if
        end do
    end function strictly_increasing

    logical function valid_partition_piece(partition, piece, branch)
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        type(eqdsk_composite_r_piece_t), intent(in) :: piece
        integer, intent(in) :: branch

        integer :: expected_left_id, expected_right_id

        valid_partition_piece = .false.
        if (branch < EQDSK_R_BRANCH_INBOARD) return
        if (branch > EQDSK_R_BRANCH_OUTBOARD) return
        if (piece%ordinal /= branch .or. piece%source_branch_id /= branch) return
        if (piece%source_graph_certificate_id /= &
                EQDSK_CUT_GRAPH_CERTIFICATE_ID) return
        if (.not. all(ieee_is_finite([piece%source_r_lo, piece%source_r_hi, &
                piece%r_lo, piece%r_hi]))) return
        if (piece%source_r_lo >= piece%source_r_hi .or. &
                piece%r_lo >= piece%r_hi) return
        if (piece%point_only .or. .not. piece%left_closed) return
        if (branch < EQDSK_R_BRANCH_OUTBOARD) then
            if (piece%right_closed) return
        else
            if (.not. piece%right_closed) return
        end if
        if (piece%r_lo /= partition%branches(branch)%r_lo .or. &
                piece%r_hi /= partition%branches(branch)%r_hi) return
        if (piece%source_r_lo /= partition%branches(branch)%source_r_lo .or. &
                piece%source_r_hi /= partition%branches(branch)%source_r_hi) return
        call canonical_partition_boundary_ids(branch, expected_left_id, &
            expected_right_id)
        if (piece%left_boundary_id /= expected_left_id .or. &
                piece%right_boundary_id /= expected_right_id) return
        valid_partition_piece = .true.
    end function valid_partition_piece

    logical function valid_interval_piece(partition, piece)
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        type(eqdsk_composite_r_piece_t), intent(in) :: piece

        integer :: branch, expected_id

        valid_interval_piece = .false.
        if (piece%ordinal /= piece%source_branch_id) return
        branch = piece%source_branch_id
        if (branch < EQDSK_R_BRANCH_INBOARD) return
        if (branch > EQDSK_R_BRANCH_OUTBOARD) return
        if (piece%source_graph_certificate_id /= &
                EQDSK_CUT_GRAPH_CERTIFICATE_ID) return
        if (.not. all(ieee_is_finite([piece%source_r_lo, piece%source_r_hi, &
                piece%r_lo, piece%r_hi]))) return
        if (piece%source_r_lo >= piece%source_r_hi) return
        if (piece%source_r_lo /= partition%branches(branch)%source_r_lo .or. &
                piece%source_r_hi /= partition%branches(branch)%source_r_hi) return
        if (piece%r_lo < partition%branches(branch)%r_lo) return
        if (piece%r_hi > partition%branches(branch)%r_hi) return
        if (piece%point_only) then
            if (piece%r_lo /= piece%r_hi) return
            if (.not. piece%left_closed) return
            if (.not. piece%right_closed) return
            if (piece%r_lo == partition%branches(branch)%r_hi) then
                if (branch < EQDSK_R_BRANCH_OUTBOARD) return
            end if
            call canonical_point_boundary_id(partition, branch, &
                piece%r_lo, expected_id)
            if (piece%left_boundary_id /= expected_id) return
            if (piece%right_boundary_id /= expected_id) return
        else
            if (piece%r_lo >= piece%r_hi) return
            if (.not. piece%left_closed) return
            if (piece%r_hi == partition%branches(branch)%r_hi) then
                if (branch < EQDSK_R_BRANCH_OUTBOARD) then
                    if (piece%right_closed) return
                else
                    if (.not. piece%right_closed) return
                end if
            else
                if (.not. piece%right_closed) return
            end if
            call expected_interval_boundary_id(partition, branch, &
                piece%r_lo, .true., expected_id)
            if (piece%left_boundary_id /= expected_id) return
            call expected_interval_boundary_id(partition, branch, &
                piece%r_hi, .false., expected_id)
            if (piece%right_boundary_id /= expected_id) return
        end if
        valid_interval_piece = .true.
    end function valid_interval_piece

    subroutine make_point_piece(partition, branch, radius, piece)
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        integer, intent(in) :: branch
        real(dp), intent(in) :: radius
        type(eqdsk_composite_r_piece_t), intent(out) :: piece

        integer :: boundary_id

        piece = partition%branches(branch)
        piece%r_lo = radius
        piece%r_hi = radius
        piece%point_only = .true.
        piece%left_closed = .true.
        piece%right_closed = .true.
        call canonical_point_boundary_id(partition, branch, radius, boundary_id)
        piece%left_boundary_id = boundary_id
        piece%right_boundary_id = boundary_id
    end subroutine make_point_piece

    subroutine expected_point_branch(partition, radius, branch, status)
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        real(dp), intent(in) :: radius
        integer, intent(out) :: branch, status

        integer :: seam_id, source_id

        call locate_eqdsk_composite_r_point(partition, radius, branch, seam_id, &
            source_id, status)
    end subroutine expected_point_branch

    subroutine canonical_partition_boundary_ids(branch, left_id, right_id)
        integer, intent(in) :: branch
        integer, intent(out) :: left_id, right_id

        left_id = 0
        right_id = 0
        select case (branch)
        case (EQDSK_R_BRANCH_INBOARD)
            left_id = EQDSK_R_LEFT_ENDPOINT_ID
            right_id = EQDSK_R_SEAM_INBOARD_AXIS_ID
        case (EQDSK_R_BRANCH_AXIS)
            left_id = EQDSK_R_SEAM_INBOARD_AXIS_ID
            right_id = EQDSK_R_SEAM_AXIS_OUTBOARD_ID
        case (EQDSK_R_BRANCH_OUTBOARD)
            left_id = EQDSK_R_SEAM_AXIS_OUTBOARD_ID
            right_id = EQDSK_R_RIGHT_ENDPOINT_ID
        end select
    end subroutine canonical_partition_boundary_ids

    subroutine canonical_point_boundary_id(partition, branch, radius, id)
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        integer, intent(in) :: branch
        real(dp), intent(in) :: radius
        integer, intent(out) :: id

        id = 0
        if (radius == partition%branches(branch)%r_lo) then
            id = partition%branches(branch)%left_boundary_id
        else if (radius == partition%branches(branch)%r_hi .and. &
                branch == EQDSK_R_BRANCH_OUTBOARD) then
            id = partition%branches(branch)%right_boundary_id
        end if
    end subroutine canonical_point_boundary_id

    subroutine expected_interval_boundary_id(partition, branch, radius, is_left, &
            id)
        type(eqdsk_composite_r_partition_t), intent(in) :: partition
        integer, intent(in) :: branch
        real(dp), intent(in) :: radius
        logical, intent(in) :: is_left
        integer, intent(out) :: id

        id = 0
        if (is_left) then
            if (radius == partition%branches(branch)%r_lo) then
                id = partition%branches(branch)%left_boundary_id
            end if
        else
            if (radius == partition%branches(branch)%r_hi) then
                id = partition%branches(branch)%right_boundary_id
            end if
        end if
    end subroutine expected_interval_boundary_id

end module neort_gc_eqdsk_composite_r_ownership
