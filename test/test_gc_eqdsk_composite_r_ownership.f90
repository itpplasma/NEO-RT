program test_gc_eqdsk_composite_r_ownership
    use, intrinsic :: ieee_arithmetic, only: ieee_quiet_nan, ieee_value
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_eqdsk_composite_r_ownership, only: &
        EQDSK_R_BRANCH_AXIS, EQDSK_R_BRANCH_INBOARD, &
        EQDSK_R_BRANCH_OUTBOARD, EQDSK_R_LEFT_ENDPOINT_ID, &
        EQDSK_R_OWNERSHIP_CERTIFICATE_ID, EQDSK_R_OWNERSHIP_INVALID_INPUT, &
        EQDSK_R_OWNERSHIP_INVALID_PARTITION, EQDSK_R_OWNERSHIP_OUT_OF_RANGE, &
        EQDSK_R_OWNERSHIP_SUCCESS, EQDSK_R_RIGHT_ENDPOINT_ID, &
        EQDSK_R_SEAM_AXIS_OUTBOARD_ID, EQDSK_R_SEAM_INBOARD_AXIS_ID, &
        build_eqdsk_composite_r_partition_from_bounds, &
        decompose_eqdsk_composite_r_interval, &
        eqdsk_composite_r_interval_t, eqdsk_composite_r_partition_t, &
        locate_eqdsk_composite_r_point, &
        validate_eqdsk_composite_r_interval, &
        validate_eqdsk_composite_r_partition
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID
    implicit none

    type(eqdsk_composite_r_partition_t) :: partition, damaged
    type(eqdsk_composite_r_interval_t) :: decomposition, damaged_decomposition
    real(dp) :: physical_bounds(4), source_bounds(2, 3)
    integer :: source_certificate_ids(3), status, branch, seam, source_id

    physical_bounds = [1.0_dp, 4.0_dp, 7.0_dp, 10.0_dp]
    source_bounds(:, 1) = [0.5_dp, 4.0_dp]
    source_bounds(:, 2) = [4.0_dp, 7.0_dp]
    source_bounds(:, 3) = [7.0_dp, 10.5_dp]
    source_certificate_ids = [EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID, EQDSK_CUT_GRAPH_CERTIFICATE_ID]

    call build_eqdsk_composite_r_partition_from_bounds(physical_bounds, &
        source_bounds, source_certificate_ids, partition, status)
    call require(status == EQDSK_R_OWNERSHIP_SUCCESS, &
        'manufactured partition did not validate')
    call require(partition%certificate_id == &
        EQDSK_R_OWNERSHIP_CERTIFICATE_ID, 'partition theorem ID changed')
    call require(partition%branches(1)%source_graph_certificate_id == &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID .and. &
        partition%branches(2)%source_graph_certificate_id == &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID .and. &
        partition%branches(3)%source_graph_certificate_id == &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
        'source graph provenance was not retained')
    call require(partition%branches(1)%left_closed .and. &
        .not. partition%branches(1)%right_closed .and. &
        partition%branches(2)%left_closed .and. &
        .not. partition%branches(2)%right_closed .and. &
        partition%branches(3)%left_closed .and. &
        partition%branches(3)%right_closed, &
        'half-open branch ownership is wrong')
    call require(partition%seams(1)%seam_id == EQDSK_R_SEAM_INBOARD_AXIS_ID &
        .and. partition%seams(2)%seam_id == EQDSK_R_SEAM_AXIS_OUTBOARD_ID, &
        'seam IDs are not stable')

    call locate_and_check(1.0_dp, EQDSK_R_BRANCH_INBOARD, &
        EQDSK_R_LEFT_ENDPOINT_ID, EQDSK_CUT_GRAPH_CERTIFICATE_ID)
    call locate_and_check(4.0_dp, EQDSK_R_BRANCH_AXIS, &
        EQDSK_R_SEAM_INBOARD_AXIS_ID, EQDSK_CUT_GRAPH_CERTIFICATE_ID)
    call locate_and_check(5.0_dp, EQDSK_R_BRANCH_AXIS, 0, &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID)
    call locate_and_check(7.0_dp, EQDSK_R_BRANCH_OUTBOARD, &
        EQDSK_R_SEAM_AXIS_OUTBOARD_ID, EQDSK_CUT_GRAPH_CERTIFICATE_ID)
    call locate_and_check(10.0_dp, EQDSK_R_BRANCH_OUTBOARD, &
        EQDSK_R_RIGHT_ENDPOINT_ID, EQDSK_CUT_GRAPH_CERTIFICATE_ID)

    call decompose_eqdsk_composite_r_interval(partition, 1.0_dp, 10.0_dp, &
        decomposition, status)
    call require(status == EQDSK_R_OWNERSHIP_SUCCESS, &
        'full physical interval did not decompose')
    call require(decomposition%piece_count == 3, &
        'full physical interval does not have three pieces')
    call require(decomposition%pieces(1)%r_lo == 1.0_dp .and. &
        decomposition%pieces(1)%r_hi == 4.0_dp .and. &
        decomposition%pieces(2)%r_lo == 4.0_dp .and. &
        decomposition%pieces(2)%r_hi == 7.0_dp .and. &
        decomposition%pieces(3)%r_lo == 7.0_dp .and. &
        decomposition%pieces(3)%r_hi == 10.0_dp, &
        'full physical interval boundaries are wrong')
    call require(decomposition%pieces(1)%source_graph_certificate_id == &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID .and. &
        decomposition%pieces(2)%source_graph_certificate_id == &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID .and. &
        decomposition%pieces(3)%source_graph_certificate_id == &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
        'decomposition lost branch provenance')

    call decompose_eqdsk_composite_r_interval(partition, 1.0_dp, 4.0_dp, &
        decomposition, status)
    call require(status == EQDSK_R_OWNERSHIP_SUCCESS .and. &
        decomposition%piece_count == 2, &
        'closed interval ending at seam was not represented')
    call require(decomposition%pieces(1)%source_branch_id == &
        EQDSK_R_BRANCH_INBOARD .and. &
        decomposition%pieces(2)%source_branch_id == EQDSK_R_BRANCH_AXIS .and. &
        decomposition%pieces(2)%point_only .and. &
        decomposition%pieces(2)%r_lo == 4.0_dp .and. &
        decomposition%pieces(2)%left_boundary_id == &
        EQDSK_R_SEAM_INBOARD_AXIS_ID, &
        'first seam point was assigned to the wrong branch')

    call decompose_eqdsk_composite_r_interval(partition, 4.5_dp, 8.0_dp, &
        decomposition, status)
    call require(status == EQDSK_R_OWNERSHIP_SUCCESS .and. &
        decomposition%piece_count == 2, 'cross-seam query was not split')
    call require(decomposition%pieces(1)%r_lo == 4.5_dp .and. &
        decomposition%pieces(1)%r_hi == 7.0_dp .and. &
        .not. decomposition%pieces(1)%right_closed .and. &
        decomposition%pieces(2)%r_lo == 7.0_dp .and. &
        decomposition%pieces(2)%r_hi == 8.0_dp .and. &
        decomposition%pieces(2)%left_closed, &
        'second seam ownership or closure is wrong')
    call validate_eqdsk_composite_r_interval(partition, decomposition, status)
    call require(status == EQDSK_R_OWNERSHIP_SUCCESS, &
        'valid manufactured decomposition failed its verifier')

    call decompose_eqdsk_composite_r_interval(partition, 4.5_dp, 7.0_dp, &
        decomposition, status)
    call require(status == EQDSK_R_OWNERSHIP_SUCCESS .and. &
        decomposition%piece_count == 2, &
        'query ending at seam2 was not represented')
    call require(decomposition%pieces(1)%source_branch_id == &
        EQDSK_R_BRANCH_AXIS .and. decomposition%pieces(1)%r_hi == 7.0_dp .and. &
        .not. decomposition%pieces(1)%right_closed .and. &
        decomposition%pieces(2)%source_branch_id == EQDSK_R_BRANCH_OUTBOARD .and. &
        decomposition%pieces(2)%point_only .and. &
        decomposition%pieces(2)%left_boundary_id == &
        EQDSK_R_SEAM_AXIS_OUTBOARD_ID, &
        'second seam endpoint was assigned to the wrong branch')
    call validate_eqdsk_composite_r_interval(partition, decomposition, status)
    call require(status == EQDSK_R_OWNERSHIP_SUCCESS, &
        'seam2 endpoint decomposition failed its verifier')

    call decompose_eqdsk_composite_r_interval(partition, 7.0_dp, 8.0_dp, &
        decomposition, status)
    call require(status == EQDSK_R_OWNERSHIP_SUCCESS .and. &
        decomposition%piece_count == 1 .and. &
        decomposition%pieces(1)%source_branch_id == EQDSK_R_BRANCH_OUTBOARD .and. &
        decomposition%pieces(1)%left_boundary_id == &
        EQDSK_R_SEAM_AXIS_OUTBOARD_ID, &
        'query starting at seam2 was assigned to the wrong branch')

    call check_closed_point(4.0_dp, EQDSK_R_BRANCH_AXIS, &
        EQDSK_R_SEAM_INBOARD_AXIS_ID)
    call check_closed_point(7.0_dp, EQDSK_R_BRANCH_OUTBOARD, &
        EQDSK_R_SEAM_AXIS_OUTBOARD_ID)
    call check_closed_point(5.0_dp, EQDSK_R_BRANCH_AXIS, 0)
    call check_closed_point(10.0_dp, EQDSK_R_BRANCH_OUTBOARD, &
        EQDSK_R_RIGHT_ENDPOINT_ID)

    damaged = partition
    damaged%branches(1)%right_closed = .true.
    call validate_eqdsk_composite_r_partition(damaged, status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_PARTITION, &
        'closed inboard branch was accepted')

    damaged = partition
    damaged%branches(2)%source_r_lo = 4.1_dp
    call validate_eqdsk_composite_r_partition(damaged, status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_PARTITION, &
        'source seam gap was accepted')

    damaged = partition
    damaged%seams(1)%seam_id = 999
    call validate_eqdsk_composite_r_partition(damaged, status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_PARTITION, &
        'forged seam certificate was accepted')

    damaged = partition
    damaged%branches(2)%left_boundary_id = 999
    call validate_eqdsk_composite_r_partition(damaged, status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_PARTITION, &
        'forged branch boundary certificate was accepted')

    damaged = partition
    damaged%branches(1)%source_graph_certificate_id = &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID - 1
    call validate_eqdsk_composite_r_partition(damaged, status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_PARTITION, &
        'forged source graph certificate was accepted')

    call decompose_eqdsk_composite_r_interval(partition, 0.0_dp, 4.0_dp, &
        decomposition, status)
    call require(status == EQDSK_R_OWNERSHIP_OUT_OF_RANGE, &
        'out-of-domain query was accepted')

    damaged = partition
    damaged%certificate_id = 0
    call validate_eqdsk_composite_r_partition(damaged, status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_PARTITION, &
        'missing partition certificate was accepted')

    source_certificate_ids(2) = EQDSK_CUT_GRAPH_CERTIFICATE_ID - 1
    call build_eqdsk_composite_r_partition_from_bounds(physical_bounds, &
        source_bounds, source_certificate_ids, damaged, status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_INPUT, &
        'missing source certificate was accepted')

    physical_bounds(3) = ieee_value(0.0_dp, ieee_quiet_nan)
    call build_eqdsk_composite_r_partition_from_bounds(physical_bounds, &
        source_bounds, [EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
        EQDSK_CUT_GRAPH_CERTIFICATE_ID, EQDSK_CUT_GRAPH_CERTIFICATE_ID], damaged, status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_INPUT, &
        'nonfinite seam was accepted')

    call decompose_eqdsk_composite_r_interval(partition, 1.0_dp, 10.0_dp, &
        damaged_decomposition, status)
    damaged_decomposition%pieces(1)%r_hi = 5.0_dp
    call validate_eqdsk_composite_r_interval(partition, damaged_decomposition, &
        status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_PARTITION, &
        'piece outside its assigned branch was accepted')

    call decompose_eqdsk_composite_r_interval(partition, 1.0_dp, 10.0_dp, &
        damaged_decomposition, status)
    damaged_decomposition%pieces(2)%source_branch_id = EQDSK_R_BRANCH_OUTBOARD
    call validate_eqdsk_composite_r_interval(partition, damaged_decomposition, &
        status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_PARTITION, &
        'forged branch coverage was accepted')

    call decompose_eqdsk_composite_r_interval(partition, 1.0_dp, 10.0_dp, &
        damaged_decomposition, status)
    damaged_decomposition%pieces(1)%right_boundary_id = 999
    call validate_eqdsk_composite_r_interval(partition, damaged_decomposition, &
        status)
    call require(status == EQDSK_R_OWNERSHIP_INVALID_PARTITION, &
        'forged interval boundary certificate was accepted')

    write (*, '(a)') 'test_gc_eqdsk_composite_r_ownership OK'

contains

    subroutine check_closed_point(radius, expected_branch, expected_boundary)
        real(dp), intent(in) :: radius
        integer, intent(in) :: expected_branch, expected_boundary

        call decompose_eqdsk_composite_r_interval(partition, radius, radius, &
            decomposition, status)
        call require(status == EQDSK_R_OWNERSHIP_SUCCESS .and. &
            decomposition%piece_count == 1 .and. &
            decomposition%pieces(1)%point_only .and. &
            decomposition%pieces(1)%source_branch_id == expected_branch .and. &
            decomposition%pieces(1)%left_boundary_id == expected_boundary .and. &
            decomposition%pieces(1)%right_boundary_id == expected_boundary, &
            'closed point ownership oracle failed')
        call validate_eqdsk_composite_r_interval(partition, decomposition, &
            status)
        call require(status == EQDSK_R_OWNERSHIP_SUCCESS, &
            'closed point failed its proof-carrying verifier')
    end subroutine check_closed_point

    subroutine locate_and_check(radius, expected_branch, expected_seam, &
            expected_source)
        real(dp), intent(in) :: radius
        integer, intent(in) :: expected_branch, expected_seam, expected_source

        call locate_eqdsk_composite_r_point(partition, radius, branch, seam, &
            source_id, status)
        call require(status == EQDSK_R_OWNERSHIP_SUCCESS .and. &
            branch == expected_branch .and. seam == expected_seam .and. &
            source_id == expected_source, 'point ownership oracle failed')
    end subroutine locate_and_check

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_composite_r_ownership
