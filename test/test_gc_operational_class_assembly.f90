program test_gc_operational_class_assembly
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_certified_interval_roots, only: gc_interval_t
    use neort_gc_operational_class_assembly, only: &
        GC_CLASS_ASSEMBLY_DUPLICATE_BOUNDARY, GC_CLASS_ASSEMBLY_SUCCESS, &
        GC_CLASS_ASSEMBLY_INVALID_MEASURE, &
        GC_CLASS_BOUNDARY_REGULAR, GC_CLASS_BOUNDARY_TURNING, &
        assemble_gc_operational_classes, &
        gc_operational_allowed_interval_t, gc_operational_class_assembly_result_t, &
        gc_operational_canonical_measure_evidence_t
    use neort_gc_operational_fixed_points, only: &
        GC_FIXED_POINT_O, GC_FIXED_POINT_SUCCESS, GC_FIXED_POINT_X, &
        gc_operational_fixed_point_t, &
        gc_operational_fixed_point_result_t
    use neort_gc_operational_partner_crossings, only: &
        GC_PARTNER_BOUNDARY_USUAL, GC_PARTNER_BOUNDARY_X, GC_PARTNER_SUCCESS, &
        gc_operational_separatrix_boundary_t, &
        gc_operational_partner_result_t
    use util_for_test, only: pass_test
    implicit none

    type(gc_operational_allowed_interval_t) :: allowed
    type(gc_operational_fixed_point_result_t) :: fixed_points
    type(gc_operational_partner_result_t) :: partners
    type(gc_operational_class_assembly_result_t) :: result
    type(gc_operational_canonical_measure_evidence_t), allocatable :: evidence(:)

    call make_manufactured_input(allowed, fixed_points, partners)
    allocate(evidence(3))
    call set_measure(evidence(1), 4.75_dp, 5.25_dp, 901)
    call set_measure(evidence(2), 7.50_dp, 8.50_dp, 902)
    call set_measure(evidence(3), 2.75_dp, 3.25_dp, 903)
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result, &
        evidence)
    call require(result%status == GC_CLASS_ASSEMBLY_SUCCESS, &
        'manufactured class assembly failed')
    call require(result%complete .and. result%nclasses == 3, &
        'class count or completion gate is wrong')
    call require(all(result%classes(1:3)%class_id == [1, 2, 3]), &
        'class identifiers are not ordered')
    call require(all(result%classes(1:3)%ifuntype == [14, 43, 32]), &
        'POTATO ifuntype map is wrong')
    call require(all(result%classes(1:3)%left_kind == [1, 4, 3]) .and. &
        all(result%classes(1:3)%right_kind == [4, 3, 2]), &
        'left/right boundary kinds are wrong')
    call require(maxval(abs([result%classes%x_lo]-[0.0_dp, 4.0_dp, 7.0_dp])) &
        < 1.0e-12_dp, 'class lower coordinates are not sorted')
    call require(maxval(abs([result%classes%x_hi]-[4.0_dp, 7.0_dp, 10.0_dp])) &
        < 1.0e-12_dp, 'class upper coordinates are not sorted')
    call require(maxval(abs([result%classes%canonical_total_variation]- &
        [5.0_dp, 8.0_dp, 3.0_dp])) < 1.0e-12_dp, &
        'certified canonical measure midpoint is wrong')
    call require(result%classes(1)%canonical_measure_certified .and. &
        result%classes(1)%canonical_measure_certificate_id == 901 .and. &
        result%classes(1)%canonical_measure_enclosure%lo == 4.75_dp .and. &
        result%classes(1)%canonical_measure_enclosure%hi == 5.25_dp, &
        'certified canonical measure evidence was not preserved')
    call require(abs(sum(result%classes%canonical_total_variation)-16.0_dp) &
        < 1.0e-12_dp, 'class total variation does not partition the fixture')

    evidence(2)%enclosure = gc_interval_t(-0.1_dp, 8.5_dp)
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result, &
        evidence)
    call require(result%status == GC_CLASS_ASSEMBLY_INVALID_MEASURE .and. &
        .not. result%complete, 'out-of-range measure evidence was accepted')
    evidence(2)%enclosure = gc_interval_t(7.50_dp, 8.50_dp)
    evidence(2)%certified = .false.
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result, &
        evidence)
    call require(result%status == GC_CLASS_ASSEMBLY_INVALID_MEASURE .and. &
        .not. result%complete, 'uncertified measure evidence was accepted')
    evidence(2)%certified = .true.
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result, &
        evidence)
    call require(result%status == GC_CLASS_ASSEMBLY_SUCCESS .and. &
        result%complete, 'valid measure evidence was not reusable')
    call require(result%classes(1)%left_boundary_id == 0 .and. &
        result%classes(1)%right_boundary_id == 2 .and. &
        result%classes(2)%left_boundary_id == 2 .and. &
        result%classes(2)%right_boundary_id == 1 .and. &
        result%classes(3)%left_boundary_id == 1 .and. &
        result%classes(3)%right_boundary_id == 0, &
        'boundary ownership identifiers are wrong')

    partners%boundaries(1)%x = partners%boundaries(2)%x
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result)
    call require(result%status == GC_CLASS_ASSEMBLY_DUPLICATE_BOUNDARY .and. &
        .not. result%complete, 'duplicate separatrix boundary was accepted')

    call make_no_separatrix_input(allowed, fixed_points, partners)
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result)
    call require(result%status == GC_CLASS_ASSEMBLY_SUCCESS .and. &
        result%nclasses == 1, 'single allowed interval was not retained')
    call require(.not. result%classes(1)%canonical_measure_certified, &
        'compatibility path falsely claimed certified measure evidence')
    call require(result%classes(1)%ifuntype == 11 .and. &
        abs(result%classes(1)%canonical_total_variation-4.0_dp) < 1.0e-12_dp, &
        'single-class canonical variation is wrong')

    call make_x_without_partners_input(allowed, fixed_points, partners)
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result)
    call require_rejected(result, 'fixed X point without partners was accepted')

    call make_pair_without_x_boundary_input(allowed, fixed_points, partners)
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result)
    call require_rejected(result, 'partner pair without its X boundary was accepted')

    call make_single_x_valid_input(allowed, fixed_points, partners)
    fixed_points%points(2)%point_id = fixed_points%points(1)%point_id
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result)
    call require_rejected(result, 'duplicate fixed-point identifier was accepted')

    call make_single_x_valid_input(allowed, fixed_points, partners)
    partners%boundaries(2)%boundary_id = partners%boundaries(1)%boundary_id
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result)
    call require_rejected(result, 'duplicate partner-boundary identifier was accepted')

    call make_endpoint_fixed_point_input(allowed, fixed_points, partners)
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result)
    call require_rejected(result, 'endpoint fixed point was accepted')

    call make_single_x_valid_input(allowed, fixed_points, partners)
    partners%boundaries(2)%canonical_residual = 1.0_dp
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result)
    call require_rejected(result, 'large partner residual was accepted')

    call make_multiple_x_input(allowed, fixed_points, partners)
    call assemble_gc_operational_classes(allowed, fixed_points, partners, result)
    call require(result%status == GC_CLASS_ASSEMBLY_SUCCESS .and. &
        result%complete .and. result%nclasses == 7, &
        'multiple-X class assembly failed')
    call require(all(result%classes%ifuntype == [14, 43, 33, 33, 34, 43, 31]), &
        'multiple-X ifuntype map is wrong')
    call require(maxval(abs([result%classes%x_lo]- &
        [0.0_dp, 4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp, 14.0_dp, 16.0_dp])) &
        < 1.0e-12_dp, 'multiple-X lower boundaries are wrong')
    call require(maxval(abs([result%classes%x_hi]- &
        [4.0_dp, 6.0_dp, 8.0_dp, 10.0_dp, 14.0_dp, 16.0_dp, 20.0_dp])) &
        < 1.0e-12_dp, 'multiple-X upper boundaries are wrong')
    call require(all(result%classes%left_boundary_id == [0, 1, 2, 3, 4, 5, 6]) .and. &
        all(result%classes%right_boundary_id == [1, 2, 3, 4, 5, 6, 0]), &
        'multiple-X boundary ownership is wrong')
    call require(maxval(abs([result%classes%canonical_total_variation]- &
        [4.0_dp, 0.0_dp, 4.0_dp, 4.0_dp, 4.0_dp, 0.0_dp, 2.0_dp])) &
        < 1.0e-12_dp, 'multiple-X canonical variation is wrong')
    call require(abs(sum(result%classes%canonical_total_variation)-18.0_dp) &
        < 1.0e-12_dp, 'multiple-X variation does not partition the fixture')

    call pass_test

contains

    subroutine make_manufactured_input(allowed, fixed_points, partners)
        type(gc_operational_allowed_interval_t), intent(out) :: allowed
        type(gc_operational_fixed_point_result_t), intent(out) :: fixed_points
        type(gc_operational_partner_result_t), intent(out) :: partners

        allowed%component_id = 9
        allowed%sigma = 1
        allowed%x_lo = 0.0_dp
        allowed%x_hi = 10.0_dp
        allowed%left_kind = GC_CLASS_BOUNDARY_REGULAR
        allowed%right_kind = GC_CLASS_BOUNDARY_TURNING
        allowed%canonical_lo = 0.0_dp
        allowed%canonical_hi = 2.0_dp

        fixed_points%status = GC_FIXED_POINT_SUCCESS
        fixed_points%complete = .true.
        fixed_points%npoints = 3
        fixed_points%n_o_points = 2
        fixed_points%n_x_points = 1
        allocate(fixed_points%points(3))
        fixed_points%points(1)%point_id = 12
        fixed_points%points(1)%kind = GC_FIXED_POINT_O
        fixed_points%points(1)%x = 5.0_dp
        fixed_points%points(1)%canonical_momentum = 1.0_dp
        fixed_points%points(2)%point_id = 10
        fixed_points%points(2)%kind = GC_FIXED_POINT_O
        fixed_points%points(2)%x = 2.0_dp
        fixed_points%points(2)%canonical_momentum = 3.0_dp
        fixed_points%points(3)%point_id = 11
        fixed_points%points(3)%kind = GC_FIXED_POINT_X
        fixed_points%points(3)%x = 4.0_dp
        fixed_points%points(3)%canonical_momentum = 5.0_dp

        partners%status = GC_PARTNER_SUCCESS
        partners%complete = .true.
        partners%npairs = 1
        partners%nboundaries = 2
        allocate(partners%pairs(1), partners%boundaries(2))
        partners%pairs(1)%pair_id = 1
        partners%pairs(1)%x_point_id = 11
        partners%pairs(1)%regular_crossing_count = 1
        partners%pairs(1)%canonical_momentum = 5.0_dp
        partners%boundaries(1)%boundary_id = 1
        partners%boundaries(1)%kind = GC_PARTNER_BOUNDARY_USUAL
        partners%boundaries(1)%pair_id = 1
        partners%boundaries(1)%x_point_id = 11
        partners%boundaries(1)%x = 7.0_dp
        partners%boundaries(1)%canonical_momentum = 5.0_dp
        partners%boundaries(2)%boundary_id = 2
        partners%boundaries(2)%kind = GC_PARTNER_BOUNDARY_X
        partners%boundaries(2)%pair_id = 1
        partners%boundaries(2)%x_point_id = 11
        partners%boundaries(2)%x = 4.0_dp
        partners%boundaries(2)%canonical_momentum = 5.0_dp
    end subroutine make_manufactured_input

    subroutine make_no_separatrix_input(allowed, fixed_points, partners)
        type(gc_operational_allowed_interval_t), intent(out) :: allowed
        type(gc_operational_fixed_point_result_t), intent(out) :: fixed_points
        type(gc_operational_partner_result_t), intent(out) :: partners

        allowed%component_id = 1
        allowed%sigma = -1
        allowed%x_lo = 0.0_dp
        allowed%x_hi = 1.0_dp
        allowed%left_kind = GC_CLASS_BOUNDARY_REGULAR
        allowed%right_kind = GC_CLASS_BOUNDARY_REGULAR
        allowed%canonical_lo = 1.0_dp
        allowed%canonical_hi = 5.0_dp

        fixed_points%status = GC_FIXED_POINT_SUCCESS
        fixed_points%complete = .true.
        fixed_points%npoints = 1
        fixed_points%n_o_points = 1
        fixed_points%n_x_points = 0
        allocate(fixed_points%points(1))
        fixed_points%points(1)%point_id = 1
        fixed_points%points(1)%kind = GC_FIXED_POINT_O
        fixed_points%points(1)%x = 0.5_dp
        fixed_points%points(1)%canonical_momentum = 1.0_dp

        partners%status = GC_PARTNER_SUCCESS
        partners%complete = .true.
        partners%npairs = 0
        partners%nboundaries = 0
        allocate(partners%pairs(0), partners%boundaries(0))
    end subroutine make_no_separatrix_input

    subroutine make_x_without_partners_input(allowed, fixed_points, partners)
        type(gc_operational_allowed_interval_t), intent(out) :: allowed
        type(gc_operational_fixed_point_result_t), intent(out) :: fixed_points
        type(gc_operational_partner_result_t), intent(out) :: partners

        call make_single_x_interval(allowed)
        fixed_points%status = GC_FIXED_POINT_SUCCESS
        fixed_points%complete = .true.
        fixed_points%npoints = 1
        fixed_points%n_o_points = 0
        fixed_points%n_x_points = 1
        allocate(fixed_points%points(1))
        fixed_points%points(1)%point_id = 11
        fixed_points%points(1)%kind = GC_FIXED_POINT_X
        fixed_points%points(1)%x = 4.0_dp
        fixed_points%points(1)%canonical_momentum = 4.0_dp

        partners%status = GC_PARTNER_SUCCESS
        partners%complete = .true.
        partners%npairs = 0
        partners%nboundaries = 0
        allocate(partners%pairs(0), partners%boundaries(0))
    end subroutine make_x_without_partners_input

    subroutine make_pair_without_x_boundary_input(allowed, fixed_points, partners)
        type(gc_operational_allowed_interval_t), intent(out) :: allowed
        type(gc_operational_fixed_point_result_t), intent(out) :: fixed_points
        type(gc_operational_partner_result_t), intent(out) :: partners

        call make_single_x_interval(allowed)
        fixed_points%status = GC_FIXED_POINT_SUCCESS
        fixed_points%complete = .true.
        fixed_points%npoints = 1
        fixed_points%n_o_points = 0
        fixed_points%n_x_points = 1
        allocate(fixed_points%points(1))
        fixed_points%points(1)%point_id = 11
        fixed_points%points(1)%kind = GC_FIXED_POINT_X
        fixed_points%points(1)%x = 4.0_dp
        fixed_points%points(1)%canonical_momentum = 4.0_dp

        partners%status = GC_PARTNER_SUCCESS
        partners%complete = .true.
        partners%npairs = 1
        partners%nboundaries = 1
        allocate(partners%pairs(1), partners%boundaries(1))
        partners%pairs(1)%pair_id = 1
        partners%pairs(1)%x_point_id = 11
        partners%pairs(1)%regular_crossing_count = 1
        partners%pairs(1)%canonical_momentum = 4.0_dp
        partners%boundaries(1)%boundary_id = 1
        partners%boundaries(1)%kind = GC_PARTNER_BOUNDARY_USUAL
        partners%boundaries(1)%pair_id = 1
        partners%boundaries(1)%x_point_id = 11
        partners%boundaries(1)%x = 8.0_dp
        partners%boundaries(1)%canonical_momentum = 4.0_dp
    end subroutine make_pair_without_x_boundary_input

    subroutine make_single_x_valid_input(allowed, fixed_points, partners)
        type(gc_operational_allowed_interval_t), intent(out) :: allowed
        type(gc_operational_fixed_point_result_t), intent(out) :: fixed_points
        type(gc_operational_partner_result_t), intent(out) :: partners

        call make_single_x_interval(allowed)
        fixed_points%status = GC_FIXED_POINT_SUCCESS
        fixed_points%complete = .true.
        fixed_points%npoints = 2
        fixed_points%n_o_points = 1
        fixed_points%n_x_points = 1
        allocate(fixed_points%points(2))
        fixed_points%points(1)%point_id = 11
        fixed_points%points(1)%kind = GC_FIXED_POINT_X
        fixed_points%points(1)%x = 4.0_dp
        fixed_points%points(1)%canonical_momentum = 4.0_dp
        fixed_points%points(2)%point_id = 12
        fixed_points%points(2)%kind = GC_FIXED_POINT_O
        fixed_points%points(2)%x = 6.0_dp
        fixed_points%points(2)%canonical_momentum = 6.0_dp

        partners%status = GC_PARTNER_SUCCESS
        partners%complete = .true.
        partners%npairs = 1
        partners%nboundaries = 2
        allocate(partners%pairs(1), partners%boundaries(2))
        partners%pairs(1)%pair_id = 1
        partners%pairs(1)%x_point_id = 11
        partners%pairs(1)%regular_crossing_count = 1
        partners%pairs(1)%canonical_momentum = 4.0_dp
        partners%boundaries(1)%boundary_id = 1
        partners%boundaries(1)%kind = GC_PARTNER_BOUNDARY_X
        partners%boundaries(1)%pair_id = 1
        partners%boundaries(1)%x_point_id = 11
        partners%boundaries(1)%x = 4.0_dp
        partners%boundaries(1)%canonical_momentum = 4.0_dp
        partners%boundaries(2)%boundary_id = 2
        partners%boundaries(2)%kind = GC_PARTNER_BOUNDARY_USUAL
        partners%boundaries(2)%pair_id = 1
        partners%boundaries(2)%x_point_id = 11
        partners%boundaries(2)%x = 8.0_dp
        partners%boundaries(2)%canonical_momentum = 4.0_dp
    end subroutine make_single_x_valid_input

    subroutine make_endpoint_fixed_point_input(allowed, fixed_points, partners)
        type(gc_operational_allowed_interval_t), intent(out) :: allowed
        type(gc_operational_fixed_point_result_t), intent(out) :: fixed_points
        type(gc_operational_partner_result_t), intent(out) :: partners

        call make_single_x_interval(allowed)
        fixed_points%status = GC_FIXED_POINT_SUCCESS
        fixed_points%complete = .true.
        fixed_points%npoints = 1
        fixed_points%n_o_points = 1
        fixed_points%n_x_points = 0
        allocate(fixed_points%points(1))
        fixed_points%points(1)%point_id = 21
        fixed_points%points(1)%kind = GC_FIXED_POINT_O
        fixed_points%points(1)%x = allowed%x_lo
        fixed_points%points(1)%canonical_momentum = allowed%canonical_lo

        partners%status = GC_PARTNER_SUCCESS
        partners%complete = .true.
        partners%npairs = 0
        partners%nboundaries = 0
        allocate(partners%pairs(0), partners%boundaries(0))
    end subroutine make_endpoint_fixed_point_input

    subroutine make_single_x_interval(allowed)
        type(gc_operational_allowed_interval_t), intent(out) :: allowed

        allowed%component_id = 3
        allowed%sigma = 1
        allowed%x_lo = 0.0_dp
        allowed%x_hi = 10.0_dp
        allowed%left_kind = GC_CLASS_BOUNDARY_REGULAR
        allowed%right_kind = GC_CLASS_BOUNDARY_REGULAR
        allowed%canonical_lo = 0.0_dp
        allowed%canonical_hi = 10.0_dp
    end subroutine make_single_x_interval

    subroutine make_multiple_x_input(allowed, fixed_points, partners)
        type(gc_operational_allowed_interval_t), intent(out) :: allowed
        type(gc_operational_fixed_point_result_t), intent(out) :: fixed_points
        type(gc_operational_partner_result_t), intent(out) :: partners

        allowed%component_id = 5
        allowed%sigma = -1
        allowed%x_lo = 0.0_dp
        allowed%x_hi = 20.0_dp
        allowed%left_kind = GC_CLASS_BOUNDARY_REGULAR
        allowed%right_kind = GC_CLASS_BOUNDARY_REGULAR
        allowed%canonical_lo = 0.0_dp
        allowed%canonical_hi = 10.0_dp

        fixed_points%status = GC_FIXED_POINT_SUCCESS
        fixed_points%complete = .true.
        fixed_points%npoints = 6
        fixed_points%n_o_points = 4
        fixed_points%n_x_points = 2
        allocate(fixed_points%points(6))
        call set_point(fixed_points%points(1), 101, GC_FIXED_POINT_O, 2.0_dp, 2.0_dp)
        call set_point(fixed_points%points(2), 102, GC_FIXED_POINT_X, 4.0_dp, 4.0_dp)
        call set_point(fixed_points%points(3), 103, GC_FIXED_POINT_O, 7.0_dp, 6.0_dp)
        call set_point(fixed_points%points(4), 104, GC_FIXED_POINT_O, 12.0_dp, 6.0_dp)
        call set_point(fixed_points%points(5), 105, GC_FIXED_POINT_X, 14.0_dp, 8.0_dp)
        call set_point(fixed_points%points(6), 106, GC_FIXED_POINT_O, 18.0_dp, 9.0_dp)

        partners%status = GC_PARTNER_SUCCESS
        partners%complete = .true.
        partners%npairs = 2
        partners%nboundaries = 6
        allocate(partners%pairs(2), partners%boundaries(6))
        partners%pairs(1)%pair_id = 1
        partners%pairs(1)%x_point_id = 102
        partners%pairs(1)%regular_crossing_count = 2
        partners%pairs(1)%canonical_momentum = 4.0_dp
        partners%pairs(2)%pair_id = 2
        partners%pairs(2)%x_point_id = 105
        partners%pairs(2)%regular_crossing_count = 2
        partners%pairs(2)%canonical_momentum = 8.0_dp

        call set_boundary(partners%boundaries(1), 1, GC_PARTNER_BOUNDARY_X, &
            1, 102, 4.0_dp, 4.0_dp)
        call set_boundary(partners%boundaries(2), 2, GC_PARTNER_BOUNDARY_USUAL, &
            1, 102, 6.0_dp, 4.0_dp)
        call set_boundary(partners%boundaries(3), 3, GC_PARTNER_BOUNDARY_USUAL, &
            2, 105, 8.0_dp, 8.0_dp)
        call set_boundary(partners%boundaries(4), 4, GC_PARTNER_BOUNDARY_USUAL, &
            1, 102, 10.0_dp, 4.0_dp)
        call set_boundary(partners%boundaries(5), 5, GC_PARTNER_BOUNDARY_X, &
            2, 105, 14.0_dp, 8.0_dp)
        call set_boundary(partners%boundaries(6), 6, GC_PARTNER_BOUNDARY_USUAL, &
            2, 105, 16.0_dp, 8.0_dp)
    end subroutine make_multiple_x_input

    subroutine set_point(point, point_id, kind, x, canonical_momentum)
        type(gc_operational_fixed_point_t), intent(out) :: point
        integer, intent(in) :: point_id, kind
        real(dp), intent(in) :: x, canonical_momentum

        point%point_id = point_id
        point%kind = kind
        point%x = x
        point%canonical_momentum = canonical_momentum
    end subroutine set_point

    subroutine set_boundary(boundary, boundary_id, kind, pair_id, x_point_id, &
            x, canonical_momentum)
        type(gc_operational_separatrix_boundary_t), intent(out) :: boundary
        integer, intent(in) :: boundary_id, kind, pair_id, x_point_id
        real(dp), intent(in) :: x, canonical_momentum

        boundary%boundary_id = boundary_id
        boundary%kind = kind
        boundary%pair_id = pair_id
        boundary%x_point_id = x_point_id
        boundary%x = x
        boundary%canonical_momentum = canonical_momentum
    end subroutine set_boundary

    subroutine set_measure(evidence, lo, hi, certificate_id)
        type(gc_operational_canonical_measure_evidence_t), intent(out) :: evidence
        real(dp), intent(in) :: lo, hi
        integer, intent(in) :: certificate_id

        evidence = gc_operational_canonical_measure_evidence_t()
        evidence%enclosure = gc_interval_t(lo, hi)
        evidence%certificate_id = certificate_id
        evidence%certified = .true.
    end subroutine set_measure

    subroutine require_rejected(result, message)
        type(gc_operational_class_assembly_result_t), intent(in) :: result
        character(len=*), intent(in) :: message

        call require(result%status /= GC_CLASS_ASSEMBLY_SUCCESS .and. &
            .not. result%complete, message)
    end subroutine require_rejected

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_operational_class_assembly
