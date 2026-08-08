program test_gc_operational_topology_chain
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_operational_class_assembly, only: &
        GC_CLASS_ASSEMBLY_SUCCESS, GC_CLASS_BOUNDARY_REGULAR, &
        assemble_gc_operational_classes, gc_operational_allowed_interval_t, &
        gc_operational_class_assembly_result_t
    use neort_gc_operational_fixed_points, only: &
        GC_FIXED_POINT_O, GC_FIXED_POINT_SUCCESS, GC_FIXED_POINT_X, &
        find_gc_operational_fixed_points, gc_operational_fixed_point_result_t
    use neort_gc_operational_partner_crossings, only: &
        GC_PARTNER_BOUNDARY_USUAL, GC_PARTNER_BOUNDARY_X, &
        GC_PARTNER_SUCCESS, find_gc_operational_partner_crossings, &
        gc_operational_partner_options_t, gc_operational_partner_result_t
    use neort_gc_operational_scalar_roots, only: &
        gc_operational_root_options_t
    use util_for_test, only: pass_test
    implicit none

    type(gc_operational_root_options_t) :: root_options
    type(gc_operational_partner_options_t) :: partner_options
    type(gc_operational_allowed_interval_t) :: allowed
    type(gc_operational_fixed_point_result_t) :: fixed_points
    type(gc_operational_partner_result_t) :: partners
    type(gc_operational_class_assembly_result_t) :: classes

    root_options%initial_intervals = 64
    root_options%maximum_iterations = 256
    root_options%maximum_roots = 16
    root_options%relative_tolerance = 1.0e-13_dp
    partner_options%root = root_options
    partner_options%x_exclusion_fraction = 1.0e-7_dp

    allowed%component_id = 1
    allowed%sigma = 1
    allowed%x_lo = 0.0_dp
    allowed%x_hi = 8.0_dp
    allowed%left_kind = GC_CLASS_BOUNDARY_REGULAR
    allowed%right_kind = GC_CLASS_BOUNDARY_REGULAR
    allowed%canonical_lo = -24.0_dp
    allowed%canonical_hi = 72.0_dp

    call find_gc_operational_fixed_points(canonical_jet, classifier, &
        allowed%x_lo, allowed%x_hi, root_options, fixed_points)
    call require(fixed_points%status == GC_FIXED_POINT_SUCCESS, &
        'full-chain fixed-point search failed')
    call require(fixed_points%complete .and. fixed_points%npoints == 2, &
        'full-chain fixed-point set is incomplete')
    call require(fixed_points%n_o_points == 1 .and. fixed_points%n_x_points == 1, &
        'full-chain O/X counts are wrong')
    call require(all(fixed_points%points%kind == [GC_FIXED_POINT_X, &
        GC_FIXED_POINT_O]), 'fixed-point classification is wrong')
    call require(maxval(abs([fixed_points%points%x] - &
        [2.0_dp, 14.0_dp/3.0_dp])) < 1.0e-9_dp, &
        'fixed-point locations are wrong')
    call require(maxval(abs([fixed_points%points%canonical_momentum] - &
        [0.0_dp, -256.0_dp/27.0_dp])) < 1.0e-8_dp, &
        'fixed-point canonical levels are wrong')
    call require(maxval(abs([fixed_points%points%first_derivative])) < 1.0e-8_dp, &
        'stationary first derivatives are wrong')
    call require(maxval(abs([fixed_points%points%second_derivative] - &
        [-8.0_dp, 8.0_dp])) < 1.0e-8_dp, &
        'stationary second derivatives are wrong')
    call require(maxval(abs([fixed_points%points%flow_discriminant] - &
        [8.0_dp, -8.0_dp])) < 1.0e-8_dp, &
        'classifier D=-P'''' sign is wrong')

    call find_gc_operational_partner_crossings(canonical_jet, fixed_points, &
        allowed%x_lo, allowed%x_hi, partner_options, partners)
    call require(partners%status == GC_PARTNER_SUCCESS, &
        'full-chain partner search failed')
    call require(partners%complete .and. partners%npairs == 1 .and. &
        partners%nboundaries == 2, 'full-chain partner set is incomplete')
    call require(partners%pairs(1)%x_point_id == &
        fixed_points%points(1)%point_id .and. &
        partners%pairs(1)%regular_crossing_count == 1, &
        'X-point partner cardinality is wrong')
    call require(all(partners%boundaries%kind == &
        [GC_PARTNER_BOUNDARY_X, GC_PARTNER_BOUNDARY_USUAL]), &
        'partner boundary roles are wrong')
    call require(maxval(abs([partners%boundaries%x] - &
        [2.0_dp, 6.0_dp])) < 1.0e-8_dp, &
        'partner boundary locations are wrong')
    call require(maxval(abs([partners%boundaries%canonical_momentum] - &
        [0.0_dp, 0.0_dp])) < 1.0e-8_dp, &
        'partner boundary levels are wrong')
    call require(maxval(abs([partners%boundaries%canonical_residual])) < 1.0e-8_dp, &
        'partner residuals are wrong')

    call assemble_gc_operational_classes(allowed, fixed_points, partners, classes)
    call require(classes%status == GC_CLASS_ASSEMBLY_SUCCESS .and. &
        classes%complete .and. classes%nclasses == 3, &
        'full-chain class assembly failed')
    call require(maxval(abs([classes%classes%x_lo] - &
        [0.0_dp, 2.0_dp, 6.0_dp])) < 1.0e-8_dp .and. &
        maxval(abs([classes%classes%x_hi] - &
        [2.0_dp, 6.0_dp, 8.0_dp])) < 1.0e-8_dp, &
        'class intervals are wrong')
    call require(all(classes%classes%ifuntype == [14, 43, 31]), &
        'POTATO ifuntype map is wrong')
    call require(maxval(abs([classes%classes%canonical_total_variation] - &
        [24.0_dp, 512.0_dp/27.0_dp, 72.0_dp])) < 1.0e-7_dp, &
        'class canonical variations are wrong')
    call require(abs(sum(classes%classes%canonical_total_variation) - &
        3104.0_dp/27.0_dp) < 1.0e-7_dp, &
        'total canonical variation is wrong')
    call require(all(classes%classes%left_boundary_id == [0, 1, 2]) .and. &
        all(classes%classes%right_boundary_id == [1, 2, 0]), &
        'class boundary ownership is wrong')

    call pass_test

contains

    subroutine canonical_jet(x, canonical_momentum, first_derivative, &
            second_derivative, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: canonical_momentum
        real(dp), intent(out) :: first_derivative, second_derivative
        integer, intent(out) :: status

        canonical_momentum = (x-2.0_dp)**2*(x-6.0_dp)
        first_derivative = (x-2.0_dp)*(3.0_dp*x-14.0_dp)
        second_derivative = 6.0_dp*x-20.0_dp
        status = 0
    end subroutine canonical_jet

    subroutine classifier(x, discriminant, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: discriminant
        integer, intent(out) :: status

        discriminant = -(6.0_dp*x-20.0_dp)
        status = 0
    end subroutine classifier

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_operational_topology_chain
