program test_gc_eqdsk_cut_graph_atlas
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: nrad, nzet, zet
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        eqdsk_cylindrical_field_t, initialize_eqdsk_cylindrical_field
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_ATLAS_SUCCESS, EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
        build_eqdsk_cut_graph_atlas, clear_eqdsk_cut_graph_atlas, &
        eqdsk_cut_graph_atlas_options_t, eqdsk_cut_graph_atlas_t, &
        map_eqdsk_cut_graph_atlas, validate_eqdsk_cut_graph_atlas
    implicit none

    type(eqdsk_cylindrical_field_t) :: field
    type(eqdsk_cut_graph_atlas_t) :: inboard_atlas, outboard_atlas, full_atlas
    type(eqdsk_cut_graph_atlas_options_t) :: options, full_options
    character(len=1024) :: path
    real(dp) :: axis_R, inboard_lo, inboard_hi, outboard_lo, outboard_hi
    integer :: status

    call get_environment_variable('EQDSK_FILE', path)
    call require(len_trim(path) > 0, 'EQDSK_FILE is required')
    call initialize_eqdsk_cylindrical_field(trim(path), 1.0_dp, field, status)
    call require(status == 0, 'failed to initialize circular EQDSK')
    call require(nrad >= 3 .and. nzet >= 3, 'EQDSK grid is too small')

    options%max_r_depth = 12
    options%max_z_depth = 12
    options%map_absolute_tolerance = 1.0e-11_dp
    options%map_relative_tolerance = 1.0e-13_dp
    axis_R = 0.5_dp*(field%domain_R_min+field%domain_R_max)
    inboard_lo = field%domain_R_min+0.10_dp*(field%domain_R_max- &
        field%domain_R_min)
    inboard_hi = axis_R-0.10_dp*(field%domain_R_max- &
        field%domain_R_min)
    outboard_lo = axis_R+0.10_dp*(field%domain_R_max- &
        field%domain_R_min)
    outboard_hi = field%domain_R_max-0.10_dp*(field%domain_R_max- &
        field%domain_R_min)
    call build_eqdsk_cut_graph_atlas(inboard_atlas, inboard_lo, inboard_hi, &
        -10.0_dp, 10.0_dp, 0.0_dp, 1.0_dp, options, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'regular inboard circular cut branch was not certified')
    call check_regular_branch(inboard_atlas, inboard_lo, inboard_hi)
    call build_eqdsk_cut_graph_atlas(outboard_atlas, outboard_lo, &
        outboard_hi, -10.0_dp, 10.0_dp, 0.0_dp, 1.0_dp, options, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'regular outboard circular cut branch was not certified')
    call check_regular_branch(outboard_atlas, outboard_lo, outboard_hi)

    ! The public POTATO fixture clamps psi outside the LCFS before applying a
    ! quintic spline.  That nonsmooth SOL transition creates a resolved pair of
    ! off-midplane Eq.13 roots with psi_hat<1.  A physical-domain completeness
    ! certificate must reject that extra topology rather than silently select
    ! the symmetry branch.
    full_options = options
    full_options%max_r_depth = 6
    full_options%max_z_depth = 6
    call build_eqdsk_cut_graph_atlas(full_atlas, field%domain_R_min, &
        field%domain_R_max, zet(1), zet(nzet), 0.0_dp, 1.0_dp, full_options, &
        status)
    call require(status /= EQDSK_CUT_ATLAS_SUCCESS, &
        'clamped-SOL extra cut branches were incorrectly certified away')
    call require(.not. full_atlas%global_completeness_certified, &
        'failed full-domain atlas retained a completeness claim')

    call clear_eqdsk_cut_graph_atlas(inboard_atlas)
    call validate_eqdsk_cut_graph_atlas(inboard_atlas, status)
    call require(status /= EQDSK_CUT_ATLAS_SUCCESS, &
        'cleared atlas remained valid')
    write (*, '(a)') 'test_gc_eqdsk_cut_graph_atlas OK'

contains

    subroutine check_regular_branch(atlas, r_lo, r_hi)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        real(dp), intent(in) :: r_lo, r_hi

        real(dp) :: radius, position(3), dposition(3), dZ_dR, dpsihat_dR
        real(dp) :: scale
        integer :: i, local_status

        call require(atlas%global_completeness_certified, &
            'regular branch did not record rectangular completeness')
        call require(atlas%certificate_id == EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
            'regular branch certificate ID is not the generated-cut contract')
        call validate_eqdsk_cut_graph_atlas(atlas, local_status)
        call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
            'fresh regular branch failed structural validation')
        scale = max(1.0_dp, maxval(abs(zet)))
        do i = 0, 8
            radius = r_lo+(r_hi-r_lo)*real(i,dp)/8.0_dp
            call map_eqdsk_cut_graph_atlas(atlas, radius, position, &
                dposition, dZ_dR, dpsihat_dR, local_status)
            call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
                'certified regular circular cut map failed')
            call require(abs(position(2)) <= 5.0e-10_dp*scale, &
                'regular circular Eq13 branch departed from the midplane')
            call require(abs(dposition(1)-1.0_dp) <= 1.0e-13_dp, &
                'cut map changed its physical R coordinate')
            call require(abs(dposition(2)) <= 5.0e-9_dp, &
                'regular circular cut tangent departed from the R direction')
            call require(dZ_dR == dposition(2), &
                'optional and vector cut slopes disagree')
            call require(dpsihat_dR == dpsihat_dR, &
                'cut flux derivative is nonfinite')
        end do
    end subroutine check_regular_branch

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_cut_graph_atlas
