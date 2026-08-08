program test_gc_eqdsk_composite_cut_atlas
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: R0
    use field_eq_mod, only: nrad, nzet, rad, zet
    use geoflux_coordinates, only: geoflux_get_axis
    use neort_gc_eqdsk_composite_cut_atlas, only: &
        EQDSK_COMPOSITE_ATLAS_SUCCESS, build_eqdsk_composite_cut_atlas, &
        eqdsk_composite_cut_atlas_options_t, eqdsk_composite_cut_atlas_t, &
        validate_eqdsk_composite_cut_atlas
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        eqdsk_cylindrical_field_t, initialize_eqdsk_cylindrical_field
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_ATLAS_SUCCESS, build_eqdsk_cut_graph_atlas, &
        eqdsk_cut_graph_atlas_options_t, eqdsk_cut_graph_atlas_t, &
        map_eqdsk_cut_graph_atlas_flux
    implicit none

    type(eqdsk_cylindrical_field_t) :: field
    type(eqdsk_cut_graph_atlas_t) :: inboard_seed_atlas
    type(eqdsk_cut_graph_atlas_t) :: outboard_seed_atlas
    type(eqdsk_cut_graph_atlas_options_t) :: seed_options
    type(eqdsk_composite_cut_atlas_options_t) :: options
    type(eqdsk_composite_cut_atlas_t) :: atlas
    character(len=1024) :: path
    real(dp), parameter :: target = 0.5_dp
    real(dp) :: axis_Z, inboard_lo, inboard_hi, outboard_lo, outboard_hi
    real(dp) :: inboard_position(3), outboard_position(3), tangent(3)
    real(dp) :: dZ_dR, dpsihat_dR, half_R, half_Z
    real(dp) :: axis_box(4), inboard_box(4), outboard_box(4)
    integer :: axis_R_index, axis_Z_index, status

    call get_environment_variable('EQDSK_FILE', path)
    call require(len_trim(path) > 0, 'EQDSK_FILE is required')
    call initialize_eqdsk_cylindrical_field(trim(path), 1.0_dp, field, status)
    call require(status == 0, 'failed to initialize circular EQDSK')
    call geoflux_get_axis(axis_box(1), axis_Z)
    call require(abs(axis_box(1)-R0) <= 1.0e-12_dp*max(1.0_dp, abs(R0)), &
        'geoflux and direct-field axis radii disagree')

    axis_R_index = minloc(abs(rad-R0), dim=1)
    axis_Z_index = minloc(abs(zet-axis_Z), dim=1)
    call require(axis_R_index > 1 .and. axis_R_index < nrad, &
        'axis R box is not interior')
    call require(axis_Z_index > 1 .and. axis_Z_index < nzet, &
        'axis Z box is not interior')
    axis_box = [rad(axis_R_index-1), rad(axis_R_index+1), &
        zet(axis_Z_index-1), zet(axis_Z_index+1)]

    inboard_lo = field%domain_R_min+0.10_dp*(field%domain_R_max- &
        field%domain_R_min)
    inboard_hi = R0-0.10_dp*(field%domain_R_max-field%domain_R_min)
    outboard_lo = R0+0.10_dp*(field%domain_R_max-field%domain_R_min)
    outboard_hi = field%domain_R_max-0.10_dp*(field%domain_R_max- &
        field%domain_R_min)
    seed_options%max_r_depth = 12
    seed_options%max_z_depth = 12
    seed_options%map_absolute_tolerance = 1.0e-12_dp
    seed_options%map_relative_tolerance = 1.0e-13_dp
    call build_eqdsk_cut_graph_atlas(inboard_seed_atlas, inboard_lo, &
        inboard_hi, -10.0_dp, 10.0_dp, 0.0_dp, 1.0_dp, seed_options, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'inboard seed atlas failed')
    call build_eqdsk_cut_graph_atlas(outboard_seed_atlas, outboard_lo, &
        outboard_hi, -10.0_dp, 10.0_dp, 0.0_dp, 1.0_dp, seed_options, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'outboard seed atlas failed')
    call map_eqdsk_cut_graph_atlas_flux(inboard_seed_atlas, target, &
        inboard_position, tangent, dZ_dR, dpsihat_dR, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'inboard endpoint seed failed')
    call map_eqdsk_cut_graph_atlas_flux(outboard_seed_atlas, target, &
        outboard_position, tangent, dZ_dR, dpsihat_dR, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'outboard endpoint seed failed')

    half_R = 0.20_dp*minval(rad(2:nrad)-rad(1:nrad-1))
    half_Z = 0.20_dp*minval(zet(2:nzet)-zet(1:nzet-1))
    inboard_box = [inboard_position(1)-half_R, &
        inboard_position(1)+half_R, inboard_position(2)-half_Z, &
        inboard_position(2)+half_Z]
    outboard_box = [outboard_position(1)-half_R, &
        outboard_position(1)+half_R, outboard_position(2)-half_Z, &
        outboard_position(2)+half_Z]
    options%graph = seed_options
    options%endpoint%max_newton_iterations = 8
    call build_eqdsk_composite_cut_atlas(axis_box, inboard_box, &
        outboard_box, inboard_position(1:2), outboard_position(1:2), target, &
        1.0_dp, -10.0_dp, 10.0_dp, options, atlas, status)
    call require(status == EQDSK_COMPOSITE_ATLAS_SUCCESS, &
        'interior circular composite cut atlas failed')
    call require(atlas%geometric_completeness_certified .and. &
        atlas%branch_connectivity_certified .and. &
        atlas%surface_intersection_pair_certified, &
        'composite atlas omitted a geometric theorem gate')
    call require(.not. atlas%orbit_crossing_multiplicity_certified, &
        'geometric atlas overclaimed orbit-crossing multiplicity')
    call require(atlas%inboard_graph%flux_monotonicity_sign == -1 .and. &
        atlas%outboard_graph%flux_monotonicity_sign == 1, &
        'composite branch flux orientations are wrong')
    call validate_eqdsk_composite_cut_atlas(atlas, status)
    call require(status == EQDSK_COMPOSITE_ATLAS_SUCCESS, &
        'fresh composite atlas failed validation')

    write (*, '(a)') 'test_gc_eqdsk_composite_cut_atlas OK'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_composite_cut_atlas
