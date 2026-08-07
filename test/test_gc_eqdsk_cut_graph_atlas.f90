program test_gc_eqdsk_cut_graph_atlas
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: nrad, nzet, rad, zet
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        eqdsk_cylindrical_field_t, initialize_eqdsk_cylindrical_field
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_ATLAS_SUCCESS, EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
        build_eqdsk_cut_graph_atlas, clear_eqdsk_cut_graph_atlas, &
        eqdsk_cut_graph_atlas_options_t, eqdsk_cut_graph_atlas_t, &
        map_eqdsk_cut_graph_atlas, validate_eqdsk_cut_graph_atlas
    implicit none

    type(eqdsk_cylindrical_field_t) :: field
    type(eqdsk_cut_graph_atlas_t) :: atlas, boundary_atlas
    type(eqdsk_cut_graph_atlas_options_t) :: options, boundary_options
    character(len=1024) :: path
    real(dp) :: radius, position(3), dposition(3), dZ_dR, dpsihat_dR
    real(dp) :: scale
    integer :: i, status

    call get_environment_variable('EQDSK_FILE', path)
    call require(len_trim(path) > 0, 'EQDSK_FILE is required')
    call initialize_eqdsk_cylindrical_field(trim(path), 1.0_dp, field, status)
    call require(status == 0, 'failed to initialize circular EQDSK')
    call require(nrad >= 3 .and. nzet >= 3, 'EQDSK grid is too small')

    options%max_r_depth = 12
    options%max_z_depth = 12
    options%map_absolute_tolerance = 1.0e-11_dp
    options%map_relative_tolerance = 1.0e-13_dp
    call build_eqdsk_cut_graph_atlas(atlas, rad(1), rad(nrad), zet(1), &
        zet(nzet), options, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'full circular cut graph was not certified')
    call require(atlas%global_completeness_certified, &
        'atlas did not record global completeness')
    call require(atlas%certificate_id == EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
        'atlas certificate ID is not the generated-cut contract')
    call validate_eqdsk_cut_graph_atlas(atlas, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'fresh circular atlas failed structural validation')

    scale = max(1.0_dp, maxval(abs(zet)))
    do i = 0, 16
        radius = rad(1)+(rad(nrad)-rad(1))*real(i,dp)/16.0_dp
        call map_eqdsk_cut_graph_atlas(atlas, radius, position, dposition, &
            dZ_dR, dpsihat_dR, status)
        call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
            'certified circular cut map failed')
        call require(abs(position(2)) <= 5.0e-10_dp*scale, &
            'circular Eq13 cut departed from the symmetry midplane')
        call require(abs(dposition(1)-1.0_dp) <= 1.0e-13_dp, &
            'cut map changed its physical R coordinate')
        call require(abs(dposition(2)) <= 5.0e-9_dp, &
            'circular cut tangent departed from the R direction')
        call require(dZ_dR == dposition(2), &
            'optional and vector cut slopes disagree')
        call require(dpsihat_dR == dpsihat_dR, &
            'cut flux derivative is nonfinite')
    end do

    boundary_options = options
    boundary_options%max_r_depth = 1
    boundary_options%max_z_depth = 1
    call build_eqdsk_cut_graph_atlas(boundary_atlas, rad(1), rad(nrad), &
        0.0_dp, zet(nzet), boundary_options, status)
    call require(status /= EQDSK_CUT_ATLAS_SUCCESS, &
        'cut on the requested Z boundary was incorrectly certified')
    call require(.not. boundary_atlas%global_completeness_certified, &
        'failed boundary atlas retained a completeness claim')

    call clear_eqdsk_cut_graph_atlas(atlas)
    call validate_eqdsk_cut_graph_atlas(atlas, status)
    call require(status /= EQDSK_CUT_ATLAS_SUCCESS, &
        'cleared atlas remained valid')
    write (*, '(a)') 'test_gc_eqdsk_cut_graph_atlas OK'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_cut_graph_atlas
