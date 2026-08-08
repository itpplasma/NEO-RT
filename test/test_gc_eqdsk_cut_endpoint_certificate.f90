program test_gc_eqdsk_cut_endpoint_certificate
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: nrad, nzet, psi_sep, rad, zet
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        eqdsk_cylindrical_field_t, initialize_eqdsk_cylindrical_field
    use neort_gc_eqdsk_cut_endpoint_certificate, only: &
        EQDSK_ENDPOINT_CERT_AXIS_REQUIRES_LIMIT, &
        EQDSK_ENDPOINT_CERTIFICATE_ID, &
        EQDSK_ENDPOINT_CERT_SUCCESS, build_eqdsk_cut_endpoint_certificate, &
        eqdsk_cut_endpoint_certificate_t, eqdsk_cut_endpoint_options_t, &
        validate_eqdsk_cut_endpoint_certificate
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_ATLAS_SUCCESS, build_eqdsk_cut_graph_atlas, &
        eqdsk_cut_graph_atlas_options_t, eqdsk_cut_graph_atlas_t, &
        map_eqdsk_cut_graph_atlas, map_eqdsk_cut_graph_atlas_flux
    use neort_gc_eqdsk_cut_jet, only: EQDSK_CUT_JET_SUCCESS, &
        eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
    use neort_gc_eqdsk_flux_profile_map, only: EQDSK_FLUX_MAP_CERTIFICATE_ID
    implicit none

    type(eqdsk_cylindrical_field_t) :: field
    type(eqdsk_cut_graph_atlas_t) :: outboard_atlas
    type(eqdsk_cut_graph_atlas_options_t) :: atlas_options
    type(eqdsk_cut_endpoint_options_t) :: endpoint_options
    type(eqdsk_cut_endpoint_certificate_t) :: certificate, rejected
    type(eqdsk_cut_jet_t) :: jet
    character(len=1024) :: path
    real(dp) :: axis_R, outboard_lo, outboard_hi
    real(dp) :: endpoint_flux(2), target, seed(3), tangent(3)
    real(dp) :: half_R, half_Z, dZ_dR, dpsihat_dR
    integer :: status, jet_status

    call get_environment_variable('EQDSK_FILE', path)
    call require(len_trim(path) > 0, 'EQDSK_FILE is required')
    call initialize_eqdsk_cylindrical_field(trim(path), 1.0_dp, field, status)
    call require(status == 0, 'failed to initialize circular EQDSK')
    call require(nrad >= 3 .and. nzet >= 3, 'EQDSK grid is too small')

    axis_R = 0.5_dp*(field%domain_R_min+field%domain_R_max)
    outboard_lo = axis_R+0.10_dp*(field%domain_R_max-field%domain_R_min)
    outboard_hi = field%domain_R_max-0.10_dp*(field%domain_R_max- &
        field%domain_R_min)
    atlas_options%max_r_depth = 12
    atlas_options%max_z_depth = 12
    atlas_options%map_absolute_tolerance = 1.0e-12_dp
    atlas_options%map_relative_tolerance = 1.0e-13_dp
    call build_eqdsk_cut_graph_atlas(outboard_atlas, outboard_lo, &
        outboard_hi, -10.0_dp, 10.0_dp, 0.0_dp, 1.0_dp, atlas_options, &
        status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'outboard seed atlas was not certified')

    call endpoint_psihat(outboard_atlas, outboard_lo, endpoint_flux(1))
    call endpoint_psihat(outboard_atlas, outboard_hi, endpoint_flux(2))
    target = 0.5_dp*sum(endpoint_flux)
    call map_eqdsk_cut_graph_atlas_flux(outboard_atlas, target, seed, &
        tangent, dZ_dR, dpsihat_dR, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'failed to construct endpoint seed')

    half_R = 0.20_dp*minval(rad(2:nrad)-rad(1:nrad-1))
    half_Z = 0.20_dp*minval(zet(2:nzet)-zet(1:nzet-1))
    endpoint_options%max_newton_iterations = 8
    call build_eqdsk_cut_endpoint_certificate(seed(1)-half_R, &
        seed(1)+half_R, seed(2)-half_Z, seed(2)+half_Z, target, 1.0_dp, &
        seed(1), seed(2), endpoint_options, certificate, status)
    call require(status == EQDSK_ENDPOINT_CERT_SUCCESS, &
        'regular circular endpoint was not Krawczyk-certified')
    call require(certificate%cell_tiles_covered >= 2, &
        'endpoint certificate did not cover the grid seam')
    call require(certificate%regular_domain_certified .and. &
        certificate%strict_inclusion_certified .and. &
        certificate%unique_endpoint_certified, &
        'endpoint certificate omitted a theorem gate')
    call require(certificate%krawczyk_Z_lo <= 0.0_dp .and. &
        certificate%krawczyk_Z_hi >= 0.0_dp, &
        'circular endpoint enclosure excluded the exact symmetry midplane')
    call validate_eqdsk_cut_endpoint_certificate(certificate, status)
    call require(status == EQDSK_ENDPOINT_CERT_SUCCESS, &
        'fresh endpoint certificate failed validation')
    call require(EQDSK_ENDPOINT_CERTIFICATE_ID /= &
        EQDSK_FLUX_MAP_CERTIFICATE_ID, &
        'endpoint and flux-map certificate IDs collide')
    call require(certificate%certificate_id == EQDSK_ENDPOINT_CERTIFICATE_ID, &
        'endpoint certificate used the wrong named ID')

    call evaluate_eqdsk_cut_jet([certificate%newton_point_R, &
        certificate%newton_point_Z, 0.0_dp], 1.0_dp, 1, &
        [0.0_dp, 0.0_dp, 0.0_dp], jet, jet_status)
    call require(jet_status == EQDSK_CUT_JET_SUCCESS, &
        'certified endpoint point evaluation failed')
    call require(abs(jet%psi_jet(1)/psi_sep-target) <= 1.0e-11_dp, &
        'endpoint Newton point missed the target flux')

    call build_eqdsk_cut_endpoint_certificate(outboard_lo, outboard_hi, &
        -half_Z, half_Z, target, -1.0_dp, seed(1), seed(2), &
        endpoint_options, rejected, status)
    call require(status /= EQDSK_ENDPOINT_CERT_SUCCESS, &
        'negative field scale passed the regular-domain gate')
    call build_eqdsk_cut_endpoint_certificate(seed(1)-half_R, &
        seed(1)+half_R, seed(2)-half_Z, seed(2)+half_Z, 0.0_dp, 1.0_dp, &
        seed(1), seed(2), endpoint_options, rejected, status)
    call require(status == EQDSK_ENDPOINT_CERT_AXIS_REQUIRES_LIMIT, &
        'magnetic axis was accepted as a regular endpoint')

    ! This box spans both the inboard and outboard intersections at target.
    ! A one-root certificate must fail rather than selecting the supplied seed.
    call build_eqdsk_cut_endpoint_certificate( &
        field%domain_R_min+0.10_dp*(field%domain_R_max-field%domain_R_min), &
        outboard_hi, -half_Z, half_Z, target, 1.0_dp, seed(1), seed(2), &
        endpoint_options, rejected, status)
    call require(status /= EQDSK_ENDPOINT_CERT_SUCCESS, &
        'two-endpoint box was incorrectly certified unique')
    call require(.not. rejected%unique_endpoint_certified, &
        'failed two-endpoint box retained a uniqueness claim')

    write (*, '(a)') 'test_gc_eqdsk_cut_endpoint_certificate OK'

contains

    subroutine endpoint_psihat(atlas, radius, psihat)
        type(eqdsk_cut_graph_atlas_t), intent(in) :: atlas
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: psihat

        real(dp) :: position(3), local_tangent(3)
        type(eqdsk_cut_jet_t) :: local_jet
        integer :: local_status

        call map_eqdsk_cut_graph_atlas(atlas, radius, position, &
            local_tangent, status=local_status)
        call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
            'failed to map branch endpoint')
        call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], local_jet, local_status)
        call require(local_status == EQDSK_CUT_JET_SUCCESS, &
            'failed to evaluate branch endpoint flux')
        psihat = local_jet%psi_jet(1)/psi_sep
    end subroutine endpoint_psihat

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_cut_endpoint_certificate
