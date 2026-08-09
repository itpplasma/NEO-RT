program test_gc_eqdsk_cut_graph_atlas
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: nrad, nzet, psi_sep, rad, zet
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        eqdsk_cylindrical_field_t, initialize_eqdsk_cylindrical_field, &
        map_eqdsk_flux_position
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_ATLAS_INVALID_CERTIFICATE, &
        EQDSK_CUT_ATLAS_SUCCESS, EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
        build_eqdsk_cut_graph_atlas, clear_eqdsk_cut_graph_atlas, &
        enclose_eqdsk_cut_graph_strip, &
        eqdsk_cut_graph_atlas_options_t, eqdsk_cut_graph_atlas_t, &
        map_eqdsk_cut_graph_atlas, map_eqdsk_cut_graph_atlas_flux, &
        validate_eqdsk_cut_graph_atlas
    use neort_gc_eqdsk_cut_jet, only: EQDSK_CUT_JET_SUCCESS, &
        eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
    use neort_gc_eqdsk_cut_interval, only: eqdsk_cut_interval_result_t
    use neort_gc_eqdsk_flux_profile_map, only: &
        EQDSK_FLUX_MAP_SUCCESS, eqdsk_flux_profile_map_t, &
        initialize_eqdsk_flux_profile_map
    use neort_gc_eqdsk_radial_interval_map, only: &
        EQDSK_RADIAL_INTERVAL_INVALID_ATLAS, &
        EQDSK_RADIAL_INTERVAL_NOT_MONOTONE, &
        EQDSK_RADIAL_INTERVAL_NORMALIZATION, &
        EQDSK_RADIAL_INTERVAL_OUT_OF_RANGE, &
        EQDSK_RADIAL_INTERVAL_SUCCESS, &
        eqdsk_rho_interval_provenance_t, &
        map_eqdsk_outboard_r_interval_to_rho_tor
    use neort_gc_outward_interval, only: gc_outward_interval_t
    use neort_gc_eqdsk_axis_certificate, only: &
        EQDSK_AXIS_CERT_SUCCESS, build_eqdsk_axis_certificate, &
        eqdsk_axis_certificate_t, validate_eqdsk_axis_certificate
    use util, only: pi
    implicit none

    type(eqdsk_cylindrical_field_t) :: field
    type(eqdsk_cut_graph_atlas_t) :: inboard_atlas, outboard_atlas, full_atlas
    type(eqdsk_cut_graph_atlas_t) :: mutated_atlas
    type(eqdsk_cut_graph_atlas_options_t) :: options, full_options
    type(eqdsk_axis_certificate_t) :: axis_certificate
    type(eqdsk_flux_profile_map_t) :: flux_map, unit_scale_flux_map
    type(eqdsk_flux_profile_map_t) :: mismatched_flux_map
    type(eqdsk_rho_interval_provenance_t) :: radial_provenance
    type(eqdsk_rho_interval_provenance_t) :: unit_scale_provenance
    type(eqdsk_rho_interval_provenance_t) :: mismatched_provenance
    type(gc_outward_interval_t) :: rho_interval
    type(gc_outward_interval_t) :: unit_scale_rho_interval
    type(gc_outward_interval_t) :: mismatched_rho_interval
    character(len=1024) :: path
    real(dp) :: axis_R, inboard_lo, inboard_hi, outboard_lo, outboard_hi
    real(dp), parameter :: branch_surface_lo = 0.25_dp
    real(dp), parameter :: branch_surface_hi = 0.64_dp
    real(dp) :: flux_position(3)
    real(dp) :: mismatched_psi_sep
    integer :: status, axis_R_index, axis_Z_index

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
    axis_R_index = minloc(abs(rad-axis_R), dim=1)
    axis_Z_index = minloc(abs(zet), dim=1)
    call require(axis_R_index > 1 .and. axis_R_index < nrad, &
        'circular magnetic axis is not interior to the R grid')
    call require(axis_Z_index > 1 .and. axis_Z_index < nzet, &
        'circular magnetic axis is not interior to the Z grid')
    call build_eqdsk_axis_certificate(rad(axis_R_index-1), &
        rad(axis_R_index+1), zet(axis_Z_index-1), zet(axis_Z_index+1), &
        axis_certificate, status)
    call require(status == EQDSK_AXIS_CERT_SUCCESS, &
        'nondegenerate circular magnetic axis was not certified')
    call validate_eqdsk_axis_certificate(axis_certificate, status)
    call require(status == EQDSK_AXIS_CERT_SUCCESS, &
        'fresh magnetic-axis certificate failed validation')
    call map_eqdsk_flux_position(branch_surface_hi, pi, 0.0_dp, &
        flux_position, status)
    call require(status == 0, 'failed to map inboard branch surface')
    inboard_lo = flux_position(1)
    call map_eqdsk_flux_position(branch_surface_lo, pi, 0.0_dp, &
        flux_position, status)
    call require(status == 0, 'failed to map inboard branch surface')
    inboard_hi = flux_position(1)
    call map_eqdsk_flux_position(branch_surface_lo, 0.0_dp, 0.0_dp, &
        flux_position, status)
    call require(status == 0, 'failed to map outboard branch surface')
    outboard_lo = flux_position(1)
    call map_eqdsk_flux_position(branch_surface_hi, 0.0_dp, 0.0_dp, &
        flux_position, status)
    call require(status == 0, 'failed to map outboard branch surface')
    outboard_hi = flux_position(1)
    call build_eqdsk_cut_graph_atlas(inboard_atlas, inboard_lo, inboard_hi, &
        -10.0_dp, 10.0_dp, 0.0_dp, 1.0_dp, options, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'regular inboard circular cut branch was not certified')
    call require(inboard_atlas%flux_monotonicity_certified .and. &
        inboard_atlas%flux_monotonicity_sign == -1, &
        'inboard cut flux was not certified strictly decreasing')
    call check_regular_branch(inboard_atlas, inboard_lo, inboard_hi)
    call build_eqdsk_cut_graph_atlas(outboard_atlas, outboard_lo, &
        outboard_hi, -10.0_dp, 10.0_dp, 0.0_dp, 1.0_dp, options, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'regular outboard circular cut branch was not certified')
    call require(outboard_atlas%flux_monotonicity_certified .and. &
        outboard_atlas%flux_monotonicity_sign == 1, &
        'outboard cut flux was not certified strictly increasing')
    call require(outboard_atlas%raw_psi_sep_valid, &
        'outboard atlas omitted raw separatrix provenance validity')
    call require(outboard_atlas%raw_psi_sep == psi_sep, &
        'outboard atlas did not retain the raw separatrix flux')
    call check_regular_branch(outboard_atlas, outboard_lo, outboard_hi)
    call initialize_eqdsk_flux_profile_map([0.0_dp, 1.0_dp], &
        [0.0_dp, 3.0_dp*psi_sep], 3.0_dp, psi_sep, flux_map, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'scaled linear flux profile was not certified')
    call check_radial_interval(outboard_atlas, flux_map, &
        outboard_lo+0.3_dp*(outboard_hi-outboard_lo), &
        outboard_lo+0.7_dp*(outboard_hi-outboard_lo), rho_interval, &
        radial_provenance)
    call initialize_eqdsk_flux_profile_map([0.0_dp, 1.0_dp], &
        [0.0_dp, psi_sep], 1.0_dp, psi_sep, unit_scale_flux_map, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'unit-scale flux profile was not certified')
    call map_eqdsk_outboard_r_interval_to_rho_tor(outboard_atlas, &
        unit_scale_flux_map, outboard_lo+0.3_dp*(outboard_hi-outboard_lo), &
        outboard_lo+0.7_dp*(outboard_hi-outboard_lo), &
        unit_scale_rho_interval, unit_scale_provenance, status)
    call require(status == EQDSK_RADIAL_INTERVAL_SUCCESS, &
        'unit-scale radial map was rejected')
    call require(abs(midpoint(rho_interval)-midpoint(unit_scale_rho_interval)) &
        <= 1.0e-10_dp*max(1.0_dp, abs(midpoint(rho_interval))), &
        'field-scale cancellation changed the normalized radial map')

    mutated_atlas = outboard_atlas
    mutated_atlas%raw_psi_sep = 1.1_dp*mutated_atlas%raw_psi_sep
    call validate_eqdsk_cut_graph_atlas(mutated_atlas, status)
    call require(status == EQDSK_CUT_ATLAS_INVALID_CERTIFICATE, &
        'mutated raw atlas separatrix provenance was accepted')
    call map_eqdsk_outboard_r_interval_to_rho_tor(mutated_atlas, flux_map, &
        outboard_lo, outboard_hi, mismatched_rho_interval, &
        mismatched_provenance, status)
    call require(status == EQDSK_RADIAL_INTERVAL_INVALID_ATLAS, &
        'radial map accepted a mutated atlas separatrix provenance')

    mismatched_psi_sep = 1.1_dp*psi_sep
    call initialize_eqdsk_flux_profile_map([0.0_dp, 1.0_dp], &
        [0.0_dp, 3.0_dp*mismatched_psi_sep], 3.0_dp, mismatched_psi_sep, &
        mismatched_flux_map, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'mismatched-separatrix profile was not independently certified')
    call map_eqdsk_outboard_r_interval_to_rho_tor(outboard_atlas, &
        mismatched_flux_map, outboard_lo, outboard_hi, &
        mismatched_rho_interval, mismatched_provenance, status)
    call require(status == EQDSK_RADIAL_INTERVAL_NORMALIZATION, &
        'raw psi_sep/profile provenance mismatch was accepted')

    call require(size(outboard_atlas%strips) > 1, &
        'fixture does not provide a multistrip outboard interval')
    call check_radial_interval(outboard_atlas, flux_map, &
        outboard_atlas%strips(1)%r_lo, &
        outboard_atlas%strips(size(outboard_atlas%strips))%r_hi, rho_interval, &
        radial_provenance)
    call require(radial_provenance%nstrips > 1, &
        'multistrip radial interval was reduced to one strip')

    call map_eqdsk_outboard_r_interval_to_rho_tor(outboard_atlas, flux_map, &
        outboard_lo-1.0e-6_dp, outboard_hi, rho_interval, radial_provenance, &
        status)
    call require(status == EQDSK_RADIAL_INTERVAL_OUT_OF_RANGE, &
        'R box below the certified domain was accepted')
    call map_eqdsk_outboard_r_interval_to_rho_tor(outboard_atlas, flux_map, &
        outboard_lo, outboard_hi+1.0e-6_dp, rho_interval, radial_provenance, &
        status)
    call require(status == EQDSK_RADIAL_INTERVAL_OUT_OF_RANGE, &
        'R box above the certified domain was accepted')

    call map_eqdsk_outboard_r_interval_to_rho_tor(inboard_atlas, flux_map, &
        inboard_lo, inboard_hi, rho_interval, radial_provenance, status)
    call require(status == EQDSK_RADIAL_INTERVAL_NOT_MONOTONE, &
        'decreasing inboard graph was accepted as the outboard rho map')

    ! The rectangular EQDSK box includes the SOL, whereas this atlas was
    ! requested for the closed-flux interval psi_hat in [0,1].  It must fail
    ! closed rather than inventing a graph through the unrequested SOL.
    full_options = options
    full_options%max_r_depth = 6
    full_options%max_z_depth = 6
    call build_eqdsk_cut_graph_atlas(full_atlas, field%domain_R_min, &
        field%domain_R_max, zet(1), zet(nzet), 0.0_dp, 1.0_dp, full_options, &
        status)
    call require(status /= EQDSK_CUT_ATLAS_SUCCESS, &
        'SOL-containing closed-flux atlas was incorrectly certified')
    call require(.not. full_atlas%global_completeness_certified, &
        'failed SOL-containing atlas retained a completeness certificate')

    call clear_eqdsk_cut_graph_atlas(inboard_atlas)
    call require(.not. inboard_atlas%raw_psi_sep_valid .and. &
        inboard_atlas%raw_psi_sep == 0.0_dp, &
        'clearing atlas did not remove raw separatrix provenance')
    call validate_eqdsk_cut_graph_atlas(inboard_atlas, status)
    call require(status /= EQDSK_CUT_ATLAS_SUCCESS, &
        'cleared atlas remained valid')
    write (*, '(a)') 'test_gc_eqdsk_cut_graph_atlas OK'

contains

    subroutine check_regular_branch(atlas, r_lo, r_hi)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        real(dp), intent(in) :: r_lo, r_hi

        real(dp) :: radius, position(3), dposition(3), dZ_dR, dpsihat_dR
        real(dp) :: endpoint_position(3), endpoint_tangent(3)
        real(dp) :: scale, endpoint_flux(2), target_flux
        type(eqdsk_cut_jet_t) :: jet
        type(eqdsk_cut_interval_result_t), allocatable :: enclosures(:)
        integer :: i, local_status, jet_status, strip_index

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
        do strip_index = 1, size(atlas%strips)
            radius = 0.5_dp*(atlas%strips(strip_index)%r_lo + &
                atlas%strips(strip_index)%r_hi)
            call enclose_eqdsk_cut_graph_strip(atlas, strip_index, radius, &
                radius, enclosures, local_status)
            call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
                'certified graph point enclosure failed')
            call require(allocated(enclosures), &
                'certified graph point enclosure was not allocated')
            call require(size(enclosures) > 0, &
                'certified graph point enclosure was empty')
            call map_eqdsk_cut_graph_atlas(atlas, radius, position, &
                dposition, status=local_status)
            call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
                'point-enclosure scalar map failed')
            call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, &
                [0.0_dp, 0.0_dp, 0.0_dp], jet, jet_status)
            call require(jet_status == EQDSK_CUT_JET_SUCCESS, &
                'point-enclosure scalar jet failed')
            call require(any_enclosure_contains(enclosures, jet), &
                'interval graph cover excluded its analytic midplane root')
        end do
        call map_eqdsk_cut_graph_atlas(atlas, r_lo, endpoint_position, &
            endpoint_tangent, status=local_status)
        call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
            'lower branch endpoint map failed')
        call evaluate_eqdsk_cut_jet(endpoint_position, 1.0_dp, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], jet, jet_status)
        call require(jet_status == EQDSK_CUT_JET_SUCCESS, &
            'lower branch endpoint flux failed')
        endpoint_flux(1) = jet%psi_jet(1)/psi_sep
        call map_eqdsk_cut_graph_atlas(atlas, r_hi, endpoint_position, &
            endpoint_tangent, status=local_status)
        call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
            'upper branch endpoint map failed')
        call evaluate_eqdsk_cut_jet(endpoint_position, 1.0_dp, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], jet, jet_status)
        call require(jet_status == EQDSK_CUT_JET_SUCCESS, &
            'upper branch endpoint flux failed')
        endpoint_flux(2) = jet%psi_jet(1)/psi_sep
        target_flux = 0.5_dp*sum(endpoint_flux)
        call map_eqdsk_cut_graph_atlas_flux(atlas, target_flux, position, &
            dposition, dZ_dR, dpsihat_dR, local_status)
        call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
            'inverse certified branch map failed')
        call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], jet, jet_status)
        call require(jet_status == EQDSK_CUT_JET_SUCCESS, &
            'inverse branch flux evaluation failed')
        call require(abs(jet%psi_jet(1)/psi_sep-target_flux) <= &
            2.0e-10_dp*max(1.0_dp, abs(target_flux)), &
            'inverse branch map missed its normalized-flux target')
    end subroutine check_regular_branch

    subroutine check_radial_interval(atlas, map, r_lo, r_hi, rho_box, &
            provenance)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        type(eqdsk_flux_profile_map_t), intent(in) :: map
        real(dp), intent(in) :: r_lo, r_hi
        type(gc_outward_interval_t), intent(out) :: rho_box
        type(eqdsk_rho_interval_provenance_t), intent(out) :: provenance

        type(eqdsk_cut_jet_t) :: jet
        real(dp) :: position(3), tangent(3), expected_rho(2)
        integer :: i, local_status, jet_status

        call map_eqdsk_outboard_r_interval_to_rho_tor(atlas, map, r_lo, &
            r_hi, rho_box, provenance, local_status)
        call require(local_status == EQDSK_RADIAL_INTERVAL_SUCCESS, &
            'regular outboard R interval did not map to rho_tor')
        call require(provenance%mapping_certified, &
            'radial interval omitted its completed certificate')
        call require(provenance%nstrips >= 1, &
            'radial interval retained no strip provenance')
        call require(provenance%n_graph_enclosures >= provenance%nstrips, &
            'radial interval lost graph enclosure provenance')
        do i = 1, 2
            call map_eqdsk_cut_graph_atlas(atlas, merge(r_lo, r_hi, i == 1), &
                position, tangent, status=local_status)
            call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
                'radial interval endpoint scalar map failed')
            call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, &
                [0.0_dp, 0.0_dp, 0.0_dp], jet, jet_status)
            call require(jet_status == EQDSK_CUT_JET_SUCCESS, &
                'radial interval endpoint scalar jet failed')
            expected_rho(i) = sqrt(jet%psi_jet(1)/psi_sep)
            call require(rho_box%lo <= expected_rho(i) .and. &
                    rho_box%hi >= expected_rho(i), &
                'radial interval excluded its scalar endpoint oracle')
            call require(provenance%psi_physical%lo <= &
                    3.0_dp*jet%psi_jet(1) .and. &
                    provenance%psi_physical%hi >= &
                    3.0_dp*jet%psi_jet(1), &
                'physical flux did not retain exactly one field scale')
        end do
    end subroutine check_radial_interval

    pure real(dp) function midpoint(interval)
        type(gc_outward_interval_t), intent(in) :: interval

        midpoint = 0.5_dp*(interval%lo+interval%hi)
    end function midpoint

    logical function any_enclosure_contains(enclosures, jet)
        type(eqdsk_cut_interval_result_t), intent(in) :: enclosures(:)
        type(eqdsk_cut_jet_t), intent(in) :: jet

        integer :: i

        any_enclosure_contains = .false.
        do i = 1, size(enclosures)
            if (enclosures(i)%psi%lo > jet%psi_jet(1)) cycle
            if (enclosures(i)%psi%hi < jet%psi_jet(1)) cycle
            if (enclosures(i)%psi_R%lo > jet%psi_jet(2)) cycle
            if (enclosures(i)%psi_R%hi < jet%psi_jet(2)) cycle
            if (enclosures(i)%psi_Z%lo > jet%psi_jet(3)) cycle
            if (enclosures(i)%psi_Z%hi < jet%psi_jet(3)) cycle
            if (enclosures(i)%numerator%lo > 0.0_dp) cycle
            if (enclosures(i)%numerator%hi < 0.0_dp) cycle
            any_enclosure_contains = .true.
            return
        end do
    end function any_enclosure_contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_cut_graph_atlas
