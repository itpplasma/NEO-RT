program test_gc_eqdsk_allowed_region_cut_box
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: nrad, nzet, psi_sep, rad, zet
    use neort_eqdsk_physical_flux_map_symbolic, only: &
        evaluate_neort_eqdsk_physical_flux_map
    use neort_gc_eqdsk_allowed_region_cut_box, only: &
        EQDSK_CUT_BOX_INVALID_ATLAS, EQDSK_CUT_BOX_INVALID_PROFILE, &
        EQDSK_CUT_BOX_NOT_OUTBOARD, EQDSK_CUT_BOX_OUT_OF_RANGE, &
        EQDSK_CUT_BOX_SUCCESS, &
        eqdsk_allowed_region_cut_provenance_t, &
        eqdsk_potential_profile_nodes_t, &
        evaluate_eqdsk_allowed_region_cut_box, &
        validate_eqdsk_potential_profile_nodes
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        eqdsk_cylindrical_field_t, initialize_eqdsk_cylindrical_field
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_ATLAS_SUCCESS, build_eqdsk_cut_graph_atlas, &
        enclose_eqdsk_cut_graph_strip, eqdsk_cut_graph_atlas_options_t, &
        eqdsk_cut_graph_atlas_t, map_eqdsk_cut_graph_atlas
    use neort_gc_eqdsk_cut_jet, only: EQDSK_CUT_JET_SUCCESS, &
        eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
    use neort_gc_eqdsk_cut_interval, only: eqdsk_cut_interval_result_t
    use neort_gc_eqdsk_allowed_region_interval, only: &
        eqdsk_allowed_interval_result_t
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_t
    implicit none

    type(eqdsk_cylindrical_field_t) :: field
    type(eqdsk_cut_graph_atlas_t) :: atlas, inboard_atlas, mutated_atlas
    type(eqdsk_cut_graph_atlas_options_t) :: options
    type(eqdsk_potential_profile_nodes_t) :: profile, bad_profile
    type(eqdsk_allowed_interval_result_t) :: result
    type(eqdsk_allowed_region_cut_provenance_t) :: provenance
    type(gc_outward_interval_t) :: query
    character(len=1024) :: path
    real(dp) :: axis_R, outboard_lo, outboard_hi, inboard_lo, inboard_hi
    real(dp) :: radius, position(3), tangent(3), scalar_psi
    real(dp) :: scalar_dpsi, query_width
    real(dp) :: profile_psi(2), profile_phi(2), profile_omega(2)
    type(eqdsk_cut_interval_result_t), allocatable :: enclosures(:)
    type(eqdsk_cut_jet_t) :: jet
    integer :: i, expected_leaf_count, status, local_status
    integer :: axis_R_index, axis_Z_index

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

    outboard_lo = axis_R+0.10_dp*(field%domain_R_max-axis_R)
    outboard_hi = field%domain_R_max-0.10_dp*( &
        field%domain_R_max-field%domain_R_min)
    inboard_lo = field%domain_R_min+0.10_dp*( &
        field%domain_R_max-field%domain_R_min)
    inboard_hi = axis_R-0.10_dp*(field%domain_R_max-axis_R)
    call build_eqdsk_cut_graph_atlas(atlas, outboard_lo, outboard_hi, &
        -10.0_dp, 10.0_dp, 0.0_dp, 1.0_dp, options, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'outboard graph atlas was not certified')
    call require(size(atlas%strips) > 1, &
        'fixture does not exercise multistrip coverage')

    profile_psi = [-psi_sep, 2.0_dp*psi_sep]
    profile_phi = [0.0_dp, 0.0_dp]
    profile_omega = [0.0_dp, 0.0_dp]
    call validate_eqdsk_potential_profile_nodes(profile_psi, profile_phi, &
        profile_omega, profile, status)
    call require(status == EQDSK_CUT_BOX_SUCCESS, &
        'potential profile nodes were not structurally validated')
    call require(profile%structurally_validated, &
        'profile-node structural validation flag was not set')

    query = gc_outward_interval(outboard_lo, outboard_hi)
    call evaluate_eqdsk_allowed_region_cut_box(atlas, query, 1.0_dp, &
        psi_sep, profile, 1.0e6_dp, 0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp, 1, &
        result, provenance, status)
    call require(status == EQDSK_CUT_BOX_SUCCESS, &
        'certified outboard cut-box evaluation failed')
    call require(provenance%certified, 'cut-box result was not certified')
    call require(provenance%nstrips > 1, &
        'cut-box provenance collapsed a multistrip request')
    call require(provenance%n_graph_enclosures >= provenance%nstrips, &
        'cut-box provenance retained no graph leaves')
    call require(size(provenance%leaf_enclosures) == &
        provenance%n_graph_enclosures, 'leaf enclosure count was not retained')
    call require(size(provenance%leaf_strip_index) == &
        provenance%n_graph_enclosures, 'leaf strip IDs were not retained')
    call require(provenance%profile_inputs_structurally_validated, &
        'profile input validation was not recorded separately')
    ! This ID identifies the generated interpolation kernel only; it does
    ! not certify the supplied profile-node values.
    call require(provenance%generated_profile_interpolation_certificate_id > 0, &
        'generated profile interpolation identity was not recorded')
    call require(len_trim( &
        provenance%generated_profile_interpolation_fingerprint) > 0, &
        'generated profile interpolation fingerprint was not recorded')
    call check_containing_strip_queries(atlas, provenance)

    expected_leaf_count = 0
    do i = 1, size(atlas%strips)
        call enclose_eqdsk_cut_graph_strip(atlas, i, &
            atlas%strips(i)%r_lo, atlas%strips(i)%r_hi, enclosures, &
            local_status)
        call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
            'independent strip enclosure failed')
        expected_leaf_count = expected_leaf_count+size(enclosures)
    end do
    call require(expected_leaf_count == provenance%n_graph_enclosures, &
        'cut-box did not retain every graph leaf')

    do i = 1, size(atlas%strips)
        radius = 0.5_dp*(atlas%strips(i)%r_lo+atlas%strips(i)%r_hi)
        call map_eqdsk_cut_graph_atlas(atlas, radius, position, tangent, &
            status=local_status)
        call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
            'subbox scalar graph map failed')
        call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], jet, local_status)
        call require(local_status == EQDSK_CUT_JET_SUCCESS, &
            'subbox scalar jet failed')
        call require(any_leaf_contains(provenance, i, jet), &
            'retained graph leaves excluded a scalar subbox oracle')
    end do

    query_width = outboard_hi-outboard_lo
    do i = 1, 2
        radius = outboard_lo+query_width*real(i-1,dp)
        call map_eqdsk_cut_graph_atlas(atlas, radius, position, tangent, &
            status=local_status)
        call require(local_status == EQDSK_CUT_ATLAS_SUCCESS, &
            'endpoint graph map failed')
        call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], jet, local_status)
        call require(local_status == EQDSK_CUT_JET_SUCCESS, &
            'endpoint scalar jet failed')
        call evaluate_neort_eqdsk_physical_flux_map(jet%psi_jet(1), &
            jet%psi_jet(2), jet%psi_jet(3), tangent(2), 1.0_dp, &
            scalar_psi, scalar_dpsi)
        call require(result%psi_physical%lo <= scalar_psi .and. &
                result%psi_physical%hi >= scalar_psi, &
            'allowed cut box excluded an independent endpoint flux oracle')
    end do

    mutated_atlas = atlas
    mutated_atlas%raw_psi_sep = 1.1_dp*mutated_atlas%raw_psi_sep
    call evaluate_eqdsk_allowed_region_cut_box(mutated_atlas, query, 1.0_dp, &
        psi_sep, profile, 1.0e6_dp, 0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp, 1, &
        result, provenance, status)
    call require(status == EQDSK_CUT_BOX_INVALID_ATLAS, &
        'mutated atlas provenance was accepted')

    call evaluate_eqdsk_allowed_region_cut_box(atlas, &
        gc_outward_interval(outboard_lo-1.0e-6_dp, outboard_hi), 1.0_dp, &
        psi_sep, profile, 1.0e6_dp, 0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp, 1, &
        result, provenance, status)
    call require(status == EQDSK_CUT_BOX_OUT_OF_RANGE, &
        'out-of-domain lower R box was accepted')

    call build_eqdsk_cut_graph_atlas(inboard_atlas, inboard_lo, inboard_hi, &
        -10.0_dp, 10.0_dp, 0.0_dp, 1.0_dp, options, status)
    call require(status == EQDSK_CUT_ATLAS_SUCCESS, &
        'inboard graph atlas was not certified')
    call evaluate_eqdsk_allowed_region_cut_box(inboard_atlas, &
        gc_outward_interval(inboard_lo, inboard_hi), 1.0_dp, psi_sep, &
        profile, 1.0e6_dp, 0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp, 1, result, &
        provenance, status)
    call require(status == EQDSK_CUT_BOX_NOT_OUTBOARD, &
        'inboard decreasing graph was accepted as outboard')

    bad_profile = profile
    bad_profile%psi(2) = bad_profile%psi(1)
    call evaluate_eqdsk_allowed_region_cut_box(atlas, query, 1.0_dp, &
        psi_sep, bad_profile, 1.0e6_dp, 0.0_dp, 2.0_dp, -1.0_dp, 1.0_dp, 1, &
        result, provenance, status)
    call require(status == EQDSK_CUT_BOX_INVALID_PROFILE, &
        'mutated potential profile was accepted')

    write (*, '(a)') 'test_gc_eqdsk_allowed_region_cut_box: PASS'

contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

    subroutine check_containing_strip_queries(atlas, provenance)
        type(eqdsk_cut_graph_atlas_t), intent(in) :: atlas
        type(eqdsk_allowed_region_cut_provenance_t), intent(in) :: provenance

        integer :: i, j, strip_position
        logical :: found

        ! These arrays describe the strip query supplied to enclose_..., not
        ! exact internal leaf boxes, which that API does not expose.
        do i = 1, provenance%n_graph_enclosures
            found = .false.
            strip_position = 0
            do j = 1, provenance%nstrips
                if (provenance%strip_index(j) == &
                        provenance%leaf_strip_index(i)) then
                    found = .true.
                    strip_position = j
                    exit
                end if
            end do
            call require(found, 'leaf strip provenance lost its strip query')
            call require(provenance%leaf_strip_r_lo(i) == &
                    provenance%strip_r_lo(strip_position), &
                'leaf provenance overclaimed its containing R query')
            call require(provenance%leaf_strip_r_hi(i) == &
                    provenance%strip_r_hi(strip_position), &
                'leaf provenance overclaimed its containing R query')
            call require(provenance%leaf_strip_z_lo(i) == &
                    atlas%strips(provenance%leaf_strip_index(i))%z_lo, &
                'leaf provenance changed its containing Z query')
            call require(provenance%leaf_strip_z_hi(i) == &
                    atlas%strips(provenance%leaf_strip_index(i))%z_hi, &
                'leaf provenance changed its containing Z query')
        end do
    end subroutine check_containing_strip_queries

    logical function any_leaf_contains(provenance, strip_index, jet)
        type(eqdsk_allowed_region_cut_provenance_t), intent(in) :: provenance
        integer, intent(in) :: strip_index
        type(eqdsk_cut_jet_t), intent(in) :: jet

        integer :: i

        any_leaf_contains = .false.
        do i = 1, provenance%n_graph_enclosures
            if (provenance%leaf_strip_index(i) /= strip_index) cycle
            if (provenance%leaf_enclosures(i)%psi%lo > &
                    jet%psi_jet(1)) cycle
            if (provenance%leaf_enclosures(i)%psi%hi < &
                    jet%psi_jet(1)) cycle
            if (provenance%leaf_enclosures(i)%psi_R%lo > &
                    jet%psi_jet(2)) cycle
            if (provenance%leaf_enclosures(i)%psi_R%hi < &
                    jet%psi_jet(2)) cycle
            if (provenance%leaf_enclosures(i)%psi_Z%lo > &
                    jet%psi_jet(3)) cycle
            if (provenance%leaf_enclosures(i)%psi_Z%hi < &
                    jet%psi_jet(3)) cycle
            any_leaf_contains = .true.
            return
        end do
    end function any_leaf_contains

end program test_gc_eqdsk_allowed_region_cut_box
