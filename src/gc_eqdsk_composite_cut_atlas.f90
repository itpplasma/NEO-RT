module neort_gc_eqdsk_composite_cut_atlas
    !! Composite geometric certificate for the complete interior Eq.13 cut.
    !!
    !! The component certificates retain distinct theorem meanings:
    !! graph atlases prove zero-set completeness, the axis certificate proves
    !! one nondegenerate stationary point and its limiting chart, and endpoint
    !! Krawczyk boxes prove transverse intersections with the requested outer
    !! flux surface.  This module proves their overlap/connectivity.  It does
    !! not claim Buchholz orbit-crossing multiplicity or monotonic B on arcs.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: psi_sep
    use neort_eqdsk_cut_axis_limit_symbolic, only: &
        evaluate_neort_eqdsk_cut_axis_limit
    use neort_eqdsk_cut_axis_rho_limit_symbolic, only: &
        evaluate_neort_eqdsk_cut_axis_rho_limit
    use neort_eqdsk_cut_flux_coordinate_symbolic, only: &
        evaluate_neort_eqdsk_cut_flux_coordinate
    use neort_eqdsk_cut_r_flux_chart_symbolic, only: &
        evaluate_neort_eqdsk_cut_r_flux_chart
    use neort_eqdsk_rho_tor_map_symbolic, only: &
        evaluate_neort_eqdsk_rho_tor_map
    use neort_gc_eqdsk_axis_certificate, only: &
        EQDSK_AXIS_CERT_SUCCESS, build_eqdsk_axis_certificate, &
        eqdsk_axis_certificate_t, validate_eqdsk_axis_certificate
    use neort_gc_eqdsk_cut_endpoint_certificate, only: &
        EQDSK_ENDPOINT_CERT_SUCCESS, build_eqdsk_cut_endpoint_certificate, &
        eqdsk_cut_endpoint_certificate_t, eqdsk_cut_endpoint_options_t, &
        validate_eqdsk_cut_endpoint_certificate
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_ATLAS_SUCCESS, build_eqdsk_cut_graph_atlas, &
        eqdsk_cut_graph_atlas_options_t, eqdsk_cut_graph_atlas_t, &
        map_eqdsk_cut_graph_atlas_flux, validate_eqdsk_cut_graph_atlas
    use neort_gc_eqdsk_cut_jet, only: &
        EQDSK_CUT_JET_SUCCESS, eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
    use neort_gc_eqdsk_flux_profile_map, only: &
        EQDSK_FLUX_MAP_SUCCESS, eqdsk_flux_profile_map_t, &
        map_eqdsk_s_tor_to_psihat
    implicit none
    private

    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_SUCCESS = 0
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_AXIS_FAILURE = 2
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_ENDPOINT_FAILURE = 3
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_GRAPH_FAILURE = 4
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_NOT_MONOTONE = 5
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_NOT_CONNECTED = 6
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_INVALID_CERTIFICATE = 7
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_MAPPING_FAILURE = 8
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_OUT_OF_RANGE = 9
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_CERTIFICATE_ID = 130016
    integer, parameter, public :: EQDSK_COMPOSITE_CUT_INBOARD = -1
    integer, parameter, public :: EQDSK_COMPOSITE_CUT_OUTBOARD = 1

    type, public :: eqdsk_composite_cut_atlas_options_t
        type(eqdsk_cut_graph_atlas_options_t) :: graph
        type(eqdsk_cut_endpoint_options_t) :: endpoint
    end type eqdsk_composite_cut_atlas_options_t

    type, public :: eqdsk_composite_cut_atlas_t
        real(dp) :: target_psihat = 0.0_dp
        real(dp) :: field_scale = 0.0_dp
        real(dp) :: requested_z_lo = 0.0_dp
        real(dp) :: requested_z_hi = 0.0_dp
        type(eqdsk_cut_endpoint_options_t) :: endpoint_options
        type(eqdsk_axis_certificate_t) :: axis
        type(eqdsk_cut_endpoint_certificate_t) :: inboard_endpoint
        type(eqdsk_cut_endpoint_certificate_t) :: outboard_endpoint
        type(eqdsk_cut_graph_atlas_t) :: inboard_graph
        type(eqdsk_cut_graph_atlas_t) :: axis_graph
        type(eqdsk_cut_graph_atlas_t) :: outboard_graph
        integer :: certificate_id = 0
        logical :: geometric_completeness_certified = .false.
        logical :: branch_connectivity_certified = .false.
        logical :: surface_intersection_pair_certified = .false.
        logical :: orbit_crossing_multiplicity_certified = .false.
    end type eqdsk_composite_cut_atlas_t

    public :: build_eqdsk_composite_cut_atlas
    public :: map_eqdsk_composite_cut_atlas_rho
    public :: validate_eqdsk_composite_cut_atlas

contains

    subroutine build_eqdsk_composite_cut_atlas(axis_box, inboard_box, &
            outboard_box, inboard_seed, outboard_seed, target_psihat, &
            field_scale, z_lo, z_hi, options, atlas, status)
        real(dp), intent(in) :: axis_box(4), inboard_box(4), outboard_box(4)
        real(dp), intent(in) :: inboard_seed(2), outboard_seed(2)
        real(dp), intent(in) :: target_psihat, field_scale, z_lo, z_hi
        type(eqdsk_composite_cut_atlas_options_t), intent(in) :: options
        type(eqdsk_composite_cut_atlas_t), intent(out) :: atlas
        integer, intent(out) :: status

        integer :: local_status

        atlas = eqdsk_composite_cut_atlas_t()
        status = EQDSK_COMPOSITE_ATLAS_INVALID_INPUT
        if (.not. all(ieee_is_finite([axis_box, inboard_box, outboard_box, &
                inboard_seed, outboard_seed, target_psihat, field_scale, &
                z_lo, z_hi]))) return
        if (target_psihat <= 0.0_dp .or. target_psihat > 1.0_dp) return
        if (field_scale <= 0.0_dp .or. z_lo >= z_hi) return
        if (.not. valid_box(axis_box) .or. .not. valid_box(inboard_box) .or. &
                .not. valid_box(outboard_box)) return
        if (.not. point_in_box(inboard_seed, inboard_box) .or. &
                .not. point_in_box(outboard_seed, outboard_box)) return
        if (inboard_box(2) >= axis_box(1) .or. &
                outboard_box(1) <= axis_box(2)) return
        if (inboard_box(1) >= outboard_box(1)) return

        call build_eqdsk_axis_certificate(axis_box(1), axis_box(2), &
            axis_box(3), axis_box(4), atlas%axis, local_status)
        if (local_status /= EQDSK_AXIS_CERT_SUCCESS) then
            status = EQDSK_COMPOSITE_ATLAS_AXIS_FAILURE
            return
        end if
        call build_eqdsk_cut_endpoint_certificate(inboard_box(1), &
            inboard_box(2), inboard_box(3), inboard_box(4), target_psihat, &
            field_scale, inboard_seed(1), inboard_seed(2), options%endpoint, &
            atlas%inboard_endpoint, local_status)
        if (local_status /= EQDSK_ENDPOINT_CERT_SUCCESS) then
            status = EQDSK_COMPOSITE_ATLAS_ENDPOINT_FAILURE
            return
        end if
        call build_eqdsk_cut_endpoint_certificate(outboard_box(1), &
            outboard_box(2), outboard_box(3), outboard_box(4), target_psihat, &
            field_scale, outboard_seed(1), outboard_seed(2), options%endpoint, &
            atlas%outboard_endpoint, local_status)
        if (local_status /= EQDSK_ENDPOINT_CERT_SUCCESS) then
            status = EQDSK_COMPOSITE_ATLAS_ENDPOINT_FAILURE
            return
        end if

        call build_eqdsk_cut_graph_atlas(atlas%inboard_graph, &
            inboard_box(1), axis_box(1), z_lo, z_hi, 0.0_dp, 1.0_dp, &
            options%graph, local_status)
        if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
            status = EQDSK_COMPOSITE_ATLAS_GRAPH_FAILURE
            return
        end if
        call build_eqdsk_cut_graph_atlas(atlas%axis_graph, axis_box(1), &
            axis_box(2), z_lo, z_hi, 0.0_dp, 1.0_dp, options%graph, &
            local_status)
        if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
            status = EQDSK_COMPOSITE_ATLAS_GRAPH_FAILURE
            return
        end if
        call build_eqdsk_cut_graph_atlas(atlas%outboard_graph, axis_box(2), &
            outboard_box(2), z_lo, z_hi, 0.0_dp, 1.0_dp, options%graph, &
            local_status)
        if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
            status = EQDSK_COMPOSITE_ATLAS_GRAPH_FAILURE
            return
        end if
        if (.not. atlas%inboard_graph%flux_monotonicity_certified .or. &
                atlas%inboard_graph%flux_monotonicity_sign /= -1) then
            status = EQDSK_COMPOSITE_ATLAS_NOT_MONOTONE
            return
        end if
        if (.not. atlas%outboard_graph%flux_monotonicity_certified .or. &
                atlas%outboard_graph%flux_monotonicity_sign /= 1) then
            status = EQDSK_COMPOSITE_ATLAS_NOT_MONOTONE
            return
        end if

        atlas%target_psihat = target_psihat
        atlas%field_scale = field_scale
        atlas%requested_z_lo = z_lo
        atlas%requested_z_hi = z_hi
        atlas%endpoint_options = options%endpoint
        atlas%geometric_completeness_certified = .true.
        atlas%branch_connectivity_certified = .true.
        atlas%surface_intersection_pair_certified = .true.
        ! Orbit multiplicity additionally needs monotonic B on both arcs and
        ! the complete-return theorem.  Geometry alone must never set it.
        atlas%orbit_crossing_multiplicity_certified = .false.
        atlas%certificate_id = EQDSK_COMPOSITE_ATLAS_CERTIFICATE_ID
        call validate_eqdsk_composite_cut_atlas(atlas, status)
    end subroutine build_eqdsk_composite_cut_atlas

    subroutine validate_eqdsk_composite_cut_atlas(atlas, status)
        type(eqdsk_composite_cut_atlas_t), intent(inout) :: atlas
        integer, intent(out) :: status

        integer :: local_status

        status = EQDSK_COMPOSITE_ATLAS_INVALID_CERTIFICATE
        if (atlas%certificate_id /= EQDSK_COMPOSITE_ATLAS_CERTIFICATE_ID) return
        if (.not. atlas%geometric_completeness_certified) return
        if (.not. atlas%branch_connectivity_certified) return
        if (.not. atlas%surface_intersection_pair_certified) return
        if (atlas%orbit_crossing_multiplicity_certified) return
        if (.not. all(ieee_is_finite([atlas%target_psihat, &
                atlas%field_scale, atlas%requested_z_lo, &
                atlas%requested_z_hi]))) return
        if (atlas%target_psihat <= 0.0_dp .or. &
                atlas%target_psihat > 1.0_dp) return
        if (atlas%field_scale <= 0.0_dp .or. &
                atlas%requested_z_lo >= atlas%requested_z_hi) return
        call validate_eqdsk_axis_certificate(atlas%axis, local_status)
        if (local_status /= EQDSK_AXIS_CERT_SUCCESS) return
        call validate_eqdsk_cut_endpoint_certificate( &
            atlas%inboard_endpoint, local_status)
        if (local_status /= EQDSK_ENDPOINT_CERT_SUCCESS) return
        call validate_eqdsk_cut_endpoint_certificate( &
            atlas%outboard_endpoint, local_status)
        if (local_status /= EQDSK_ENDPOINT_CERT_SUCCESS) return
        call validate_eqdsk_cut_graph_atlas(atlas%inboard_graph, local_status)
        if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) return
        call validate_eqdsk_cut_graph_atlas(atlas%axis_graph, local_status)
        if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) return
        call validate_eqdsk_cut_graph_atlas(atlas%outboard_graph, local_status)
        if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) return
        if (.not. atlas%inboard_graph%flux_monotonicity_certified .or. &
                atlas%inboard_graph%flux_monotonicity_sign /= -1) return
        if (.not. atlas%outboard_graph%flux_monotonicity_certified .or. &
                atlas%outboard_graph%flux_monotonicity_sign /= 1) return
        if (atlas%inboard_endpoint%target_psihat /= atlas%target_psihat .or. &
                atlas%outboard_endpoint%target_psihat /= &
                atlas%target_psihat) return
        if (atlas%inboard_graph%requested_r_hi /= atlas%axis%r_lo) return
        if (atlas%axis_graph%requested_r_lo /= atlas%axis%r_lo .or. &
                atlas%axis_graph%requested_r_hi /= atlas%axis%r_hi) return
        if (atlas%outboard_graph%requested_r_lo /= atlas%axis%r_hi) return
        if (atlas%inboard_graph%requested_r_lo /= &
                atlas%inboard_endpoint%r_lo) return
        if (atlas%outboard_graph%requested_r_hi /= &
                atlas%outboard_endpoint%r_hi) return
        status = EQDSK_COMPOSITE_ATLAS_SUCCESS
    end subroutine validate_eqdsk_composite_cut_atlas

    subroutine map_eqdsk_composite_cut_atlas_rho(atlas, flux_map, rho_tor, &
            branch_sign, position, dposition_drho_tor, status)
        type(eqdsk_composite_cut_atlas_t), intent(inout) :: atlas
        type(eqdsk_flux_profile_map_t), intent(in) :: flux_map
        real(dp), intent(in) :: rho_tor
        integer, intent(in) :: branch_sign
        real(dp), intent(out) :: position(3), dposition_drho_tor(3)
        integer, intent(out) :: status

        type(eqdsk_cut_jet_t) :: axis_jet
        real(dp) :: s_tor, dstor_drho_tor, psihat, dpsihat_dstor
        real(dp) :: dposition_dR(3), dZ_dR, dpsihat_dR
        real(dp) :: axis_position(3), axis_slope
        real(dp) :: unused_psihat, axis_curvature
        integer :: local_status

        position = 0.0_dp
        dposition_drho_tor = 0.0_dp
        status = EQDSK_COMPOSITE_ATLAS_INVALID_INPUT
        if (.not. ieee_is_finite(rho_tor)) return
        if (rho_tor < 0.0_dp .or. rho_tor > 1.0_dp) then
            status = EQDSK_COMPOSITE_ATLAS_OUT_OF_RANGE
            return
        end if
        if (abs(branch_sign) /= 1) return
        call validate_eqdsk_composite_cut_atlas(atlas, local_status)
        if (local_status /= EQDSK_COMPOSITE_ATLAS_SUCCESS) then
            status = EQDSK_COMPOSITE_ATLAS_INVALID_CERTIFICATE
            return
        end if
        call evaluate_neort_eqdsk_rho_tor_map(rho_tor, s_tor, &
            dstor_drho_tor)
        call map_eqdsk_s_tor_to_psihat(flux_map, s_tor, psihat, &
            dpsihat_dstor, local_status)
        if (local_status /= EQDSK_FLUX_MAP_SUCCESS) then
            status = EQDSK_COMPOSITE_ATLAS_OUT_OF_RANGE
            return
        end if
        if (psihat > atlas%target_psihat) then
            status = EQDSK_COMPOSITE_ATLAS_OUT_OF_RANGE
            return
        end if

        call evaluate_axis_chart(atlas, axis_position, axis_slope, axis_jet, &
            local_status)
        if (local_status /= EQDSK_COMPOSITE_ATLAS_SUCCESS) then
            status = local_status
            return
        end if
        if (rho_tor == 0.0_dp) then
            call evaluate_neort_eqdsk_cut_axis_rho_limit(axis_slope, &
                axis_jet%psi_jet(4), axis_jet%psi_jet(5), &
                axis_jet%psi_jet(6), psi_sep, dpsihat_dstor, &
                real(branch_sign, dp), axis_curvature, &
                dposition_drho_tor(1), dposition_drho_tor(2))
            position = axis_position
            if (.not. all(ieee_is_finite([position, &
                    dposition_drho_tor, axis_curvature])) .or. &
                    axis_curvature <= 0.0_dp) then
                position = 0.0_dp
                dposition_drho_tor = 0.0_dp
                status = EQDSK_COMPOSITE_ATLAS_MAPPING_FAILURE
                return
            end if
            status = EQDSK_COMPOSITE_ATLAS_SUCCESS
            return
        end if

        if (branch_sign == EQDSK_COMPOSITE_CUT_INBOARD) then
            call map_eqdsk_cut_graph_atlas_flux(atlas%inboard_graph, &
                psihat, position, dposition_dR, dZ_dR, dpsihat_dR, &
                local_status)
        else
            call map_eqdsk_cut_graph_atlas_flux(atlas%outboard_graph, &
                psihat, position, dposition_dR, dZ_dR, dpsihat_dR, &
                local_status)
        end if
        if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
            call map_axis_neighborhood(atlas, axis_position, axis_slope, &
                axis_jet, psihat, branch_sign, position, dZ_dR, &
                dpsihat_dR, local_status)
            if (local_status /= EQDSK_COMPOSITE_ATLAS_SUCCESS) then
                status = local_status
                return
            end if
        end if
        call evaluate_neort_eqdsk_cut_flux_coordinate(dstor_drho_tor, &
            dpsihat_dstor, dpsihat_dR, dZ_dR, unused_psihat, &
            dposition_drho_tor(1), dposition_drho_tor(2))
        if (.not. all(ieee_is_finite([position, dposition_drho_tor, &
                unused_psihat])) .or. dpsihat_dR*real(branch_sign, dp) <= &
                0.0_dp) then
            position = 0.0_dp
            dposition_drho_tor = 0.0_dp
            status = EQDSK_COMPOSITE_ATLAS_MAPPING_FAILURE
            return
        end if
        status = EQDSK_COMPOSITE_ATLAS_SUCCESS
    end subroutine map_eqdsk_composite_cut_atlas_rho

    subroutine evaluate_axis_chart(atlas, position, slope, jet, status)
        type(eqdsk_composite_cut_atlas_t), intent(in) :: atlas
        real(dp), intent(out) :: position(3), slope
        type(eqdsk_cut_jet_t), intent(out) :: jet
        integer, intent(out) :: status

        real(dp) :: unused_dpsihat_dR
        integer :: local_status

        position = [atlas%axis%axis_newton_R, atlas%axis%axis_newton_Z, &
            0.0_dp]
        slope = 0.0_dp
        jet = eqdsk_cut_jet_t()
        call evaluate_eqdsk_cut_jet(position, atlas%field_scale, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], jet, local_status)
        if (local_status /= EQDSK_CUT_JET_SUCCESS) then
            status = EQDSK_COMPOSITE_ATLAS_MAPPING_FAILURE
            return
        end if
        call evaluate_neort_eqdsk_cut_r_flux_chart( &
            jet%d_cut_numerator_d_R, jet%d_cut_numerator_d_Z, &
            jet%psi_jet(2), jet%psi_jet(3), psi_sep, slope, &
            unused_dpsihat_dR)
        if (.not. all(ieee_is_finite([position, slope, &
                unused_dpsihat_dR]))) then
            status = EQDSK_COMPOSITE_ATLAS_MAPPING_FAILURE
            return
        end if
        status = EQDSK_COMPOSITE_ATLAS_SUCCESS
    end subroutine evaluate_axis_chart

    subroutine map_axis_neighborhood(atlas, axis_position, axis_slope, &
            axis_jet, psihat, branch_sign, position, dZ_dR, dpsihat_dR, &
            status)
        type(eqdsk_composite_cut_atlas_t), intent(in) :: atlas
        real(dp), intent(in) :: axis_position(3), axis_slope, psihat
        type(eqdsk_cut_jet_t), intent(in) :: axis_jet
        integer, intent(in) :: branch_sign
        real(dp), intent(out) :: position(3), dZ_dR, dpsihat_dR
        integer, intent(out) :: status

        type(eqdsk_cut_endpoint_certificate_t) :: endpoint
        type(eqdsk_cut_jet_t) :: jet
        real(dp) :: axis_curvature, absolute_dpsihat_dR, delta_R
        real(dp) :: seed(2), half_R, half_Z, box(4)
        integer :: local_status

        position = 0.0_dp
        dZ_dR = 0.0_dp
        dpsihat_dR = 0.0_dp
        status = EQDSK_COMPOSITE_ATLAS_MAPPING_FAILURE
        call evaluate_neort_eqdsk_cut_axis_limit(axis_slope, &
            axis_jet%psi_jet(4), axis_jet%psi_jet(5), &
            axis_jet%psi_jet(6), psi_sep, psihat, real(branch_sign, dp), &
            axis_curvature, absolute_dpsihat_dR, delta_R)
        if (.not. all(ieee_is_finite([axis_curvature, &
                absolute_dpsihat_dR, delta_R])) .or. &
                axis_curvature <= 0.0_dp .or. abs(delta_R) <= &
                tiny(delta_R)) return
        seed = axis_position(1:2)+[delta_R, axis_slope*delta_R]
        half_R = 0.45_dp*abs(delta_R)
        half_Z = half_R*max(1.0_dp, abs(axis_slope))
        box = [seed(1)-half_R, seed(1)+half_R, seed(2)-half_Z, &
            seed(2)+half_Z]
        if (.not. valid_box(box)) return
        if (box(1) < atlas%axis%r_lo .or. box(2) > atlas%axis%r_hi .or. &
                box(3) < atlas%axis%z_lo .or. &
                box(4) > atlas%axis%z_hi) return
        call build_eqdsk_cut_endpoint_certificate(box(1), box(2), box(3), &
            box(4), psihat, atlas%field_scale, seed(1), seed(2), &
            atlas%endpoint_options, endpoint, local_status)
        if (local_status /= EQDSK_ENDPOINT_CERT_SUCCESS) return
        position = [endpoint%newton_point_R, endpoint%newton_point_Z, 0.0_dp]
        call evaluate_eqdsk_cut_jet(position, atlas%field_scale, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], jet, local_status)
        if (local_status /= EQDSK_CUT_JET_SUCCESS) return
        call evaluate_neort_eqdsk_cut_r_flux_chart( &
            jet%d_cut_numerator_d_R, jet%d_cut_numerator_d_Z, &
            jet%psi_jet(2), jet%psi_jet(3), psi_sep, dZ_dR, dpsihat_dR)
        if (.not. all(ieee_is_finite([position, dZ_dR, &
                dpsihat_dR]))) return
        status = EQDSK_COMPOSITE_ATLAS_SUCCESS
    end subroutine map_axis_neighborhood

    pure logical function valid_box(box)
        real(dp), intent(in) :: box(4)

        valid_box = box(1) > 0.0_dp .and. box(1) < box(2) .and. &
            box(3) < box(4)
    end function valid_box

    pure logical function point_in_box(point, box)
        real(dp), intent(in) :: point(2), box(4)

        point_in_box = point(1) > box(1) .and. point(1) < box(2) .and. &
            point(2) > box(3) .and. point(2) < box(4)
    end function point_in_box

end module neort_gc_eqdsk_composite_cut_atlas
