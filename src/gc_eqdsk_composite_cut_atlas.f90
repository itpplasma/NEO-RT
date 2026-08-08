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
        validate_eqdsk_cut_graph_atlas
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
    integer, parameter, public :: EQDSK_COMPOSITE_ATLAS_CERTIFICATE_ID = 130016

    type, public :: eqdsk_composite_cut_atlas_options_t
        type(eqdsk_cut_graph_atlas_options_t) :: graph
        type(eqdsk_cut_endpoint_options_t) :: endpoint
    end type eqdsk_composite_cut_atlas_options_t

    type, public :: eqdsk_composite_cut_atlas_t
        real(dp) :: target_psihat = 0.0_dp
        real(dp) :: field_scale = 0.0_dp
        real(dp) :: requested_z_lo = 0.0_dp
        real(dp) :: requested_z_hi = 0.0_dp
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
