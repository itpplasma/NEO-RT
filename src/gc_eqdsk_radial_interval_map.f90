module neort_gc_eqdsk_radial_interval_map
    !! Certified regular-outboard map from physical R boxes to rho_tor boxes.
    !!
    !! This module owns only certificate validation, strip partitioning,
    !! interval intersection, and hulling.  Fortsym-generated kernels own
    !! every physical and coordinate formula in the chain
    !!
    !!   raw EQDSK psi -> physical psi -> psi_hat -> s_tor -> rho_tor.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: psi_sep
    use neort_eqdsk_physical_flux_map_interval_symbolic, only: &
        evaluate_neort_eqdsk_physical_flux_map_interval
    use neort_eqdsk_physical_flux_normalization_interval_symbolic, only: &
        evaluate_neort_eqdsk_physical_flux_normalization_interval
    use neort_eqdsk_s_tor_to_rho_interval_symbolic, only: &
        evaluate_neort_eqdsk_s_tor_to_rho_interval
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_ATLAS_SUCCESS, EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
        enclose_eqdsk_cut_graph_strip, eqdsk_cut_graph_atlas_t, &
        validate_eqdsk_cut_graph_atlas
    use neort_gc_eqdsk_cut_interval, only: eqdsk_cut_interval_result_t
    use neort_gc_eqdsk_flux_profile_map, only: &
        EQDSK_FLUX_MAP_CERTIFICATE_ID, EQDSK_FLUX_MAP_SUCCESS, &
        eqdsk_flux_profile_map_t, map_eqdsk_psihat_interval_to_s_tor, &
        validate_eqdsk_flux_profile_map
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_is_valid, gc_outward_interval_t
    implicit none
    private

    integer, parameter, public :: EQDSK_RADIAL_INTERVAL_SUCCESS = 0
    integer, parameter, public :: EQDSK_RADIAL_INTERVAL_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_RADIAL_INTERVAL_INVALID_ATLAS = 2
    integer, parameter, public :: EQDSK_RADIAL_INTERVAL_INVALID_PROFILE = 3
    integer, parameter, public :: EQDSK_RADIAL_INTERVAL_OUT_OF_RANGE = 4
    integer, parameter, public :: EQDSK_RADIAL_INTERVAL_GRAPH_FAILURE = 5
    integer, parameter, public :: EQDSK_RADIAL_INTERVAL_NONFINITE = 6
    integer, parameter, public :: EQDSK_RADIAL_INTERVAL_NOT_MONOTONE = 7
    integer, parameter, public :: EQDSK_RADIAL_INTERVAL_NORMALIZATION = 8

    type, public :: eqdsk_rho_interval_provenance_t
        type(gc_outward_interval_t) :: radius
        type(gc_outward_interval_t) :: psi_physical
        type(gc_outward_interval_t) :: psihat
        type(gc_outward_interval_t) :: s_tor
        type(gc_outward_interval_t) :: rho_tor
        integer :: graph_certificate_id = 0
        integer :: profile_certificate_id = 0
        integer :: nstrips = 0
        integer :: n_graph_enclosures = 0
        integer :: profile_segments_covered = 0
        integer, allocatable :: strip_index(:)
        real(dp), allocatable :: strip_r_lo(:)
        real(dp), allocatable :: strip_r_hi(:)
        logical :: axis_limit_used = .false.
        logical :: mapping_certified = .false.
    end type eqdsk_rho_interval_provenance_t

    public :: map_eqdsk_outboard_r_interval_to_rho_tor

contains

    subroutine map_eqdsk_outboard_r_interval_to_rho_tor(atlas, flux_map, &
            r_lo, r_hi, rho_tor, provenance, status)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        type(eqdsk_flux_profile_map_t), intent(in) :: flux_map
        real(dp), intent(in) :: r_lo, r_hi
        type(gc_outward_interval_t), intent(out) :: rho_tor
        type(eqdsk_rho_interval_provenance_t), intent(out) :: provenance
        integer, intent(out) :: status

        type(eqdsk_cut_interval_result_t), allocatable :: enclosures(:)
        type(gc_outward_interval_t) :: physical_psi, physical_derivative
        type(gc_outward_interval_t) :: normalized_psi, s_tor
        type(gc_outward_interval_t) :: point_scale, point_psi_sep
        real(dp) :: sub_lo, sub_hi
        integer :: i, j, local_status, strip_count, profile_segments
        logical :: first_enclosure

        rho_tor = gc_outward_interval(0.0_dp, 0.0_dp)
        provenance = eqdsk_rho_interval_provenance_t()
        status = EQDSK_RADIAL_INTERVAL_INVALID_INPUT
        if (.not. all(ieee_is_finite([r_lo, r_hi]))) return
        if (r_lo <= 0.0_dp .or. r_hi < r_lo) return
        call validate_eqdsk_cut_graph_atlas(atlas, local_status)
        if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
            status = EQDSK_RADIAL_INTERVAL_INVALID_ATLAS
            return
        end if
        if (.not. atlas%flux_monotonicity_certified) then
            status = EQDSK_RADIAL_INTERVAL_NOT_MONOTONE
            return
        end if
        if (atlas%flux_monotonicity_sign /= 1) then
            status = EQDSK_RADIAL_INTERVAL_NOT_MONOTONE
            return
        end if
        call validate_eqdsk_flux_profile_map(flux_map, local_status)
        if (local_status /= EQDSK_FLUX_MAP_SUCCESS) then
            status = EQDSK_RADIAL_INTERVAL_INVALID_PROFILE
            return
        end if
        if (flux_map%psi_sep /= psi_sep) then
            status = EQDSK_RADIAL_INTERVAL_NORMALIZATION
            return
        end if
        if (r_lo < atlas%requested_r_lo .or. &
                r_hi > atlas%requested_r_hi) then
            status = EQDSK_RADIAL_INTERVAL_OUT_OF_RANGE
            return
        end if

        strip_count = 0
        do i = 1, size(atlas%strips)
            sub_lo = max(r_lo, atlas%strips(i)%r_lo)
            sub_hi = min(r_hi, atlas%strips(i)%r_hi)
            if (sub_hi < sub_lo) cycle
            strip_count = strip_count+1
        end do
        if (strip_count < 1) then
            status = EQDSK_RADIAL_INTERVAL_OUT_OF_RANGE
            return
        end if
        allocate(provenance%strip_index(strip_count), &
            provenance%strip_r_lo(strip_count), &
            provenance%strip_r_hi(strip_count))
        provenance%radius = gc_outward_interval(r_lo, r_hi)
        provenance%nstrips = strip_count
        provenance%graph_certificate_id = EQDSK_CUT_GRAPH_CERTIFICATE_ID
        provenance%profile_certificate_id = EQDSK_FLUX_MAP_CERTIFICATE_ID
        point_scale = gc_outward_interval(flux_map%field_scale, &
            flux_map%field_scale)
        point_psi_sep = gc_outward_interval(flux_map%psi_sep, &
            flux_map%psi_sep)

        strip_count = 0
        first_enclosure = .true.
        do i = 1, size(atlas%strips)
            sub_lo = max(r_lo, atlas%strips(i)%r_lo)
            sub_hi = min(r_hi, atlas%strips(i)%r_hi)
            if (sub_hi < sub_lo) cycle
            strip_count = strip_count+1
            provenance%strip_index(strip_count) = i
            provenance%strip_r_lo(strip_count) = sub_lo
            provenance%strip_r_hi(strip_count) = sub_hi
            call enclose_eqdsk_cut_graph_strip(atlas, i, sub_lo, sub_hi, &
                enclosures, local_status)
            if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
                status = EQDSK_RADIAL_INTERVAL_GRAPH_FAILURE
                return
            end if
            if (.not. allocated(enclosures)) then
                status = EQDSK_RADIAL_INTERVAL_GRAPH_FAILURE
                return
            end if
            if (size(enclosures) < 1) then
                status = EQDSK_RADIAL_INTERVAL_GRAPH_FAILURE
                return
            end if
            do j = 1, size(enclosures)
                if (.not. enclosures(j)%r_chart_certified) then
                    status = EQDSK_RADIAL_INTERVAL_INVALID_ATLAS
                    return
                end if
                if (enclosures(j)%dpsihat_dR%lo <= 0.0_dp) then
                    status = EQDSK_RADIAL_INTERVAL_NOT_MONOTONE
                    return
                end if
                call evaluate_neort_eqdsk_physical_flux_map_interval( &
                    enclosures(j)%psi, enclosures(j)%psi_R, &
                    enclosures(j)%psi_Z, enclosures(j)%dZ_dR, point_scale, &
                    physical_psi, physical_derivative)
                if (.not. gc_outward_interval_is_valid(physical_psi)) then
                    status = EQDSK_RADIAL_INTERVAL_NONFINITE
                    return
                end if
                if (.not. gc_outward_interval_is_valid( &
                        physical_derivative)) then
                    status = EQDSK_RADIAL_INTERVAL_NONFINITE
                    return
                end if
                if (physical_derivative%lo <= 0.0_dp) then
                    status = EQDSK_RADIAL_INTERVAL_NOT_MONOTONE
                    return
                end if
                call evaluate_neort_eqdsk_physical_flux_normalization_interval( &
                    physical_psi, point_scale, point_psi_sep, normalized_psi)
                if (.not. gc_outward_interval_is_valid(normalized_psi)) then
                    status = EQDSK_RADIAL_INTERVAL_NONFINITE
                    return
                end if
                normalized_psi%lo = max(normalized_psi%lo, &
                    enclosures(j)%psi_hat%lo)
                normalized_psi%hi = min(normalized_psi%hi, &
                    enclosures(j)%psi_hat%hi)
                if (normalized_psi%hi < normalized_psi%lo) then
                    status = EQDSK_RADIAL_INTERVAL_NORMALIZATION
                    return
                end if
                if (first_enclosure) then
                    provenance%psi_physical = physical_psi
                    provenance%psihat = normalized_psi
                    first_enclosure = .false.
                else
                    call include_interval(provenance%psi_physical, &
                        physical_psi)
                    call include_interval(provenance%psihat, normalized_psi)
                end if
                provenance%n_graph_enclosures = &
                    provenance%n_graph_enclosures+1
            end do
        end do
        if (first_enclosure) then
            status = EQDSK_RADIAL_INTERVAL_GRAPH_FAILURE
            return
        end if
        if (provenance%psihat%lo < 0.0_dp .or. &
                provenance%psihat%hi > 1.0_dp) then
            status = EQDSK_RADIAL_INTERVAL_OUT_OF_RANGE
            return
        end if
        call map_eqdsk_psihat_interval_to_s_tor(flux_map, &
            provenance%psihat, s_tor, profile_segments, local_status)
        if (local_status /= EQDSK_FLUX_MAP_SUCCESS) then
            status = EQDSK_RADIAL_INTERVAL_INVALID_PROFILE
            return
        end if
        call evaluate_neort_eqdsk_s_tor_to_rho_interval(s_tor, rho_tor)
        if (.not. gc_outward_interval_is_valid(rho_tor)) then
            status = EQDSK_RADIAL_INTERVAL_NONFINITE
            return
        end if
        if (rho_tor%lo < 0.0_dp .or. rho_tor%hi > 1.0_dp) then
            status = EQDSK_RADIAL_INTERVAL_OUT_OF_RANGE
            return
        end if
        provenance%s_tor = s_tor
        provenance%rho_tor = rho_tor
        provenance%profile_segments_covered = profile_segments
        provenance%mapping_certified = .true.
        status = EQDSK_RADIAL_INTERVAL_SUCCESS
    end subroutine map_eqdsk_outboard_r_interval_to_rho_tor

    pure subroutine include_interval(hull, value)
        type(gc_outward_interval_t), intent(inout) :: hull
        type(gc_outward_interval_t), intent(in) :: value

        hull%lo = min(hull%lo, value%lo)
        hull%hi = max(hull%hi, value%hi)
    end subroutine include_interval

end module neort_gc_eqdsk_radial_interval_map
