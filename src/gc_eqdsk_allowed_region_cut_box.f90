module neort_gc_eqdsk_allowed_region_cut_box
    !! Certified runtime orchestration for one physical-R graph box.
    !!
    !! This module owns no magnetic, electric, energy, or canonical formula.
    !! It validates source certificates and explicit profile-node structure,
    !! partitions the request at certified graph strips, preserves every
    !! returned leaf, and forwards those leaves to the Fortsym-backed
    !! allowed-region evaluator.  Profile data validation is not a data
    !! certificate; the generated interpolation-kernel identity is recorded
    !! separately in the returned provenance.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_generated_certificate_registry, only: certificate_fingerprint, &
        certificate_index
    use neort_gc_eqdsk_allowed_region_interval, only: &
        EQDSK_ALLOWED_INTERVAL_SUCCESS, &
        eqdsk_allowed_interval_result_t, &
        evaluate_eqdsk_allowed_region_interval
    use neort_gc_eqdsk_cut_graph_atlas, only: &
        EQDSK_CUT_ATLAS_SUCCESS, EQDSK_CUT_GRAPH_CERTIFICATE_ID, &
        enclose_eqdsk_cut_graph_strip, enclose_eqdsk_cut_graph_strip_tight, &
        eqdsk_cut_graph_atlas_t, validate_eqdsk_cut_graph_atlas
    use neort_gc_eqdsk_cut_interval, only: eqdsk_cut_interval_result_t
    use neort_gc_outward_interval, only: gc_outward_interval_is_valid, &
        gc_outward_interval_t
    implicit none
    private

    integer, parameter, public :: EQDSK_CUT_BOX_SUCCESS = 0
    integer, parameter, public :: EQDSK_CUT_BOX_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_CUT_BOX_INVALID_ATLAS = 2
    integer, parameter, public :: EQDSK_CUT_BOX_INVALID_PROFILE = 3
    integer, parameter, public :: EQDSK_CUT_BOX_OUT_OF_RANGE = 4
    integer, parameter, public :: EQDSK_CUT_BOX_NOT_MONOTONE = 5
    ! Compatibility spelling retained for callers of the former outboard-only
    ! API.  The accepted atlas orientation is now explicitly recorded below.
    integer, parameter, public :: EQDSK_CUT_BOX_NOT_OUTBOARD = &
        EQDSK_CUT_BOX_NOT_MONOTONE
    integer, parameter, public :: EQDSK_CUT_BOX_GRAPH_FAILURE = 6
    integer, parameter, public :: EQDSK_CUT_BOX_ALLOWED_FAILURE = 7
    integer, parameter, public :: EQDSK_CUT_BOX_NONFINITE = 8

    type, public :: eqdsk_potential_profile_nodes_t
        real(dp), allocatable :: psi(:)
        real(dp), allocatable :: phi(:)
        real(dp), allocatable :: omega(:)
        logical :: structurally_validated = .false.
    end type eqdsk_potential_profile_nodes_t

    type, public :: eqdsk_allowed_region_cut_provenance_t
        type(gc_outward_interval_t) :: query_radius
        real(dp) :: field_scale = 0.0_dp
        real(dp) :: raw_psi_sep = 0.0_dp
        integer :: graph_certificate_id = 0
        integer :: graph_flux_orientation_sign = 0
        logical :: profile_inputs_structurally_validated = .false.
        integer :: generated_profile_interpolation_certificate_id = 0
        integer :: generated_energy_certificate_id = 0
        integer :: generated_canonical_certificate_id = 0
        character(len=64) :: generated_profile_interpolation_fingerprint = ''
        character(len=64) :: generated_energy_fingerprint = ''
        character(len=64) :: generated_canonical_fingerprint = ''
        integer :: nstrips = 0
        integer :: n_graph_enclosures = 0
        integer, allocatable :: strip_index(:)
        real(dp), allocatable :: strip_r_lo(:)
        real(dp), allocatable :: strip_r_hi(:)
        integer, allocatable :: leaf_strip_index(:)
        ! These are the containing strip query bounds for each returned
        ! leaf, not exact internal leaf-R or leaf-Z boxes.
        real(dp), allocatable :: leaf_strip_r_lo(:)
        real(dp), allocatable :: leaf_strip_r_hi(:)
        real(dp), allocatable :: leaf_strip_z_lo(:)
        real(dp), allocatable :: leaf_strip_z_hi(:)
        type(eqdsk_cut_interval_result_t), allocatable :: leaf_enclosures(:)
        logical :: certified = .false.
    end type eqdsk_allowed_region_cut_provenance_t

    public :: validate_eqdsk_potential_profile_nodes
    public :: clear_eqdsk_potential_profile_nodes
    public :: evaluate_eqdsk_allowed_region_cut_box

contains

    subroutine validate_eqdsk_potential_profile_nodes(psi, phi, omega, &
            nodes, status)
        real(dp), intent(in) :: psi(:), phi(:), omega(:)
        type(eqdsk_potential_profile_nodes_t), intent(out) :: nodes
        integer, intent(out) :: status

        integer :: i, count

        call clear_eqdsk_potential_profile_nodes(nodes)
        status = EQDSK_CUT_BOX_INVALID_PROFILE
        count = size(psi)
        if (count < 2) return
        if (size(phi) /= count) return
        if (size(omega) /= count) return
        if (.not. all(ieee_is_finite(psi))) then
            status = EQDSK_CUT_BOX_NONFINITE
            return
        end if
        if (.not. all(ieee_is_finite(phi))) then
            status = EQDSK_CUT_BOX_NONFINITE
            return
        end if
        if (.not. all(ieee_is_finite(omega))) then
            status = EQDSK_CUT_BOX_NONFINITE
            return
        end if
        do i = 1, count-1
            if (psi(i+1) <= psi(i)) return
        end do

        allocate(nodes%psi(count), nodes%phi(count), nodes%omega(count))
        nodes%psi = psi
        nodes%phi = phi
        nodes%omega = omega
        nodes%structurally_validated = .true.
        status = EQDSK_CUT_BOX_SUCCESS
    end subroutine validate_eqdsk_potential_profile_nodes

    subroutine clear_eqdsk_potential_profile_nodes(nodes)
        type(eqdsk_potential_profile_nodes_t), intent(inout) :: nodes

        if (allocated(nodes%psi)) deallocate(nodes%psi)
        if (allocated(nodes%phi)) deallocate(nodes%phi)
        if (allocated(nodes%omega)) deallocate(nodes%omega)
        nodes%structurally_validated = .false.
    end subroutine clear_eqdsk_potential_profile_nodes

    subroutine evaluate_eqdsk_allowed_region_cut_box(atlas, radius, &
            field_scale, raw_psi_sep, profile, h0, j_k, mass, charge, &
            c_light, sigma, result, provenance, status, tighten_z_depth)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        type(gc_outward_interval_t), intent(in) :: radius
        real(dp), intent(in) :: field_scale, raw_psi_sep
        type(eqdsk_potential_profile_nodes_t), intent(in) :: profile
        real(dp), intent(in) :: h0, j_k, mass, charge, c_light
        integer, intent(in) :: sigma
        type(eqdsk_allowed_interval_result_t), intent(out) :: result
        type(eqdsk_allowed_region_cut_provenance_t), intent(out) :: &
            provenance
        integer, intent(out) :: status
        integer, intent(in), optional :: tighten_z_depth

        type(eqdsk_cut_interval_result_t), allocatable :: enclosures(:)
        type(eqdsk_cut_interval_result_t), allocatable :: collected(:)
        integer, allocatable :: leaf_strip(:)
        real(dp), allocatable :: leaf_strip_r_lo(:), leaf_strip_r_hi(:)
        real(dp), allocatable :: leaf_strip_z_lo(:), leaf_strip_z_hi(:)
        real(dp) :: sub_lo, sub_hi
        integer :: i, j, strip_count, strip_position
        integer :: local_status, atlas_status
        integer :: z_refinement_depth
        logical :: valid

        result = eqdsk_allowed_interval_result_t()
        provenance = eqdsk_allowed_region_cut_provenance_t()
        status = EQDSK_CUT_BOX_INVALID_INPUT
        z_refinement_depth = 0
        if (present(tighten_z_depth)) z_refinement_depth = tighten_z_depth
        if (z_refinement_depth < 0) return

        if (.not. gc_outward_interval_is_valid(radius)) return
        if (radius%lo <= 0.0_dp .or. radius%hi < radius%lo) return
        if (.not. all(ieee_is_finite([field_scale, raw_psi_sep, h0, &
                j_k, mass, charge, c_light]))) then
            status = EQDSK_CUT_BOX_NONFINITE
            return
        end if
        if (field_scale <= 0.0_dp .or. raw_psi_sep <= 0.0_dp) return
        if (mass <= 0.0_dp .or. c_light <= 0.0_dp) return
        if (j_k < 0.0_dp .or. abs(charge) <= tiny(charge)) return
        if (abs(sigma) /= 1) return

        call validate_eqdsk_cut_graph_atlas(atlas, atlas_status)
        if (atlas_status /= EQDSK_CUT_ATLAS_SUCCESS) then
            status = EQDSK_CUT_BOX_INVALID_ATLAS
            return
        end if
        if (atlas%certificate_id /= EQDSK_CUT_GRAPH_CERTIFICATE_ID) then
            status = EQDSK_CUT_BOX_INVALID_ATLAS
            return
        end if
        if (raw_psi_sep /= atlas%raw_psi_sep) then
            status = EQDSK_CUT_BOX_INVALID_ATLAS
            return
        end if
        ! Allowed energy and canonical momentum are evaluated in physical R.
        ! They do not invert the flux coordinate, so the complete axis graph
        ! (where d psi/dR changes sign) is a valid source box as well.  Retain
        ! the graph orientation as provenance when one exists, but do not make
        ! flux monotonicity an unstated physics precondition.
        if (radius%lo < atlas%requested_r_lo .or. &
                radius%hi > atlas%requested_r_hi) then
            status = EQDSK_CUT_BOX_OUT_OF_RANGE
            return
        end if

        ! This validates only array shape, finiteness, and strict node order;
        ! it does not certify the supplied physical profile data.
        valid = valid_profile_nodes(profile)
        if (.not. valid) then
            status = EQDSK_CUT_BOX_INVALID_PROFILE
            return
        end if
        if (.not. generated_identities_available(provenance)) then
            status = EQDSK_CUT_BOX_INVALID_INPUT
            return
        end if

        strip_count = 0
        do i = 1, size(atlas%strips)
            sub_lo = max(radius%lo, atlas%strips(i)%r_lo)
            sub_hi = min(radius%hi, atlas%strips(i)%r_hi)
            if (sub_hi < sub_lo) cycle
            strip_count = strip_count+1
        end do
        if (strip_count < 1) then
            status = EQDSK_CUT_BOX_OUT_OF_RANGE
            return
        end if

        allocate(provenance%strip_index(strip_count), &
            provenance%strip_r_lo(strip_count), &
            provenance%strip_r_hi(strip_count))
        provenance%query_radius = radius
        provenance%field_scale = field_scale
        provenance%raw_psi_sep = raw_psi_sep
        provenance%graph_certificate_id = atlas%certificate_id
        provenance%graph_flux_orientation_sign = &
            atlas%flux_monotonicity_sign
        provenance%profile_inputs_structurally_validated = .true.
        provenance%nstrips = strip_count

        strip_position = 0
        do i = 1, size(atlas%strips)
            sub_lo = max(radius%lo, atlas%strips(i)%r_lo)
            sub_hi = min(radius%hi, atlas%strips(i)%r_hi)
            if (sub_hi < sub_lo) cycle
            strip_position = strip_position+1
            provenance%strip_index(strip_position) = i
            provenance%strip_r_lo(strip_position) = sub_lo
            provenance%strip_r_hi(strip_position) = sub_hi
            if (z_refinement_depth > 0) then
                call enclose_eqdsk_cut_graph_strip_tight(atlas, i, sub_lo, &
                    sub_hi, z_refinement_depth, enclosures, local_status)
            else
                call enclose_eqdsk_cut_graph_strip(atlas, i, sub_lo, sub_hi, &
                    enclosures, local_status)
            end if
            if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
                status = EQDSK_CUT_BOX_GRAPH_FAILURE
                return
            end if
            if (.not. allocated(enclosures)) then
                status = EQDSK_CUT_BOX_GRAPH_FAILURE
                return
            end if
            if (size(enclosures) < 1) then
                status = EQDSK_CUT_BOX_GRAPH_FAILURE
                return
            end if
            do j = 1, size(enclosures)
                call append_leaf(enclosures(j), i, sub_lo, sub_hi, &
                    atlas%strips(i)%z_lo, atlas%strips(i)%z_hi, collected, &
                    leaf_strip, leaf_strip_r_lo, leaf_strip_r_hi, &
                    leaf_strip_z_lo, leaf_strip_z_hi)
            end do
        end do

        if (strip_position /= provenance%nstrips) then
            status = EQDSK_CUT_BOX_GRAPH_FAILURE
            return
        end if
        if (.not. allocated(collected)) then
            status = EQDSK_CUT_BOX_GRAPH_FAILURE
            return
        end if
        if (.not. allocated(leaf_strip)) then
            status = EQDSK_CUT_BOX_GRAPH_FAILURE
            return
        end if
        if (size(leaf_strip) < 1) then
            status = EQDSK_CUT_BOX_GRAPH_FAILURE
            return
        end if

        call evaluate_eqdsk_allowed_region_interval(radius, field_scale, &
            raw_psi_sep, collected, profile%psi, profile%phi, &
            profile%omega, h0, j_k, mass, charge, c_light, sigma, &
            result, local_status)
        if (local_status /= EQDSK_ALLOWED_INTERVAL_SUCCESS) then
            status = EQDSK_CUT_BOX_ALLOWED_FAILURE
            return
        end if

        call move_alloc(leaf_strip, provenance%leaf_strip_index)
        call move_alloc(leaf_strip_r_lo, provenance%leaf_strip_r_lo)
        call move_alloc(leaf_strip_r_hi, provenance%leaf_strip_r_hi)
        call move_alloc(leaf_strip_z_lo, provenance%leaf_strip_z_lo)
        call move_alloc(leaf_strip_z_hi, provenance%leaf_strip_z_hi)
        call move_alloc(collected, provenance%leaf_enclosures)
        provenance%n_graph_enclosures = size(provenance%leaf_enclosures)
        ! The output certificate rests on the validated graph leaves, the
        ! generated interval kernels called by the delegated evaluator, and
        ! the explicit inputs checked above; no profile-data certificate is
        ! inferred from the generated interpolation-kernel identity.
        provenance%certified = .true.
        status = EQDSK_CUT_BOX_SUCCESS
    end subroutine evaluate_eqdsk_allowed_region_cut_box

    logical function valid_profile_nodes(profile)
        type(eqdsk_potential_profile_nodes_t), intent(in) :: profile

        integer :: i, count

        valid_profile_nodes = .false.
        if (.not. profile%structurally_validated) return
        if (.not. allocated(profile%psi)) return
        if (.not. allocated(profile%phi)) return
        if (.not. allocated(profile%omega)) return
        count = size(profile%psi)
        if (count < 2) return
        if (size(profile%phi) /= count) return
        if (size(profile%omega) /= count) return
        if (.not. all(ieee_is_finite(profile%psi))) return
        if (.not. all(ieee_is_finite(profile%phi))) return
        if (.not. all(ieee_is_finite(profile%omega))) return
        do i = 1, count-1
            if (profile%psi(i+1) <= profile%psi(i)) return
        end do
        valid_profile_nodes = .true.
    end function valid_profile_nodes

    logical function generated_identities_available(provenance)
        ! These IDs identify generated kernels used by the evaluator.  They
        ! are not certificates of the caller-supplied profile arrays.
        type(eqdsk_allowed_region_cut_provenance_t), intent(inout) :: &
            provenance
        integer :: profile_id, energy_id, canonical_id

        generated_identities_available = .false.
        profile_id = certificate_index('profile_potential_map')
        energy_id = certificate_index('eqdsk_allowed_energy')
        canonical_id = certificate_index('eqdsk_canonical_cut')
        if (profile_id <= 0) return
        if (energy_id <= 0) return
        if (canonical_id <= 0) return
        provenance%generated_profile_interpolation_certificate_id = &
            profile_id
        provenance%generated_energy_certificate_id = energy_id
        provenance%generated_canonical_certificate_id = canonical_id
        provenance%generated_profile_interpolation_fingerprint = &
            certificate_fingerprint(profile_id)
        provenance%generated_energy_fingerprint = &
            certificate_fingerprint(energy_id)
        provenance%generated_canonical_fingerprint = &
            certificate_fingerprint(canonical_id)
        generated_identities_available = .true.
    end function generated_identities_available

    subroutine append_leaf(leaf, strip_index, r_lo, r_hi, z_lo, z_hi, &
            collected, strip_ids, strip_r_los, strip_r_his, strip_z_los, &
            strip_z_his)
        type(eqdsk_cut_interval_result_t), intent(in) :: leaf
        integer, intent(in) :: strip_index
        real(dp), intent(in) :: r_lo, r_hi, z_lo, z_hi
        type(eqdsk_cut_interval_result_t), allocatable, intent(inout) :: &
            collected(:)
        integer, allocatable, intent(inout) :: strip_ids(:)
        real(dp), allocatable, intent(inout) :: strip_r_los(:), strip_r_his(:)
        real(dp), allocatable, intent(inout) :: strip_z_los(:), strip_z_his(:)

        type(eqdsk_cut_interval_result_t), allocatable :: new_collected(:)
        integer, allocatable :: new_strip_ids(:)
        real(dp), allocatable :: new_strip_r_los(:), new_strip_r_his(:)
        real(dp), allocatable :: new_strip_z_los(:), new_strip_z_his(:)
        integer :: count

        if (allocated(collected)) then
            count = size(collected)
        else
            count = 0
        end if
        allocate(new_collected(count+1), new_strip_ids(count+1), &
            new_strip_r_los(count+1), new_strip_r_his(count+1), &
            new_strip_z_los(count+1), new_strip_z_his(count+1))
        if (count > 0) then
            new_collected(1:count) = collected
            new_strip_ids(1:count) = strip_ids
            new_strip_r_los(1:count) = strip_r_los
            new_strip_r_his(1:count) = strip_r_his
            new_strip_z_los(1:count) = strip_z_los
            new_strip_z_his(1:count) = strip_z_his
        end if
        new_collected(count+1) = leaf
        new_strip_ids(count+1) = strip_index
        new_strip_r_los(count+1) = r_lo
        new_strip_r_his(count+1) = r_hi
        new_strip_z_los(count+1) = z_lo
        new_strip_z_his(count+1) = z_hi
        call move_alloc(new_collected, collected)
        call move_alloc(new_strip_ids, strip_ids)
        call move_alloc(new_strip_r_los, strip_r_los)
        call move_alloc(new_strip_r_his, strip_r_his)
        call move_alloc(new_strip_z_los, strip_z_los)
        call move_alloc(new_strip_z_his, strip_z_his)
    end subroutine append_leaf

end module neort_gc_eqdsk_allowed_region_cut_box
