module neort_gc_eqdsk_cut_graph_atlas
    !! Fail-closed certificate for the regular R-graph of the Buchholz cut.
    !!
    !! This module owns only grid ownership, subdivision, interval sign
    !! classification, array management, bisection, and certificate checks.
    !! The interval evaluator owns every physical enclosure.  Numerical
    !! mapping is permitted only after the complete certificate has closed and
    !! uses the scalar Fortsym-generated cut jet and R-flux chart.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: nrad, nzet, rad, zet, psi_sep
    use neort_gc_eqdsk_cut_interval, only: &
        EQDSK_CUT_INTERVAL_SUCCESS, eqdsk_cut_interval_result_t, &
        evaluate_eqdsk_cut_interval_box
    use neort_gc_eqdsk_cut_jet, only: &
        EQDSK_CUT_JET_SUCCESS, eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
    use neort_eqdsk_cut_r_flux_chart_symbolic, only: &
        evaluate_neort_eqdsk_cut_r_flux_chart
    implicit none
    private

    integer, parameter, public :: EQDSK_CUT_ATLAS_SUCCESS = 0
    integer, parameter, public :: EQDSK_CUT_ATLAS_INVALID_INPUT = 1
    integer, parameter, public :: EQDSK_CUT_ATLAS_UNINITIALIZED = 2
    integer, parameter, public :: EQDSK_CUT_ATLAS_GRID_INVALID = 3
    integer, parameter, public :: EQDSK_CUT_ATLAS_INTERVAL_FAILURE = 4
    integer, parameter, public :: EQDSK_CUT_ATLAS_UNRESOLVED = 5
    integer, parameter, public :: EQDSK_CUT_ATLAS_MULTIPLE_BANDS = 6
    integer, parameter, public :: EQDSK_CUT_ATLAS_INCOMPLETE = 7
    integer, parameter, public :: EQDSK_CUT_ATLAS_INVALID_CERTIFICATE = 8
    integer, parameter, public :: EQDSK_CUT_ATLAS_OUT_OF_RANGE = 9
    integer, parameter, public :: EQDSK_CUT_ATLAS_NO_BRACKET = 10
    integer, parameter, public :: EQDSK_CUT_ATLAS_MAPPING_FAILURE = 11
    integer, parameter, public :: EQDSK_CUT_ATLAS_NONFINITE = 12
    integer, parameter, public :: EQDSK_CUT_GRAPH_CERTIFICATE_ID = 130013

    type, public :: eqdsk_cut_graph_atlas_options_t
        integer :: max_r_depth = 24
        integer :: max_z_depth = 24
        integer :: max_bisection_iterations = 128
        real(dp) :: minimum_r_width = 0.0_dp
        real(dp) :: minimum_z_width = 0.0_dp
        real(dp) :: map_absolute_tolerance = 1.0e-12_dp
        real(dp) :: map_relative_tolerance = 1.0e-12_dp
    end type eqdsk_cut_graph_atlas_options_t

    type, public :: eqdsk_cut_graph_atlas_strip_t
        real(dp) :: r_lo = 0.0_dp
        real(dp) :: r_hi = 0.0_dp
        real(dp) :: z_lo = 0.0_dp
        real(dp) :: z_hi = 0.0_dp
        real(dp) :: dZ_dR_lo = 0.0_dp
        real(dp) :: dZ_dR_hi = 0.0_dp
        real(dp) :: dpsihat_dR_lo = 0.0_dp
        real(dp) :: dpsihat_dR_hi = 0.0_dp
        integer :: cell_R = 0
        integer :: normal_sign = 0
        integer :: candidate_leaf_count = 0
    end type eqdsk_cut_graph_atlas_strip_t

    type, public :: eqdsk_cut_graph_atlas_t
        real(dp) :: requested_r_lo = 0.0_dp
        real(dp) :: requested_r_hi = 0.0_dp
        real(dp) :: requested_z_lo = 0.0_dp
        real(dp) :: requested_z_hi = 0.0_dp
        real(dp) :: requested_psihat_lo = 0.0_dp
        real(dp) :: requested_psihat_hi = 0.0_dp
        type(eqdsk_cut_graph_atlas_options_t) :: options
        type(eqdsk_cut_graph_atlas_strip_t), allocatable :: strips(:)
        integer :: certificate_id = 0
        real(dp) :: raw_psi_sep = 0.0_dp
        logical :: raw_psi_sep_valid = .false.
        logical :: global_completeness_certified = .false.
        integer :: flux_monotonicity_sign = 0
        logical :: flux_monotonicity_certified = .false.
        integer :: failure_cell_R = 0
        integer :: failure_cell_Z = 0
        integer :: failure_r_depth = 0
        integer :: failure_z_depth = 0
        integer :: failure_stage = 0
        real(dp) :: failure_r_lo = 0.0_dp
        real(dp) :: failure_r_hi = 0.0_dp
        real(dp) :: failure_z_lo = 0.0_dp
        real(dp) :: failure_z_hi = 0.0_dp
    end type eqdsk_cut_graph_atlas_t

    public :: build_eqdsk_cut_graph_atlas
    public :: clear_eqdsk_cut_graph_atlas
    public :: validate_eqdsk_cut_graph_atlas
    public :: map_eqdsk_cut_graph_atlas
    public :: map_eqdsk_cut_graph_atlas_flux
    public :: enclose_eqdsk_cut_graph_strip

    integer, parameter :: Z_COVER_SUCCESS = 0
    integer, parameter :: Z_COVER_UNRESOLVED = 1
    integer, parameter :: Z_COVER_MULTIPLE = 2
    integer, parameter :: Z_COVER_FATAL = 3

    type :: candidate_leaf_t
        real(dp) :: z_lo = 0.0_dp
        real(dp) :: z_hi = 0.0_dp
        type(eqdsk_cut_interval_result_t) :: interval
    end type candidate_leaf_t

contains

    subroutine build_eqdsk_cut_graph_atlas(atlas, r_lo, r_hi, z_lo, z_hi, &
            psihat_lo, psihat_hi, options, status)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        real(dp), intent(in) :: r_lo, r_hi, z_lo, z_hi
        real(dp), intent(in) :: psihat_lo, psihat_hi
        type(eqdsk_cut_graph_atlas_options_t), intent(in) :: options
        integer, intent(out) :: status

        integer :: cell_R, local_status
        real(dp) :: slab_lo, slab_hi

        call clear_eqdsk_cut_graph_atlas(atlas)
        status = EQDSK_CUT_ATLAS_INVALID_INPUT
        if (.not. valid_options(options)) return
        if (.not. all(ieee_is_finite([r_lo, r_hi, z_lo, z_hi, &
                psihat_lo, psihat_hi]))) then
            status = EQDSK_CUT_ATLAS_NONFINITE
            return
        end if
        if (r_lo >= r_hi .or. z_lo >= z_hi .or. &
                psihat_lo >= psihat_hi) return
        status = validate_grid()
        if (status /= EQDSK_CUT_ATLAS_SUCCESS) return
        if (.not. ieee_is_finite(psi_sep)) then
            status = EQDSK_CUT_ATLAS_NONFINITE
            return
        end if
        if (psi_sep <= 0.0_dp) return
        if (r_lo < rad(1) .or. r_hi > rad(nrad) .or. &
                z_lo < zet(1) .or. z_hi > zet(nzet)) then
            status = EQDSK_CUT_ATLAS_OUT_OF_RANGE
            return
        end if

        atlas%requested_r_lo = r_lo
        atlas%requested_r_hi = r_hi
        atlas%requested_z_lo = z_lo
        atlas%requested_z_hi = z_hi
        atlas%requested_psihat_lo = psihat_lo
        atlas%requested_psihat_hi = psihat_hi
        atlas%options = options
        atlas%raw_psi_sep = psi_sep
        atlas%raw_psi_sep_valid = .true.

        do cell_R = 1, nrad-1
            slab_lo = max(r_lo, rad(cell_R))
            slab_hi = min(r_hi, rad(cell_R+1))
            if (slab_hi <= slab_lo) cycle
            call certify_r_slab(atlas, cell_R, slab_lo, slab_hi, 0, &
                local_status)
            if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
                call invalidate_eqdsk_cut_graph_atlas(atlas)
                status = local_status
                return
            end if
        end do

        status = validate_atlas_structure(atlas, .false.)
        if (status /= EQDSK_CUT_ATLAS_SUCCESS) then
            call clear_eqdsk_cut_graph_atlas(atlas)
            return
        end if
        call certify_flux_monotonicity(atlas)
        atlas%global_completeness_certified = .true.
        atlas%certificate_id = EQDSK_CUT_GRAPH_CERTIFICATE_ID
        status = validate_atlas_structure(atlas, .true.)
        if (status /= EQDSK_CUT_ATLAS_SUCCESS) then
            call clear_eqdsk_cut_graph_atlas(atlas)
            return
        end if
    end subroutine build_eqdsk_cut_graph_atlas

    subroutine clear_eqdsk_cut_graph_atlas(atlas)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas

        if (allocated(atlas%strips)) deallocate(atlas%strips)
        atlas%requested_r_lo = 0.0_dp
        atlas%requested_r_hi = 0.0_dp
        atlas%requested_z_lo = 0.0_dp
        atlas%requested_z_hi = 0.0_dp
        atlas%requested_psihat_lo = 0.0_dp
        atlas%requested_psihat_hi = 0.0_dp
        atlas%options = eqdsk_cut_graph_atlas_options_t()
        atlas%certificate_id = 0
        atlas%raw_psi_sep = 0.0_dp
        atlas%raw_psi_sep_valid = .false.
        atlas%global_completeness_certified = .false.
        atlas%flux_monotonicity_sign = 0
        atlas%flux_monotonicity_certified = .false.
        atlas%failure_cell_R = 0
        atlas%failure_cell_Z = 0
        atlas%failure_r_depth = 0
        atlas%failure_z_depth = 0
        atlas%failure_stage = 0
        atlas%failure_r_lo = 0.0_dp
        atlas%failure_r_hi = 0.0_dp
        atlas%failure_z_lo = 0.0_dp
        atlas%failure_z_hi = 0.0_dp
    end subroutine clear_eqdsk_cut_graph_atlas

    subroutine validate_eqdsk_cut_graph_atlas(atlas, status)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        integer, intent(out) :: status

        status = validate_atlas_structure(atlas, .true.)
    end subroutine validate_eqdsk_cut_graph_atlas

    subroutine enclose_eqdsk_cut_graph_strip(atlas, strip_index, r_lo, r_hi, &
            enclosures, status)
        type(eqdsk_cut_graph_atlas_t), intent(in) :: atlas
        integer, intent(in) :: strip_index
        real(dp), intent(in) :: r_lo, r_hi
        type(eqdsk_cut_interval_result_t), allocatable, intent(out) :: &
            enclosures(:)
        integer, intent(out) :: status

        call enclose_eqdsk_cut_graph_strip_depth(atlas, strip_index, r_lo, &
            r_hi, 0, enclosures, status)
    end subroutine enclose_eqdsk_cut_graph_strip

    recursive subroutine enclose_eqdsk_cut_graph_strip_depth(atlas, strip_index, &
            r_lo, r_hi, depth, enclosures, status)
        type(eqdsk_cut_graph_atlas_t), intent(in) :: atlas
        integer, intent(in) :: strip_index, depth
        real(dp), intent(in) :: r_lo, r_hi
        type(eqdsk_cut_interval_result_t), allocatable, intent(out) :: &
            enclosures(:)
        integer, intent(out) :: status

        type(eqdsk_cut_graph_atlas_t) :: working
        type(eqdsk_cut_graph_atlas_strip_t) :: certified_strip
        type(candidate_leaf_t), allocatable :: candidates(:)
        type(eqdsk_cut_interval_result_t), allocatable :: left(:), right(:)
        integer :: i, n_candidates, outcome, local_status
        real(dp) :: midpoint

        if (allocated(enclosures)) deallocate(enclosures)
        status = EQDSK_CUT_ATLAS_INVALID_INPUT
        if (.not. all(ieee_is_finite([r_lo, r_hi]))) return
        if (r_hi < r_lo) return
        if (.not. allocated(atlas%strips)) return
        if (strip_index < 1 .or. strip_index > size(atlas%strips)) return
        if (r_lo < atlas%strips(strip_index)%r_lo .or. &
                r_hi > atlas%strips(strip_index)%r_hi) return
        if (atlas%certificate_id /= EQDSK_CUT_GRAPH_CERTIFICATE_ID .or. &
                .not. atlas%global_completeness_certified) then
            status = EQDSK_CUT_ATLAS_INVALID_CERTIFICATE
            return
        end if

        working = atlas
        call collect_slab_candidates(working, &
            atlas%strips(strip_index)%cell_R, r_lo, r_hi, candidates, &
            n_candidates, outcome, local_status)
        if (outcome /= Z_COVER_SUCCESS) then
            status = local_status
            if (status == EQDSK_CUT_ATLAS_SUCCESS) then
                status = EQDSK_CUT_ATLAS_UNRESOLVED
            end if
            return
        end if
        do i = 1, n_candidates
            if (candidates(i)%interval%psi_hat%lo < &
                    atlas%requested_psihat_lo-1.0e-12_dp .or. &
                    candidates(i)%interval%psi_hat%hi > &
                    atlas%requested_psihat_hi+1.0e-12_dp) then
                if (depth >= atlas%options%max_r_depth .or. &
                        r_hi-r_lo <= atlas%options%minimum_r_width) then
                    status = EQDSK_CUT_ATLAS_UNRESOLVED
                    return
                end if
                midpoint = 0.5_dp*(r_lo+r_hi)
                if (midpoint <= r_lo .or. midpoint >= r_hi) then
                    status = EQDSK_CUT_ATLAS_UNRESOLVED
                    return
                end if
                call enclose_eqdsk_cut_graph_strip_depth(atlas, strip_index, &
                    r_lo, midpoint, depth+1, left, local_status)
                if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
                    status = local_status
                    return
                end if
                call enclose_eqdsk_cut_graph_strip_depth(atlas, strip_index, &
                    midpoint, r_hi, depth+1, right, local_status)
                if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
                    status = local_status
                    return
                end if
                allocate(enclosures(size(left)+size(right)))
                if (size(left) > 0) enclosures(1:size(left)) = left
                if (size(right) > 0) enclosures(size(left)+1:) = right
                return
            end if
        end do
        call assemble_candidate_band(working, &
            atlas%strips(strip_index)%cell_R, r_lo, r_hi, candidates, &
            n_candidates, certified_strip, outcome, local_status)
        if (outcome /= Z_COVER_SUCCESS) then
            status = local_status
            if (status == EQDSK_CUT_ATLAS_SUCCESS) then
                status = EQDSK_CUT_ATLAS_UNRESOLVED
            end if
            return
        end if
        allocate(enclosures(n_candidates))
        do i = 1, n_candidates
            enclosures(i) = candidates(i)%interval
        end do
        status = EQDSK_CUT_ATLAS_SUCCESS
    end subroutine enclose_eqdsk_cut_graph_strip_depth

    subroutine map_eqdsk_cut_graph_atlas(atlas, radius, position, &
            dposition_dR, dZ_dR, dpsihat_dR, status)
        type(eqdsk_cut_graph_atlas_t), intent(in) :: atlas
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: position(3)
        real(dp), intent(out) :: dposition_dR(3)
        real(dp), intent(out), optional :: dZ_dR, dpsihat_dR
        integer, intent(out) :: status

        integer :: strip_index, iteration, jet_status
        real(dp) :: left, right, midpoint, f_left, f_right, f_midpoint
        real(dp) :: z_root, tolerance, dZ_root, dpsihat_root
        type(eqdsk_cut_jet_t) :: jet
        logical :: converged, have_midpoint

        position = 0.0_dp
        dposition_dR = 0.0_dp
        if (present(dZ_dR)) dZ_dR = 0.0_dp
        if (present(dpsihat_dR)) dpsihat_dR = 0.0_dp
        status = validate_atlas_structure(atlas, .true.)
        if (status /= EQDSK_CUT_ATLAS_SUCCESS) return
        if (.not. ieee_is_finite(radius)) then
            status = EQDSK_CUT_ATLAS_NONFINITE
            return
        end if

        strip_index = locate_strip(atlas, radius)
        if (strip_index < 1) then
            status = EQDSK_CUT_ATLAS_OUT_OF_RANGE
            return
        end if
        left = atlas%strips(strip_index)%z_lo
        right = atlas%strips(strip_index)%z_hi
        call evaluate_scalar_numerator(radius, left, f_left, jet_status)
        if (jet_status /= EQDSK_CUT_JET_SUCCESS) then
            status = EQDSK_CUT_ATLAS_MAPPING_FAILURE
            return
        end if
        call evaluate_scalar_numerator(radius, right, f_right, jet_status)
        if (jet_status /= EQDSK_CUT_JET_SUCCESS) then
            status = EQDSK_CUT_ATLAS_MAPPING_FAILURE
            return
        end if
        if (.not. strict_bracket(atlas%strips(strip_index)%normal_sign, &
                f_left, f_right)) then
            status = EQDSK_CUT_ATLAS_NO_BRACKET
            return
        end if

        midpoint = 0.5_dp*(left+right)
        f_midpoint = 0.0_dp
        converged = .false.
        have_midpoint = .false.
        do iteration = 1, atlas%options%max_bisection_iterations
            midpoint = 0.5_dp*(left+right)
            call evaluate_scalar_numerator(radius, midpoint, f_midpoint, &
                jet_status)
            if (jet_status /= EQDSK_CUT_JET_SUCCESS) then
                status = EQDSK_CUT_ATLAS_MAPPING_FAILURE
                return
            end if
            have_midpoint = .true.
            tolerance = max(atlas%options%map_absolute_tolerance, &
                atlas%options%map_relative_tolerance*max(1.0_dp, &
                abs(midpoint)))
            if (abs(right-left) <= tolerance .or. f_midpoint == 0.0_dp) then
                converged = .true.
                exit
            end if
            if (atlas%strips(strip_index)%normal_sign > 0) then
                if (f_midpoint < 0.0_dp) then
                    left = midpoint
                else
                    right = midpoint
                end if
            else
                if (f_midpoint > 0.0_dp) then
                    left = midpoint
                else
                    right = midpoint
                end if
            end if
        end do
        if (.not. have_midpoint .or. .not. converged) then
            status = EQDSK_CUT_ATLAS_MAPPING_FAILURE
            return
        end if
        z_root = midpoint
        call evaluate_scalar_jet(radius, z_root, jet, jet_status)
        if (jet_status /= EQDSK_CUT_JET_SUCCESS) then
            status = EQDSK_CUT_ATLAS_MAPPING_FAILURE
            return
        end if
        call evaluate_neort_eqdsk_cut_r_flux_chart( &
            jet%d_cut_numerator_d_R, jet%d_cut_numerator_d_Z, &
            jet%psi_jet(2), jet%psi_jet(3), atlas%raw_psi_sep, dZ_root, &
            dpsihat_root)
        if (.not. all(ieee_is_finite([dZ_root, dpsihat_root]))) then
            status = EQDSK_CUT_ATLAS_MAPPING_FAILURE
            return
        end if
        position = [radius, z_root, 0.0_dp]
        dposition_dR = [1.0_dp, dZ_root, 0.0_dp]
        if (present(dZ_dR)) dZ_dR = dZ_root
        if (present(dpsihat_dR)) dpsihat_dR = dpsihat_root
        status = EQDSK_CUT_ATLAS_SUCCESS
    end subroutine map_eqdsk_cut_graph_atlas

    subroutine map_eqdsk_cut_graph_atlas_flux(atlas, target_psihat, position, &
            dposition_dR, dZ_dR, dpsihat_dR, status)
        type(eqdsk_cut_graph_atlas_t), intent(in) :: atlas
        real(dp), intent(in) :: target_psihat
        real(dp), intent(out) :: position(3), dposition_dR(3)
        real(dp), intent(out) :: dZ_dR, dpsihat_dR
        integer, intent(out) :: status

        real(dp) :: left, right, midpoint, flux_left, flux_right, flux_midpoint
        real(dp) :: tolerance
        integer :: iteration, local_status
        logical :: converged

        position = 0.0_dp
        dposition_dR = 0.0_dp
        dZ_dR = 0.0_dp
        dpsihat_dR = 0.0_dp
        status = validate_atlas_structure(atlas, .true.)
        if (status /= EQDSK_CUT_ATLAS_SUCCESS) return
        if (.not. atlas%flux_monotonicity_certified) then
            status = EQDSK_CUT_ATLAS_INVALID_CERTIFICATE
            return
        end if
        if (.not. ieee_is_finite(target_psihat)) then
            status = EQDSK_CUT_ATLAS_NONFINITE
            return
        end if

        left = atlas%requested_r_lo
        right = atlas%requested_r_hi
        call evaluate_graph_flux(atlas, left, flux_left, local_status)
        if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
            status = local_status
            return
        end if
        call evaluate_graph_flux(atlas, right, flux_right, local_status)
        if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
            status = local_status
            return
        end if
        if (target_psihat < min(flux_left, flux_right) .or. &
                target_psihat > max(flux_left, flux_right)) then
            status = EQDSK_CUT_ATLAS_OUT_OF_RANGE
            return
        end if

        converged = .false.
        midpoint = left
        do iteration = 1, atlas%options%max_bisection_iterations
            midpoint = left+0.5_dp*(right-left)
            call evaluate_graph_flux(atlas, midpoint, flux_midpoint, &
                local_status)
            if (local_status /= EQDSK_CUT_ATLAS_SUCCESS) then
                status = local_status
                return
            end if
            tolerance = max(atlas%options%map_absolute_tolerance, &
                atlas%options%map_relative_tolerance*max(1.0_dp, &
                abs(midpoint)))
            if (abs(right-left) <= tolerance .or. &
                    flux_midpoint == target_psihat) then
                converged = .true.
                exit
            end if
            if ((flux_midpoint < target_psihat) .eqv. &
                    (atlas%flux_monotonicity_sign > 0)) then
                left = midpoint
            else
                right = midpoint
            end if
        end do
        if (.not. converged) then
            status = EQDSK_CUT_ATLAS_MAPPING_FAILURE
            return
        end if
        call map_eqdsk_cut_graph_atlas(atlas, midpoint, position, &
            dposition_dR, dZ_dR, dpsihat_dR, status)
    end subroutine map_eqdsk_cut_graph_atlas_flux

    subroutine evaluate_graph_flux(atlas, radius, psihat, status)
        type(eqdsk_cut_graph_atlas_t), intent(in) :: atlas
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: psihat
        integer, intent(out) :: status

        type(eqdsk_cut_jet_t) :: jet
        real(dp) :: position(3), tangent(3)
        integer :: jet_status

        psihat = 0.0_dp
        call map_eqdsk_cut_graph_atlas(atlas, radius, position, tangent, &
            status=status)
        if (status /= EQDSK_CUT_ATLAS_SUCCESS) return
        call evaluate_scalar_jet(position(1), position(2), jet, jet_status)
        if (jet_status /= EQDSK_CUT_JET_SUCCESS) then
            status = EQDSK_CUT_ATLAS_MAPPING_FAILURE
            return
        end if
        psihat = jet%psi_jet(1)/atlas%raw_psi_sep
        if (.not. ieee_is_finite(psihat)) then
            psihat = 0.0_dp
            status = EQDSK_CUT_ATLAS_NONFINITE
            return
        end if
        status = EQDSK_CUT_ATLAS_SUCCESS
    end subroutine evaluate_graph_flux

    recursive subroutine certify_r_slab(atlas, cell_R, r_lo, r_hi, depth, &
            status)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        integer, intent(in) :: cell_R, depth
        real(dp), intent(in) :: r_lo, r_hi
        integer, intent(out) :: status

        integer :: outcome
        real(dp) :: midpoint
        type(eqdsk_cut_graph_atlas_strip_t) :: strip

        call certify_slab_contents(atlas, cell_R, r_lo, r_hi, outcome, strip, &
            status)
        if (outcome == Z_COVER_FATAL) return
        if (outcome == Z_COVER_SUCCESS) then
            call append_strip(atlas, strip)
            status = EQDSK_CUT_ATLAS_SUCCESS
            return
        end if
        if (depth >= atlas%options%max_r_depth .or. &
                r_hi-r_lo <= atlas%options%minimum_r_width) then
            atlas%failure_cell_R = cell_R
            atlas%failure_r_depth = depth
            atlas%failure_r_lo = r_lo
            atlas%failure_r_hi = r_hi
            if (outcome == Z_COVER_MULTIPLE) then
                status = EQDSK_CUT_ATLAS_MULTIPLE_BANDS
            else
                status = EQDSK_CUT_ATLAS_UNRESOLVED
            end if
            return
        end if
        midpoint = 0.5_dp*(r_lo+r_hi)
        if (midpoint <= r_lo .or. midpoint >= r_hi) then
            status = EQDSK_CUT_ATLAS_UNRESOLVED
            return
        end if
        call certify_r_slab(atlas, cell_R, r_lo, midpoint, depth+1, status)
        if (status /= EQDSK_CUT_ATLAS_SUCCESS) return
        call certify_r_slab(atlas, cell_R, midpoint, r_hi, depth+1, status)
    end subroutine certify_r_slab

    subroutine certify_slab_contents(atlas, cell_R, r_lo, r_hi, outcome, &
            strip, status)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        integer, intent(in) :: cell_R
        real(dp), intent(in) :: r_lo, r_hi
        integer, intent(out) :: outcome
        type(eqdsk_cut_graph_atlas_strip_t), intent(out) :: strip
        integer, intent(out) :: status

        integer :: n_candidates
        type(candidate_leaf_t), allocatable :: candidates(:)

        strip = eqdsk_cut_graph_atlas_strip_t()
        outcome = Z_COVER_UNRESOLVED
        status = EQDSK_CUT_ATLAS_SUCCESS
        call collect_slab_candidates(atlas, cell_R, r_lo, r_hi, candidates, &
            n_candidates, outcome, status)
        if (outcome /= Z_COVER_SUCCESS) return
        call assemble_candidate_band(atlas, cell_R, r_lo, r_hi, candidates, &
            n_candidates, strip, outcome, status)
        if (allocated(candidates)) deallocate(candidates)
    end subroutine certify_slab_contents

    subroutine collect_slab_candidates(atlas, cell_R, r_lo, r_hi, candidates, &
            n_candidates, outcome, status)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        integer, intent(in) :: cell_R
        real(dp), intent(in) :: r_lo, r_hi
        type(candidate_leaf_t), allocatable, intent(out) :: candidates(:)
        integer, intent(out) :: n_candidates, outcome, status

        integer :: cell_Z, local_status, local_outcome
        real(dp) :: box_lo, box_hi

        if (allocated(candidates)) deallocate(candidates)
        n_candidates = 0
        outcome = Z_COVER_UNRESOLVED
        status = EQDSK_CUT_ATLAS_SUCCESS
        do cell_Z = 1, nzet-1
            box_lo = max(atlas%requested_z_lo, zet(cell_Z))
            box_hi = min(atlas%requested_z_hi, zet(cell_Z+1))
            if (box_hi <= box_lo) cycle
            call cover_z_box(atlas, cell_R, cell_Z, r_lo, r_hi, box_lo, &
                box_hi, 0, candidates, n_candidates, local_outcome, &
                local_status)
            if (local_outcome == Z_COVER_FATAL) then
                outcome = Z_COVER_FATAL
                status = local_status
                return
            end if
            if (local_outcome /= Z_COVER_SUCCESS) then
                outcome = local_outcome
                return
            end if
        end do
        if (n_candidates <= 0) then
            atlas%failure_stage = 2
            return
        end if
        outcome = Z_COVER_SUCCESS
    end subroutine collect_slab_candidates

    recursive subroutine cover_z_box(atlas, cell_R, cell_Z, r_lo, r_hi, z_lo, &
            z_hi, depth, candidates, n_candidates, outcome, status)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        integer, intent(in) :: cell_R, cell_Z, depth
        real(dp), intent(in) :: r_lo, r_hi, z_lo, z_hi
        type(candidate_leaf_t), allocatable, intent(inout) :: candidates(:)
        integer, intent(inout) :: n_candidates
        integer, intent(out) :: outcome, status

        integer :: interval_status
        real(dp) :: midpoint
        type(eqdsk_cut_interval_result_t) :: interval
        type(candidate_leaf_t) :: leaf

        outcome = Z_COVER_UNRESOLVED
        status = EQDSK_CUT_ATLAS_SUCCESS
        call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, r_lo, r_hi, &
            z_lo, z_hi, interval, interval_status)
        if (interval_status /= EQDSK_CUT_INTERVAL_SUCCESS) then
            outcome = Z_COVER_FATAL
            status = EQDSK_CUT_ATLAS_INTERVAL_FAILURE
            return
        end if
        ! Buchholz's cut is the segment inside an explicit boundary flux
        ! surface, not the complete zero set in the rectangular EQDSK spline
        ! box.  Boxes certified outside that physical flux interval cannot
        ! contribute a cut component.  Boundary-straddling boxes remain
        ! unresolved until N or the flux enclosure separates them.
        if (interval%psi_hat%hi < atlas%requested_psihat_lo .or. &
                interval%psi_hat%lo > atlas%requested_psihat_hi) then
            outcome = Z_COVER_SUCCESS
            return
        end if
        if (.not. interval_contains_zero(interval%numerator)) then
            outcome = Z_COVER_SUCCESS
            return
        end if
        if (regular_candidate(interval)) then
            leaf%z_lo = z_lo
            leaf%z_hi = z_hi
            leaf%interval = interval
            call append_candidate(candidates, n_candidates, leaf)
            outcome = Z_COVER_SUCCESS
            return
        end if
        if (depth >= atlas%options%max_z_depth .or. &
                z_hi-z_lo <= atlas%options%minimum_z_width) then
            atlas%failure_cell_Z = cell_Z
            atlas%failure_z_depth = depth
            atlas%failure_stage = 1
            atlas%failure_z_lo = z_lo
            atlas%failure_z_hi = z_hi
            outcome = Z_COVER_UNRESOLVED
            return
        end if
        midpoint = 0.5_dp*(z_lo+z_hi)
        if (midpoint <= z_lo .or. midpoint >= z_hi) then
            outcome = Z_COVER_UNRESOLVED
            return
        end if
        call cover_z_box(atlas, cell_R, cell_Z, r_lo, r_hi, z_lo, midpoint, &
            depth+1, candidates, n_candidates, outcome, status)
        if (outcome /= Z_COVER_SUCCESS) return
        call cover_z_box(atlas, cell_R, cell_Z, r_lo, r_hi, midpoint, z_hi, &
            depth+1, candidates, n_candidates, outcome, status)
    end subroutine cover_z_box

    subroutine assemble_candidate_band(atlas, cell_R, r_lo, r_hi, candidates, &
            n_candidates, strip, outcome, status)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        integer, intent(in) :: cell_R, n_candidates
        real(dp), intent(in) :: r_lo, r_hi
        type(candidate_leaf_t), intent(in) :: candidates(:)
        type(eqdsk_cut_graph_atlas_strip_t), intent(out) :: strip
        integer, intent(out) :: outcome, status

        integer :: i, normal_sign, face_cell, face_status
        real(dp) :: band_lo, band_hi
        type(eqdsk_cut_interval_result_t) :: lower_face, upper_face

        strip = eqdsk_cut_graph_atlas_strip_t()
        outcome = Z_COVER_UNRESOLVED
        status = EQDSK_CUT_ATLAS_SUCCESS
        band_lo = candidates(1)%z_lo
        band_hi = candidates(n_candidates)%z_hi
        do i = 2, n_candidates
            if (candidates(i)%z_lo /= candidates(i-1)%z_hi) then
                atlas%failure_stage = 3
                outcome = Z_COVER_MULTIPLE
                return
            end if
        end do
        if (candidates(1)%interval%numerator_Z%lo > 0.0_dp) then
            normal_sign = 1
        else if (candidates(1)%interval%numerator_Z%hi < 0.0_dp) then
            normal_sign = -1
        else
            atlas%failure_stage = 4
            outcome = Z_COVER_UNRESOLVED
            return
        end if
        do i = 1, n_candidates
            if (.not. regular_candidate(candidates(i)%interval)) then
                outcome = Z_COVER_UNRESOLVED
                return
            end if
            if (normal_sign > 0) then
                if (candidates(i)%interval%numerator_Z%lo <= 0.0_dp) then
                    atlas%failure_stage = 4
                    outcome = Z_COVER_UNRESOLVED
                    return
                end if
            else
                if (candidates(i)%interval%numerator_Z%hi >= 0.0_dp) then
                    atlas%failure_stage = 4
                    outcome = Z_COVER_UNRESOLVED
                    return
                end if
            end if
        end do
        call locate_z_cell(band_lo, face_cell)
        if (face_cell < 1) then
            status = EQDSK_CUT_ATLAS_GRID_INVALID
            outcome = Z_COVER_FATAL
            return
        end if
        call evaluate_eqdsk_cut_interval_box(cell_R, face_cell, r_lo, r_hi, &
            band_lo, band_lo, lower_face, face_status)
        if (face_status /= EQDSK_CUT_INTERVAL_SUCCESS) then
            status = EQDSK_CUT_ATLAS_INTERVAL_FAILURE
            outcome = Z_COVER_FATAL
            return
        end if
        call locate_z_cell(band_hi, face_cell)
        if (face_cell < 1) then
            status = EQDSK_CUT_ATLAS_GRID_INVALID
            outcome = Z_COVER_FATAL
            return
        end if
        call evaluate_eqdsk_cut_interval_box(cell_R, face_cell, r_lo, r_hi, &
            band_hi, band_hi, upper_face, face_status)
        if (face_status /= EQDSK_CUT_INTERVAL_SUCCESS) then
            status = EQDSK_CUT_ATLAS_INTERVAL_FAILURE
            outcome = Z_COVER_FATAL
            return
        end if
        if (.not. opposite_face_signs(normal_sign, lower_face%numerator, &
                upper_face%numerator)) then
            atlas%failure_stage = 5
            atlas%failure_z_lo = band_lo
            atlas%failure_z_hi = band_hi
            outcome = Z_COVER_UNRESOLVED
            return
        end if

        strip%r_lo = r_lo
        strip%r_hi = r_hi
        strip%z_lo = band_lo
        strip%z_hi = band_hi
        strip%cell_R = cell_R
        strip%normal_sign = normal_sign
        strip%candidate_leaf_count = n_candidates
        strip%dZ_dR_lo = candidates(1)%interval%dZ_dR%lo
        strip%dZ_dR_hi = candidates(1)%interval%dZ_dR%hi
        strip%dpsihat_dR_lo = candidates(1)%interval%dpsihat_dR%lo
        strip%dpsihat_dR_hi = candidates(1)%interval%dpsihat_dR%hi
        do i = 2, n_candidates
            strip%dZ_dR_lo = min(strip%dZ_dR_lo, &
                candidates(i)%interval%dZ_dR%lo)
            strip%dZ_dR_hi = max(strip%dZ_dR_hi, &
                candidates(i)%interval%dZ_dR%hi)
            strip%dpsihat_dR_lo = min(strip%dpsihat_dR_lo, &
                candidates(i)%interval%dpsihat_dR%lo)
            strip%dpsihat_dR_hi = max(strip%dpsihat_dR_hi, &
                candidates(i)%interval%dpsihat_dR%hi)
        end do
        outcome = Z_COVER_SUCCESS
    end subroutine assemble_candidate_band

    logical function regular_candidate(interval)
        type(eqdsk_cut_interval_result_t), intent(in) :: interval

        regular_candidate = interval%positive_denominator%lo > 0.0_dp
        if (.not. regular_candidate) return
        if (.not. interval%denominator_positive_certified) then
            regular_candidate = .false.
            return
        end if
        regular_candidate = interval%r_chart_certified
        if (.not. regular_candidate) return
        if (.not. strict_interval_sign(interval%numerator_Z)) then
            regular_candidate = .false.
        end if
    end function regular_candidate

    logical function strict_interval_sign(interval)
        use neort_gc_outward_interval, only: gc_outward_interval_t
        type(gc_outward_interval_t), intent(in) :: interval

        strict_interval_sign = interval%lo > 0.0_dp .or. &
            interval%hi < 0.0_dp
    end function strict_interval_sign

    logical function interval_contains_zero(interval)
        use neort_gc_outward_interval, only: gc_outward_interval_t
        type(gc_outward_interval_t), intent(in) :: interval

        interval_contains_zero = interval%lo <= 0.0_dp .and. &
            interval%hi >= 0.0_dp
    end function interval_contains_zero

    logical function opposite_face_signs(normal_sign, lower, upper)
        use neort_gc_outward_interval, only: gc_outward_interval_t
        integer, intent(in) :: normal_sign
        type(gc_outward_interval_t), intent(in) :: lower, upper

        opposite_face_signs = .false.
        if (normal_sign > 0) then
            if (lower%hi >= 0.0_dp) return
            if (upper%lo <= 0.0_dp) return
        else
            if (lower%lo <= 0.0_dp) return
            if (upper%hi >= 0.0_dp) return
        end if
        opposite_face_signs = .true.
    end function opposite_face_signs

    logical function strict_bracket(normal_sign, lower, upper)
        integer, intent(in) :: normal_sign
        real(dp), intent(in) :: lower, upper

        strict_bracket = .false.
        if (normal_sign > 0) then
            if (lower >= 0.0_dp .or. upper <= 0.0_dp) return
        else
            if (lower <= 0.0_dp .or. upper >= 0.0_dp) return
        end if
        strict_bracket = .true.
    end function strict_bracket

    subroutine evaluate_scalar_numerator(radius, height, numerator, status)
        real(dp), intent(in) :: radius, height
        real(dp), intent(out) :: numerator
        integer, intent(out) :: status

        type(eqdsk_cut_jet_t) :: jet

        call evaluate_scalar_jet(radius, height, jet, status)
        numerator = jet%cut_numerator
    end subroutine evaluate_scalar_numerator

    subroutine evaluate_scalar_jet(radius, height, jet, status)
        real(dp), intent(in) :: radius, height
        type(eqdsk_cut_jet_t), intent(out) :: jet
        integer, intent(out) :: status

        call evaluate_eqdsk_cut_jet([radius, height, 0.0_dp], 1.0_dp, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], jet, status)
    end subroutine evaluate_scalar_jet

    integer function locate_strip(atlas, radius)
        type(eqdsk_cut_graph_atlas_t), intent(in) :: atlas
        real(dp), intent(in) :: radius
        integer :: i

        locate_strip = 0
        if (.not. allocated(atlas%strips)) return
        do i = 1, size(atlas%strips)
            if (radius >= atlas%strips(i)%r_lo .and. &
                    radius <= atlas%strips(i)%r_hi) then
                locate_strip = i
                return
            end if
        end do
    end function locate_strip

    subroutine locate_z_cell(height, cell_Z)
        real(dp), intent(in) :: height
        integer, intent(out) :: cell_Z
        integer :: i

        cell_Z = 0
        do i = 1, nzet-1
            if (height >= zet(i) .and. height <= zet(i+1)) then
                cell_Z = i
                return
            end if
        end do
    end subroutine locate_z_cell

    subroutine append_candidate(candidates, count, candidate)
        type(candidate_leaf_t), allocatable, intent(inout) :: candidates(:)
        integer, intent(inout) :: count
        type(candidate_leaf_t), intent(in) :: candidate

        type(candidate_leaf_t), allocatable :: extended(:)

        allocate(extended(count+1))
        if (count > 0) extended(1:count) = candidates
        extended(count+1) = candidate
        call move_alloc(extended, candidates)
        count = count+1
    end subroutine append_candidate

    subroutine append_strip(atlas, strip)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas
        type(eqdsk_cut_graph_atlas_strip_t), intent(in) :: strip

        type(eqdsk_cut_graph_atlas_strip_t), allocatable :: extended(:)
        integer :: count

        if (allocated(atlas%strips)) then
            count = size(atlas%strips)
        else
            count = 0
        end if
        allocate(extended(count+1))
        if (count > 0) extended(1:count) = atlas%strips
        extended(count+1) = strip
        call move_alloc(extended, atlas%strips)
    end subroutine append_strip

    subroutine invalidate_eqdsk_cut_graph_atlas(atlas)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas

        if (allocated(atlas%strips)) deallocate(atlas%strips)
        atlas%certificate_id = 0
        atlas%raw_psi_sep = 0.0_dp
        atlas%raw_psi_sep_valid = .false.
        atlas%global_completeness_certified = .false.
        atlas%flux_monotonicity_sign = 0
        atlas%flux_monotonicity_certified = .false.
    end subroutine invalidate_eqdsk_cut_graph_atlas

    subroutine certify_flux_monotonicity(atlas)
        type(eqdsk_cut_graph_atlas_t), intent(inout) :: atlas

        logical :: all_negative, all_positive
        integer :: i

        atlas%flux_monotonicity_sign = 0
        atlas%flux_monotonicity_certified = .false.
        if (.not. allocated(atlas%strips)) return
        if (size(atlas%strips) < 1) return
        all_negative = .true.
        all_positive = .true.
        do i = 1, size(atlas%strips)
            all_negative = all_negative .and. &
                atlas%strips(i)%dpsihat_dR_hi < 0.0_dp
            all_positive = all_positive .and. &
                atlas%strips(i)%dpsihat_dR_lo > 0.0_dp
        end do
        if (all_negative) then
            atlas%flux_monotonicity_sign = -1
            atlas%flux_monotonicity_certified = .true.
        else if (all_positive) then
            atlas%flux_monotonicity_sign = 1
            atlas%flux_monotonicity_certified = .true.
        end if
    end subroutine certify_flux_monotonicity

    integer function validate_grid()
        integer :: i

        validate_grid = EQDSK_CUT_ATLAS_GRID_INVALID
        if (nrad < 2 .or. nzet < 2) return
        if (.not. allocated(rad)) return
        if (.not. allocated(zet)) return
        if (size(rad) /= nrad .or. size(zet) /= nzet) return
        if (.not. all(ieee_is_finite(rad))) return
        if (.not. all(ieee_is_finite(zet))) return
        do i = 1, nrad-1
            if (rad(i+1) <= rad(i)) return
        end do
        do i = 1, nzet-1
            if (zet(i+1) <= zet(i)) return
        end do
        validate_grid = EQDSK_CUT_ATLAS_SUCCESS
    end function validate_grid

    logical function valid_options(options)
        type(eqdsk_cut_graph_atlas_options_t), intent(in) :: options

        valid_options = options%max_r_depth >= 0 .and. &
            options%max_z_depth >= 0 .and. &
            options%max_bisection_iterations > 0 .and. &
            ieee_is_finite(options%minimum_r_width) .and. &
            ieee_is_finite(options%minimum_z_width) .and. &
            ieee_is_finite(options%map_absolute_tolerance) .and. &
            ieee_is_finite(options%map_relative_tolerance) .and. &
            options%minimum_r_width >= 0.0_dp .and. &
            options%minimum_z_width >= 0.0_dp .and. &
            options%map_absolute_tolerance > 0.0_dp .and. &
            options%map_relative_tolerance >= 0.0_dp
    end function valid_options

    integer function validate_atlas_structure(atlas, require_global)
        type(eqdsk_cut_graph_atlas_t), intent(in) :: atlas
        logical, intent(in) :: require_global
        integer :: i
        real(dp) :: previous_r

        validate_atlas_structure = EQDSK_CUT_ATLAS_INVALID_CERTIFICATE
        if (validate_grid() /= EQDSK_CUT_ATLAS_SUCCESS) return
        if (.not. valid_options(atlas%options)) return
        validate_atlas_structure = validate_source_provenance(atlas)
        if (validate_atlas_structure /= EQDSK_CUT_ATLAS_SUCCESS) return
        if (.not. all(ieee_is_finite([atlas%requested_r_lo, &
                atlas%requested_r_hi, atlas%requested_z_lo, &
                atlas%requested_z_hi, atlas%requested_psihat_lo, &
                atlas%requested_psihat_hi]))) then
            validate_atlas_structure = EQDSK_CUT_ATLAS_NONFINITE
            return
        end if
        if (atlas%requested_r_lo >= atlas%requested_r_hi .or. &
                atlas%requested_z_lo >= atlas%requested_z_hi .or. &
                atlas%requested_psihat_lo >= atlas%requested_psihat_hi) return
        if (atlas%requested_r_lo < rad(1) .or. &
                atlas%requested_r_hi > rad(nrad) .or. &
                atlas%requested_z_lo < zet(1) .or. &
                atlas%requested_z_hi > zet(nzet)) return
        if (require_global .and. .not. &
                atlas%global_completeness_certified) then
            validate_atlas_structure = EQDSK_CUT_ATLAS_INCOMPLETE
            return
        end if
        if (require_global .and. &
                atlas%certificate_id /= EQDSK_CUT_GRAPH_CERTIFICATE_ID) return
        if (atlas%flux_monotonicity_certified) then
            if (abs(atlas%flux_monotonicity_sign) /= 1) return
        else
            if (atlas%flux_monotonicity_sign /= 0) return
        end if
        if (.not. allocated(atlas%strips)) return
        if (size(atlas%strips) < 1) return
        previous_r = atlas%requested_r_lo
        do i = 1, size(atlas%strips)
            if (.not. valid_strip(atlas%strips(i))) return
            if (atlas%flux_monotonicity_certified) then
                if (atlas%flux_monotonicity_sign < 0 .and. &
                        atlas%strips(i)%dpsihat_dR_hi >= 0.0_dp) return
                if (atlas%flux_monotonicity_sign > 0 .and. &
                        atlas%strips(i)%dpsihat_dR_lo <= 0.0_dp) return
            end if
            if (atlas%strips(i)%r_lo /= previous_r) return
            previous_r = atlas%strips(i)%r_hi
        end do
        if (previous_r /= atlas%requested_r_hi) return
        validate_atlas_structure = EQDSK_CUT_ATLAS_SUCCESS
    end function validate_atlas_structure

    integer function validate_source_provenance(atlas)
        type(eqdsk_cut_graph_atlas_t), intent(in) :: atlas

        validate_source_provenance = EQDSK_CUT_ATLAS_INVALID_CERTIFICATE
        if (.not. atlas%raw_psi_sep_valid) return
        if (.not. ieee_is_finite(atlas%raw_psi_sep)) then
            validate_source_provenance = EQDSK_CUT_ATLAS_NONFINITE
            return
        end if
        if (atlas%raw_psi_sep <= 0.0_dp) return
        if (.not. ieee_is_finite(psi_sep)) then
            validate_source_provenance = EQDSK_CUT_ATLAS_NONFINITE
            return
        end if
        if (psi_sep <= 0.0_dp) return
        if (atlas%raw_psi_sep /= psi_sep) return
        validate_source_provenance = EQDSK_CUT_ATLAS_SUCCESS
    end function validate_source_provenance

    logical function valid_strip(strip)
        type(eqdsk_cut_graph_atlas_strip_t), intent(in) :: strip

        valid_strip = .false.
        if (strip%cell_R < 1 .or. strip%cell_R >= nrad) return
        if (strip%normal_sign /= 1 .and. strip%normal_sign /= -1) return
        if (strip%candidate_leaf_count < 1) return
        if (strip%r_lo >= strip%r_hi .or. strip%z_lo >= strip%z_hi) return
        if (strip%r_lo < rad(strip%cell_R) .or. &
                strip%r_hi > rad(strip%cell_R+1)) return
        if (strip%z_lo < zet(1) .or. strip%z_hi > zet(nzet)) return
        if (.not. all(ieee_is_finite([strip%r_lo, strip%r_hi, &
                strip%z_lo, strip%z_hi, strip%dZ_dR_lo, strip%dZ_dR_hi, &
                strip%dpsihat_dR_lo, strip%dpsihat_dR_hi]))) return
        if (strip%dZ_dR_lo > strip%dZ_dR_hi) return
        if (strip%dpsihat_dR_lo > strip%dpsihat_dR_hi) return
        valid_strip = .true.
    end function valid_strip

end module neort_gc_eqdsk_cut_graph_atlas
