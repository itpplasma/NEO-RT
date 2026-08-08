module gc_poincare_cut_atlas
    !! Contracts for the Buchholz Poincare-cut construction.
    !!
    !! This module is deliberately a validation boundary.  It does not turn
    !! caller-set flags or integer labels into proofs.  Complete branches,
    !! interval connectedness, Cdot simplicity, the exactly-two theorem, state
    !! return tolerances, and homoclinic asymptotics must be supplied by the
    !! corresponding Fortsym-generated providers.  Missing providers fail
    !! closed with PCA_PROVIDER_UNAVAILABLE.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    private

    integer, parameter, public :: PCA_SUCCESS = 0
    integer, parameter, public :: PCA_INVALID_CONTRACT = 1
    integer, parameter, public :: PCA_CALLBACK_FAILURE = 2
    integer, parameter, public :: PCA_UNRESOLVED_TANGENCY = 3
    integer, parameter, public :: PCA_UNKNOWN_MULTIPLICITY = 4
    integer, parameter, public :: PCA_INCOMPLETE_ATLAS = 5
    integer, parameter, public :: PCA_PROVIDER_UNAVAILABLE = 6
    integer, parameter, public :: PCA_INVALID_STATE = 7
    integer, parameter, public :: PCA_INVALID_CERTIFICATE = 8

    integer, parameter, public :: PCA_WALL_NONE = 0
    integer, parameter, public :: PCA_WALL_HIT = 1
    integer, parameter, public :: PCA_WALL_UNRESOLVED = 2

    integer, parameter, public :: PCA_DOMAIN_NONE = 0
    integer, parameter, public :: PCA_DOMAIN_EXIT = 1
    integer, parameter, public :: PCA_DOMAIN_UNRESOLVED = 2

    integer, parameter, public :: PCA_BOUNDARY_REGULAR = 1
    integer, parameter, public :: PCA_BOUNDARY_REFLECTING = 2
    integer, parameter, public :: PCA_BOUNDARY_HOMOCLINIC_X = 3
    integer, parameter, public :: PCA_BOUNDARY_HOMOCLINIC_USUAL = 4

    type, public :: poincare_cut_branch_t
        integer :: branch_id = 0
        integer :: component_identity_seed = 0
        integer :: field_input_id = 0
        integer :: equilibrium_input_id = 0
        integer :: orientation = 0
        character(len=32) :: parameter_units = ''
        real(dp) :: parameter_lower = 0.0_dp
        real(dp) :: parameter_upper = 0.0_dp
        real(dp) :: endpoint_start(2) = 0.0_dp
        real(dp) :: endpoint_end(2) = 0.0_dp
        real(dp) :: cut_residual_tolerance = 0.0_dp
        logical :: endpoint_start_known = .false.
        logical :: endpoint_end_known = .false.
        logical :: physical_r_ordered = .false.
        logical :: continuous_inboard_axis_outboard = .false.
        logical :: outboard_only = .false.
    end type poincare_cut_branch_t

    type, public :: poincare_cut_atlas_t
        type(poincare_cut_branch_t), allocatable :: branches(:)
        integer :: expected_branch_count = 0
        logical :: symmetric_midplane_specialization = .false.
    end type poincare_cut_atlas_t

    type, public :: poincare_validated_cut_atlas_t
        private
        type(poincare_cut_atlas_t) :: atlas
        logical :: issued = .false.
    end type poincare_validated_cut_atlas_t

    type, public :: poincare_phase_space_context_t
        !! Invariant context required by the Buchholz Eq. 14 measure.
        real(dp) :: H0 = 0.0_dp
        real(dp) :: JK = 0.0_dp
        real(dp) :: sigma = 0.0_dp
        real(dp) :: psi_star = 0.0_dp
        real(dp) :: dpsi_star_dRc = 0.0_dp
        integer :: field_input_id = 0
        integer :: equilibrium_input_id = 0
        integer :: profile_input_id = 0
    end type poincare_phase_space_context_t

    type, public :: poincare_state_evidence_t
        !! State data are evidence for a provider, not a proof flag.
        real(dp) :: time = 0.0_dp
        real(dp) :: state(5) = 0.0_dp
        type(poincare_phase_space_context_t) :: context
    end type poincare_state_evidence_t

    type, public :: poincare_allowed_interval_t
        !! One connected allowed interval for one launch sigma and one
        !! psi-star parameterization copy.
        integer :: branch_id = 0
        integer :: sigma = 0
        integer :: local_component_id = 0
        integer :: psi_star_parameterization_copy_id = 0
        integer :: lower_boundary_kind = 0
        integer :: upper_boundary_kind = 0
        integer :: lower_homoclinic_pair_id = 0
        integer :: upper_homoclinic_pair_id = 0
        real(dp) :: lower = 0.0_dp
        real(dp) :: upper = 0.0_dp
    end type poincare_allowed_interval_t

    type, public :: poincare_homoclinic_pair_t
        !! Exactly one X and one usual crossing of one separatrix.
        integer :: pair_id = 0
        integer :: separatrix_identity = 0
        integer :: x_endpoint_count = 0
        integer :: usual_endpoint_count = 0
        real(dp) :: x_coefficient = 0.0_dp
        real(dp) :: usual_coefficient = 0.0_dp
    end type poincare_homoclinic_pair_t

    type, public :: poincare_global_component_t
        integer :: global_component_id = 0
        integer :: component_identity_seed = 0
        integer :: branch_id = 0
        integer :: sigma = 0
        integer :: local_component_id = 0
        integer :: psi_star_parameterization_copy_id = 0
        integer :: lower_boundary_kind = 0
        integer :: upper_boundary_kind = 0
        integer :: lower_homoclinic_pair_id = 0
        integer :: upper_homoclinic_pair_id = 0
        real(dp) :: lower = 0.0_dp
        real(dp) :: upper = 0.0_dp
    end type poincare_global_component_t

    type, public :: poincare_crossing_evidence_t
        integer :: branch_id = 0
        real(dp) :: parameter = 0.0_dp
        real(dp) :: position(2) = 0.0_dp
        type(poincare_phase_space_context_t) :: context
        real(dp) :: time = 0.0_dp
    end type poincare_crossing_evidence_t

    type, public :: poincare_orbit_diagnostics_t
        !! Diagnostics only; never used to define a numerical class.
        logical :: trapped_passing_label = .false.
        logical :: parallel_velocity_sign_change = .false.
        logical :: axis_encirclement = .false.
        logical :: winding_available = .false.
        real(dp) :: winding_value = 0.0_dp
    end type poincare_orbit_diagnostics_t

    type, public :: poincare_return_cycle_evidence_t
        integer :: launch_branch_id = 0
        real(dp) :: launch_parameter = 0.0_dp
        real(dp) :: launch_position(2) = 0.0_dp
        type(poincare_state_evidence_t) :: launch_state
        type(poincare_state_evidence_t) :: return_state
        type(poincare_crossing_evidence_t), allocatable :: crossings(:)
        real(dp) :: state_tolerance = 0.0_dp
        integer :: wall_outcome = PCA_WALL_UNRESOLVED
        integer :: domain_outcome = PCA_DOMAIN_UNRESOLVED
        type(poincare_orbit_diagnostics_t) :: diagnostics
    end type poincare_return_cycle_evidence_t

    type, public :: poincare_exactly_two_certificate_t
        !! Theorem evidence supplied independently of numerical integration.
        !!
        !! crossing_count is the multiplicity proved by the provider; it is
        !! deliberately not inferred from evidence%crossings.  certificate_id
        !! identifies the theorem/provider instance and must be positive.
        integer :: certificate_id = 0
        integer :: crossing_count = 0
        logical :: exactly_two_proved = .false.
    end type poincare_exactly_two_certificate_t

    abstract interface
        subroutine poincare_cut_branch_evaluator(branch, parameter, position, &
                parameterization_value, parameterization_derivative, status)
            import :: dp, poincare_cut_branch_t
            type(poincare_cut_branch_t), intent(in) :: branch
            real(dp), intent(in) :: parameter
            real(dp), intent(out) :: position(2)
            real(dp), intent(out) :: parameterization_value
            real(dp), intent(out) :: parameterization_derivative
            integer, intent(out) :: status
        end subroutine poincare_cut_branch_evaluator

        subroutine poincare_fortsym_eq13_cut_kernel(branch, parameter, position, &
                residual, residual_derivative, status)
            import :: dp, poincare_cut_branch_t
            type(poincare_cut_branch_t), intent(in) :: branch
            real(dp), intent(in) :: parameter
            real(dp), intent(in) :: position(2)
            real(dp), intent(out) :: residual
            real(dp), intent(out) :: residual_derivative
            integer, intent(out) :: status
        end subroutine poincare_fortsym_eq13_cut_kernel

        subroutine poincare_branch_disjointness(branch_a, branch_b, disjoint, status)
            import :: poincare_cut_branch_t
            type(poincare_cut_branch_t), intent(in) :: branch_a
            type(poincare_cut_branch_t), intent(in) :: branch_b
            logical, intent(out) :: disjoint
            integer, intent(out) :: status
        end subroutine poincare_branch_disjointness

        subroutine poincare_fortsym_branch_certificate(branch, status)
            import :: poincare_cut_branch_t
            type(poincare_cut_branch_t), intent(in) :: branch
            integer, intent(out) :: status
        end subroutine poincare_fortsym_branch_certificate

        subroutine poincare_fortsym_allowed_interval_certificate(interval, status)
            import :: poincare_allowed_interval_t
            type(poincare_allowed_interval_t), intent(in) :: interval
            integer, intent(out) :: status
        end subroutine poincare_fortsym_allowed_interval_certificate

        subroutine poincare_fortsym_homoclinic_pair_certificate(pair, status)
            import :: poincare_homoclinic_pair_t
            type(poincare_homoclinic_pair_t), intent(in) :: pair
            integer, intent(out) :: status
        end subroutine poincare_fortsym_homoclinic_pair_certificate

        subroutine poincare_fortsym_cdot_enclosure_verifier(branch, parameter, &
                position, context, lower, upper, status)
            import :: dp, poincare_cut_branch_t, poincare_phase_space_context_t
            type(poincare_cut_branch_t), intent(in) :: branch
            real(dp), intent(in) :: parameter
            real(dp), intent(in) :: position(2)
            type(poincare_phase_space_context_t), intent(in) :: context
            real(dp), intent(out) :: lower
            real(dp), intent(out) :: upper
            integer, intent(out) :: status
        end subroutine poincare_fortsym_cdot_enclosure_verifier

        subroutine poincare_fortsym_phase_space_density_kernel(branch, parameter, &
                position, context, signed_jacobian, positive_density, status)
            import :: dp, poincare_cut_branch_t, poincare_phase_space_context_t
            type(poincare_cut_branch_t), intent(in) :: branch
            real(dp), intent(in) :: parameter
            real(dp), intent(in) :: position(2)
            type(poincare_phase_space_context_t), intent(in) :: context
            real(dp), intent(out) :: signed_jacobian
            real(dp), intent(out) :: positive_density
            integer, intent(out) :: status
        end subroutine poincare_fortsym_phase_space_density_kernel

        subroutine poincare_fortsym_context_binding_certificate(launch_context, &
                candidate_context, status)
            import :: poincare_phase_space_context_t
            type(poincare_phase_space_context_t), intent(in) :: launch_context
            type(poincare_phase_space_context_t), intent(in) :: candidate_context
            integer, intent(out) :: status
        end subroutine poincare_fortsym_context_binding_certificate

        subroutine poincare_fortsym_exactly_two_provider(atlas, evidence, &
                certificate, status)
            import :: poincare_validated_cut_atlas_t
            import :: poincare_return_cycle_evidence_t
            import :: poincare_exactly_two_certificate_t
            type(poincare_validated_cut_atlas_t), intent(in) :: atlas
            type(poincare_return_cycle_evidence_t), intent(in) :: evidence
            type(poincare_exactly_two_certificate_t), intent(out) :: certificate
            integer, intent(out) :: status
        end subroutine poincare_fortsym_exactly_two_provider

        subroutine poincare_fortsym_return_state_certificate(launch_state, &
                return_state, tolerance, status)
            import :: dp, poincare_state_evidence_t
            type(poincare_state_evidence_t), intent(in) :: launch_state
            type(poincare_state_evidence_t), intent(in) :: return_state
            real(dp), intent(in) :: tolerance
            integer, intent(out) :: status
        end subroutine poincare_fortsym_return_state_certificate
    end interface

    public :: poincare_cut_branch_evaluator
    public :: poincare_fortsym_eq13_cut_kernel
    public :: poincare_branch_disjointness
    public :: poincare_fortsym_branch_certificate
    public :: poincare_fortsym_allowed_interval_certificate
    public :: poincare_fortsym_homoclinic_pair_certificate
    public :: poincare_fortsym_cdot_enclosure_verifier
    public :: poincare_fortsym_phase_space_density_kernel
    public :: poincare_fortsym_context_binding_certificate
    public :: poincare_fortsym_exactly_two_provider
    public :: poincare_fortsym_return_state_certificate
    public :: make_symmetric_midplane_atlas
    public :: validate_poincare_cut_atlas
    public :: combine_allowed_intervals
    public :: validate_return_cycle_evidence

contains

    subroutine make_symmetric_midplane_atlas(branch, atlas, status, message)
        type(poincare_cut_branch_t), intent(in) :: branch
        type(poincare_cut_atlas_t), intent(out) :: atlas
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        allocate(atlas%branches(1))
        atlas%branches(1) = branch
        atlas%expected_branch_count = 1
        atlas%symmetric_midplane_specialization = .true.
        status = PCA_INCOMPLETE_ATLAS
        message = 'atlas constructed but Fortsym completeness provider is absent'
    end subroutine make_symmetric_midplane_atlas

    subroutine validate_poincare_cut_atlas(atlas, evaluate, cut_kernel, disjointness, &
            branch_certificate, validated_atlas, status, message)
        type(poincare_cut_atlas_t), intent(in) :: atlas
        procedure(poincare_cut_branch_evaluator) :: evaluate
        procedure(poincare_fortsym_eq13_cut_kernel) :: cut_kernel
        procedure(poincare_branch_disjointness) :: disjointness
        procedure(poincare_fortsym_branch_certificate) :: branch_certificate
        type(poincare_validated_cut_atlas_t), intent(out) :: validated_atlas
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: i, j, k, local_status, branch_count
        real(dp) :: parameter, previous_value, value, derivative
        real(dp) :: position(2), residual, residual_derivative
        logical :: disjoint, have_previous

        validated_atlas%issued = .false.
        status = PCA_SUCCESS
        message = 'complete Poincare cut atlas'
        if (.not. allocated(atlas%branches)) then
            status = PCA_INCOMPLETE_ATLAS
            message = 'atlas has no branches'
            return
        end if
        branch_count = size(atlas%branches)
        if (branch_count < 1) then
            status = PCA_INCOMPLETE_ATLAS
            message = 'atlas has no branches'
            return
        end if
        if (atlas%expected_branch_count /= branch_count) then
            status = PCA_INVALID_CONTRACT
            message = 'branch count disagrees with atlas contract'
            return
        end if
        if (atlas%symmetric_midplane_specialization .and. branch_count /= 1) then
            status = PCA_INVALID_CONTRACT
            message = 'symmetric midplane atlas must have one branch'
            return
        end if

        !! Establish the identity domain and uniqueness before validating any
        !! branch content.  The array position is the canonical branch
        !! identity; neither a caller label nor a provider callback may
        !! substitute for it.
        do i = 1, branch_count
            if (atlas%branches(i)%branch_id < 1 .or. &
                    atlas%branches(i)%branch_id > branch_count) then
                status = PCA_INVALID_CONTRACT
                message = 'branch identity is outside the atlas domain'
                return
            end if
        end do
        do i = 1, branch_count
            do j = 1, i - 1
                if (atlas%branches(i)%branch_id == atlas%branches(j)%branch_id) then
                    status = PCA_INVALID_CONTRACT
                    message = 'branch identity is duplicated'
                    return
                end if
            end do
        end do
        !! Complete all cross-branch structural checks before invoking a
        !! provider for any one branch.  A fail-closed callback on branch one
        !! must not mask an invalid later branch.
        do i = 1, branch_count
            if (atlas%branches(i)%branch_id /= i) then
                status = PCA_INVALID_CONTRACT
                message = 'branch identity is not the canonical array position'
                return
            end if
            if (atlas%branches(i)%component_identity_seed <= 0) then
                status = PCA_INVALID_CONTRACT
                message = 'branch has no positive component identity seed'
                return
            end if
            if (atlas%branches(i)%field_input_id <= 0 .or. &
                    atlas%branches(i)%equilibrium_input_id <= 0) then
                status = PCA_INVALID_CONTRACT
                message = 'branch has no positive field/equilibrium identity'
                return
            end if
            if (i > 1) then
                if (atlas%branches(i)%field_input_id /= &
                        atlas%branches(1)%field_input_id .or. &
                        atlas%branches(i)%equilibrium_input_id /= &
                        atlas%branches(1)%equilibrium_input_id) then
                    status = PCA_INVALID_CONTRACT
                    message = 'atlas branches mix field/equilibrium identities'
                    return
                end if
            end if
            do j = 1, i - 1
                if (atlas%branches(i)%component_identity_seed == &
                        atlas%branches(j)%component_identity_seed) then
                    status = PCA_INVALID_CONTRACT
                    message = 'component identity seed is duplicated'
                    return
                end if
            end do
            if (abs(atlas%branches(i)%orientation) /= 1) then
                status = PCA_INVALID_CONTRACT
                message = 'branch orientation must be plus or minus one'
                return
            end if
            if (.not. ieee_is_finite(atlas%branches(i)%parameter_lower) .or. &
                    .not. ieee_is_finite(atlas%branches(i)%parameter_upper)) then
                status = PCA_INVALID_CONTRACT
                message = 'branch parameter interval is not finite'
                return
            end if
            if (atlas%branches(i)%parameter_upper <= &
                    atlas%branches(i)%parameter_lower) then
                status = PCA_INVALID_CONTRACT
                message = 'branch parameter interval is not ordered'
                return
            end if
            if (atlas%branches(i)%cut_residual_tolerance < 0.0_dp .or. &
                    .not. ieee_is_finite(atlas%branches(i)%cut_residual_tolerance)) then
                status = PCA_INVALID_CONTRACT
                message = 'cut residual tolerance is invalid'
                return
            end if
            if (.not. atlas%branches(i)%endpoint_start_known .or. &
                    .not. atlas%branches(i)%endpoint_end_known) then
                status = PCA_INCOMPLETE_ATLAS
                message = 'branch endpoints are not supplied'
                return
            end if
            if (any(.not. ieee_is_finite(atlas%branches(i)%endpoint_start)) .or. &
                    any(.not. ieee_is_finite(atlas%branches(i)%endpoint_end))) then
                status = PCA_INVALID_CONTRACT
                message = 'branch endpoint is not finite'
                return
            end if
            if (atlas%symmetric_midplane_specialization) then
                if (.not. atlas%branches(i)%physical_r_ordered .or. &
                        .not. atlas%branches(i)%continuous_inboard_axis_outboard .or. &
                        atlas%branches(i)%outboard_only) then
                    status = PCA_INVALID_CONTRACT
                    message = 'symmetric specialization is not a complete physical branch'
                    return
                end if
                if (trim(atlas%branches(i)%parameter_units) == 's_tor') then
                    status = PCA_INVALID_CONTRACT
                    message = 'symmetric specialization cannot use s_tor'
                    return
                end if
            end if

            !! This is the non-negotiable global branch proof.  A missing
            !! Fortsym provider must not be replaced by flags or samples.
            call branch_certificate(atlas%branches(i), local_status)
            if (local_status == PCA_PROVIDER_UNAVAILABLE) then
                status = PCA_PROVIDER_UNAVAILABLE
                message = 'Fortsym branch completeness provider is unavailable'
                return
            end if
            if (local_status /= PCA_SUCCESS) then
                status = PCA_CALLBACK_FAILURE
                message = 'Fortsym branch completeness certificate failed'
                return
            end if

            have_previous = .false.
            do k = 0, 4
                parameter = atlas%branches(i)%parameter_lower + real(k, dp) * &
                    (atlas%branches(i)%parameter_upper - &
                    atlas%branches(i)%parameter_lower) / 4.0_dp
                call evaluate(atlas%branches(i), parameter, position, value, &
                    derivative, local_status)
                if (local_status /= PCA_SUCCESS) then
                    status = PCA_CALLBACK_FAILURE
                    message = 'branch evaluator failed'
                    return
                end if
                if (any(.not. ieee_is_finite(position)) .or. &
                        .not. ieee_is_finite(value) .or. &
                        .not. ieee_is_finite(derivative)) then
                    status = PCA_INVALID_CONTRACT
                    message = 'branch sample is not finite'
                    return
                end if
                if (atlas%branches(i)%orientation == 1 .and. derivative <= 0.0_dp) then
                    status = PCA_INVALID_CONTRACT
                    message = 'branch sample violates orientation'
                    return
                end if
                if (atlas%branches(i)%orientation == -1 .and. derivative >= 0.0_dp) then
                    status = PCA_INVALID_CONTRACT
                    message = 'branch sample violates orientation'
                    return
                end if
                if (have_previous) then
                    if (atlas%branches(i)%orientation == 1 .and. value <= previous_value) then
                        status = PCA_INVALID_CONTRACT
                        message = 'branch parameterization reverses direction'
                        return
                    end if
                    if (atlas%branches(i)%orientation == -1 .and. value >= previous_value) then
                        status = PCA_INVALID_CONTRACT
                        message = 'branch parameterization reverses direction'
                        return
                    end if
                end if
                previous_value = value
                have_previous = .true.
                call cut_kernel(atlas%branches(i), parameter, position, residual, &
                    residual_derivative, local_status)
                if (local_status == PCA_PROVIDER_UNAVAILABLE) then
                    status = PCA_PROVIDER_UNAVAILABLE
                    message = 'Fortsym Eq.13 cut provider is unavailable'
                    return
                end if
                if (local_status /= PCA_SUCCESS) then
                    status = PCA_CALLBACK_FAILURE
                    message = 'Eq.13 cut kernel failed'
                    return
                end if
                if (.not. ieee_is_finite(residual) .or. &
                        .not. ieee_is_finite(residual_derivative)) then
                    status = PCA_INVALID_CONTRACT
                    message = 'Eq.13 cut sample is not finite'
                    return
                end if
                if (abs(residual) > atlas%branches(i)%cut_residual_tolerance) then
                    status = PCA_INVALID_CONTRACT
                    message = 'branch does not satisfy Eq.13 cut contract'
                    return
                end if
            end do
        end do

        do i = 1, branch_count
            do j = 1, i - 1
                call disjointness(atlas%branches(i), atlas%branches(j), disjoint, &
                    local_status)
                if (local_status == PCA_PROVIDER_UNAVAILABLE) then
                    status = PCA_PROVIDER_UNAVAILABLE
                    message = 'Fortsym branch disjointness provider is unavailable'
                    return
                end if
                if (local_status /= PCA_SUCCESS) then
                    status = PCA_CALLBACK_FAILURE
                    message = 'branch disjointness callback failed'
                    return
                end if
                if (.not. disjoint) then
                    status = PCA_INVALID_CONTRACT
                    message = 'cut branches overlap'
                    return
                end if
            end do
        end do

        validated_atlas%atlas = atlas
        validated_atlas%issued = .true.
    end subroutine validate_poincare_cut_atlas

    subroutine combine_allowed_intervals(validated_atlas, local_intervals, pairs, &
            interval_certificate, pair_certificate, components, status, message)
        type(poincare_validated_cut_atlas_t), intent(in) :: validated_atlas
        type(poincare_allowed_interval_t), intent(in) :: local_intervals(:)
        type(poincare_homoclinic_pair_t), intent(in) :: pairs(:)
        procedure(poincare_fortsym_allowed_interval_certificate) :: &
            interval_certificate
        procedure(poincare_fortsym_homoclinic_pair_certificate) :: pair_certificate
        type(poincare_global_component_t), allocatable, intent(out) :: components(:)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        type(poincare_allowed_interval_t), allocatable :: ordered(:)
        integer :: i, j, branch_index, local_status
        logical :: found

        status = PCA_SUCCESS
        message = 'globally unique component identities'
        if (.not. validated_atlas%issued .or. &
                .not. allocated(validated_atlas%atlas%branches)) then
            status = PCA_INCOMPLETE_ATLAS
            message = 'intervals require a validated atlas token'
            return
        end if
        if (size(local_intervals) < 1) then
            status = PCA_INCOMPLETE_ATLAS
            message = 'no allowed intervals were supplied'
            return
        end if
        allocate(ordered(size(local_intervals)))
        ordered = local_intervals
        call canonical_sort_intervals(ordered)

        do i = 1, size(ordered)
            found = .false.
            branch_index = 0
            do j = 1, size(validated_atlas%atlas%branches)
                if (validated_atlas%atlas%branches(j)%branch_id == ordered(i)%branch_id) then
                    found = .true.
                    branch_index = j
                    exit
                end if
            end do
            if (.not. found) then
                status = PCA_INVALID_CONTRACT
                message = 'interval references an unknown branch'
                return
            end if
            if (abs(ordered(i)%sigma) /= 1 .or. ordered(i)%local_component_id <= 0 .or. &
                    ordered(i)%psi_star_parameterization_copy_id <= 0) then
                status = PCA_INVALID_CONTRACT
                message = 'interval identity is invalid'
                return
            end if
            if (.not. valid_boundary_kind(ordered(i)%lower_boundary_kind) .or. &
                    .not. valid_boundary_kind(ordered(i)%upper_boundary_kind)) then
                status = PCA_INVALID_CONTRACT
                message = 'interval boundary kind is unknown'
                return
            end if
            if (.not. ieee_is_finite(ordered(i)%lower) .or. &
                    .not. ieee_is_finite(ordered(i)%upper) .or. &
                    ordered(i)%upper <= ordered(i)%lower) then
                status = PCA_INVALID_CONTRACT
                message = 'interval endpoints are invalid'
                return
            end if
            if (.not. valid_boundary_pair_reference(ordered(i)%lower_boundary_kind, &
                    ordered(i)%lower_homoclinic_pair_id) .or. &
                    .not. valid_boundary_pair_reference(ordered(i)%upper_boundary_kind, &
                    ordered(i)%upper_homoclinic_pair_id)) then
                status = PCA_INVALID_CONTRACT
                message = 'homoclinic boundary reference is invalid'
                return
            end if
            call interval_certificate(ordered(i), local_status)
            if (local_status == PCA_PROVIDER_UNAVAILABLE) then
                status = PCA_PROVIDER_UNAVAILABLE
                message = 'Fortsym allowed-interval provider is unavailable'
                return
            end if
            if (local_status /= PCA_SUCCESS) then
                status = PCA_CALLBACK_FAILURE
                message = 'Fortsym allowed-interval certificate failed'
                return
            end if
            do j = 1, i - 1
                if (same_interval_identity(ordered(j), ordered(i))) then
                    status = PCA_INVALID_CONTRACT
                    message = 'local component identity is duplicated'
                    return
                end if
                if (same_parameterization(ordered(j), ordered(i))) then
                    if (max(ordered(j)%lower, ordered(i)%lower) < &
                            min(ordered(j)%upper, ordered(i)%upper)) then
                        status = PCA_INVALID_CONTRACT
                        message = 'allowed intervals overlap'
                        return
                    end if
                end if
            end do
            if (branch_index < 1) then
                status = PCA_INVALID_CONTRACT
                message = 'internal branch lookup failure'
                return
            end if
        end do

        call validate_homoclinic_pairs(ordered, pairs, pair_certificate, status, message)
        if (status /= PCA_SUCCESS) return
        allocate(components(size(ordered)))
        do i = 1, size(ordered)
            components(i)%global_component_id = i
            do j = 1, size(validated_atlas%atlas%branches)
                if (validated_atlas%atlas%branches(j)%branch_id == ordered(i)%branch_id) then
                    components(i)%component_identity_seed = &
                        validated_atlas%atlas%branches(j)%component_identity_seed
                    exit
                end if
            end do
            components(i)%branch_id = ordered(i)%branch_id
            components(i)%sigma = ordered(i)%sigma
            components(i)%local_component_id = ordered(i)%local_component_id
            components(i)%psi_star_parameterization_copy_id = &
                ordered(i)%psi_star_parameterization_copy_id
            components(i)%lower_boundary_kind = ordered(i)%lower_boundary_kind
            components(i)%upper_boundary_kind = ordered(i)%upper_boundary_kind
            components(i)%lower_homoclinic_pair_id = ordered(i)%lower_homoclinic_pair_id
            components(i)%upper_homoclinic_pair_id = ordered(i)%upper_homoclinic_pair_id
            components(i)%lower = ordered(i)%lower
            components(i)%upper = ordered(i)%upper
        end do
    end subroutine combine_allowed_intervals

    subroutine validate_return_cycle_evidence(validated_atlas, evidence, &
            cdot_enclosure, density_kernel, context_binding_certificate, &
            exactly_two_provider, return_state_certificate, status, message)
        type(poincare_validated_cut_atlas_t), intent(in) :: validated_atlas
        type(poincare_return_cycle_evidence_t), intent(in) :: evidence
        procedure(poincare_fortsym_cdot_enclosure_verifier) :: cdot_enclosure
        procedure(poincare_fortsym_phase_space_density_kernel) :: density_kernel
        procedure(poincare_fortsym_context_binding_certificate) :: &
            context_binding_certificate
        procedure(poincare_fortsym_exactly_two_provider) :: exactly_two_provider
        procedure(poincare_fortsym_return_state_certificate) :: &
            return_state_certificate
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        integer :: i, branch_index, local_status
        integer :: launch_sign, first_sign, second_sign, crossing_sign
        real(dp) :: signed_jacobian, positive_density
        logical :: found
        type(poincare_exactly_two_certificate_t) :: exactly_two_certificate

        status = PCA_SUCCESS
        message = 'certified two-intersection full return'
        if (.not. validated_atlas%issued .or. &
                .not. allocated(validated_atlas%atlas%branches)) then
            status = PCA_INCOMPLETE_ATLAS
            message = 'return evidence requires a validated complete atlas token'
            return
        end if
        if (.not. allocated(evidence%crossings)) then
            status = PCA_UNKNOWN_MULTIPLICITY
            message = 'return evidence has no crossing sequence'
            return
        end if
        if (size(evidence%crossings) /= 2) then
            status = PCA_UNKNOWN_MULTIPLICITY
            message = 'return evidence does not contain exactly two crossings'
            return
        end if
        if (.not. valid_context(evidence%launch_state%context) .or. &
                .not. valid_context(evidence%return_state%context)) then
            status = PCA_INVALID_STATE
            message = 'launch/return invariant context is invalid'
            return
        end if
        do i = 1, 2
            if (.not. valid_context(evidence%crossings(i)%context)) then
                status = PCA_INVALID_STATE
                message = 'crossing invariant context is invalid'
                return
            end if
        end do
        if (.not. same_context_ids(evidence%launch_state%context, &
                evidence%return_state%context)) then
            status = PCA_INVALID_CONTRACT
            message = 'return context input identities differ from launch'
            return
        end if
        do i = 1, 2
            if (.not. same_context_ids(evidence%launch_state%context, &
                    evidence%crossings(i)%context)) then
                status = PCA_INVALID_CONTRACT
                message = 'crossing context input identities differ from launch'
                return
            end if
        end do
        call validate_context_branch_bindings(validated_atlas, evidence, status, message)
        if (status /= PCA_SUCCESS) return
        if (.not. ieee_is_finite(evidence%state_tolerance) .or. &
                evidence%state_tolerance <= 0.0_dp) then
            status = PCA_INVALID_STATE
            message = 'return state tolerance is invalid'
            return
        end if
        if (evidence%wall_outcome /= PCA_WALL_NONE .or. &
                evidence%domain_outcome /= PCA_DOMAIN_NONE) then
            status = PCA_INVALID_CONTRACT
            message = 'return cycle has a wall or domain exit'
            return
        end if
        if (.not. ieee_is_finite(evidence%launch_state%time) .or. &
                .not. ieee_is_finite(evidence%return_state%time) .or. &
                .not. ieee_is_finite(evidence%crossings(1)%time) .or. &
                .not. ieee_is_finite(evidence%crossings(2)%time) .or. &
                evidence%launch_state%time >= evidence%crossings(1)%time .or. &
                evidence%crossings(1)%time >= evidence%crossings(2)%time .or. &
                evidence%crossings(2)%time > evidence%return_state%time) then
            status = PCA_INVALID_STATE
            message = 'crossing events are not ordered after launch'
            return
        end if
        !! Input identities are exact above.  H0, JK, and psi_star are
        !! floating invariants and are bound here by generated, unit-aware
        !! tolerances.  Sigma and dpsi_star_dRc are point-local and may differ.
        call context_binding_certificate(evidence%launch_state%context, &
            evidence%return_state%context, local_status)
        if (local_status == PCA_PROVIDER_UNAVAILABLE) then
            status = PCA_PROVIDER_UNAVAILABLE
            message = 'Fortsym invariant-context binding provider is unavailable'
            return
        end if
        if (local_status /= PCA_SUCCESS) then
            status = PCA_INVALID_STATE
            message = 'return invariants do not match launch within generated tolerance'
            return
        end if
        do i = 1, 2
            call context_binding_certificate(evidence%launch_state%context, &
                evidence%crossings(i)%context, local_status)
            if (local_status == PCA_PROVIDER_UNAVAILABLE) then
                status = PCA_PROVIDER_UNAVAILABLE
                message = 'Fortsym invariant-context binding provider is unavailable'
                return
            end if
            if (local_status /= PCA_SUCCESS) then
                status = PCA_INVALID_STATE
                message = 'crossing invariants differ from launch beyond tolerance'
                return
            end if
        end do
        call return_state_certificate(evidence%launch_state, evidence%return_state, &
            evidence%state_tolerance, local_status)
        if (local_status == PCA_PROVIDER_UNAVAILABLE) then
            status = PCA_PROVIDER_UNAVAILABLE
            message = 'Fortsym return-state certificate is unavailable'
            return
        end if
        if (local_status /= PCA_SUCCESS) then
            status = PCA_INVALID_STATE
            message = 'launch/return state certificate failed'
            return
        end if
        !! The observed crossing array is evidence only.  Its length cannot
        !! prove that no additional crossings were missed.  Only the separate
        !! theorem/provider may issue the typed multiplicity certificate.
        exactly_two_certificate = poincare_exactly_two_certificate_t()
        call exactly_two_provider(validated_atlas, evidence, exactly_two_certificate, &
            local_status)
        if (local_status == PCA_PROVIDER_UNAVAILABLE) then
            status = PCA_PROVIDER_UNAVAILABLE
            message = 'Fortsym exactly-two theorem provider is unavailable'
            return
        end if
        if (local_status /= PCA_SUCCESS) then
            status = PCA_INVALID_CERTIFICATE
            message = 'Fortsym exactly-two theorem provider failed'
            return
        end if
        if (exactly_two_certificate%certificate_id <= 0 .or. &
                exactly_two_certificate%crossing_count /= 2 .or. &
                .not. exactly_two_certificate%exactly_two_proved) then
            status = PCA_INVALID_CERTIFICATE
            message = 'exactly-two provider returned no valid theorem certificate'
            return
        end if

        found = .false.
        do branch_index = 1, size(validated_atlas%atlas%branches)
            if (validated_atlas%atlas%branches(branch_index)%branch_id == &
                    evidence%launch_branch_id) then
                found = .true.
                exit
            end if
        end do
        if (.not. found) then
            status = PCA_INVALID_CONTRACT
            message = 'launch branch is not in the atlas'
            return
        end if
        if (.not. ieee_is_finite(evidence%launch_parameter) .or. &
                any(.not. ieee_is_finite(evidence%launch_position)) .or. &
                evidence%launch_parameter < &
                validated_atlas%atlas%branches(branch_index)%parameter_lower .or. &
                evidence%launch_parameter > &
                validated_atlas%atlas%branches(branch_index)%parameter_upper) then
            status = PCA_INVALID_CONTRACT
            message = 'launch point lies outside its certified branch interval'
            return
        end if
        call certified_cdot_sign(validated_atlas%atlas%branches(branch_index), &
            evidence%launch_parameter, evidence%launch_position, &
            evidence%launch_state%context, cdot_enclosure, launch_sign, status, message)
        if (status /= PCA_SUCCESS) return

        do i = 1, 2
            found = .false.
            do branch_index = 1, size(validated_atlas%atlas%branches)
                if (validated_atlas%atlas%branches(branch_index)%branch_id == &
                        evidence%crossings(i)%branch_id) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                status = PCA_INVALID_CONTRACT
                message = 'crossing references an unknown branch'
                return
            end if
            if (.not. ieee_is_finite(evidence%crossings(i)%parameter) .or. &
                    any(.not. ieee_is_finite(evidence%crossings(i)%position)) .or. &
                    evidence%crossings(i)%parameter < &
                    validated_atlas%atlas%branches(branch_index)%parameter_lower .or. &
                    evidence%crossings(i)%parameter > &
                    validated_atlas%atlas%branches(branch_index)%parameter_upper) then
                status = PCA_INVALID_CONTRACT
                message = 'crossing lies outside its certified branch interval'
                return
            end if
            call certified_cdot_sign(validated_atlas%atlas%branches(branch_index), &
                evidence%crossings(i)%parameter, evidence%crossings(i)%position, &
                evidence%crossings(i)%context, cdot_enclosure, crossing_sign, status, &
                message)
            if (status /= PCA_SUCCESS) return
            call density_kernel(validated_atlas%atlas%branches(branch_index), &
                evidence%crossings(i)%parameter, evidence%crossings(i)%position, &
                evidence%crossings(i)%context, signed_jacobian, positive_density, &
                local_status)
            if (local_status == PCA_PROVIDER_UNAVAILABLE) then
                status = PCA_PROVIDER_UNAVAILABLE
                message = 'Fortsym phase-space density provider is unavailable'
                return
            end if
            if (local_status /= PCA_SUCCESS) then
                status = PCA_CALLBACK_FAILURE
                message = 'Fortsym phase-space density provider failed'
                return
            end if
            if (.not. ieee_is_finite(signed_jacobian) .or. &
                    .not. ieee_is_finite(positive_density) .or. &
                    signed_jacobian == 0.0_dp .or. positive_density <= 0.0_dp) then
                status = PCA_INVALID_CONTRACT
                message = 'generated phase-space measure is invalid'
                return
            end if
            if (i == 1) first_sign = crossing_sign
            if (i == 2) second_sign = crossing_sign
        end do
        if (first_sign == launch_sign .or. second_sign /= launch_sign) then
            status = PCA_INVALID_CONTRACT
            message = 'crossing orientations do not form opposite then same return'
            return
        end if
    end subroutine validate_return_cycle_evidence

    subroutine validate_context_branch_bindings(validated_atlas, evidence, status, message)
        type(poincare_validated_cut_atlas_t), intent(in) :: validated_atlas
        type(poincare_return_cycle_evidence_t), intent(in) :: evidence
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer :: i, branch_index
        logical :: found

        status = PCA_SUCCESS
        message = 'context identities match atlas branches'
        found = .false.
        do branch_index = 1, size(validated_atlas%atlas%branches)
            if (validated_atlas%atlas%branches(branch_index)%branch_id == &
                    evidence%launch_branch_id) then
                found = .true.
                exit
            end if
        end do
        if (.not. found) then
            status = PCA_INVALID_CONTRACT
            message = 'launch branch is not in the atlas'
            return
        end if
        if (evidence%launch_state%context%field_input_id /= &
                validated_atlas%atlas%branches(branch_index)%field_input_id .or. &
                evidence%launch_state%context%equilibrium_input_id /= &
                validated_atlas%atlas%branches(branch_index)%equilibrium_input_id) then
            status = PCA_INVALID_CONTRACT
            message = 'launch context does not match branch field/equilibrium identity'
            return
        end if
        do i = 1, 2
            found = .false.
            do branch_index = 1, size(validated_atlas%atlas%branches)
                if (validated_atlas%atlas%branches(branch_index)%branch_id == &
                        evidence%crossings(i)%branch_id) then
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                status = PCA_INVALID_CONTRACT
                message = 'crossing branch is not in the atlas'
                return
            end if
            if (evidence%crossings(i)%context%field_input_id /= &
                    validated_atlas%atlas%branches(branch_index)%field_input_id .or. &
                    evidence%crossings(i)%context%equilibrium_input_id /= &
                    validated_atlas%atlas%branches(branch_index)%equilibrium_input_id) then
                status = PCA_INVALID_CONTRACT
                message = 'crossing context does not match branch field identity'
                return
            end if
        end do
    end subroutine validate_context_branch_bindings

    subroutine certified_cdot_sign(branch, parameter, position, context, verifier, &
            sign, status, message)
        type(poincare_cut_branch_t), intent(in) :: branch
        real(dp), intent(in) :: parameter, position(2)
        type(poincare_phase_space_context_t), intent(in) :: context
        procedure(poincare_fortsym_cdot_enclosure_verifier) :: verifier
        integer, intent(out) :: sign, status
        character(len=*), intent(out) :: message
        real(dp) :: lower, upper
        call verifier(branch, parameter, position, context, lower, upper, status)
        if (status == PCA_PROVIDER_UNAVAILABLE) then
            message = 'Fortsym Cdot enclosure provider is unavailable'
            return
        end if
        if (status /= PCA_SUCCESS) then
            message = 'Fortsym Cdot enclosure provider failed'
            return
        end if
        if (.not. ieee_is_finite(lower) .or. .not. ieee_is_finite(upper) .or. &
                upper <= lower) then
            status = PCA_UNRESOLVED_TANGENCY
            message = 'Cdot enclosure is invalid'
            return
        end if
        if (lower <= 0.0_dp .and. upper >= 0.0_dp) then
            status = PCA_UNRESOLVED_TANGENCY
            message = 'Cdot enclosure does not exclude zero'
            return
        end if
        sign = 1
        if (upper < 0.0_dp) sign = -1
        status = PCA_SUCCESS
        message = 'simple Cdot crossing certified'
    end subroutine certified_cdot_sign

    subroutine validate_homoclinic_pairs(intervals, pairs, pair_certificate, status, message)
        type(poincare_allowed_interval_t), intent(in) :: intervals(:)
        type(poincare_homoclinic_pair_t), intent(in) :: pairs(:)
        procedure(poincare_fortsym_homoclinic_pair_certificate) :: pair_certificate
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer :: i, j, x_count, usual_count, local_status
        logical :: referenced

        status = PCA_SUCCESS
        message = 'homoclinic pair identities are complete'
        do i = 1, size(pairs)
            if (pairs(i)%pair_id <= 0 .or. pairs(i)%separatrix_identity <= 0 .or. &
                    pairs(i)%x_endpoint_count /= 1 .or. &
                    pairs(i)%usual_endpoint_count /= 1) then
                status = PCA_INVALID_CONTRACT
                message = 'homoclinic pair must contain exactly one X and one usual endpoint'
                return
            end if
            if (.not. ieee_is_finite(pairs(i)%x_coefficient) .or. &
                    .not. ieee_is_finite(pairs(i)%usual_coefficient) .or. &
                    pairs(i)%x_coefficient <= 0.0_dp .or. &
                    pairs(i)%usual_coefficient <= 0.0_dp) then
                status = PCA_INVALID_CONTRACT
                message = 'homoclinic coefficients are invalid'
                return
            end if
            do j = 1, i - 1
                if (pairs(j)%pair_id == pairs(i)%pair_id .or. &
                        pairs(j)%separatrix_identity == pairs(i)%separatrix_identity) then
                    status = PCA_INVALID_CONTRACT
                    message = 'homoclinic pair identity is duplicated'
                    return
                end if
            end do
            call pair_certificate(pairs(i), local_status)
            if (local_status == PCA_PROVIDER_UNAVAILABLE) then
                status = PCA_PROVIDER_UNAVAILABLE
                message = 'Fortsym homoclinic-pair provider is unavailable'
                return
            end if
            if (local_status /= PCA_SUCCESS) then
                status = PCA_CALLBACK_FAILURE
                message = 'Fortsym homoclinic-pair certificate failed'
                return
            end if
            x_count = 0
            usual_count = 0
            do j = 1, size(intervals)
                if (intervals(j)%lower_homoclinic_pair_id == pairs(i)%pair_id) then
                    if (intervals(j)%lower_boundary_kind == PCA_BOUNDARY_HOMOCLINIC_X) &
                        x_count = x_count + 1
                    if (intervals(j)%lower_boundary_kind == PCA_BOUNDARY_HOMOCLINIC_USUAL) &
                        usual_count = usual_count + 1
                end if
                if (intervals(j)%upper_homoclinic_pair_id == pairs(i)%pair_id) then
                    if (intervals(j)%upper_boundary_kind == PCA_BOUNDARY_HOMOCLINIC_X) &
                        x_count = x_count + 1
                    if (intervals(j)%upper_boundary_kind == PCA_BOUNDARY_HOMOCLINIC_USUAL) &
                        usual_count = usual_count + 1
                end if
            end do
            if (x_count /= 1 .or. usual_count /= 1) then
                status = PCA_INVALID_CONTRACT
                message = 'homoclinic pair does not occur exactly once as X and usual'
                return
            end if
        end do
        do i = 1, size(intervals)
            referenced = intervals(i)%lower_homoclinic_pair_id > 0 .or. &
                intervals(i)%upper_homoclinic_pair_id > 0
            if (referenced .and. size(pairs) == 0) then
                status = PCA_INCOMPLETE_ATLAS
                message = 'homoclinic interval has no pair record'
                return
            end if
            if (referenced) then
                if (intervals(i)%lower_homoclinic_pair_id > 0) then
                    if (.not. pair_id_is_present( &
                            intervals(i)%lower_homoclinic_pair_id, pairs)) then
                        status = PCA_INVALID_CONTRACT
                        message = 'lower endpoint references an unknown homoclinic pair'
                        return
                    end if
                end if
                if (intervals(i)%upper_homoclinic_pair_id > 0) then
                    if (.not. pair_id_is_present( &
                            intervals(i)%upper_homoclinic_pair_id, pairs)) then
                        status = PCA_INVALID_CONTRACT
                        message = 'upper endpoint references an unknown homoclinic pair'
                        return
                    end if
                end if
            end if
        end do
    end subroutine validate_homoclinic_pairs

    pure logical function valid_boundary_kind(kind)
        integer, intent(in) :: kind
        valid_boundary_kind = kind >= PCA_BOUNDARY_REGULAR .and. &
            kind <= PCA_BOUNDARY_HOMOCLINIC_USUAL
    end function valid_boundary_kind

    pure logical function valid_boundary_pair_reference(kind, pair_id)
        integer, intent(in) :: kind, pair_id
        if (kind == PCA_BOUNDARY_HOMOCLINIC_X .or. &
                kind == PCA_BOUNDARY_HOMOCLINIC_USUAL) then
            valid_boundary_pair_reference = pair_id > 0
        else
            valid_boundary_pair_reference = pair_id == 0
        end if
    end function valid_boundary_pair_reference

    pure logical function same_interval_identity(first, second)
        type(poincare_allowed_interval_t), intent(in) :: first, second
        same_interval_identity = first%branch_id == second%branch_id .and. &
            first%sigma == second%sigma .and. &
            first%local_component_id == second%local_component_id .and. &
            first%psi_star_parameterization_copy_id == &
            second%psi_star_parameterization_copy_id
    end function same_interval_identity

    pure logical function same_parameterization(first, second)
        type(poincare_allowed_interval_t), intent(in) :: first, second
        same_parameterization = first%branch_id == second%branch_id .and. &
            first%sigma == second%sigma .and. &
            first%psi_star_parameterization_copy_id == &
            second%psi_star_parameterization_copy_id
    end function same_parameterization

    pure logical function pair_id_is_present(pair_id, pairs)
        integer, intent(in) :: pair_id
        type(poincare_homoclinic_pair_t), intent(in) :: pairs(:)
        integer :: i
        pair_id_is_present = .false.
        do i = 1, size(pairs)
            if (pairs(i)%pair_id == pair_id) then
                pair_id_is_present = .true.
                return
            end if
        end do
    end function pair_id_is_present

    pure logical function valid_context(context)
        type(poincare_phase_space_context_t), intent(in) :: context
        valid_context = ieee_is_finite(context%H0) .and. &
            ieee_is_finite(context%JK) .and. ieee_is_finite(context%sigma) .and. &
            ieee_is_finite(context%psi_star) .and. &
            ieee_is_finite(context%dpsi_star_dRc) .and. &
            abs(context%sigma) == 1.0_dp .and. &
            context%dpsi_star_dRc /= 0.0_dp .and. &
            context%field_input_id > 0 .and. &
            context%equilibrium_input_id > 0 .and. &
            context%profile_input_id > 0
    end function valid_context

    pure logical function same_context_ids(first, second)
        type(poincare_phase_space_context_t), intent(in) :: first, second
        same_context_ids = first%field_input_id == second%field_input_id .and. &
            first%equilibrium_input_id == second%equilibrium_input_id .and. &
            first%profile_input_id == second%profile_input_id
    end function same_context_ids

    subroutine canonical_sort_intervals(intervals)
        type(poincare_allowed_interval_t), intent(inout) :: intervals(:)
        type(poincare_allowed_interval_t) :: held
        integer :: i, j
        do i = 2, size(intervals)
            held = intervals(i)
            j = i - 1
            do while (j >= 1)
                if (.not. canonical_interval_less(held, intervals(j))) exit
                intervals(j + 1) = intervals(j)
                j = j - 1
            end do
            intervals(j + 1) = held
        end do
    end subroutine canonical_sort_intervals

    pure logical function canonical_interval_less(first, second)
        type(poincare_allowed_interval_t), intent(in) :: first, second
        canonical_interval_less = first%branch_id < second%branch_id
        if (first%branch_id == second%branch_id) then
            canonical_interval_less = first%sigma < second%sigma
            if (first%sigma == second%sigma) then
                canonical_interval_less = first%local_component_id < &
                    second%local_component_id
                if (first%local_component_id == second%local_component_id) then
                    canonical_interval_less = first%psi_star_parameterization_copy_id < &
                        second%psi_star_parameterization_copy_id
                    if (first%psi_star_parameterization_copy_id == &
                            second%psi_star_parameterization_copy_id) then
                        canonical_interval_less = first%lower < second%lower
                        if (first%lower == second%lower) then
                            canonical_interval_less = first%upper < second%upper
                        end if
                    end if
                end if
            end if
        end if
    end function canonical_interval_less

end module gc_poincare_cut_atlas
