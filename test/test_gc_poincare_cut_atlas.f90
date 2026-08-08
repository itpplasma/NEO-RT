program test_gc_poincare_cut_atlas
    !! The success callbacks below are independent test oracles only.  They
    !! are not production Fortsym providers; production wiring must return
    !! PCA_PROVIDER_UNAVAILABLE until those generated providers exist.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gc_poincare_cut_atlas, only: PCA_CALLBACK_FAILURE, PCA_DOMAIN_NONE, &
        PCA_INCOMPLETE_ATLAS, PCA_INVALID_CERTIFICATE, PCA_INVALID_CONTRACT, &
        PCA_INVALID_STATE, &
        PCA_UNKNOWN_MULTIPLICITY, &
        PCA_PROVIDER_UNAVAILABLE, PCA_SUCCESS, PCA_UNRESOLVED_TANGENCY, &
        PCA_WALL_NONE, PCA_BOUNDARY_REGULAR, PCA_BOUNDARY_REFLECTING, &
        PCA_BOUNDARY_HOMOCLINIC_USUAL, PCA_BOUNDARY_HOMOCLINIC_X, &
        combine_allowed_intervals, make_symmetric_midplane_atlas, &
        poincare_allowed_interval_t, poincare_branch_disjointness, &
        poincare_cut_atlas_t, poincare_cut_branch_evaluator, poincare_cut_branch_t, &
        poincare_fortsym_allowed_interval_certificate, &
        poincare_fortsym_branch_certificate, poincare_fortsym_cdot_enclosure_verifier, &
        poincare_fortsym_context_binding_certificate, &
        poincare_fortsym_eq13_cut_kernel, poincare_fortsym_exactly_two_provider, &
        poincare_fortsym_homoclinic_pair_certificate, &
        poincare_fortsym_phase_space_density_kernel, &
        poincare_fortsym_return_state_certificate, poincare_global_component_t, &
        poincare_homoclinic_pair_t, poincare_phase_space_context_t, &
        poincare_return_cycle_evidence_t, poincare_state_evidence_t, &
        poincare_exactly_two_certificate_t, &
        poincare_validated_cut_atlas_t, &
        validate_poincare_cut_atlas, validate_return_cycle_evidence
    implicit none

    call test_constructor_fails_closed()
    call test_missing_branch_provider_fails_closed()
    call test_structural_branch_identity_rejected()
    call test_missing_outboard_branch_rejected()
    call test_validated_token_is_issued_only_after_provider()
    call test_unissued_token_rejected()
    call test_unordered_events_rejected()
    call test_cdot_tangency_rejected()
    call test_orientation_is_derived_from_cdot()
    call test_invalid_exactly_two_certificate_rejected()
    call test_crossing_input_identity_mismatch_rejected()
    call test_return_invariant_mismatch_rejected()
    call test_missing_context_provider_fails_closed()
    call test_homoclinic_pair_cardinality_rejected()
    call test_homoclinic_pair_provider_is_required()
    call test_homoclinic_pair_contract()
    call test_one_sided_homoclinic_endpoints()
    call test_two_different_pair_endpoints()
    write (*, '(a)') 'PASS: gc_poincare_cut_atlas adversarial contract tests'

contains

    subroutine test_constructor_fails_closed()
        type(poincare_cut_branch_t) :: branch
        type(poincare_cut_atlas_t) :: atlas
        integer :: status
        character(len=256) :: message
        branch = test_branch()
        call make_symmetric_midplane_atlas(branch, atlas, status, message)
        call require(status == PCA_INCOMPLETE_ATLAS, &
            'unproved constructor reported success')
    end subroutine test_constructor_fails_closed

    subroutine test_missing_branch_provider_fails_closed()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, provider_unavailable, token, status, message)
        call require(status == PCA_PROVIDER_UNAVAILABLE, &
            'missing Fortsym branch provider was accepted')
    end subroutine test_missing_branch_provider_fails_closed

    subroutine test_structural_branch_identity_rejected()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        integer :: status
        character(len=256) :: message
        atlas = test_atlas(2)
        atlas%branches(2)%branch_id = 1
        call require(size(atlas%branches) == 2, &
            'duplicate-identity fixture did not contain two branches')
        call require(atlas%branches(1)%branch_id == 1, &
            'duplicate-identity fixture changed branch one')
        call require(atlas%branches(2)%branch_id == 1, &
            'duplicate-identity fixture did not change branch two')
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, provider_unavailable, token, status, message)
        call require(status == PCA_INVALID_CONTRACT, &
            'duplicate branch identity was accepted')
    end subroutine test_structural_branch_identity_rejected

    subroutine test_missing_outboard_branch_rejected()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        atlas%branches(1)%outboard_only = .true.
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, provider_unavailable, token, status, message)
        call require(status == PCA_INVALID_CONTRACT, &
            'outboard-only symmetric branch was accepted')
    end subroutine test_missing_outboard_branch_rejected

    subroutine test_validated_token_is_issued_only_after_provider()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test provider did not issue token')
    end subroutine test_validated_token_is_issued_only_after_provider

    subroutine test_unissued_token_rejected()
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_return_cycle_evidence_t) :: evidence
        integer :: status
        character(len=256) :: message
        call validate_return_cycle_evidence(token, evidence, certify_cdot, &
            certify_density, certify_context_binding, certify_exactly_two, &
            certify_return_state, status, message)
        call require(status == PCA_INCOMPLETE_ATLAS, &
            'raw or unissued atlas token was accepted')
    end subroutine test_unissued_token_rejected

    subroutine test_unordered_events_rejected()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_return_cycle_evidence_t) :: evidence
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        evidence = test_evidence()
        evidence%crossings(1)%time = 2.0_dp
        evidence%crossings(2)%time = 1.0_dp
        call validate_return_cycle_evidence(token, evidence, certify_cdot, &
            certify_density, certify_context_binding, certify_exactly_two, &
            certify_return_state, status, message)
        call require(status == PCA_INVALID_STATE, 'unordered events were accepted')
    end subroutine test_unordered_events_rejected

    subroutine test_cdot_tangency_rejected()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_return_cycle_evidence_t) :: evidence
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        evidence = test_evidence()
        call validate_return_cycle_evidence(token, evidence, certify_tangent_cdot, &
            certify_density, certify_context_binding, certify_exactly_two, &
            certify_return_state, status, message)
        call require(status == PCA_UNRESOLVED_TANGENCY, &
            'zero-containing Cdot enclosure was accepted')
    end subroutine test_cdot_tangency_rejected

    subroutine test_orientation_is_derived_from_cdot()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_return_cycle_evidence_t) :: evidence
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        evidence = test_evidence()
        call validate_return_cycle_evidence(token, evidence, certify_cdot, &
            certify_density, certify_context_binding, certify_exactly_two, &
            certify_return_state, status, message)
        call require(status == PCA_SUCCESS, &
            'Cdot-derived opposite/same event ordering was rejected')
    end subroutine test_orientation_is_derived_from_cdot

    subroutine test_invalid_exactly_two_certificate_rejected()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_return_cycle_evidence_t) :: evidence
        integer :: status
        character(len=256) :: message

        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        evidence = test_evidence()
        call validate_return_cycle_evidence(token, evidence, certify_cdot, &
            certify_density, certify_context_binding, malformed_exactly_two, &
            certify_return_state, status, message)
        call require(status == PCA_INVALID_CERTIFICATE, &
            'untyped two-event observation was accepted as a theorem')
    end subroutine test_invalid_exactly_two_certificate_rejected

    subroutine test_crossing_input_identity_mismatch_rejected()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_return_cycle_evidence_t) :: evidence
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        evidence = test_evidence()
        evidence%crossings(1)%context%field_input_id = 999
        call validate_return_cycle_evidence(token, evidence, certify_cdot, &
            certify_density, certify_context_binding, certify_exactly_two, &
            certify_return_state, status, message)
        call require(status == PCA_INVALID_CONTRACT, &
            'crossing field identity mismatch was accepted')
        evidence = test_evidence()
        evidence%return_state%context%profile_input_id = 999
        call validate_return_cycle_evidence(token, evidence, certify_cdot, &
            certify_density, certify_context_binding, certify_exactly_two, &
            certify_return_state, status, message)
        call require(status == PCA_INVALID_CONTRACT, &
            'return profile identity mismatch was accepted')
    end subroutine test_crossing_input_identity_mismatch_rejected

    subroutine test_return_invariant_mismatch_rejected()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_return_cycle_evidence_t) :: evidence
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        evidence = test_evidence()
        evidence%return_state%context%H0 = &
            evidence%launch_state%context%H0 + 1.0e-3_dp
        call validate_return_cycle_evidence(token, evidence, certify_cdot, &
            certify_density, certify_context_binding, certify_exactly_two, &
            certify_return_state, status, message)
        call require(status == PCA_INVALID_STATE, &
            'return invariant mismatch was accepted')
    end subroutine test_return_invariant_mismatch_rejected

    subroutine test_missing_context_provider_fails_closed()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_return_cycle_evidence_t) :: evidence
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        evidence = test_evidence()
        call validate_return_cycle_evidence(token, evidence, certify_cdot, &
            certify_density, provider_unavailable_context, certify_exactly_two, &
            certify_return_state, status, message)
        call require(status == PCA_PROVIDER_UNAVAILABLE, &
            'missing invariant-context provider was accepted')
    end subroutine test_missing_context_provider_fails_closed

    subroutine test_homoclinic_pair_cardinality_rejected()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_allowed_interval_t) :: intervals(2)
        type(poincare_homoclinic_pair_t) :: pairs(1)
        type(poincare_global_component_t), allocatable :: components(:)
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        intervals(1) = test_homoclinic_interval()
        intervals(2) = test_homoclinic_interval()
        intervals(2)%local_component_id = 2
        intervals(2)%lower = 2.0_dp
        intervals(2)%upper = 3.0_dp
        pairs(1) = test_pair()
        call combine_allowed_intervals(token, intervals, pairs, certify_interval, &
            certify_pair, components, status, message)
        call require(status == PCA_INVALID_CONTRACT, &
            'homoclinic pair with two X endpoints was accepted')
    end subroutine test_homoclinic_pair_cardinality_rejected

    subroutine test_homoclinic_pair_provider_is_required()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_allowed_interval_t) :: intervals(1)
        type(poincare_homoclinic_pair_t) :: pairs(1)
        type(poincare_global_component_t), allocatable :: components(:)
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        intervals(1) = test_homoclinic_interval()
        pairs(1) = test_pair()
        call combine_allowed_intervals(token, intervals, pairs, certify_interval, &
            provider_unavailable_pair, components, status, message)
        call require(status == PCA_PROVIDER_UNAVAILABLE, &
            'missing homoclinic Fortsym provider was accepted')
    end subroutine test_homoclinic_pair_provider_is_required

    subroutine test_homoclinic_pair_contract()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_allowed_interval_t) :: intervals(1)
        type(poincare_homoclinic_pair_t) :: pairs(1)
        type(poincare_global_component_t), allocatable :: components(:)
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        intervals(1) = test_homoclinic_interval()
        pairs(1) = test_pair()
        call combine_allowed_intervals(token, intervals, pairs, certify_interval, &
            certify_pair, components, status, message)
        call require(status == PCA_SUCCESS, 'valid homoclinic pair was rejected')
        call require(size(components) == 1, 'component count changed')
    end subroutine test_homoclinic_pair_contract

    subroutine test_one_sided_homoclinic_endpoints()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_allowed_interval_t) :: intervals(2)
        type(poincare_homoclinic_pair_t) :: pairs(1)
        type(poincare_global_component_t), allocatable :: components(:)
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        intervals(1) = test_homoclinic_interval()
        intervals(1)%lower_boundary_kind = PCA_BOUNDARY_REGULAR
        intervals(1)%lower_homoclinic_pair_id = 0
        intervals(1)%upper_boundary_kind = PCA_BOUNDARY_HOMOCLINIC_X
        intervals(1)%upper_homoclinic_pair_id = 7
        intervals(2) = test_homoclinic_interval()
        intervals(2)%local_component_id = 2
        intervals(2)%lower_boundary_kind = PCA_BOUNDARY_HOMOCLINIC_USUAL
        intervals(2)%lower_homoclinic_pair_id = 7
        intervals(2)%upper_boundary_kind = PCA_BOUNDARY_REFLECTING
        intervals(2)%upper_homoclinic_pair_id = 0
        intervals(2)%lower = 1.0_dp
        intervals(2)%upper = 2.0_dp
        pairs(1) = test_pair()
        call combine_allowed_intervals(token, intervals, pairs, certify_interval, &
            certify_pair, components, status, message)
        call require(status == PCA_SUCCESS, &
            'valid one-sided homoclinic endpoint was rejected')
    end subroutine test_one_sided_homoclinic_endpoints

    subroutine test_two_different_pair_endpoints()
        type(poincare_cut_atlas_t) :: atlas
        type(poincare_validated_cut_atlas_t) :: token
        type(poincare_allowed_interval_t) :: intervals(3)
        type(poincare_homoclinic_pair_t) :: pairs(2)
        type(poincare_global_component_t), allocatable :: components(:)
        integer :: status
        character(len=256) :: message
        atlas = test_atlas()
        call validate_poincare_cut_atlas(atlas, evaluate_line, evaluate_cut, &
            certify_disjoint, certify_branch, token, status, message)
        call require(status == PCA_SUCCESS, 'test token setup failed')
        intervals(1) = test_homoclinic_interval()
        intervals(1)%lower_boundary_kind = PCA_BOUNDARY_HOMOCLINIC_X
        intervals(1)%lower_homoclinic_pair_id = 7
        intervals(1)%upper_boundary_kind = PCA_BOUNDARY_REGULAR
        intervals(1)%upper_homoclinic_pair_id = 0
        intervals(2) = test_homoclinic_interval()
        intervals(2)%local_component_id = 2
        intervals(2)%lower_boundary_kind = PCA_BOUNDARY_HOMOCLINIC_USUAL
        intervals(2)%lower_homoclinic_pair_id = 7
        intervals(2)%upper_boundary_kind = PCA_BOUNDARY_HOMOCLINIC_X
        intervals(2)%upper_homoclinic_pair_id = 8
        intervals(2)%lower = 1.0_dp
        intervals(2)%upper = 2.0_dp
        intervals(3) = test_homoclinic_interval()
        intervals(3)%local_component_id = 3
        intervals(3)%lower_boundary_kind = PCA_BOUNDARY_HOMOCLINIC_USUAL
        intervals(3)%lower_homoclinic_pair_id = 8
        intervals(3)%upper_boundary_kind = PCA_BOUNDARY_REFLECTING
        intervals(3)%upper_homoclinic_pair_id = 0
        intervals(3)%lower = 2.0_dp
        intervals(3)%upper = 3.0_dp
        pairs(1) = test_pair()
        pairs(2) = test_pair(8)
        call combine_allowed_intervals(token, intervals, pairs, certify_interval, &
            certify_pair, components, status, message)
        call require(status == PCA_SUCCESS, &
            'interval with two different homoclinic pairs was rejected')
    end subroutine test_two_different_pair_endpoints

    function test_branch() result(branch)
        type(poincare_cut_branch_t) :: branch
        branch%branch_id = 1
        branch%component_identity_seed = 101
        branch%field_input_id = 11
        branch%equilibrium_input_id = 13
        branch%orientation = 1
        branch%parameter_units = 'physical_R'
        branch%parameter_lower = -1.0_dp
        branch%parameter_upper = 1.0_dp
        branch%endpoint_start = [1.0_dp, 0.0_dp]
        branch%endpoint_end = [3.0_dp, 0.0_dp]
        branch%endpoint_start_known = .true.
        branch%endpoint_end_known = .true.
        branch%cut_residual_tolerance = 1.0e-12_dp
        branch%physical_r_ordered = .true.
        branch%continuous_inboard_axis_outboard = .true.
    end function test_branch

    function test_atlas(n) result(atlas)
        integer, intent(in), optional :: n
        type(poincare_cut_atlas_t) :: atlas
        integer :: count, i
        count = 1
        if (present(n)) count = n
        allocate(atlas%branches(count))
        !! Initialize each element explicitly.  This keeps the malformed
        !! duplicate-identity fixture independent of compiler handling of
        !! scalar expansion into an allocatable derived-type array.
        do i = 1, count
            atlas%branches(i) = test_branch()
        end do
        atlas%expected_branch_count = count
        atlas%symmetric_midplane_specialization = count == 1
        if (count == 2) then
            atlas%branches(2)%branch_id = 2
            atlas%branches(2)%component_identity_seed = 202
            atlas%branches(2)%physical_r_ordered = .false.
            atlas%branches(2)%continuous_inboard_axis_outboard = .false.
        end if
    end function test_atlas

    function test_evidence() result(evidence)
        type(poincare_return_cycle_evidence_t) :: evidence
        evidence%launch_branch_id = 1
        evidence%launch_parameter = -0.75_dp
        evidence%launch_position = [1.25_dp, 0.0_dp]
        evidence%launch_state%time = 0.0_dp
        evidence%launch_state%state = [1.25_dp, 0.0_dp, 0.0_dp, 1.0_dp, 0.1_dp]
        evidence%return_state = evidence%launch_state
        evidence%return_state%time = 2.0_dp
        evidence%state_tolerance = 1.0e-10_dp
        evidence%wall_outcome = PCA_WALL_NONE
        evidence%domain_outcome = PCA_DOMAIN_NONE
        evidence%launch_state%context = test_context()
        evidence%return_state%context = test_context()
        allocate(evidence%crossings(2))
        evidence%crossings(1)%branch_id = 1
        evidence%crossings(1)%parameter = -0.5_dp
        evidence%crossings(1)%position = [1.5_dp, 0.0_dp]
        evidence%crossings(1)%context = test_context()
        evidence%crossings(1)%time = 1.0_dp
        evidence%crossings(2)%branch_id = 1
        evidence%crossings(2)%parameter = 0.5_dp
        evidence%crossings(2)%position = [2.5_dp, 0.0_dp]
        evidence%crossings(2)%context = test_context()
        evidence%crossings(2)%time = 2.0_dp
    end function test_evidence

    function test_context() result(context)
        type(poincare_phase_space_context_t) :: context
        context%H0 = 1.0_dp
        context%JK = 0.2_dp
        context%sigma = 1.0_dp
        context%psi_star = 0.3_dp
        context%dpsi_star_dRc = 0.4_dp
        context%field_input_id = 11
        context%equilibrium_input_id = 13
        context%profile_input_id = 12
    end function test_context

    function test_homoclinic_interval() result(interval)
        type(poincare_allowed_interval_t) :: interval
        interval%branch_id = 1
        interval%sigma = 1
        interval%local_component_id = 1
        interval%psi_star_parameterization_copy_id = 1
        interval%lower_boundary_kind = PCA_BOUNDARY_HOMOCLINIC_X
        interval%upper_boundary_kind = PCA_BOUNDARY_HOMOCLINIC_USUAL
        interval%lower_homoclinic_pair_id = 7
        interval%upper_homoclinic_pair_id = 7
        interval%lower = 0.0_dp
        interval%upper = 1.0_dp
    end function test_homoclinic_interval

    function test_pair(requested_id) result(pair)
        integer, intent(in), optional :: requested_id
        type(poincare_homoclinic_pair_t) :: pair
        pair%pair_id = 7
        if (present(requested_id)) pair%pair_id = requested_id
        pair%separatrix_identity = 10 * pair%pair_id
        pair%x_endpoint_count = 1
        pair%usual_endpoint_count = 1
        pair%x_coefficient = 2.0_dp
        pair%usual_coefficient = 1.0_dp
    end function test_pair

    subroutine evaluate_line(branch, parameter, position, value, derivative, status)
        type(poincare_cut_branch_t), intent(in) :: branch
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: position(2), value, derivative
        integer, intent(out) :: status
        position = [2.0_dp + parameter, 0.0_dp]
        value = parameter
        derivative = real(branch%orientation, dp)
        status = PCA_SUCCESS
    end subroutine evaluate_line

    subroutine evaluate_cut(branch, parameter, position, residual, derivative, status)
        type(poincare_cut_branch_t), intent(in) :: branch
        real(dp), intent(in) :: parameter, position(2)
        real(dp), intent(out) :: residual, derivative
        integer, intent(out) :: status
        residual = 0.0_dp
        derivative = 1.0_dp
        status = PCA_SUCCESS
    end subroutine evaluate_cut

    subroutine certify_disjoint(branch_a, branch_b, disjoint, status)
        type(poincare_cut_branch_t), intent(in) :: branch_a, branch_b
        logical, intent(out) :: disjoint
        integer, intent(out) :: status
        disjoint = branch_a%branch_id /= branch_b%branch_id
        status = PCA_SUCCESS
    end subroutine certify_disjoint

    subroutine certify_branch(branch, status)
        type(poincare_cut_branch_t), intent(in) :: branch
        integer, intent(out) :: status
        if (branch%branch_id /= 1 .or. branch%component_identity_seed /= 101 .or. &
                branch%field_input_id /= 11 .or. &
                branch%equilibrium_input_id /= 13 .or. &
                branch%parameter_lower /= -1.0_dp .or. &
                branch%parameter_upper /= 1.0_dp .or. &
                .not. branch%endpoint_start_known .or. &
                .not. branch%endpoint_end_known .or. &
                .not. branch%physical_r_ordered .or. &
                .not. branch%continuous_inboard_axis_outboard .or. &
                branch%outboard_only) then
            status = PCA_INVALID_CONTRACT
            return
        end if
        status = PCA_SUCCESS
    end subroutine certify_branch

    subroutine provider_unavailable(branch, status)
        type(poincare_cut_branch_t), intent(in) :: branch
        integer, intent(out) :: status
        status = PCA_PROVIDER_UNAVAILABLE
    end subroutine provider_unavailable

    subroutine certify_interval(interval, status)
        type(poincare_allowed_interval_t), intent(in) :: interval
        integer, intent(out) :: status
        if (interval%branch_id /= 1 .or. interval%sigma /= 1 .or. &
                interval%local_component_id <= 0 .or. &
                interval%psi_star_parameterization_copy_id /= 1 .or. &
                interval%lower < 0.0_dp .or. interval%upper > 3.0_dp .or. &
                interval%upper <= interval%lower) then
            status = PCA_INVALID_CONTRACT
            return
        end if
        status = PCA_SUCCESS
    end subroutine certify_interval

    subroutine certify_pair(pair, status)
        type(poincare_homoclinic_pair_t), intent(in) :: pair
        integer, intent(out) :: status
        if (pair%pair_id <= 0 .or. &
                pair%separatrix_identity /= 10 * pair%pair_id .or. &
                pair%x_endpoint_count /= 1 .or. pair%usual_endpoint_count /= 1 .or. &
                pair%x_coefficient /= 2.0_dp .or. pair%usual_coefficient /= 1.0_dp) then
            status = PCA_INVALID_CONTRACT
            return
        end if
        status = PCA_SUCCESS
    end subroutine certify_pair

    subroutine provider_unavailable_pair(pair, status)
        type(poincare_homoclinic_pair_t), intent(in) :: pair
        integer, intent(out) :: status
        status = PCA_PROVIDER_UNAVAILABLE
    end subroutine provider_unavailable_pair

    subroutine certify_cdot(branch, parameter, position, context, lower, upper, status)
        type(poincare_cut_branch_t), intent(in) :: branch
        real(dp), intent(in) :: parameter, position(2)
        type(poincare_phase_space_context_t), intent(in) :: context
        real(dp), intent(out) :: lower, upper
        integer, intent(out) :: status
        lower = -1.0_dp
        upper = 1.0_dp
        if (parameter < -0.6_dp) then
            lower = 0.5_dp
            upper = 1.5_dp
        else if (parameter > 0.0_dp) then
            lower = 0.5_dp
            upper = 1.5_dp
        end if
        if (parameter > -0.6_dp .and. parameter < 0.0_dp) then
            lower = -1.5_dp
            upper = -0.5_dp
        end if
        status = PCA_SUCCESS
    end subroutine certify_cdot

    subroutine certify_tangent_cdot(branch, parameter, position, context, lower, upper, status)
        type(poincare_cut_branch_t), intent(in) :: branch
        real(dp), intent(in) :: parameter, position(2)
        type(poincare_phase_space_context_t), intent(in) :: context
        real(dp), intent(out) :: lower, upper
        integer, intent(out) :: status
        lower = -0.5_dp
        upper = 0.5_dp
        status = PCA_SUCCESS
    end subroutine certify_tangent_cdot

    subroutine certify_density(branch, parameter, position, context, signed_jacobian, &
            positive_density, status)
        type(poincare_cut_branch_t), intent(in) :: branch
        real(dp), intent(in) :: parameter, position(2)
        type(poincare_phase_space_context_t), intent(in) :: context
        real(dp), intent(out) :: signed_jacobian, positive_density
        integer, intent(out) :: status
        signed_jacobian = context%dpsi_star_dRc
        positive_density = 1.0_dp
        status = PCA_SUCCESS
    end subroutine certify_density

    subroutine certify_context_binding(launch_context, candidate_context, status)
        type(poincare_phase_space_context_t), intent(in) :: launch_context
        type(poincare_phase_space_context_t), intent(in) :: candidate_context
        integer, intent(out) :: status
        real(dp), parameter :: invariant_tolerance = 1.0e-12_dp
        if (abs(candidate_context%H0 - launch_context%H0) > invariant_tolerance .or. &
                abs(candidate_context%JK - launch_context%JK) > &
                invariant_tolerance .or. &
                abs(candidate_context%psi_star - launch_context%psi_star) > &
                invariant_tolerance) then
            status = PCA_INVALID_STATE
            return
        end if
        status = PCA_SUCCESS
    end subroutine certify_context_binding

    subroutine provider_unavailable_context(launch_context, candidate_context, status)
        type(poincare_phase_space_context_t), intent(in) :: launch_context
        type(poincare_phase_space_context_t), intent(in) :: candidate_context
        integer, intent(out) :: status
        status = PCA_PROVIDER_UNAVAILABLE
    end subroutine provider_unavailable_context

    subroutine certify_exactly_two(atlas, evidence, certificate, status)
        type(poincare_validated_cut_atlas_t), intent(in) :: atlas
        type(poincare_return_cycle_evidence_t), intent(in) :: evidence
        type(poincare_exactly_two_certificate_t), intent(out) :: certificate
        integer, intent(out) :: status
        certificate = poincare_exactly_two_certificate_t()
        if (.not. allocated(evidence%crossings)) then
            status = PCA_UNKNOWN_MULTIPLICITY
            return
        end if
        if (size(evidence%crossings) /= 2) then
            status = PCA_UNKNOWN_MULTIPLICITY
            return
        end if
        if (evidence%launch_branch_id /= 1) then
            status = PCA_INVALID_CONTRACT
            return
        end if
        certificate%certificate_id = 2401
        certificate%crossing_count = 2
        certificate%exactly_two_proved = .true.
        status = PCA_SUCCESS
    end subroutine certify_exactly_two

    subroutine malformed_exactly_two(atlas, evidence, certificate, status)
        type(poincare_validated_cut_atlas_t), intent(in) :: atlas
        type(poincare_return_cycle_evidence_t), intent(in) :: evidence
        type(poincare_exactly_two_certificate_t), intent(out) :: certificate
        integer, intent(out) :: status

        certificate = poincare_exactly_two_certificate_t()
        status = PCA_SUCCESS
    end subroutine malformed_exactly_two

    subroutine certify_return_state(launch_state, return_state, tolerance, status)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        type(poincare_state_evidence_t), intent(in) :: launch_state, return_state
        real(dp), intent(in) :: tolerance
        integer, intent(out) :: status
        if (.not. all(ieee_is_finite(launch_state%state)) .or. &
                .not. all(ieee_is_finite(return_state%state))) then
            status = PCA_INVALID_STATE
            return
        end if
        if (launch_state%context%field_input_id <= 0 .or. &
                launch_state%context%equilibrium_input_id <= 0 .or. &
                launch_state%context%profile_input_id <= 0) then
            status = PCA_INVALID_STATE
            return
        end if
        status = PCA_SUCCESS
    end subroutine certify_return_state

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) error stop message
    end subroutine require

end program test_gc_poincare_cut_atlas
