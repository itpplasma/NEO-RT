module neort_gc_model2_transport_dispatch
    !! Model-2-only seam between a physical provider and runtime delivery.
    !!
    !! The caller owns construction of the direct-EQDSK factory and supplies
    !! its certificates.  This module owns neither field evaluation nor a
    !! local eta reduction.  It only runs the certified nonlocal integral,
    !! preserves the provider's signed force slots, and exposes observed
    !! node coverage for the executable-owned runtime record.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_full_fow_runtime_delivery, only: &
        emit_gc_full_fow_runtime_surface_record
    use neort_gc_full_fow_runtime_metadata, only: &
        GC_FULL_FOW_CONJUGATE_POLICY, GC_FULL_FOW_PREFACTOR_CONVENTION, &
        GC_FULL_FOW_REAL_FIELD_AMPLITUDE_CONVENTION, &
        GC_FULL_FOW_ACTION_CONVENTION, GC_FULL_FOW_BOUND_METHOD, &
        format_gc_full_fow_frequency_convention, &
        gc_full_fow_runtime_backend_state_t
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_SUCCESS, gc_nonlocal_resonance_options_t
    use neort_gc_nonlocal_transport_integral, only: &
        integrate_gc_nonlocal_transport
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_MAX_FORCE_VALUES, gc_nonlocal_transport_options_t, &
        gc_nonlocal_transport_observed_evidence_t, &
        gc_nonlocal_transport_provider_t, gc_nonlocal_transport_result_t

    implicit none
    private

    integer, parameter, public :: GC_MODEL2_DISPATCH_SUCCESS = 0
    integer, parameter, public :: GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE = 7101
    integer, parameter, public :: GC_MODEL2_DISPATCH_BACKEND_UNCERTIFIED = 7102
    integer, parameter, public :: GC_MODEL2_DISPATCH_INTEGRAL_FAILED = 7103
    integer, parameter, public :: GC_MODEL2_DISPATCH_FORCE_LAYOUT_INVALID = 7104
    integer, parameter, public :: GC_MODEL2_DISPATCH_NONFINITE = 7105
    integer, parameter, public :: GC_MODEL2_DISPATCH_NOT_CERTIFIED = 7106

    type, public :: gc_model2_force_layout_t
        integer :: d11_slot = 0
        integer :: d12_slot = 0
        integer :: torque_slot = 0
    end type gc_model2_force_layout_t

    type, public :: gc_model2_transport_options_t
        type(gc_nonlocal_transport_options_t) :: integral = &
            gc_nonlocal_transport_options_t()
        type(gc_model2_force_layout_t) :: force_layout = &
            gc_model2_force_layout_t()
        integer, allocatable :: poloidal_harmonics(:)
        integer :: toroidal_harmonic = 0
    end type gc_model2_transport_options_t

    type, public :: gc_model2_backend_evidence_t
        character(len=32) :: backend = 'cylindrical_full_fow'
        character(len=16) :: coordinates = 'R,Z,phi'
        character(len=1024) :: wall_actual_path = ''
        character(len=16) :: wall_units = ''
        character(len=64) :: wall_sha256 = ''
        logical :: factory_available = .false.
        logical :: field_certified = .false.
        logical :: profile_certified = .false.
        logical :: perturbation_certified = .false.
        logical :: topology_certified = .false.
        logical :: wall_certified = .false.
        logical :: canonical_measure_certified = .false.
        logical :: component_identity_certified = .false.
        character(len=64) :: perturbation_amplitude_convention = ''
        character(len=1024) :: perturbation_input_path = ''
        logical :: perturbation_provenance_certified = .false.
        character(len=64) :: phase_space_bound_method = ''
        logical :: orbit_step_refinement_certified = .false.
        real(dp) :: orbit_base_step = 0.0_dp
        real(dp) :: orbit_refined_step = 0.0_dp
        real(dp) :: orbit_period_refinement_error = 0.0_dp
        real(dp) :: orbit_delta_phi_refinement_error = 0.0_dp
        real(dp) :: orbit_omega_b_refinement_error = 0.0_dp
        real(dp) :: orbit_omega_phi_refinement_error = 0.0_dp
        real(dp) :: orbit_h_m_refinement_error = 0.0_dp
        real(dp) :: orbit_shell_refinement_error = 0.0_dp
        integer :: legacy_backend_entries = 0
        integer :: chart_fallback_entries = 0
    end type gc_model2_backend_evidence_t

    type, public :: gc_model2_observed_evidence_t
        !! Counters copied from the physical provider after this execution.
        !! The dispatcher never infers these categories from node status.
        integer :: physical_return_attempts = 0
        integer :: invariant_successes = 0
        integer :: invariant_failures = 0
        integer :: return_successes = 0
        integer :: return_no_return = 0
        integer :: return_radial_domain = 0
        integer :: return_wall_loss = 0
        integer :: return_errors = 0
        integer :: wall_not_hit = 0
        integer :: wall_loss = 0
        integer :: wall_errors = 0
        integer :: topology_certification_attempts = 0
        integer :: topology_certification_successes = 0
        logical :: return_accounting_complete = .false.
    end type gc_model2_observed_evidence_t

    type, public :: gc_model2_harmonic_execution_t
        integer :: poloidal_harmonic = 0
        type(gc_nonlocal_transport_result_t) :: integral
        real(dp) :: d11(3) = 0.0_dp
        real(dp) :: d12(3) = 0.0_dp
        real(dp) :: torque(3) = 0.0_dp
    end type gc_model2_harmonic_execution_t

    type, public :: gc_model2_transport_execution_t
        type(gc_nonlocal_transport_result_t) :: integral
        type(gc_model2_harmonic_execution_t), allocatable :: harmonics(:)
        type(gc_model2_observed_evidence_t) :: observed
        type(gc_full_fow_runtime_backend_state_t) :: runtime
        real(dp) :: d11 = 0.0_dp
        real(dp) :: d12 = 0.0_dp
        real(dp) :: torque = 0.0_dp
        real(dp) :: d11_class(3) = 0.0_dp
        real(dp) :: d12_class(3) = 0.0_dp
        real(dp) :: torque_class(3) = 0.0_dp
        integer :: integral_status = GC_NONLOCAL_SUCCESS
        integer :: attempted_nodes = 0
        integer :: certified_nodes = 0
        integer :: unresolved_nodes = 0
        logical :: force_slots_accepted = .false.
        logical :: certified = .false.
        character(len=256) :: invariant_status_coverage = ''
        character(len=256) :: return_status_coverage = ''
        character(len=256) :: wall_status_coverage = ''
    end type gc_model2_transport_execution_t

    public :: emit_gc_model2_runtime_record
    public :: execute_gc_model2_transport
    public :: finalize_gc_model2_transport_execution
    public :: gc_model2_dispatch_required

contains

    pure logical function gc_model2_dispatch_required(frequency_model)
        integer, intent(in) :: frequency_model

        gc_model2_dispatch_required = frequency_model == 2
    end function gc_model2_dispatch_required

    subroutine execute_gc_model2_transport(provider, harmonic_m, harmonic_n, &
            backend, options, execution, status)
        class(gc_nonlocal_transport_provider_t), intent(inout) :: provider
        integer, intent(in) :: harmonic_m, harmonic_n
        type(gc_model2_backend_evidence_t), intent(in) :: backend
        type(gc_model2_transport_options_t), intent(in) :: options
        type(gc_model2_transport_execution_t), intent(out) :: execution
        integer, intent(out) :: status

        type(gc_nonlocal_transport_result_t) :: harmonic_result
        type(gc_nonlocal_transport_observed_evidence_t) :: provider_evidence
        integer, allocatable :: poloidal_harmonics(:)
        integer :: local_status, harmonic_index, harmonic_count, active_n

        execution = gc_model2_transport_execution_t()
        status = GC_MODEL2_DISPATCH_BACKEND_UNCERTIFIED
        local_status = validate_backend(backend, .false.)
        if (local_status /= GC_MODEL2_DISPATCH_SUCCESS) then
            status = local_status
            return
        end if
        local_status = validate_force_layout(options)
        if (local_status /= GC_MODEL2_DISPATCH_SUCCESS) then
            status = local_status
            return
        end if
        active_n = harmonic_n
        if (options%toroidal_harmonic /= 0) active_n = options%toroidal_harmonic
        if (active_n == 0) then
            status = GC_MODEL2_DISPATCH_FORCE_LAYOUT_INVALID
            return
        end if
        if (allocated(options%poloidal_harmonics)) then
            harmonic_count = size(options%poloidal_harmonics)
            allocate(poloidal_harmonics(harmonic_count))
            poloidal_harmonics = options%poloidal_harmonics
        else
            ! The direct perturbation representation does not authorize an
            ! implicit harmonic range.  The executable must pass the exact
            ! runtime-selected set derived from its input/configuration.
            status = GC_MODEL2_DISPATCH_FORCE_LAYOUT_INVALID
            return
        end if
        if (harmonic_count < 1) then
            status = GC_MODEL2_DISPATCH_FORCE_LAYOUT_INVALID
            return
        end if
        call provider%begin_execution(local_status)
        if (local_status /= GC_NONLOCAL_SUCCESS) then
            status = GC_MODEL2_DISPATCH_BACKEND_UNCERTIFIED
            return
        end if
        allocate(execution%harmonics(harmonic_count))
        execution%integral = gc_nonlocal_transport_result_t()
        do harmonic_index = 1, harmonic_count
            call integrate_gc_nonlocal_transport(provider, poloidal_harmonics(harmonic_index), &
                active_n, options%integral, harmonic_result, local_status)
            execution%integral_status = local_status
            execution%harmonics(harmonic_index)%poloidal_harmonic = &
                poloidal_harmonics(harmonic_index)
            execution%harmonics(harmonic_index)%integral = harmonic_result
            if (local_status /= GC_NONLOCAL_SUCCESS .or. &
                    .not. harmonic_result%certified) then
                execution%integral = harmonic_result
                call record_observed_execution(execution, backend, options, &
                    poloidal_harmonics, active_n)
                status = GC_MODEL2_DISPATCH_INTEGRAL_FAILED
                return
            end if
            if (harmonic_index == 1) then
                execution%integral = harmonic_result
            else
                call add_integral_result(execution%integral, harmonic_result, &
                    local_status)
                if (local_status /= GC_NONLOCAL_SUCCESS) then
                    status = GC_MODEL2_DISPATCH_INTEGRAL_FAILED
                    return
                end if
            end if
            call accept_harmonic_slots(execution%harmonics(harmonic_index), &
                options%force_layout, local_status)
            if (local_status /= GC_MODEL2_DISPATCH_SUCCESS) then
                status = local_status
                return
            end if
            execution%d11_class = execution%d11_class + &
                execution%harmonics(harmonic_index)%d11
            execution%d12_class = execution%d12_class + &
                execution%harmonics(harmonic_index)%d12
            execution%torque_class = execution%torque_class + &
                execution%harmonics(harmonic_index)%torque
        end do
        execution%d11 = sum(execution%d11_class)
        execution%d12 = sum(execution%d12_class)
        execution%torque = sum(execution%torque_class)
        call provider%get_execution_evidence(provider_evidence, local_status)
        if (local_status /= GC_NONLOCAL_SUCCESS) then
            status = GC_MODEL2_DISPATCH_BACKEND_UNCERTIFIED
            return
        end if
        execution%observed%physical_return_attempts = &
            provider_evidence%physical_return_attempts
        execution%observed%invariant_successes = provider_evidence%invariant_successes
        execution%observed%invariant_failures = provider_evidence%invariant_failures
        execution%observed%return_successes = provider_evidence%return_successes
        execution%observed%return_no_return = provider_evidence%return_no_return
        execution%observed%return_radial_domain = provider_evidence%return_radial_domain
        execution%observed%return_wall_loss = provider_evidence%return_wall_loss
        execution%observed%return_errors = provider_evidence%return_errors
        execution%observed%wall_not_hit = provider_evidence%wall_not_hit
        execution%observed%wall_loss = provider_evidence%wall_loss
        execution%observed%wall_errors = provider_evidence%wall_errors
        execution%observed%topology_certification_attempts = &
            provider_evidence%topology_certification_attempts
        execution%observed%topology_certification_successes = &
            provider_evidence%topology_certification_successes
        execution%observed%return_accounting_complete = &
            provider_evidence%return_accounting_complete
        call record_observed_execution(execution, backend, options, &
            poloidal_harmonics, active_n)
        if (execution%certified_nodes /= execution%attempted_nodes .or. &
            execution%unresolved_nodes /= 0) then
            status = GC_MODEL2_DISPATCH_INTEGRAL_FAILED
            return
        end if
        execution%force_slots_accepted = .true.
        status = GC_MODEL2_DISPATCH_SUCCESS
    end subroutine execute_gc_model2_transport

    subroutine add_integral_result(total, addend, status)
        type(gc_nonlocal_transport_result_t), intent(inout) :: total
        type(gc_nonlocal_transport_result_t), intent(in) :: addend
        integer, intent(out) :: status

        status = GC_MODEL2_DISPATCH_INTEGRAL_FAILED
        if (.not. addend%certified) return
        if (total%nforce /= addend%nforce) return
        if (size(total%contribution) /= size(addend%contribution)) return
        total%contribution = total%contribution + addend%contribution
        total%class_contribution = total%class_contribution + &
            addend%class_contribution
        total%weighted_nodes = total%weighted_nodes + addend%weighted_nodes
        total%certified_nodes = total%certified_nodes + addend%certified_nodes
        total%ncomponents = total%ncomponents + addend%ncomponents
        total%nroots = total%nroots + addend%nroots
        status = GC_NONLOCAL_SUCCESS
    end subroutine add_integral_result

    subroutine accept_harmonic_slots(harmonic, layout, status)
        type(gc_model2_harmonic_execution_t), intent(inout) :: harmonic
        type(gc_model2_force_layout_t), intent(in) :: layout
        integer, intent(out) :: status

        integer :: i
        status = GC_MODEL2_DISPATCH_FORCE_LAYOUT_INVALID
        do i = 1, 3
            harmonic%d11(i) = harmonic%integral%class_contribution(i, layout%d11_slot)
            harmonic%d12(i) = harmonic%integral%class_contribution(i, layout%d12_slot)
            harmonic%torque(i) = harmonic%integral%class_contribution(i, layout%torque_slot)
        end do
        if (.not. all(ieee_is_finite([harmonic%d11, harmonic%d12, &
            harmonic%torque]))) then
            status = GC_MODEL2_DISPATCH_NONFINITE
            return
        end if
        status = GC_MODEL2_DISPATCH_SUCCESS
    end subroutine accept_harmonic_slots

    subroutine record_observed_execution(execution, backend, options, &
            poloidal_harmonics, active_n)
        type(gc_model2_transport_execution_t), intent(inout) :: execution
        type(gc_model2_backend_evidence_t), intent(in) :: backend
        type(gc_model2_transport_options_t), intent(in) :: options
        integer, intent(in) :: poloidal_harmonics(:), active_n

        execution%attempted_nodes = execution%integral%weighted_nodes
        execution%certified_nodes = execution%integral%certified_nodes
        execution%unresolved_nodes = execution%integral%unresolved_nodes
        execution%runtime%backend = backend%backend
        execution%runtime%coordinates = backend%coordinates
        execution%runtime%model = 2
        execution%runtime%frequency_model = 2
        execution%runtime%wall_actual_path = backend%wall_actual_path
        execution%runtime%wall_units = backend%wall_units
        execution%runtime%wall_sha256 = backend%wall_sha256
        execution%runtime%runtime_execution_complete = .false.
        execution%runtime%orbit_backend_certified = .false.
        execution%runtime%wall_certified = backend%wall_certified
        execution%runtime%canonical_measure_certified = .false.
        execution%runtime%component_identity_certified = .false.
        execution%runtime%cylindrical_backend_entries = &
            execution%attempted_nodes
        execution%runtime%legacy_backend_entries = backend%legacy_backend_entries
        execution%runtime%chart_fallback_entries = backend%chart_fallback_entries
        execution%runtime%real_field_amplitude_convention = &
            backend%perturbation_amplitude_convention
        if (len_trim(execution%runtime%real_field_amplitude_convention) == 0) &
            execution%runtime%real_field_amplitude_convention = &
                GC_FULL_FOW_REAL_FIELD_AMPLITUDE_CONVENTION
        execution%runtime%conjugate_policy = GC_FULL_FOW_CONJUGATE_POLICY
        execution%runtime%prefactor_convention = GC_FULL_FOW_PREFACTOR_CONVENTION
        execution%runtime%action_convention = GC_FULL_FOW_ACTION_CONVENTION
        execution%runtime%phase_space_bound_method = backend%phase_space_bound_method
        execution%runtime%orbit_step_refinement_certified = &
            backend%orbit_step_refinement_certified
        execution%runtime%orbit_base_step = backend%orbit_base_step
        execution%runtime%orbit_refined_step = backend%orbit_refined_step
        execution%runtime%orbit_period_refinement_error = &
            backend%orbit_period_refinement_error
        execution%runtime%orbit_delta_phi_refinement_error = &
            backend%orbit_delta_phi_refinement_error
        execution%runtime%orbit_omega_b_refinement_error = &
            backend%orbit_omega_b_refinement_error
        execution%runtime%orbit_omega_phi_refinement_error = &
            backend%orbit_omega_phi_refinement_error
        execution%runtime%orbit_h_m_refinement_error = backend%orbit_h_m_refinement_error
        execution%runtime%orbit_shell_refinement_error = backend%orbit_shell_refinement_error
        execution%runtime%perturbation_input_path = backend%perturbation_input_path
        execution%runtime%perturbation_provenance_certified = &
            backend%perturbation_provenance_certified
        execution%runtime%quadrature_base_h0_order = options%integral%h0_order
        execution%runtime%quadrature_base_jk_order = options%integral%jk_order
        execution%runtime%quadrature_refined_h0_order = &
            2*options%integral%h0_order
        execution%runtime%quadrature_refined_jk_order = &
            2*options%integral%jk_order
        execution%runtime%quadrature_relative_tolerance = &
            options%integral%quadrature_relative_tolerance
        execution%runtime%quadrature_absolute_tolerance = &
            options%integral%quadrature_absolute_tolerance
        execution%runtime%poloidal_harmonic_min = minval(poloidal_harmonics)
        execution%runtime%poloidal_harmonic_max = maxval(poloidal_harmonics)
        execution%runtime%poloidal_harmonic_count = size(poloidal_harmonics)
        execution%runtime%executed_harmonic_count = size(execution%harmonics)
        execution%runtime%toroidal_harmonic = active_n
        call format_gc_full_fow_frequency_convention(active_n, &
            execution%runtime%frequency_convention)
        execution%runtime%quadrature_convergence_certified = &
            execution%integral%quadrature%converged
        execution%runtime%harmonic_batch_certified = &
            valid_requested_harmonic_batch(poloidal_harmonics, execution%harmonics)
        execution%runtime%class_reconstruction_certified = &
            class_reconstruction_matches(execution%integral, options%force_layout)
    end subroutine record_observed_execution

    logical function valid_requested_harmonic_batch(requested, executed)
        integer, intent(in) :: requested(:)
        type(gc_model2_harmonic_execution_t), intent(in) :: executed(:)
        integer :: i, j

        valid_requested_harmonic_batch = size(requested) > 0 .and. &
            size(executed) == size(requested)
        if (.not. valid_requested_harmonic_batch) return
        do i = 1, size(requested)
            if (count(requested == requested(i)) /= 1) then
                valid_requested_harmonic_batch = .false.
                return
            end if
            j = 1
            do while (j <= size(executed))
                if (executed(j)%poloidal_harmonic == requested(i)) exit
                j = j + 1
            end do
            if (j > size(executed)) then
                valid_requested_harmonic_batch = .false.
                return
            end if
        end do
    end function valid_requested_harmonic_batch

    logical function class_reconstruction_matches(integral, layout)
        type(gc_nonlocal_transport_result_t), intent(in) :: integral
        type(gc_model2_force_layout_t), intent(in) :: layout
        real(dp) :: difference(GC_NONLOCAL_MAX_FORCE_VALUES), scale

        class_reconstruction_matches = .false.
        difference = sum(integral%class_contribution, dim=1) - &
            integral%contribution
        scale = maxval(abs([sum(integral%class_contribution, dim=1), &
            integral%contribution]))
        if (maxval(abs(difference)) > 2.0e-11_dp*max(tiny(1.0_dp), scale)) return
        if (layout%d11_slot < 1 .or. layout%d12_slot < 1 .or. &
                layout%torque_slot < 1) return
        class_reconstruction_matches = .true.
    end function class_reconstruction_matches

    subroutine accept_force_slots(execution, layout, status)
        type(gc_model2_transport_execution_t), intent(inout) :: execution
        type(gc_model2_force_layout_t), intent(in) :: layout
        integer, intent(out) :: status

        execution%d11 = execution%integral%contribution(layout%d11_slot)
        execution%d12 = execution%integral%contribution(layout%d12_slot)
        execution%torque = execution%integral%contribution(layout%torque_slot)
        if (.not. all(ieee_is_finite([execution%d11, execution%d12, &
            execution%torque]))) then
            execution%runtime%nonlocal_transport_certified = .false.
            execution%force_slots_accepted = .false.
            status = GC_MODEL2_DISPATCH_NONFINITE
            return
        end if
        execution%force_slots_accepted = .true.
        status = GC_MODEL2_DISPATCH_SUCCESS
    end subroutine accept_force_slots

    subroutine finalize_gc_model2_transport_execution(execution, backend, observed, &
            status)
        type(gc_model2_transport_execution_t), intent(inout) :: execution
        type(gc_model2_backend_evidence_t), intent(in) :: backend
        type(gc_model2_observed_evidence_t), intent(in) :: observed
        integer, intent(out) :: status

        integer :: local_status

        status = GC_MODEL2_DISPATCH_NOT_CERTIFIED
        local_status = validate_backend(backend, .true.)
        if (local_status /= GC_MODEL2_DISPATCH_SUCCESS) then
            status = local_status
            return
        end if
        if (execution%integral_status /= GC_NONLOCAL_SUCCESS) return
        if (.not. execution%integral%certified) return
        if (.not. execution%force_slots_accepted) return
        if (execution%certified_nodes /= execution%attempted_nodes) return
        if (execution%unresolved_nodes /= 0) return
        if (.not. valid_observed_evidence(observed)) return

        call build_observed_coverage(execution, observed)
        execution%runtime%runtime_execution_complete = .true.
        execution%runtime%orbit_backend_certified = &
            backend%field_certified .and. backend%profile_certified .and. &
            backend%perturbation_certified .and. backend%topology_certified
        execution%runtime%wall_certified = backend%wall_certified
        execution%runtime%canonical_measure_certified = &
            backend%canonical_measure_certified
        execution%runtime%component_identity_certified = &
            backend%component_identity_certified
        execution%runtime%nonlocal_transport_certified = .true.
        execution%certified = .true.
        status = GC_MODEL2_DISPATCH_SUCCESS
    end subroutine finalize_gc_model2_transport_execution

    logical function valid_observed_evidence(observed)
        type(gc_model2_observed_evidence_t), intent(in) :: observed
        integer :: return_total, invariant_total, wall_total

        valid_observed_evidence = .false.
        if (.not. observed%return_accounting_complete) return
        if (observed%physical_return_attempts <= 0) return
        if (observed%topology_certification_attempts <= 0) return
        if (observed%topology_certification_successes /= &
            observed%topology_certification_attempts) return
        if (any([observed%invariant_successes, observed%invariant_failures, &
            observed%return_successes, observed%return_no_return, &
            observed%return_radial_domain, observed%return_wall_loss, &
            observed%return_errors, observed%wall_not_hit, observed%wall_loss, &
            observed%wall_errors, observed%topology_certification_attempts, &
            observed%topology_certification_successes] < 0)) return
        invariant_total = observed%invariant_successes + &
            observed%invariant_failures
        return_total = observed%return_successes + observed%return_no_return + &
            observed%return_radial_domain + observed%return_wall_loss + &
            observed%return_errors
        wall_total = observed%wall_not_hit + observed%wall_loss + &
            observed%wall_errors
        if (invariant_total /= observed%physical_return_attempts) return
        if (return_total /= observed%physical_return_attempts) return
        if (wall_total /= observed%physical_return_attempts) return
        valid_observed_evidence = .true.
    end function valid_observed_evidence

    subroutine build_observed_coverage(execution, observed)
        type(gc_model2_transport_execution_t), intent(inout) :: execution
        type(gc_model2_observed_evidence_t), intent(in) :: observed

        write (execution%invariant_status_coverage, '(A,I0,A,I0)') &
            'success:', observed%invariant_successes, ',failure:', &
            observed%invariant_failures
        write (execution%return_status_coverage, &
            '(A,I0,A,I0,A,I0,A,I0,A,I0)') 'success:', observed%return_successes, &
            ',no_return:', observed%return_no_return, ',radial_domain:', &
            observed%return_radial_domain, ',wall_loss:', observed%return_wall_loss, &
            ',error:', observed%return_errors
        write (execution%wall_status_coverage, '(A,I0,A,I0,A,I0)') 'not_hit:', &
            observed%wall_not_hit, ',wall_loss:', observed%wall_loss, ',error:', &
            observed%wall_errors
    end subroutine build_observed_coverage

    integer function validate_backend(backend, require_execution_certificates) result(status)
        type(gc_model2_backend_evidence_t), intent(in) :: backend
        logical, intent(in) :: require_execution_certificates
        logical :: wall_exists
        integer :: io_status

        status = GC_MODEL2_DISPATCH_BACKEND_UNCERTIFIED
        if (.not. backend%factory_available) then
            status = GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE
            return
        end if
        if (trim(backend%backend) /= 'cylindrical_full_fow') return
        if (trim(backend%coordinates) /= 'R,Z,phi') return
        if (.not. backend%field_certified) return
        if (.not. backend%profile_certified) return
        if (.not. backend%perturbation_certified) return
        if (.not. backend%perturbation_provenance_certified) return
        if (trim(backend%perturbation_amplitude_convention) /= &
                GC_FULL_FOW_REAL_FIELD_AMPLITUDE_CONVENTION) return
        if (len_trim(backend%perturbation_input_path) == 0) return
        if (trim(backend%phase_space_bound_method) /= GC_FULL_FOW_BOUND_METHOD) return
        if (.not. backend%wall_certified) return
        if (require_execution_certificates) then
            if (.not. backend%topology_certified) return
            if (.not. backend%canonical_measure_certified) return
            if (.not. backend%component_identity_certified) return
        end if
        if (backend%legacy_backend_entries /= 0) return
        if (backend%chart_fallback_entries /= 0) return
        if (len_trim(backend%wall_actual_path) == 0) return
        if (trim(backend%wall_units) /= 'm' .and. &
            trim(backend%wall_units) /= 'cm') return
        if (.not. valid_sha256(backend%wall_sha256)) return
        inquire (file=trim(backend%wall_actual_path), exist=wall_exists, &
            iostat=io_status)
        if (io_status /= 0 .or. .not. wall_exists) return
        status = GC_MODEL2_DISPATCH_SUCCESS
    end function validate_backend

    integer function validate_force_layout(options) result(status)
        type(gc_model2_transport_options_t), intent(in) :: options
        integer :: largest_slot
        type(gc_nonlocal_resonance_options_t) :: resonance_options

        status = GC_MODEL2_DISPATCH_FORCE_LAYOUT_INVALID
        resonance_options = options%integral%resonance_options
        largest_slot = max(options%force_layout%d11_slot, &
            options%force_layout%d12_slot, options%force_layout%torque_slot)
        if (options%force_layout%d11_slot < 1) return
        if (options%force_layout%d12_slot < 1) return
        if (options%force_layout%torque_slot < 1) return
        if (options%force_layout%d11_slot == options%force_layout%d12_slot) return
        if (options%force_layout%d11_slot == options%force_layout%torque_slot) return
        if (options%force_layout%d12_slot == options%force_layout%torque_slot) return
        if (largest_slot > GC_NONLOCAL_MAX_FORCE_VALUES) return
        if (resonance_options%force_count < largest_slot) return
        status = GC_MODEL2_DISPATCH_SUCCESS
    end function validate_force_layout

    subroutine emit_gc_model2_runtime_record(phase, surface_key, surface, &
            execution, status, message, output_path)
        character(len=*), intent(in) :: phase, surface_key
        real(dp), intent(in) :: surface
        type(gc_model2_transport_execution_t), intent(in) :: execution
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=*), intent(in), optional :: output_path

        if (.not. execution%certified) then
            status = GC_MODEL2_DISPATCH_NOT_CERTIFIED
            message = 'uncertified model-2 execution cannot be delivered'
            return
        end if
        if (present(output_path)) then
            call emit_gc_full_fow_runtime_surface_record(phase, surface_key, &
                surface, execution%runtime, execution%invariant_status_coverage, &
                execution%return_status_coverage, execution%wall_status_coverage, &
                status, message, output_path)
        else
            call emit_gc_full_fow_runtime_surface_record(phase, surface_key, &
                surface, execution%runtime, execution%invariant_status_coverage, &
                execution%return_status_coverage, execution%wall_status_coverage, &
                status, message)
        end if
    end subroutine emit_gc_model2_runtime_record

    logical function valid_sha256(value)
        character(len=*), intent(in) :: value
        integer :: i, code

        valid_sha256 = len_trim(value) == 64
        if (.not. valid_sha256) return
        do i = 1, 64
            code = iachar(value(i:i))
            if (.not. ((code >= iachar('0') .and. code <= iachar('9')) .or. &
                (code >= iachar('a') .and. code <= iachar('f')) .or. &
                (code >= iachar('A') .and. code <= iachar('F')))) then
                valid_sha256 = .false.
                return
            end if
        end do
    end function valid_sha256

end module neort_gc_model2_transport_dispatch
