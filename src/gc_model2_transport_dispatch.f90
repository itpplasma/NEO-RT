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
        gc_full_fow_runtime_backend_state_t
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_SUCCESS, gc_nonlocal_resonance_options_t
    use neort_gc_nonlocal_transport_integral, only: &
        integrate_gc_nonlocal_transport
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_MAX_FORCE_VALUES, gc_nonlocal_transport_options_t, &
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
        integer :: legacy_backend_entries = 0
        integer :: chart_fallback_entries = 0
    end type gc_model2_backend_evidence_t

    type, public :: gc_model2_transport_execution_t
        type(gc_nonlocal_transport_result_t) :: integral
        type(gc_full_fow_runtime_backend_state_t) :: runtime
        real(dp) :: d11 = 0.0_dp
        real(dp) :: d12 = 0.0_dp
        real(dp) :: torque = 0.0_dp
        integer :: integral_status = GC_NONLOCAL_SUCCESS
        integer :: attempted_nodes = 0
        integer :: certified_nodes = 0
        integer :: unresolved_nodes = 0
        logical :: certified = .false.
        character(len=256) :: invariant_status_coverage = ''
        character(len=256) :: return_status_coverage = ''
        character(len=256) :: wall_status_coverage = ''
    end type gc_model2_transport_execution_t

    public :: emit_gc_model2_runtime_record
    public :: execute_gc_model2_transport
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

        integer :: local_status

        execution = gc_model2_transport_execution_t()
        status = GC_MODEL2_DISPATCH_BACKEND_UNCERTIFIED
        local_status = validate_backend(backend)
        if (local_status /= GC_MODEL2_DISPATCH_SUCCESS) then
            status = local_status
            return
        end if
        local_status = validate_force_layout(options)
        if (local_status /= GC_MODEL2_DISPATCH_SUCCESS) then
            status = local_status
            return
        end if
        call integrate_gc_nonlocal_transport(provider, harmonic_m, harmonic_n, &
            options%integral, execution%integral, local_status)
        execution%integral_status = local_status
        call record_observed_execution(execution, backend)
        if (local_status /= GC_NONLOCAL_SUCCESS) then
            status = GC_MODEL2_DISPATCH_INTEGRAL_FAILED
            return
        end if
        if (.not. execution%integral%certified) then
            status = GC_MODEL2_DISPATCH_INTEGRAL_FAILED
            return
        end if
        if (execution%certified_nodes /= execution%attempted_nodes .or. &
            execution%unresolved_nodes /= 0) then
            status = GC_MODEL2_DISPATCH_INTEGRAL_FAILED
            return
        end if
        call accept_force_slots(execution, options%force_layout, status)
    end subroutine execute_gc_model2_transport

    subroutine record_observed_execution(execution, backend)
        type(gc_model2_transport_execution_t), intent(inout) :: execution
        type(gc_model2_backend_evidence_t), intent(in) :: backend

        execution%attempted_nodes = execution%integral%weighted_nodes
        execution%certified_nodes = execution%integral%certified_nodes
        execution%unresolved_nodes = execution%integral%unresolved_nodes
        call build_coverage(execution)
        execution%runtime%backend = backend%backend
        execution%runtime%coordinates = backend%coordinates
        execution%runtime%model = 2
        execution%runtime%frequency_model = 2
        execution%runtime%wall_actual_path = backend%wall_actual_path
        execution%runtime%wall_units = backend%wall_units
        execution%runtime%wall_sha256 = backend%wall_sha256
        execution%runtime%runtime_execution_complete = .true.
        execution%runtime%orbit_backend_certified = &
            backend%field_certified .and. backend%topology_certified
        execution%runtime%wall_certified = backend%wall_certified
        execution%runtime%canonical_measure_certified = &
            backend%canonical_measure_certified
        execution%runtime%component_identity_certified = &
            backend%component_identity_certified
        execution%runtime%cylindrical_backend_entries = &
            execution%attempted_nodes
        execution%runtime%legacy_backend_entries = backend%legacy_backend_entries
        execution%runtime%chart_fallback_entries = backend%chart_fallback_entries
    end subroutine record_observed_execution

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
            status = GC_MODEL2_DISPATCH_NONFINITE
            return
        end if
        execution%certified = .true.
        execution%runtime%nonlocal_transport_certified = .true.
        status = GC_MODEL2_DISPATCH_SUCCESS
    end subroutine accept_force_slots

    subroutine build_coverage(execution)
        type(gc_model2_transport_execution_t), intent(inout) :: execution

        integer :: success_count, failure_count

        success_count = 0
        if (allocated(execution%integral%node_status)) then
            success_count = count(execution%integral%node_status == &
                GC_NONLOCAL_SUCCESS)
        end if
        failure_count = max(0, execution%attempted_nodes - success_count)
        write (execution%invariant_status_coverage, '("success:",I0,",failure:",I0)') &
            success_count, failure_count
        write (execution%return_status_coverage, &
            '("success:",I0,",no_return:0,radial_domain:0,wall_loss:0,error:",I0)') &
            success_count, failure_count
        write (execution%wall_status_coverage, &
            '("not_hit:",I0,",wall_loss:0,error:",I0)') success_count, failure_count
    end subroutine build_coverage

    integer function validate_backend(backend) result(status)
        type(gc_model2_backend_evidence_t), intent(in) :: backend
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
        if (.not. backend%topology_certified) return
        if (.not. backend%wall_certified) return
        if (.not. backend%canonical_measure_certified) return
        if (.not. backend%component_identity_certified) return
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
