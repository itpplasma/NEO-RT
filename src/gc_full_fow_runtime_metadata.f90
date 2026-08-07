module neort_gc_full_fow_runtime_metadata
    !! Emit the executable-owned model-2 runtime contract.
    !!
    !! The campaign launcher supplies only the destination and the relative
    !! provenance name through the two documented environment variables.  All
    !! other values are observations of the running backend, or certificates
    !! returned by the running transport.  In particular, this module never
    !! upgrades an orbit-only diagnostic to a nonlocal-transport result.
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private

    character(len=*), parameter, public :: &
        GC_FULL_FOW_METADATA_SCHEMA = 'neort-full-fow-runtime-metadata-v2'
    character(len=*), parameter, public :: &
        GC_FULL_FOW_METADATA_FILENAME = 'FULL_FOW_RUNTIME_METADATA'

    integer, parameter, public :: GC_FULL_FOW_METADATA_SUCCESS = 0
    integer, parameter, public :: GC_FULL_FOW_METADATA_ENVIRONMENT_ERROR = 1
    integer, parameter, public :: GC_FULL_FOW_METADATA_INVALID_INPUT = 2
    integer, parameter, public :: GC_FULL_FOW_METADATA_NOT_RUNTIME = 3
    integer, parameter, public :: GC_FULL_FOW_METADATA_BACKEND_MISMATCH = 4
    integer, parameter, public :: GC_FULL_FOW_METADATA_WALL_MISMATCH = 5
    integer, parameter, public :: GC_FULL_FOW_METADATA_INVALID_HASH = 6
    integer, parameter, public :: GC_FULL_FOW_METADATA_FALLBACK = 7
    integer, parameter, public :: GC_FULL_FOW_METADATA_STATUS_INVALID = 8
    integer, parameter, public :: GC_FULL_FOW_METADATA_POLICY_INVALID = 9
    integer, parameter, public :: GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH = 10
    integer, parameter, public :: GC_FULL_FOW_METADATA_IO_ERROR = 11
    integer, parameter, public :: GC_FULL_FOW_METADATA_TARGET_EXISTS = 12
    integer, parameter, public :: GC_FULL_FOW_METADATA_RENAME_ERROR = 13

    character(len=*), parameter :: RUNTIME_SCHEMA = &
        'neort-full-fow-runtime-metadata-v2'
    character(len=*), parameter :: EXPECTED_BACKEND = 'cylindrical_full_fow'
    character(len=*), parameter :: EXPECTED_COORDINATES = 'R,Z,phi'
    character(len=*), parameter, public :: &
        GC_FULL_FOW_REAL_FIELD_AMPLITUDE_CONVENTION = &
        'real_field_amplitude_one_signed_n'
    character(len=*), parameter, public :: &
        GC_FULL_FOW_CONJUGATE_POLICY = 'conjugate_implicit'
    character(len=*), parameter, public :: &
        GC_FULL_FOW_PREFACTOR_CONVENTION = 'eq17_pi32_over_4_real_field'
    character(len=*), parameter, public :: &
        GC_FULL_FOW_ACTION_CONVENTION = 'J_K=mass*c*mu_phys/abs(q)'
    character(len=*), parameter, public :: &
        GC_FULL_FOW_BOUND_METHOD = 'axisymmetric_geqdsk_fpol_fluxfunction'

    interface
        function c_link(old_path, new_path) bind(C, name='link') result(code)
            import c_char, c_int
            character(kind=c_char), intent(in) :: old_path(*)
            character(kind=c_char), intent(in) :: new_path(*)
            integer(c_int) :: code
        end function c_link

        function c_unlink(path) bind(C, name='unlink') result(code)
            import c_char, c_int
            character(kind=c_char), intent(in) :: path(*)
            integer(c_int) :: code
        end function c_unlink
    end interface

    type, public :: gc_full_fow_runtime_backend_state_t
        !! Values captured after the executable has entered the model-2 path.
        character(len=32) :: backend = ''
        character(len=16) :: coordinates = ''
        integer :: model = 0
        integer :: frequency_model = 0
        character(len=1024) :: wall_actual_path = ''
        character(len=16) :: wall_units = ''
        character(len=64) :: wall_sha256 = ''
        logical :: runtime_execution_complete = .false.
        logical :: orbit_backend_certified = .false.
        logical :: wall_certified = .false.
        logical :: canonical_measure_certified = .false.
        logical :: component_identity_certified = .false.
        logical :: nonlocal_transport_certified = .false.
        integer :: cylindrical_backend_entries = 0
        integer :: legacy_backend_entries = 0
        integer :: chart_fallback_entries = 0
        character(len=64) :: real_field_amplitude_convention = ''
        character(len=32) :: conjugate_policy = ''
        character(len=64) :: prefactor_convention = ''
        character(len=64) :: action_convention = ''
        character(len=64) :: phase_space_bound_method = ''
        character(len=64) :: frequency_convention = ''
        character(len=1024) :: perturbation_input_path = ''
        logical :: perturbation_provenance_certified = .false.
        integer :: quadrature_base_h0_order = 0
        integer :: quadrature_base_jk_order = 0
        integer :: quadrature_refined_h0_order = 0
        integer :: quadrature_refined_jk_order = 0
        real(dp) :: quadrature_relative_tolerance = 0.0_dp
        real(dp) :: quadrature_absolute_tolerance = 0.0_dp
        integer :: poloidal_harmonic_min = 0
        integer :: poloidal_harmonic_max = 0
        integer :: poloidal_harmonic_count = 0
        integer :: executed_harmonic_count = 0
        integer :: toroidal_harmonic = 0
        logical :: quadrature_convergence_certified = .false.
        logical :: harmonic_batch_certified = .false.
        logical :: class_reconstruction_certified = .false.
        logical :: orbit_step_refinement_certified = .false.
        real(dp) :: orbit_base_step = 0.0_dp
        real(dp) :: orbit_refined_step = 0.0_dp
        real(dp) :: orbit_period_refinement_error = 0.0_dp
        real(dp) :: orbit_delta_phi_refinement_error = 0.0_dp
        real(dp) :: orbit_omega_b_refinement_error = 0.0_dp
        real(dp) :: orbit_omega_phi_refinement_error = 0.0_dp
        real(dp) :: orbit_h_m_refinement_error = 0.0_dp
        real(dp) :: orbit_shell_refinement_error = 0.0_dp
    end type gc_full_fow_runtime_backend_state_t

    public :: emit_gc_full_fow_runtime_metadata
    public :: format_gc_full_fow_frequency_convention

contains

    subroutine format_gc_full_fow_frequency_convention(toroidal_harmonic, text)
        integer, intent(in) :: toroidal_harmonic
        character(len=*), intent(out) :: text

        text = ''
        write (text, '("m*omega_b+",I0,"*omega_phi")') toroidal_harmonic
    end subroutine format_gc_full_fow_frequency_convention

    subroutine emit_gc_full_fow_runtime_metadata(phase, lane_kind, state, &
            invariant_status_coverage, return_status_coverage, &
            wall_status_coverage, frequency_phase_policy, diagnostic_count, &
            phase_independent_evidence, status, message)
        character(len=*), intent(in) :: phase, lane_kind
        type(gc_full_fow_runtime_backend_state_t), intent(in) :: state
        character(len=*), intent(in) :: invariant_status_coverage
        character(len=*), intent(in) :: return_status_coverage
        character(len=*), intent(in) :: wall_status_coverage
        character(len=*), intent(in) :: frequency_phase_policy
        integer, intent(in) :: diagnostic_count
        logical, intent(in) :: phase_independent_evidence
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        character(len=4096) :: output_path, temporary_path
        character(len=1024) :: relative_wall_path
        character(len=256) :: invariant_normalized
        character(len=256) :: return_normalized
        character(len=256) :: wall_normalized
        integer :: invariant_counts(5), return_counts(5), wall_counts(5)
        integer :: unit, io_status, close_status, flush_status
        logical :: target_exists, renamed

        call set_result(status, message, GC_FULL_FOW_METADATA_INVALID_INPUT, '')

        call read_runtime_environment('NEORT_FULL_FOW_METADATA_PATH', &
            output_path, status, message)
        if (status /= GC_FULL_FOW_METADATA_SUCCESS) return
        call read_runtime_environment('NEORT_FULL_FOW_WALL_RELATIVE_PATH', &
            relative_wall_path, status, message)
        if (status /= GC_FULL_FOW_METADATA_SUCCESS) return
        call validate_output_path(output_path, status, message)
        if (status /= GC_FULL_FOW_METADATA_SUCCESS) return
        call validate_relative_wall_path(relative_wall_path, status, message)
        if (status /= GC_FULL_FOW_METADATA_SUCCESS) return

        call validate_runtime_state(state, phase, lane_kind, status, message)
        if (status /= GC_FULL_FOW_METADATA_SUCCESS) return
        call validate_phase_policy(frequency_phase_policy, diagnostic_count, &
            phase_independent_evidence, status, message)
        if (status /= GC_FULL_FOW_METADATA_SUCCESS) return
        call parse_status_coverage(invariant_status_coverage, 'invariant', &
            invariant_counts, invariant_normalized, status, message)
        if (status /= GC_FULL_FOW_METADATA_SUCCESS) return
        call parse_status_coverage(return_status_coverage, 'return', &
            return_counts, return_normalized, status, message)
        if (status /= GC_FULL_FOW_METADATA_SUCCESS) return
        call parse_status_coverage(wall_status_coverage, 'wall', wall_counts, &
            wall_normalized, status, message)
        if (status /= GC_FULL_FOW_METADATA_SUCCESS) return

        ! The executable must not silently replace an existing accepted sidecar.
        inquire (file=trim(output_path), exist=target_exists, iostat=io_status)
        if (io_status /= 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_IO_ERROR, &
                'could not inspect metadata destination')
            return
        end if
        if (target_exists) then
            call set_result(status, message, GC_FULL_FOW_METADATA_TARGET_EXISTS, &
                'metadata destination already exists')
            return
        end if

        call make_temporary_path(output_path, temporary_path, status, message)
        if (status /= GC_FULL_FOW_METADATA_SUCCESS) return
        open (newunit=unit, file=trim(temporary_path), status='new', &
            action='write', form='formatted', iostat=io_status)
        if (io_status /= 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_IO_ERROR, &
                'could not create temporary metadata file')
            return
        end if

        call write_pair(unit, 'schema', RUNTIME_SCHEMA, io_status)
        if (io_status == 0) call write_pair(unit, 'metadata_origin', 'runtime', io_status)
        if (io_status == 0) call write_pair(unit, 'runtime_emitter', 'neo_rt', io_status)
        if (io_status == 0) call write_pair(unit, 'phase', trim(phase), io_status)
        if (io_status == 0) call write_pair(unit, 'lane_kind', trim(lane_kind), io_status)
        if (io_status == 0) call write_integer_pair(unit, 'model', state%model, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'frequency_model', &
            state%frequency_model, io_status)
        if (io_status == 0) call write_pair(unit, 'backend', trim(state%backend), io_status)
        if (io_status == 0) call write_pair(unit, 'coordinates', trim(state%coordinates), io_status)
        if (io_status == 0) call write_pair(unit, 'wall_path', trim(relative_wall_path), io_status)
        if (io_status == 0) call write_pair(unit, 'wall_sha256', trim(state%wall_sha256), io_status)
        if (io_status == 0) call write_pair(unit, 'wall_units', trim(state%wall_units), io_status)
        if (io_status == 0) call write_logical_pair(unit, 'wall_certified', &
            state%wall_certified, io_status)
        if (io_status == 0) call write_logical_pair(unit, 'canonical_measure_certified', &
            state%canonical_measure_certified, io_status)
        if (io_status == 0) call write_logical_pair(unit, 'component_identity_certified', &
            state%component_identity_certified, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'cylindrical_backend_entries', &
            state%cylindrical_backend_entries, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'legacy_backend_entries', &
            state%legacy_backend_entries, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'chart_fallback_entries', &
            state%chart_fallback_entries, io_status)
        if (io_status == 0) call write_pair(unit, 'real_field_amplitude_convention', &
            trim(state%real_field_amplitude_convention), io_status)
        if (io_status == 0) call write_pair(unit, 'conjugate_policy', &
            trim(state%conjugate_policy), io_status)
        if (io_status == 0) call write_pair(unit, 'prefactor_convention', &
            trim(state%prefactor_convention), io_status)
        if (io_status == 0) call write_pair(unit, 'action_convention', &
            trim(state%action_convention), io_status)
        if (io_status == 0) call write_pair(unit, 'phase_space_bound_method', &
            trim(state%phase_space_bound_method), io_status)
        if (io_status == 0) call write_pair(unit, 'frequency_convention', &
            trim(state%frequency_convention), io_status)
        if (io_status == 0) call write_pair(unit, 'perturbation_input_path', &
            trim(state%perturbation_input_path), io_status)
        if (io_status == 0) call write_logical_pair(unit, &
            'perturbation_provenance_certified', &
            state%perturbation_provenance_certified, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'quadrature_base_h0_order', &
            state%quadrature_base_h0_order, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'quadrature_base_jk_order', &
            state%quadrature_base_jk_order, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'quadrature_refined_h0_order', &
            state%quadrature_refined_h0_order, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'quadrature_refined_jk_order', &
            state%quadrature_refined_jk_order, io_status)
        if (io_status == 0) call write_real_pair(unit, 'quadrature_relative_tolerance', &
            state%quadrature_relative_tolerance, io_status)
        if (io_status == 0) call write_real_pair(unit, 'quadrature_absolute_tolerance', &
            state%quadrature_absolute_tolerance, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'poloidal_harmonic_min', &
            state%poloidal_harmonic_min, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'poloidal_harmonic_max', &
            state%poloidal_harmonic_max, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'poloidal_harmonic_count', &
            state%poloidal_harmonic_count, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'executed_harmonic_count', &
            state%executed_harmonic_count, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'toroidal_harmonic', &
            state%toroidal_harmonic, io_status)
        if (io_status == 0) call write_logical_pair(unit, &
            'quadrature_convergence_certified', &
            state%quadrature_convergence_certified, io_status)
        if (io_status == 0) call write_logical_pair(unit, 'harmonic_batch_certified', &
            state%harmonic_batch_certified, io_status)
        if (io_status == 0) call write_logical_pair(unit, &
            'class_reconstruction_certified', state%class_reconstruction_certified, &
            io_status)
        if (io_status == 0) call write_logical_pair(unit, &
            'orbit_step_refinement_certified', state%orbit_step_refinement_certified, &
            io_status)
        if (io_status == 0) call write_real_pair(unit, 'orbit_base_step', &
            state%orbit_base_step, io_status)
        if (io_status == 0) call write_real_pair(unit, 'orbit_refined_step', &
            state%orbit_refined_step, io_status)
        if (io_status == 0) call write_real_pair(unit, &
            'orbit_period_refinement_error', state%orbit_period_refinement_error, io_status)
        if (io_status == 0) call write_real_pair(unit, &
            'orbit_delta_phi_refinement_error', state%orbit_delta_phi_refinement_error, io_status)
        if (io_status == 0) call write_real_pair(unit, &
            'orbit_omega_b_refinement_error', state%orbit_omega_b_refinement_error, io_status)
        if (io_status == 0) call write_real_pair(unit, &
            'orbit_omega_phi_refinement_error', state%orbit_omega_phi_refinement_error, io_status)
        if (io_status == 0) call write_real_pair(unit, 'orbit_h_m_refinement_error', &
            state%orbit_h_m_refinement_error, io_status)
        if (io_status == 0) call write_real_pair(unit, 'orbit_shell_refinement_error', &
            state%orbit_shell_refinement_error, io_status)
        if (io_status == 0) call write_pair(unit, 'invariant_status_coverage', &
            trim(invariant_normalized), io_status)
        if (io_status == 0) call write_pair(unit, 'return_status_coverage', &
            trim(return_normalized), io_status)
        if (io_status == 0) call write_pair(unit, 'wall_status_coverage', &
            trim(wall_normalized), io_status)
        if (io_status == 0) call write_pair(unit, 'frequency_phase_policy', &
            trim(frequency_phase_policy), io_status)
        if (io_status == 0) call write_integer_pair(unit, 'diagnostic_count', &
            diagnostic_count, io_status)
        if (io_status == 0) call write_logical_pair(unit, 'phase_independent_evidence', &
            phase_independent_evidence, io_status)

        flush_status = 0
        if (io_status == 0) flush (unit, iostat=flush_status)
        close_status = 0
        close (unit, iostat=close_status)
        if (io_status /= 0 .or. flush_status /= 0 .or. close_status /= 0) then
            call unlink_file(temporary_path)
            call set_result(status, message, GC_FULL_FOW_METADATA_IO_ERROR, &
                'write of temporary metadata file failed')
            return
        end if

        call publish_file(temporary_path, output_path, renamed)
        if (.not. renamed) then
            call unlink_file(temporary_path)
            call set_result(status, message, GC_FULL_FOW_METADATA_RENAME_ERROR, &
                'atomic metadata publish failed')
            return
        end if
        call unlink_file(temporary_path)

        call set_result(status, message, GC_FULL_FOW_METADATA_SUCCESS, 'ok')
    end subroutine emit_gc_full_fow_runtime_metadata

    subroutine validate_runtime_state(state, phase, lane_kind, status, message)
        type(gc_full_fow_runtime_backend_state_t), intent(in) :: state
        character(len=*), intent(in) :: phase, lane_kind
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        logical :: wall_exists
        integer :: io_status

        call set_result(status, message, GC_FULL_FOW_METADATA_SUCCESS, 'ok')
        if (.not. state%runtime_execution_complete) then
            call set_result(status, message, GC_FULL_FOW_METADATA_NOT_RUNTIME, &
                'runtime metadata requires completed executable state')
            return
        end if
        if (trim(phase) /= 'phiI000' .and. trim(phase) /= 'phiI010') then
            call set_result(status, message, GC_FULL_FOW_METADATA_INVALID_INPUT, &
                'phase is not an actual ITER phase')
            return
        end if
        if (trim(lane_kind) /= 'frequency' .and. trim(lane_kind) /= 'torque') then
            call set_result(status, message, GC_FULL_FOW_METADATA_INVALID_INPUT, &
                'lane kind is not frequency or torque')
            return
        end if
        if (state%model /= 2 .or. state%frequency_model /= 2) then
            call set_result(status, message, GC_FULL_FOW_METADATA_BACKEND_MISMATCH, &
                'runtime state is not model-2')
            return
        end if
        if (trim(state%backend) /= EXPECTED_BACKEND) then
            call set_result(status, message, GC_FULL_FOW_METADATA_BACKEND_MISMATCH, &
                'runtime backend is not cylindrical full-FOW')
            return
        end if
        if (trim(state%coordinates) /= EXPECTED_COORDINATES) then
            call set_result(status, message, GC_FULL_FOW_METADATA_BACKEND_MISMATCH, &
                'runtime coordinates are not R,Z,phi')
            return
        end if
        if (.not. state%orbit_backend_certified) then
            call set_result(status, message, GC_FULL_FOW_METADATA_BACKEND_MISMATCH, &
                'runtime orbit/backend certificate is absent')
            return
        end if
        if (state%cylindrical_backend_entries <= 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_BACKEND_MISMATCH, &
                'cylindrical backend was not entered')
            return
        end if
        if (state%legacy_backend_entries < 0 .or. &
            state%chart_fallback_entries < 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_INVALID_INPUT, &
                'backend counters cannot be negative')
            return
        end if
        if (state%legacy_backend_entries /= 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_FALLBACK, &
                'legacy backend fallback was entered')
            return
        end if
        if (state%chart_fallback_entries /= 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_FALLBACK, &
                'chart fallback was entered')
            return
        end if
        if (len_trim(state%wall_actual_path) == 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_WALL_MISMATCH, &
                'runtime wall path is absent')
            return
        end if
        if (has_control_character(state%wall_actual_path)) then
            call set_result(status, message, GC_FULL_FOW_METADATA_WALL_MISMATCH, &
                'runtime wall path contains a control character')
            return
        end if
        inquire (file=trim(state%wall_actual_path), exist=wall_exists, &
            iostat=io_status)
        if (io_status /= 0 .or. .not. wall_exists) then
            call set_result(status, message, GC_FULL_FOW_METADATA_WALL_MISMATCH, &
                'runtime wall path is not readable')
            return
        end if
        if (trim(state%wall_units) /= 'm' .and. trim(state%wall_units) /= 'cm') then
            call set_result(status, message, GC_FULL_FOW_METADATA_WALL_MISMATCH, &
                'runtime wall units are not m or cm')
            return
        end if
        if (.not. valid_sha256(state%wall_sha256)) then
            call set_result(status, message, GC_FULL_FOW_METADATA_INVALID_HASH, &
                'runtime wall hash is not a complete SHA-256')
            return
        end if
        if (.not. state%wall_certified) then
            call set_result(status, message, GC_FULL_FOW_METADATA_WALL_MISMATCH, &
                'runtime wall certificate is absent')
            return
        end if
        if (.not. state%canonical_measure_certified) then
            call set_result(status, message, &
                GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH, &
                'runtime canonical-measure certificate is absent')
            return
        end if
        if (.not. state%component_identity_certified) then
            call set_result(status, message, &
                GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH, &
                'runtime component-identity certificate is absent')
            return
        end if
        if (trim(lane_kind) == 'torque') then
            if (.not. state%nonlocal_transport_certified) then
                call set_result(status, message, &
                    GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH, &
                    'torque lane lacks nonlocal transport certification')
                return
            end if
        end if
        if (state%nonlocal_transport_certified) then
            call validate_transport_provenance(state, status, message)
        end if
    end subroutine validate_runtime_state

    subroutine validate_transport_provenance(state, status, message)
        type(gc_full_fow_runtime_backend_state_t), intent(in) :: state
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=64) :: expected_frequency

        call set_result(status, message, GC_FULL_FOW_METADATA_SUCCESS, 'ok')
        call format_gc_full_fow_frequency_convention(state%toroidal_harmonic, &
            expected_frequency)
        if (trim(state%real_field_amplitude_convention) /= &
                GC_FULL_FOW_REAL_FIELD_AMPLITUDE_CONVENTION .or. &
                trim(state%conjugate_policy) /= GC_FULL_FOW_CONJUGATE_POLICY .or. &
                trim(state%prefactor_convention) /= GC_FULL_FOW_PREFACTOR_CONVENTION .or. &
                trim(state%action_convention) /= GC_FULL_FOW_ACTION_CONVENTION .or. &
                trim(state%phase_space_bound_method) /= GC_FULL_FOW_BOUND_METHOD .or. &
                trim(state%frequency_convention) /= trim(expected_frequency) .or. &
                .not. state%perturbation_provenance_certified .or. &
                len_trim(state%perturbation_input_path) == 0) then
            call set_result(status, message, &
                GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH, &
                'transport Fourier or prefactor convention is not exact')
            return
        end if
        if (state%quadrature_base_h0_order < 2 .or. &
                state%quadrature_base_jk_order < 2 .or. &
                state%quadrature_refined_h0_order /= &
                2*state%quadrature_base_h0_order .or. &
                state%quadrature_refined_jk_order /= &
                2*state%quadrature_base_jk_order) then
            call set_result(status, message, &
                GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH, &
                'transport quadrature orders are not an N-vs-2N pair')
            return
        end if
        if (.not. ieee_is_finite(state%quadrature_relative_tolerance) .or. &
                .not. ieee_is_finite(state%quadrature_absolute_tolerance) .or. &
                state%quadrature_relative_tolerance <= 0.0_dp .or. &
                state%quadrature_absolute_tolerance <= 0.0_dp) then
            call set_result(status, message, &
                GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH, &
                'transport quadrature tolerances are invalid')
            return
        end if
        if (state%poloidal_harmonic_count <= 0 .or. &
                state%executed_harmonic_count /= state%poloidal_harmonic_count .or. &
                state%poloidal_harmonic_min > state%poloidal_harmonic_max .or. &
                state%toroidal_harmonic == 0 .or. &
                .not. state%quadrature_convergence_certified .or. &
                .not. state%harmonic_batch_certified .or. &
                .not. state%class_reconstruction_certified .or. &
                .not. state%orbit_step_refinement_certified .or. &
                .not. ieee_is_finite(state%orbit_base_step) .or. &
                .not. ieee_is_finite(state%orbit_refined_step) .or. &
                state%orbit_base_step <= 0.0_dp .or. &
                state%orbit_refined_step <= 0.0_dp .or. &
                state%orbit_refined_step >= state%orbit_base_step .or. &
                .not. all(ieee_is_finite([state%orbit_period_refinement_error, &
                    state%orbit_delta_phi_refinement_error, &
                    state%orbit_omega_b_refinement_error, &
                    state%orbit_omega_phi_refinement_error, &
                    state%orbit_h_m_refinement_error, &
                    state%orbit_shell_refinement_error])) .or. &
                any([state%orbit_period_refinement_error, &
                    state%orbit_delta_phi_refinement_error, &
                    state%orbit_omega_b_refinement_error, &
                    state%orbit_omega_phi_refinement_error, &
                    state%orbit_h_m_refinement_error, &
                    state%orbit_shell_refinement_error] < 0.0_dp)) then
            call set_result(status, message, &
                GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH, &
                'transport batch, convergence, or class reconstruction is uncertified')
        end if
    end subroutine validate_transport_provenance

    subroutine validate_phase_policy(policy, diagnostic_count, evidence, status, message)
        character(len=*), intent(in) :: policy
        integer, intent(in) :: diagnostic_count
        logical, intent(in) :: evidence
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        call set_result(status, message, GC_FULL_FOW_METADATA_SUCCESS, 'ok')
        if (trim(policy) /= 'phase_independent' .and. &
            trim(policy) /= 'four_phase_model_diagnostics') then
            call set_result(status, message, GC_FULL_FOW_METADATA_POLICY_INVALID, &
                'unknown frequency phase policy')
            return
        end if
        if (trim(policy) == 'phase_independent') then
            if (diagnostic_count /= 2 .or. .not. evidence) then
                call set_result(status, message, GC_FULL_FOW_METADATA_POLICY_INVALID, &
                    'phase-independent policy requires two diagnostics and evidence')
                return
            end if
        else
            if (diagnostic_count /= 4 .or. evidence) then
                call set_result(status, message, GC_FULL_FOW_METADATA_POLICY_INVALID, &
                    'four-phase policy requires four diagnostics and no independent evidence')
                return
            end if
        end if
    end subroutine validate_phase_policy

    subroutine parse_status_coverage(text, kind, counts, normalized, status, message)
        character(len=*), intent(in) :: text, kind
        integer, intent(out) :: counts(:)
        character(len=*), intent(out) :: normalized
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        character(len=16) :: names(5)
        character(len=512) :: item, item_name, count_text
        integer :: required_count, seen(5), expected_index
        integer :: cursor, limit, separator, comma_position, semicolon_position
        integer :: item_end, item_length, colon, ios, value
        logical :: has_separator

        counts = 0
        normalized = ''
        names = ''
        required_count = 0
        select case (trim(kind))
        case ('invariant')
            required_count = 2
            names(1) = 'success'
            names(2) = 'failure'
        case ('return')
            required_count = 5
            names(1) = 'success'
            names(2) = 'no_return'
            names(3) = 'radial_domain'
            names(4) = 'wall_loss'
            names(5) = 'error'
        case ('wall')
            required_count = 3
            names(1) = 'not_hit'
            names(2) = 'wall_loss'
            names(3) = 'error'
        case default
            call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                'unknown status coverage kind')
            return
        end select

        seen = 0
        call set_result(status, message, GC_FULL_FOW_METADATA_SUCCESS, 'ok')
        limit = len_trim(text)
        if (limit == 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                'empty status coverage')
            return
        end if
        cursor = 1
        do
            if (cursor > limit) then
                call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                    'status coverage ends with a separator')
                return
            end if
            comma_position = index(text(cursor:limit), ',')
            semicolon_position = index(text(cursor:limit), ';')
            if (comma_position == 0) then
                separator = semicolon_position
            else if (semicolon_position == 0) then
                separator = comma_position
            else
                separator = min(comma_position, semicolon_position)
            end if
            has_separator = separator /= 0
            if (has_separator) then
                item_end = cursor + separator - 2
            else
                item_end = limit
            end if
            if (item_end < cursor) then
                call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                    'empty status coverage item')
                return
            end if
            item = trim(text(cursor:item_end))
            item_length = len_trim(item)
            if (item_length == 0) then
                call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                    'empty status coverage item')
                return
            end if
            colon = index(item(:item_length), ':')
            if (colon <= 1 .or. colon >= item_length) then
                call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                    'malformed status coverage item')
                return
            end if
            if (index(item(colon + 1:item_length), ':') /= 0) then
                call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                    'status coverage item has multiple counts')
                return
            end if
            item_name = trim(item(:colon - 1))
            count_text = trim(item(colon + 1:item_length))
            expected_index = 0
            if (required_count > 0) then
                expected_index = find_status_name(item_name, names, required_count)
            end if
            if (expected_index == 0) then
                call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                    'unknown runtime status')
                return
            end if
            if (seen(expected_index) /= 0) then
                call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                    'duplicate runtime status')
                return
            end if
            if (.not. valid_integer_text(count_text)) then
                call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                    'non-integer runtime status count')
                return
            end if
            read (count_text, *, iostat=ios) value
            if (ios /= 0 .or. value < 0) then
                call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                    'invalid runtime status count')
                return
            end if
            counts(expected_index) = value
            seen(expected_index) = 1
            if (.not. has_separator) exit
            cursor = cursor + separator
        end do

        if (any(seen(1:required_count) == 0)) then
            call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                'runtime status coverage is incomplete')
            return
        end if
        if (sum(counts(1:required_count)) <= 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                'runtime status coverage is empty')
            return
        end if
        call format_status_coverage(kind, counts, normalized, status, message)
    end subroutine parse_status_coverage

    integer function find_status_name(name, names, count) result(index_value)
        character(len=*), intent(in) :: name
        character(len=*), intent(in) :: names(:)
        integer, intent(in) :: count
        integer :: i

        index_value = 0
        do i = 1, count
            if (trim(name) == trim(names(i))) then
                index_value = i
                return
            end if
        end do
    end function find_status_name

    subroutine format_status_coverage(kind, counts, normalized, status, message)
        character(len=*), intent(in) :: kind
        integer, intent(in) :: counts(:)
        character(len=*), intent(out) :: normalized
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer :: ios

        normalized = ''
        status = GC_FULL_FOW_METADATA_SUCCESS
        message = 'ok'
        select case (trim(kind))
        case ('invariant')
            write (normalized, '(A,I0,",",A,I0)', iostat=ios) &
                'success:', counts(1), 'failure:', counts(2)
        case ('return')
            write (normalized, '(A,I0,",",A,I0,",",A,I0,",",A,I0,",",A,I0)', &
                iostat=ios) 'success:', counts(1), 'no_return:', counts(2), &
                'radial_domain:', counts(3), 'wall_loss:', counts(4), &
                'error:', counts(5)
        case ('wall')
            write (normalized, '(A,I0,",",A,I0,",",A,I0)', iostat=ios) &
                'not_hit:', counts(1), 'wall_loss:', counts(2), 'error:', counts(3)
        case default
            ios = 1
        end select
        if (ios /= 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_STATUS_INVALID, &
                'could not normalize status coverage')
        end if
    end subroutine format_status_coverage

    logical function valid_integer_text(text) result(valid)
        character(len=*), intent(in) :: text
        integer :: i, first, length_value
        character :: symbol

        valid = .false.
        length_value = len_trim(text)
        if (length_value == 0) return
        first = 1
        symbol = text(first:first)
        if (symbol == '+' .or. symbol == '-') first = first + 1
        if (first > length_value) return
        do i = first, length_value
            symbol = text(i:i)
            if (symbol < '0' .or. symbol > '9') return
        end do
        valid = .true.
    end function valid_integer_text

    logical function valid_sha256(hash) result(valid)
        character(len=*), intent(in) :: hash
        integer :: i, code
        character :: symbol

        valid = .false.
        if (len_trim(hash) /= 64) return
        do i = 1, 64
            symbol = hash(i:i)
            code = iachar(symbol)
            if (.not. ((code >= iachar('0') .and. code <= iachar('9')) .or. &
                (code >= iachar('a') .and. code <= iachar('f')) .or. &
                (code >= iachar('A') .and. code <= iachar('F')))) return
        end do
        valid = .true.
    end function valid_sha256

    logical function has_control_character(text) result(found)
        character(len=*), intent(in) :: text
        integer :: i, code

        found = .false.
        do i = 1, len_trim(text)
            code = iachar(text(i:i))
            if (code < 32 .or. code == 127) then
                found = .true.
                return
            end if
        end do
    end function has_control_character

    subroutine read_runtime_environment(name, value, status, message)
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: value
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=4096) :: buffer
        integer :: actual_length, environment_status

        value = ''
        call set_result(status, message, GC_FULL_FOW_METADATA_ENVIRONMENT_ERROR, &
            'required runtime environment variable is missing')
        buffer = ''
        call get_environment_variable(trim(name), value=buffer, &
            length=actual_length, status=environment_status)
        if (environment_status /= 0) return
        if (actual_length <= 0) return
        if (actual_length > len(value)) then
            call set_result(status, message, GC_FULL_FOW_METADATA_ENVIRONMENT_ERROR, &
                'runtime environment variable is truncated')
            return
        end if
        value = buffer(:actual_length)
        status = GC_FULL_FOW_METADATA_SUCCESS
        message = 'ok'
    end subroutine read_runtime_environment

    subroutine validate_output_path(path, status, message)
        character(len=*), intent(in) :: path
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        call set_result(status, message, GC_FULL_FOW_METADATA_SUCCESS, 'ok')
        if (len_trim(path) == 0 .or. has_control_character(path)) then
            call set_result(status, message, GC_FULL_FOW_METADATA_INVALID_INPUT, &
                'metadata output path is empty or unsafe')
            return
        end if
        if (len_trim(path) > 4000) then
            call set_result(status, message, GC_FULL_FOW_METADATA_INVALID_INPUT, &
                'metadata output path is too long')
        end if
    end subroutine validate_output_path

    subroutine validate_relative_wall_path(path, status, message)
        character(len=*), intent(in) :: path
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        call set_result(status, message, GC_FULL_FOW_METADATA_SUCCESS, 'ok')
        if (len_trim(path) == 0 .or. has_control_character(path)) then
            call set_result(status, message, GC_FULL_FOW_METADATA_WALL_MISMATCH, &
                'relative wall provenance path is empty or unsafe')
            return
        end if
        if (path(1:1) == '/') then
            call set_result(status, message, GC_FULL_FOW_METADATA_WALL_MISMATCH, &
                'wall provenance path must be relative')
            return
        end if
        if (index('/'//trim(path)//'/', '/../') /= 0) then
            call set_result(status, message, GC_FULL_FOW_METADATA_WALL_MISMATCH, &
                'wall provenance path escapes its campaign root')
            return
        end if
    end subroutine validate_relative_wall_path

    subroutine make_temporary_path(output_path, temporary_path, status, message)
        character(len=*), intent(in) :: output_path
        character(len=*), intent(out) :: temporary_path
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer(int64) :: clock_count
        character(len=64) :: suffix

        temporary_path = ''
        call set_result(status, message, GC_FULL_FOW_METADATA_INVALID_INPUT, '')
        if (len_trim(output_path) + 32 > len(temporary_path)) then
            call set_result(status, message, GC_FULL_FOW_METADATA_INVALID_INPUT, &
                'temporary metadata path is too long')
            return
        end if
        call system_clock(count=clock_count)
        write (suffix, '(I0)') clock_count
        temporary_path = trim(output_path)//'.tmp.'//trim(suffix)
        status = GC_FULL_FOW_METADATA_SUCCESS
        message = 'ok'
    end subroutine make_temporary_path

    subroutine write_pair(unit, key, value, status)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: key, value
        integer, intent(out) :: status

        write (unit, '(A,"=",A)', iostat=status) trim(key), trim(value)
    end subroutine write_pair

    subroutine write_integer_pair(unit, key, value, status)
        integer, intent(in) :: unit, value
        character(len=*), intent(in) :: key
        integer, intent(out) :: status
        character(len=32) :: text

        write (text, '(I0)') value
        call write_pair(unit, key, text, status)
    end subroutine write_integer_pair

    subroutine write_real_pair(unit, key, value, status)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: key
        real(dp), intent(in) :: value
        integer, intent(out) :: status
        character(len=64) :: text

        if (.not. ieee_is_finite(value)) then
            status = 1
            return
        end if
        write (text, '(ES24.16E3)') value
        call write_pair(unit, key, text, status)
    end subroutine write_real_pair

    subroutine write_logical_pair(unit, key, value, status)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: key
        logical, intent(in) :: value
        integer, intent(out) :: status
        character(len=5) :: text

        if (value) then
            text = 'true'
        else
            text = 'false'
        end if
        call write_pair(unit, key, text, status)
    end subroutine write_logical_pair

    subroutine publish_file(old_path, new_path, success)
        character(len=*), intent(in) :: old_path, new_path
        logical, intent(out) :: success
        character(kind=c_char), allocatable :: old_c(:), new_c(:)
        integer(c_int) :: code

        call make_c_string(old_path, old_c)
        call make_c_string(new_path, new_c)
        ! link(2) creates the destination atomically and refuses to replace an
        ! existing accepted sidecar.  The temporary name is unlinked only
        ! after the destination link succeeds.
        code = c_link(old_c, new_c)
        success = code == 0_c_int
        deallocate (old_c, new_c)
    end subroutine publish_file

    subroutine unlink_file(path)
        character(len=*), intent(in) :: path
        character(kind=c_char), allocatable :: path_c(:)
        integer(c_int) :: code

        call make_c_string(path, path_c)
        code = c_unlink(path_c)
        deallocate (path_c)
    end subroutine unlink_file

    subroutine make_c_string(path, c_path)
        character(len=*), intent(in) :: path
        character(kind=c_char), allocatable, intent(out) :: c_path(:)
        integer :: i, length_path

        length_path = len_trim(path)
        allocate (c_path(length_path + 1))
        do i = 1, length_path
            c_path(i) = path(i:i)
        end do
        c_path(length_path + 1) = c_null_char
    end subroutine make_c_string

    subroutine set_result(status, message, new_status, new_message)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer, intent(in) :: new_status
        character(len=*), intent(in) :: new_message

        status = new_status
        message = new_message
    end subroutine set_result

end module neort_gc_full_fow_runtime_metadata
