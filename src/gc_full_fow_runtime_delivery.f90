module neort_gc_full_fow_runtime_delivery
    !! Deliver one immutable executable-owned record for one model-2 surface.
    !!
    !! The old runtime metadata module remains the frequency-diagnostic
    !! contract.  This module is deliberately narrower: it is the torque
    !! lane's per-process delivery seam.  No launcher aggregate is accepted,
    !! and the record is published with link(temporary, destination), so an
    !! existing record can never be replaced.
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_gc_full_fow_runtime_metadata, only: &
        GC_FULL_FOW_CONJUGATE_POLICY, &
        GC_FULL_FOW_ACTION_CONVENTION, GC_FULL_FOW_BOUND_METHOD, &
        GC_FULL_FOW_PREFACTOR_CONVENTION, &
        GC_FULL_FOW_REAL_FIELD_AMPLITUDE_CONVENTION, &
        format_gc_full_fow_frequency_convention, gc_full_fow_runtime_backend_state_t
    implicit none
    private

    character(len=*), parameter, public :: &
        GC_FULL_FOW_DELIVERY_SCHEMA = 'neort-full-fow-runtime-metadata-v2'
    character(len=*), parameter, public :: &
        GC_FULL_FOW_DELIVERY_SUCCESS = 'surface record emitted'

    integer, parameter, public :: GC_FULL_FOW_DELIVERY_OK = 0
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_NOT_REQUESTED = 1
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_ENVIRONMENT_ERROR = 2
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_INVALID_INPUT = 3
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_NOT_RUNTIME = 4
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_CERTIFICATION_REQUIRED = 5
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_WALL_ERROR = 6
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_BACKEND_ERROR = 7
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_STATUS_ERROR = 8
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_IO_ERROR = 9
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_TARGET_EXISTS = 10
    integer, parameter, public :: GC_FULL_FOW_DELIVERY_PUBLISH_ERROR = 11

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

    public :: emit_gc_full_fow_runtime_surface_record

contains

    subroutine emit_gc_full_fow_runtime_surface_record(phase, surface_key, surface, &
            state, invariant_status_coverage, return_status_coverage, &
            wall_status_coverage, status, message, output_path)
        character(len=*), intent(in) :: phase, surface_key
        real(dp), intent(in) :: surface
        type(gc_full_fow_runtime_backend_state_t), intent(in) :: state
        character(len=*), intent(in) :: invariant_status_coverage
        character(len=*), intent(in) :: return_status_coverage
        character(len=*), intent(in) :: wall_status_coverage
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=*), intent(in), optional :: output_path

        character(len=4096) :: destination, temporary, wall_relative
        character(len=256) :: invariant_normalized, return_normalized
        character(len=256) :: wall_normalized
        character(len=16) :: invariant_names(2), return_names(5), wall_names(3)
        integer :: invariant_counts(2), return_counts(5), wall_counts(3)
        integer :: unit, io_status, flush_status, close_status
        logical :: exists, published

        invariant_names = ''
        invariant_names(1) = 'success'
        invariant_names(2) = 'failure'
        return_names = ''
        return_names(1) = 'success'
        return_names(2) = 'no_return'
        return_names(3) = 'radial_domain'
        return_names(4) = 'wall_loss'
        return_names(5) = 'error'
        wall_names = ''
        wall_names(1) = 'not_hit'
        wall_names(2) = 'wall_loss'
        wall_names(3) = 'error'

        call set_result(status, message, GC_FULL_FOW_DELIVERY_INVALID_INPUT, '')
        call read_lane_kind(status, message)
        if (status /= GC_FULL_FOW_DELIVERY_OK) return

        if (present(output_path)) then
            destination = output_path
        else
            call read_environment('NEORT_FULL_FOW_METADATA_PATH', destination, &
                status, message)
            if (status /= GC_FULL_FOW_DELIVERY_OK) return
        end if
        call validate_destination(destination, status, message)
        if (status /= GC_FULL_FOW_DELIVERY_OK) return

        call read_environment('NEORT_FULL_FOW_WALL_RELATIVE_PATH', wall_relative, &
            status, message)
        if (status /= GC_FULL_FOW_DELIVERY_OK) return
        call validate_relative_path(wall_relative, status, message)
        if (status /= GC_FULL_FOW_DELIVERY_OK) return

        call validate_record_identity(phase, surface_key, surface, state, status, &
            message)
        if (status /= GC_FULL_FOW_DELIVERY_OK) return
        call normalize_coverage(invariant_status_coverage, &
            invariant_names, invariant_counts, invariant_normalized, &
            status, message)
        if (status /= GC_FULL_FOW_DELIVERY_OK) return
        call normalize_coverage(return_status_coverage, &
            return_names, &
            return_counts, return_normalized, status, message)
        if (status /= GC_FULL_FOW_DELIVERY_OK) return
        call normalize_coverage(wall_status_coverage, &
            wall_names, wall_counts, wall_normalized, &
            status, message)
        if (status /= GC_FULL_FOW_DELIVERY_OK) return

        inquire (file=trim(destination), exist=exists, iostat=io_status)
        if (io_status /= 0) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_IO_ERROR, &
                'could not inspect surface metadata destination')
            return
        end if
        if (exists) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_TARGET_EXISTS, &
                'surface metadata destination already exists')
            return
        end if

        call make_temporary_path(destination, surface_key, temporary, status, &
            message)
        if (status /= GC_FULL_FOW_DELIVERY_OK) return
        open (newunit=unit, file=trim(temporary), status='new', action='write', &
            form='formatted', iostat=io_status)
        if (io_status /= 0) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_IO_ERROR, &
                'could not create temporary surface metadata record')
            return
        end if

        call write_pair(unit, 'schema', GC_FULL_FOW_DELIVERY_SCHEMA, io_status)
        if (io_status == 0) call write_pair(unit, 'record_kind', 'surface', io_status)
        if (io_status == 0) call write_pair(unit, 'metadata_origin', 'runtime', io_status)
        if (io_status == 0) call write_pair(unit, 'runtime_emitter', 'neo_rt', io_status)
        if (io_status == 0) call write_pair(unit, 'phase', trim(phase), io_status)
        if (io_status == 0) call write_pair(unit, 'lane_kind', 'torque', io_status)
        if (io_status == 0) call write_pair(unit, 'surface_key', trim(surface_key), io_status)
        if (io_status == 0) call write_real_pair(unit, 'surface', surface, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'model', state%model, io_status)
        if (io_status == 0) call write_integer_pair(unit, 'frequency_model', &
            state%frequency_model, io_status)
        if (io_status == 0) call write_pair(unit, 'backend', trim(state%backend), io_status)
        if (io_status == 0) call write_pair(unit, 'coordinates', trim(state%coordinates), io_status)
        if (io_status == 0) call write_pair(unit, 'wall_path', trim(wall_relative), io_status)
        if (io_status == 0) call write_pair(unit, 'wall_sha256', trim(state%wall_sha256), io_status)
        if (io_status == 0) call write_pair(unit, 'wall_units', trim(state%wall_units), io_status)
        if (io_status == 0) call write_logical_pair(unit, 'wall_certified', &
            state%wall_certified, io_status)
        if (io_status == 0) call write_logical_pair(unit, 'canonical_measure_certified', &
            state%canonical_measure_certified, io_status)
        if (io_status == 0) call write_logical_pair(unit, 'component_identity_certified', &
            state%component_identity_certified, io_status)
        if (io_status == 0) call write_logical_pair(unit, 'nonlocal_transport_certified', &
            state%nonlocal_transport_certified, io_status)
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
        if (io_status == 0) call write_pair(unit, 'invariant_status_coverage', &
            trim(invariant_normalized), io_status)
        if (io_status == 0) call write_pair(unit, 'return_status_coverage', &
            trim(return_normalized), io_status)
        if (io_status == 0) call write_pair(unit, 'wall_status_coverage', &
            trim(wall_normalized), io_status)

        flush_status = 0
        if (io_status == 0) flush (unit, iostat=flush_status)
        close_status = 0
        close (unit, iostat=close_status)
        if (io_status /= 0 .or. flush_status /= 0 .or. close_status /= 0) then
            call unlink_file(temporary)
            call set_result(status, message, GC_FULL_FOW_DELIVERY_IO_ERROR, &
                'surface metadata record write failed')
            return
        end if

        call publish_file(temporary, destination, published)
        if (.not. published) then
            call unlink_file(temporary)
            call set_result(status, message, GC_FULL_FOW_DELIVERY_PUBLISH_ERROR, &
                'immutable surface metadata publication failed')
            return
        end if
        call unlink_file(temporary)
        call set_result(status, message, GC_FULL_FOW_DELIVERY_OK, &
            GC_FULL_FOW_DELIVERY_SUCCESS)
    end subroutine emit_gc_full_fow_runtime_surface_record

    subroutine read_lane_kind(status, message)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=64) :: lane

        call read_environment('NEORT_FULL_FOW_LANE_KIND', lane, status, message)
        if (status /= GC_FULL_FOW_DELIVERY_OK) then
            status = GC_FULL_FOW_DELIVERY_NOT_REQUESTED
            message = 'surface runtime delivery was not requested'
            return
        end if
        if (trim(lane) /= 'torque') then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_NOT_REQUESTED, &
                'surface runtime delivery is only for the torque lane')
        else
            call set_result(status, message, GC_FULL_FOW_DELIVERY_OK, 'ok')
        end if
    end subroutine read_lane_kind

    subroutine validate_record_identity(phase, surface_key, surface, state, status, &
            message)
        character(len=*), intent(in) :: phase, surface_key
        real(dp), intent(in) :: surface
        type(gc_full_fow_runtime_backend_state_t), intent(in) :: state
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        logical :: wall_exists
        integer :: io_status
        character(len=64) :: expected_frequency

        call set_result(status, message, GC_FULL_FOW_DELIVERY_OK, 'ok')
        call format_gc_full_fow_frequency_convention(state%toroidal_harmonic, &
            expected_frequency)
        if (trim(phase) /= 'phiI000' .and. trim(phase) /= 'phiI010') then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_INVALID_INPUT, &
                'surface record phase is not an ITER phase')
            return
        end if
        if (.not. valid_token(surface_key)) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_INVALID_INPUT, &
                'surface record key is not a safe token')
            return
        end if
        if (surface <= 0.0_dp .or. surface >= 1.0_dp) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_INVALID_INPUT, &
                'surface is outside the native flux interval')
            return
        end if
        if (.not. state%runtime_execution_complete) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_NOT_RUNTIME, &
                'surface record is not from completed executable state')
            return
        end if
        if (state%model /= 2 .or. state%frequency_model /= 2 .or. &
                trim(state%backend) /= 'cylindrical_full_fow' .or. &
                trim(state%coordinates) /= 'R,Z,phi') then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_BACKEND_ERROR, &
                'surface record backend identity is not model-2 cylindrical full-FOW')
            return
        end if
        if (.not. state%orbit_backend_certified .or. &
                state%cylindrical_backend_entries <= 0 .or. &
                state%legacy_backend_entries /= 0 .or. &
                state%chart_fallback_entries /= 0) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_BACKEND_ERROR, &
                'surface record backend certificate is incomplete or fell back')
            return
        end if
        if (.not. state%nonlocal_transport_certified) then
            call set_result(status, message, &
                GC_FULL_FOW_DELIVERY_CERTIFICATION_REQUIRED, &
                'torque record requires certified full nonlocal transport')
            return
        end if
        if (len_trim(state%wall_actual_path) == 0 .or. &
                has_control_character(state%wall_actual_path)) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_WALL_ERROR, &
                'surface record wall path is absent or malformed')
            return
        end if
        inquire (file=trim(state%wall_actual_path), exist=wall_exists, &
            iostat=io_status)
        if (io_status /= 0 .or. .not. wall_exists .or. &
                trim(state%wall_units) /= 'm' .and. trim(state%wall_units) /= 'cm' .or. &
                .not. valid_sha256(state%wall_sha256) .or. .not. state%wall_certified) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_WALL_ERROR, &
                'surface record wall certificate is incomplete')
            return
        end if
        if (.not. state%canonical_measure_certified .or. &
                .not. state%component_identity_certified) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_CERTIFICATION_REQUIRED, &
                'surface record canonical/component certificates are incomplete')
            return
        end if
        if (trim(state%real_field_amplitude_convention) /= &
                GC_FULL_FOW_REAL_FIELD_AMPLITUDE_CONVENTION .or. &
                trim(state%conjugate_policy) /= GC_FULL_FOW_CONJUGATE_POLICY .or. &
                trim(state%prefactor_convention) /= GC_FULL_FOW_PREFACTOR_CONVENTION .or. &
                trim(state%action_convention) /= GC_FULL_FOW_ACTION_CONVENTION .or. &
                trim(state%phase_space_bound_method) /= GC_FULL_FOW_BOUND_METHOD .or. &
                trim(state%frequency_convention) /= trim(expected_frequency) .or. &
                len_trim(state%perturbation_input_path) == 0 .or. &
                .not. state%perturbation_provenance_certified .or. &
                state%quadrature_base_h0_order < 2 .or. &
                state%quadrature_base_jk_order < 2 .or. &
                state%quadrature_refined_h0_order /= &
                2*state%quadrature_base_h0_order .or. &
                state%quadrature_refined_jk_order /= &
                2*state%quadrature_base_jk_order .or. &
                .not. ieee_is_finite(state%quadrature_relative_tolerance) .or. &
                .not. ieee_is_finite(state%quadrature_absolute_tolerance) .or. &
                state%quadrature_relative_tolerance <= 0.0_dp .or. &
                state%quadrature_absolute_tolerance <= 0.0_dp .or. &
                state%poloidal_harmonic_count <= 0 .or. &
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
                GC_FULL_FOW_DELIVERY_CERTIFICATION_REQUIRED, &
                'surface transport provenance is incomplete or not exact')
        end if
    end subroutine validate_record_identity

    subroutine normalize_coverage(text, names, counts, normalized, status, message)
        character(len=*), intent(in) :: text
        character(len=*), intent(in) :: names(:)
        integer, intent(out) :: counts(:)
        character(len=*), intent(out) :: normalized
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=256) :: item, name, number
        integer :: seen(size(names)), i, start, finish, colon, ios, value
        logical :: found

        counts = 0
        seen = 0
        normalized = ''
        start = 1
        do
            if (start > len_trim(text)) exit
            finish = index(text(start:), ',')
            if (finish == 0) then
                finish = len_trim(text) - start + 1
            else
                finish = finish - 1
            end if
            item = trim(text(start:start + finish - 1))
            colon = index(item, ':')
            if (colon <= 1) then
                call set_result(status, message, GC_FULL_FOW_DELIVERY_STATUS_ERROR, &
                    'malformed runtime status coverage')
                return
            end if
            name = trim(item(:colon - 1))
            number = trim(item(colon + 1:))
            read (number, *, iostat=ios) value
            if (ios /= 0 .or. value < 0) then
                call set_result(status, message, GC_FULL_FOW_DELIVERY_STATUS_ERROR, &
                    'runtime status coverage contains a nonnegative integer error')
                return
            end if
            found = .false.
            do i = 1, size(names)
                if (name == names(i)) then
                    if (seen(i) /= 0) then
                        call set_result(status, message, &
                            GC_FULL_FOW_DELIVERY_STATUS_ERROR, &
                            'runtime status coverage contains a duplicate')
                        return
                    end if
                    seen(i) = 1
                    counts(i) = value
                    found = .true.
                    exit
                end if
            end do
            if (.not. found) then
                call set_result(status, message, GC_FULL_FOW_DELIVERY_STATUS_ERROR, &
                    'runtime status coverage contains an unknown status')
                return
            end if
            if (start + finish > len_trim(text)) exit
            start = start + finish + 1
        end do
        if (any(seen == 0) .or. sum(counts) <= 0) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_STATUS_ERROR, &
                'runtime status coverage is incomplete or empty')
            return
        end if
        do i = 1, size(names)
            if (i > 1) normalized = trim(normalized)//','
            write (number, '(I0)') counts(i)
            normalized = trim(normalized)//trim(names(i))//':'//trim(number)
        end do
        call set_result(status, message, GC_FULL_FOW_DELIVERY_OK, 'ok')
    end subroutine normalize_coverage

    subroutine read_environment(name, value, status, message)
        character(len=*), intent(in) :: name
        character(len=*), intent(out) :: value
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer :: length_value, environment_status

        value = ''
        call get_environment_variable(name, value=value, length=length_value, &
            status=environment_status)
        if (environment_status /= 0 .or. length_value <= 0) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_ENVIRONMENT_ERROR, &
                'required runtime delivery environment is missing')
            return
        end if
        if (length_value > len(value)) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_ENVIRONMENT_ERROR, &
                'runtime delivery environment value is truncated')
            return
        end if
        value = value(:length_value)
        call set_result(status, message, GC_FULL_FOW_DELIVERY_OK, 'ok')
    end subroutine read_environment

    subroutine validate_destination(path, status, message)
        character(len=*), intent(in) :: path
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        logical :: exists
        integer :: io_status

        call set_result(status, message, GC_FULL_FOW_DELIVERY_OK, 'ok')
        if (len_trim(path) == 0 .or. has_control_character(path) .or. &
                index(trim(path), 'FULL_FOW_RUNTIME_METADATA.') == 0) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_INVALID_INPUT, &
                'surface metadata destination is not a per-surface record path')
            return
        end if
        inquire (file=trim(path), exist=exists, iostat=io_status)
        if (io_status /= 0) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_IO_ERROR, &
                'could not inspect surface metadata destination')
        end if
    end subroutine validate_destination

    subroutine validate_relative_path(path, status, message)
        character(len=*), intent(in) :: path
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        character(len=1024) :: copy

        copy = trim(path)
        call set_result(status, message, GC_FULL_FOW_DELIVERY_OK, 'ok')
        if (len_trim(copy) == 0 .or. copy(1:1) == '/' .or. &
                has_control_character(copy) .or. index(copy, '..') /= 0) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_WALL_ERROR, &
                'wall path is not a safe relative path')
        end if
    end subroutine validate_relative_path

    subroutine make_temporary_path(destination, surface_key, temporary, status, message)
        character(len=*), intent(in) :: destination, surface_key
        character(len=*), intent(out) :: temporary
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer(int64) :: ticks
        character(len=64) :: suffix

        call system_clock(count=ticks)
        write (suffix, '(I0)') ticks
        temporary = trim(destination)//'.tmp.'//trim(surface_key)//'.'//trim(suffix)
        if (len_trim(temporary) >= len(temporary)) then
            call set_result(status, message, GC_FULL_FOW_DELIVERY_INVALID_INPUT, &
                'temporary surface metadata path is too long')
            return
        end if
        call set_result(status, message, GC_FULL_FOW_DELIVERY_OK, 'ok')
    end subroutine make_temporary_path

    subroutine write_pair(unit, key, value, status)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: key, value
        integer, intent(out) :: status

        write (unit, '(A,"=",A)', iostat=status) trim(key), trim(value)
    end subroutine write_pair

    subroutine write_real_pair(unit, key, value, status)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: key
        real(dp), intent(in) :: value
        integer, intent(out) :: status

        write (unit, '(A,"=",ES24.16E3)', iostat=status) trim(key), value
    end subroutine write_real_pair

    subroutine write_integer_pair(unit, key, value, status)
        integer, intent(in) :: unit, value
        character(len=*), intent(in) :: key
        integer, intent(out) :: status
        character(len=32) :: text

        write (text, '(I0)') value
        call write_pair(unit, key, text, status)
    end subroutine write_integer_pair

    subroutine write_logical_pair(unit, key, value, status)
        integer, intent(in) :: unit
        character(len=*), intent(in) :: key
        logical, intent(in) :: value
        integer, intent(out) :: status

        if (value) then
            call write_pair(unit, key, 'true', status)
        else
            call write_pair(unit, key, 'false', status)
        end if
    end subroutine write_logical_pair

    subroutine publish_file(old_path, new_path, success)
        character(len=*), intent(in) :: old_path, new_path
        logical, intent(out) :: success
        character(kind=c_char), allocatable :: old_c(:), new_c(:)
        integer(c_int) :: code

        call make_c_string(old_path, old_c)
        call make_c_string(new_path, new_c)
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

    logical function valid_token(value)
        character(len=*), intent(in) :: value
        integer :: i, code

        valid_token = len_trim(value) > 0
        if (.not. valid_token) return
        do i = 1, len_trim(value)
            code = iachar(value(i:i))
            if (.not. (code >= iachar('0') .and. code <= iachar('9')) .and. &
                    .not. (code >= iachar('A') .and. code <= iachar('Z')) .and. &
                    .not. (code >= iachar('a') .and. code <= iachar('z')) .and. &
                    value(i:i) /= '_' .and. value(i:i) /= '-') then
                valid_token = .false.
                return
            end if
        end do
    end function valid_token

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

    logical function has_control_character(value)
        character(len=*), intent(in) :: value
        integer :: i, code

        has_control_character = .false.
        do i = 1, len_trim(value)
            code = iachar(value(i:i))
            if (code < 32 .or. code == 127) then
                has_control_character = .true.
                return
            end if
        end do
    end function has_control_character

    subroutine set_result(status, message, new_status, new_message)
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        integer, intent(in) :: new_status
        character(len=*), intent(in) :: new_message

        status = new_status
        message = ''
        message = trim(new_message)
    end subroutine set_result

end module neort_gc_full_fow_runtime_delivery
