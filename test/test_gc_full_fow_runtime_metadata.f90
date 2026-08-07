program test_gc_full_fow_runtime_metadata
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use neort_gc_full_fow_runtime_metadata, only: &
        GC_FULL_FOW_METADATA_BACKEND_MISMATCH, &
        GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH, &
        GC_FULL_FOW_METADATA_ENVIRONMENT_ERROR, &
        GC_FULL_FOW_METADATA_FALLBACK, &
        GC_FULL_FOW_METADATA_INVALID_HASH, &
        GC_FULL_FOW_METADATA_NOT_RUNTIME, &
        GC_FULL_FOW_METADATA_POLICY_INVALID, &
        GC_FULL_FOW_METADATA_SCHEMA, &
        GC_FULL_FOW_METADATA_STATUS_INVALID, &
        GC_FULL_FOW_METADATA_SUCCESS, &
        GC_FULL_FOW_METADATA_TARGET_EXISTS, &
        GC_FULL_FOW_METADATA_WALL_MISMATCH, &
        emit_gc_full_fow_runtime_metadata, &
        gc_full_fow_runtime_backend_state_t
    implicit none

    interface
        function c_setenv(name, value, overwrite) bind(C, name='setenv') result(code)
            import c_char, c_int
            character(kind=c_char), intent(in) :: name(*)
            character(kind=c_char), intent(in) :: value(*)
            integer(c_int), value :: overwrite
            integer(c_int) :: code
        end function c_setenv
    end interface

    type(gc_full_fow_runtime_backend_state_t) :: state
    character(len=4096) :: base_path, wall_path, output_path
    character(len=256) :: message
    integer :: status, parse_status
    integer(int64) :: clock_count
    logical :: exists

    call system_clock(count=clock_count)
    write (base_path, '(A,I0)') '/var/tmp/ert/full-fow-runtime-metadata-test-', &
        clock_count
    wall_path = trim(base_path)//'.wall'
    output_path = trim(base_path)//'.valid'
    call write_wall(wall_path)
    call set_runtime_environment(output_path, status)
    call require(status == 0, 'could not set runtime environment')

    state = gc_full_fow_runtime_backend_state_t()
    state%backend = 'cylindrical_full_fow'
    state%coordinates = 'R,Z,phi'
    state%model = 2
    state%frequency_model = 2
    state%wall_actual_path = trim(wall_path)
    state%wall_units = 'm'
    state%wall_sha256 = '9bf116333b88191a1085078cec12932ff0b4f3fbbd41c9e441e298d316f01fe2'
    state%runtime_execution_complete = .true.
    state%orbit_backend_certified = .true.
    state%wall_certified = .true.
    state%canonical_measure_certified = .true.
    state%component_identity_certified = .true.
    state%nonlocal_transport_certified = .false.
    state%cylindrical_backend_entries = 2
    state%legacy_backend_entries = 0
    state%chart_fallback_entries = 0
    state%real_field_amplitude_convention = &
        'real_field_amplitude_one_signed_n'
    state%conjugate_policy = 'conjugate_implicit'
    state%prefactor_convention = 'eq17_pi32_over_4_real_field'
    state%quadrature_base_h0_order = 2
    state%quadrature_base_jk_order = 2
    state%quadrature_refined_h0_order = 4
    state%quadrature_refined_jk_order = 4
    state%quadrature_relative_tolerance = 1.0e-7_dp
    state%quadrature_absolute_tolerance = 1.0e-12_dp
    state%poloidal_harmonic_min = -8
    state%poloidal_harmonic_max = 8
    state%toroidal_harmonic = 3
    state%quadrature_convergence_certified = .true.
    state%harmonic_batch_certified = .true.
    state%class_reconstruction_certified = .true.

    call emit(state, output_path, status, message)
    call require(status == GC_FULL_FOW_METADATA_SUCCESS, &
        'valid frequency-only metadata did not emit')
    call parse_emitted(output_path, parse_status)
    call require(parse_status == 0, 'independent parser rejected valid metadata')

    ! These certificates describe the orbit measure and component identity;
    ! they are required by the campaign verifier even for a frequency lane.
    call require(metadata_value(output_path, 'canonical_measure_certified') == 'true', &
        'frequency lane did not retain canonical-measure certification')
    call require(metadata_value(output_path, 'component_identity_certified') == 'true', &
        'frequency lane did not retain component certification')

    ! A frequency lane does not need a nonlocal transport certificate.
    call require(metadata_value(output_path, 'lane_kind') == 'frequency', &
        'valid frequency lane was mislabeled')
    call require(metadata_value(output_path, 'real_field_amplitude_convention') == &
        'real_field_amplitude_one_signed_n', &
        'amplitude convention was not emitted exactly')
    call require(metadata_value(output_path, 'conjugate_policy') == &
        'conjugate_implicit', 'conjugate policy was not emitted exactly')
    call require(metadata_value(output_path, 'prefactor_convention') == &
        'eq17_pi32_over_4_real_field', &
        'Eq. 17 prefactor convention was not emitted exactly')
    call require(metadata_value(output_path, 'quadrature_refined_h0_order') == '4', &
        'refined H0 order was not emitted')
    call require(metadata_value(output_path, 'poloidal_harmonic_min') == '-8' .and. &
        metadata_value(output_path, 'poloidal_harmonic_max') == '8', &
        'harmonic range was not emitted')

    ! Publication is no-replace: a second emission cannot overwrite accepted
    ! runtime evidence.
    call emit(state, output_path, status, message)
    call require(status == GC_FULL_FOW_METADATA_TARGET_EXISTS, &
        'existing metadata destination was overwritten')

    output_path = trim(base_path)//'.torque-without-nonlocal'
    call emit_lane(state, output_path, 'torque', status, message)
    call require(status == GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH, &
        'torque lane accepted without nonlocal certification')
    state%nonlocal_transport_certified = .true.
    state%action_convention = 'J_K=mass*c*mu_phys/abs(q)'
    state%phase_space_bound_method = 'axisymmetric_geqdsk_fpol_fluxfunction'
    state%frequency_convention = 'm*omega_b+3*omega_phi'
    state%perturbation_input_path = 'inputs/perturbation/bmod_n.dat'
    state%perturbation_provenance_certified = .true.
    state%poloidal_harmonic_count = 17
    state%executed_harmonic_count = 17
    state%orbit_step_refinement_certified = .true.
    state%orbit_base_step = 2.0e-6_dp
    state%orbit_refined_step = 1.0e-6_dp
    state%orbit_period_refinement_error = 1.0e-10_dp
    state%orbit_delta_phi_refinement_error = 1.0e-10_dp
    state%orbit_omega_b_refinement_error = 1.0e-10_dp
    state%orbit_omega_phi_refinement_error = 1.0e-10_dp
    state%orbit_h_m_refinement_error = 1.0e-10_dp
    state%orbit_shell_refinement_error = 1.0e-10_dp
    output_path = trim(base_path)//'.torque'
    call emit_lane(state, output_path, 'torque', status, message)
    call require(status == GC_FULL_FOW_METADATA_SUCCESS, &
        'certified torque lane did not emit')
    call require(metadata_value(output_path, 'lane_kind') == 'torque', &
        'certified torque lane was mislabeled')
    state%nonlocal_transport_certified = .false.

    ! Failed runtime wall certification is fail-closed: no sidecar is
    ! published for the campaign verifier to mistake for accepted evidence.
    output_path = trim(base_path)//'.wall-false'
    call set_metadata_path(output_path, status)
    state%wall_certified = .false.
    call emit(state, output_path, status, message)
    call require(status == GC_FULL_FOW_METADATA_WALL_MISMATCH, &
        'false wall certificate was accepted')
    inquire (file=trim(output_path), exist=exists)
    call require(.not. exists, 'failed wall certificate left a metadata sidecar')
    state%wall_certified = .true.

    ! Certification claims without a nonlocal transport certificate are a
    ! structural mismatch, not an accepted approximation.
    output_path = trim(base_path)//'.certification-mismatch'
    call set_metadata_path(output_path, status)
    state%component_identity_certified = .false.
    call emit(state, output_path, status, message)
    call require(status == GC_FULL_FOW_METADATA_CERTIFICATION_MISMATCH, &
        'missing component certification was accepted')
    state%component_identity_certified = .true.

    output_path = trim(base_path)//'.launcher-only'
    call set_metadata_path(output_path, status)
    state%runtime_execution_complete = .false.
    call emit(state, output_path, status, message)
    call require(status == GC_FULL_FOW_METADATA_NOT_RUNTIME, &
        'launcher-only state was accepted')
    inquire (file=trim(output_path), exist=exists)
    call require(.not. exists, 'launcher-only state left a metadata sidecar')
    state%runtime_execution_complete = .true.

    output_path = trim(base_path)//'.truncated-hash'
    call set_metadata_path(output_path, status)
    state%wall_sha256 = repeat('a', 63)
    call emit(state, output_path, status, message)
    call require(status == GC_FULL_FOW_METADATA_INVALID_HASH, &
        'truncated wall hash was accepted')
    state%wall_sha256 = '9bf116333b88191a1085078cec12932ff0b4f3fbbd41c9e441e298d316f01fe2'

    output_path = trim(base_path)//'.unknown-status'
    call set_metadata_path(output_path, status)
    call emit_gc_full_fow_runtime_metadata('phiI000', 'frequency', state, &
        'success:1,failure:0,other:0', &
        'success:1,no_return:0,radial_domain:0,wall_loss:0,error:0', &
        'not_hit:1,wall_loss:0,error:0', 'phase_independent', 2, .true., &
        status, message)
    call require(status == GC_FULL_FOW_METADATA_STATUS_INVALID, &
        'unknown runtime status was accepted')

    output_path = trim(base_path)//'.duplicate-status'
    call set_metadata_path(output_path, status)
    call emit_gc_full_fow_runtime_metadata('phiI000', 'frequency', state, &
        'success:1,failure:0,success:0', &
        'success:1,no_return:0,radial_domain:0,wall_loss:0,error:0', &
        'not_hit:1,wall_loss:0,error:0', 'phase_independent', 2, .true., &
        status, message)
    call require(status == GC_FULL_FOW_METADATA_STATUS_INVALID, &
        'duplicate runtime status was accepted')

    output_path = trim(base_path)//'.fallback'
    call set_metadata_path(output_path, status)
    state%legacy_backend_entries = 1
    call emit(state, output_path, status, message)
    call require(status == GC_FULL_FOW_METADATA_FALLBACK, &
        'legacy backend fallback was accepted')
    state%legacy_backend_entries = 0

    output_path = trim(base_path)//'.bad-policy'
    call set_metadata_path(output_path, status)
    call emit_gc_full_fow_runtime_metadata('phiI000', 'frequency', state, &
        'success:1,failure:0', &
        'success:1,no_return:0,radial_domain:0,wall_loss:0,error:0', &
        'not_hit:1,wall_loss:0,error:0', 'phase_independent', 4, .true., &
        status, message)
    call require(status == GC_FULL_FOW_METADATA_POLICY_INVALID, &
        'phase policy mismatch was accepted')

    output_path = trim(base_path)//'.bad-backend'
    call set_metadata_path(output_path, status)
    state%backend = 'realspace_full'
    call emit(state, output_path, status, message)
    call require(status == GC_FULL_FOW_METADATA_BACKEND_MISMATCH, &
        'legacy real-space backend claim was accepted')
    state%backend = 'cylindrical_full_fow'

    call set_environment_missing(status)
    call emit_gc_full_fow_runtime_metadata('phiI000', 'frequency', state, &
        'success:1,failure:0', &
        'success:1,no_return:0,radial_domain:0,wall_loss:0,error:0', &
        'not_hit:1,wall_loss:0,error:0', 'phase_independent', 2, .true., &
        status, message)
    call require(status == GC_FULL_FOW_METADATA_ENVIRONMENT_ERROR, &
        'missing metadata environment was accepted')

    print '(A)', 'test_gc_full_fow_runtime_metadata OK'

contains

    subroutine emit(current_state, path, local_status, local_message)
        type(gc_full_fow_runtime_backend_state_t), intent(in) :: current_state
        character(len=*), intent(in) :: path
        integer, intent(out) :: local_status
        character(len=*), intent(out) :: local_message

        call emit_lane(current_state, path, 'frequency', local_status, local_message)
    end subroutine emit

    subroutine emit_lane(current_state, path, lane, local_status, local_message)
        type(gc_full_fow_runtime_backend_state_t), intent(in) :: current_state
        character(len=*), intent(in) :: path, lane
        integer, intent(out) :: local_status
        character(len=*), intent(out) :: local_message

        call set_metadata_path(path, local_status)
        call emit_gc_full_fow_runtime_metadata('phiI000', lane, &
            current_state, 'success:1,failure:0', &
            'success:1,no_return:0,radial_domain:0,wall_loss:0,error:0', &
            'not_hit:1,wall_loss:0,error:0', 'phase_independent', 2, .true., &
            local_status, local_message)
    end subroutine emit_lane

    subroutine write_wall(path)
        character(len=*), intent(in) :: path
        integer :: unit, io_status

        open (newunit=unit, file=trim(path), status='replace', action='write', &
            iostat=io_status)
        call require(io_status == 0, 'could not create test wall')
        write (unit, '(A)', iostat=io_status) '1.0 0.0'
        call require(io_status == 0, 'could not write test wall')
        write (unit, '(A)', iostat=io_status) '1.1 0.0'
        call require(io_status == 0, 'could not write test wall')
        write (unit, '(A)', iostat=io_status) '1.0 0.1'
        call require(io_status == 0, 'could not write test wall')
        close (unit, iostat=io_status)
        call require(io_status == 0, 'could not close test wall')
    end subroutine write_wall

    subroutine set_runtime_environment(path, local_status)
        character(len=*), intent(in) :: path
        integer, intent(out) :: local_status

        call setenv_value('NEORT_FULL_FOW_METADATA_PATH', path, local_status)
        if (local_status /= 0) return
        call setenv_value('NEORT_FULL_FOW_WALL_RELATIVE_PATH', &
            'inputs/wall/full_fow_wall.dat', local_status)
    end subroutine set_runtime_environment

    subroutine set_metadata_path(path, local_status)
        character(len=*), intent(in) :: path
        integer, intent(out) :: local_status

        call setenv_value('NEORT_FULL_FOW_METADATA_PATH', path, local_status)
    end subroutine set_metadata_path

    subroutine set_environment_missing(local_status)
        integer, intent(out) :: local_status

        call setenv_value('NEORT_FULL_FOW_METADATA_PATH', '', local_status)
        if (local_status /= 0) return
        call setenv_value('NEORT_FULL_FOW_WALL_RELATIVE_PATH', '', local_status)
    end subroutine set_environment_missing

    subroutine setenv_value(name, value, local_status)
        character(len=*), intent(in) :: name, value
        integer, intent(out) :: local_status
        character(kind=c_char), allocatable :: name_c(:), value_c(:)
        integer(c_int) :: code
        integer :: i, name_length, value_length

        name_length = len_trim(name)
        value_length = len_trim(value)
        allocate (name_c(name_length + 1), value_c(value_length + 1))
        do i = 1, name_length
            name_c(i) = name(i:i)
        end do
        do i = 1, value_length
            value_c(i) = value(i:i)
        end do
        name_c(name_length + 1) = c_null_char
        value_c(value_length + 1) = c_null_char
        code = c_setenv(name_c, value_c, 1_c_int)
        local_status = int(code)
        deallocate (name_c, value_c)
    end subroutine setenv_value

    subroutine parse_emitted(path, local_status)
        character(len=*), intent(in) :: path
        integer, intent(out) :: local_status
        character(len=4096) :: line, key, value
        character(len=40) :: required(55)
        character(len=4096) :: values(55)
        logical :: seen(55)
        integer :: unit, io_status, equal_position, nkeys, key_index
        logical :: found, duplicate

        required = ''
        required(1) = 'schema'
        required(2) = 'metadata_origin'
        required(3) = 'runtime_emitter'
        required(4) = 'phase'
        required(5) = 'lane_kind'
        required(6) = 'model'
        required(7) = 'frequency_model'
        required(8) = 'backend'
        required(9) = 'coordinates'
        required(10) = 'wall_path'
        required(11) = 'wall_sha256'
        required(12) = 'wall_units'
        required(13) = 'wall_certified'
        required(14) = 'canonical_measure_certified'
        required(15) = 'component_identity_certified'
        required(16) = 'cylindrical_backend_entries'
        required(17) = 'legacy_backend_entries'
        required(18) = 'chart_fallback_entries'
        required(19) = 'real_field_amplitude_convention'
        required(20) = 'conjugate_policy'
        required(21) = 'prefactor_convention'
        required(22) = 'action_convention'
        required(23) = 'phase_space_bound_method'
        required(24) = 'frequency_convention'
        required(25) = 'perturbation_input_path'
        required(26) = 'perturbation_provenance_certified'
        required(27) = 'quadrature_base_h0_order'
        required(28) = 'quadrature_base_jk_order'
        required(29) = 'quadrature_refined_h0_order'
        required(30) = 'quadrature_refined_jk_order'
        required(31) = 'quadrature_relative_tolerance'
        required(32) = 'quadrature_absolute_tolerance'
        required(33) = 'poloidal_harmonic_min'
        required(34) = 'poloidal_harmonic_max'
        required(35) = 'poloidal_harmonic_count'
        required(36) = 'executed_harmonic_count'
        required(37) = 'toroidal_harmonic'
        required(38) = 'quadrature_convergence_certified'
        required(39) = 'harmonic_batch_certified'
        required(40) = 'class_reconstruction_certified'
        required(41) = 'orbit_step_refinement_certified'
        required(42) = 'orbit_base_step'
        required(43) = 'orbit_refined_step'
        required(44) = 'orbit_period_refinement_error'
        required(45) = 'orbit_delta_phi_refinement_error'
        required(46) = 'orbit_omega_b_refinement_error'
        required(47) = 'orbit_omega_phi_refinement_error'
        required(48) = 'orbit_h_m_refinement_error'
        required(49) = 'orbit_shell_refinement_error'
        required(50) = 'invariant_status_coverage'
        required(51) = 'return_status_coverage'
        required(52) = 'wall_status_coverage'
        required(53) = 'frequency_phase_policy'
        required(54) = 'diagnostic_count'
        required(55) = 'phase_independent_evidence'
        values = ''
        seen = .false.
        nkeys = 0
        local_status = 1
        open (newunit=unit, file=trim(path), status='old', action='read', &
            iostat=io_status)
        if (io_status /= 0) return
        do
            read (unit, '(A)', iostat=io_status) line
            if (io_status < 0) exit
            if (io_status /= 0) then
                close (unit)
                return
            end if
            if (len_trim(line) == 0) cycle
            if (line(1:1) == '#') cycle
            equal_position = index(line, '=')
            if (equal_position <= 1) then
                close (unit)
                return
            end if
            key = trim(line(:equal_position - 1))
            value = trim(line(equal_position + 1:))
            key_index = find_required_key(key, required)
            if (key_index == 0) then
                close (unit)
                return
            end if
            duplicate = seen(key_index)
            if (duplicate) then
                close (unit)
                return
            end if
            nkeys = nkeys + 1
            if (nkeys > size(required)) then
                close (unit)
                return
            end if
            if (len_trim(value) > len(values)) then
                close (unit)
                return
            end if
            values(key_index) = value
            seen(key_index) = .true.
        end do
        close (unit, iostat=io_status)
        if (nkeys /= size(required)) return
        do key_index = 1, size(required)
            found = seen(key_index)
            if (.not. found) return
        end do
        if (metadata_value_from_arrays('schema', required, values) /= &
            GC_FULL_FOW_METADATA_SCHEMA) return
        if (metadata_value_from_arrays('metadata_origin', required, values) /= &
            'runtime') return
        if (metadata_value_from_arrays('runtime_emitter', required, values) /= &
            'neo_rt') return
        if (metadata_value_from_arrays('phase', required, values) /= 'phiI000') return
        if (metadata_value_from_arrays('lane_kind', required, values) /= 'frequency') return
        if (metadata_value_from_arrays('wall_certified', required, values) &
            /= 'true') return
        if (metadata_value_from_arrays('canonical_measure_certified', required, values) &
            /= 'true') return
        if (metadata_value_from_arrays('component_identity_certified', required, values) &
            /= 'true') return
        local_status = 0
    end subroutine parse_emitted

    integer function find_required_key(key, required) result(key_index)
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: required(:)
        integer :: i

        key_index = 0
        do i = 1, size(required)
            if (trim(key) == trim(required(i))) then
                key_index = i
                return
            end if
        end do
    end function find_required_key

    character(len=4096) function metadata_value(path, key) result(value)
        character(len=*), intent(in) :: path, key
        character(len=4096) :: line, line_key, line_value
        integer :: unit, io_status, equal_position

        value = ''
        open (newunit=unit, file=trim(path), status='old', action='read', &
            iostat=io_status)
        if (io_status /= 0) return
        do
            read (unit, '(A)', iostat=io_status) line
            if (io_status /= 0) exit
            equal_position = index(line, '=')
            if (equal_position <= 1) cycle
            line_key = trim(line(:equal_position - 1))
            if (trim(line_key) /= trim(key)) cycle
            line_value = trim(line(equal_position + 1:))
            value = line_value
            exit
        end do
        close (unit)
    end function metadata_value

    character(len=4096) function metadata_value_from_arrays(key, required, values) &
            result(value)
        character(len=*), intent(in) :: key
        character(len=*), intent(in) :: required(:), values(:)
        integer :: key_index

        value = ''
        key_index = find_required_key(key, required)
        if (key_index == 0) return
        value = values(key_index)
    end function metadata_value_from_arrays

    subroutine require(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        if (.not. condition) then
            write (0, '(A)') trim(description)
            error stop 1
        end if
    end subroutine require

end program test_gc_full_fow_runtime_metadata
