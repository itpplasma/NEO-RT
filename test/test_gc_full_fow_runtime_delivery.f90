program test_gc_full_fow_runtime_delivery
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use neort_gc_full_fow_runtime_delivery, only: &
        GC_FULL_FOW_DELIVERY_CERTIFICATION_REQUIRED, &
        GC_FULL_FOW_DELIVERY_INVALID_INPUT, GC_FULL_FOW_DELIVERY_OK, &
        GC_FULL_FOW_DELIVERY_TARGET_EXISTS, &
        emit_gc_full_fow_runtime_surface_record
    use neort_gc_full_fow_runtime_metadata, only: &
        GC_FULL_FOW_ACTION_CONVENTION, GC_FULL_FOW_BOUND_METHOD, &
        format_gc_full_fow_frequency_convention, &
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
    character(len=64) :: surface_key
    character(len=64) :: expected_frequency
    character(len=256) :: message
    character(len=256) :: messages(20)
    integer :: statuses(20), status, i, io_status
    integer(int64) :: clock_count
    logical :: exists

    call system_clock(count=clock_count)
    write (base_path, '(A,I0)') '/var/tmp/ert/full-fow-runtime-delivery-test-', &
        clock_count
    call execute_command_line('mkdir -p ' // trim(base_path), exitstat=status)
    call require(status == 0, 'could not create delivery test directory')
    wall_path = trim(base_path)//'/wall.dat'
    call write_wall(wall_path)
    call setenv_value('NEORT_FULL_FOW_LANE_KIND', 'torque', status)
    call require(status == 0, 'could not set delivery lane')
    call setenv_value('NEORT_FULL_FOW_WALL_RELATIVE_PATH', &
        'inputs/wall/full_fow_wall.dat', status)
    call require(status == 0, 'could not set wall relative path')

    state = gc_full_fow_runtime_backend_state_t()
    state%backend = 'cylindrical_full_fow'
    state%coordinates = 'R,Z,phi'
    state%model = 2
    state%frequency_model = 2
    state%wall_actual_path = trim(wall_path)
    state%wall_units = 'm'
    state%wall_sha256 = repeat('a', 64)
    state%runtime_execution_complete = .true.
    state%orbit_backend_certified = .true.
    state%wall_certified = .true.
    state%canonical_measure_certified = .true.
    state%component_identity_certified = .true.
    state%nonlocal_transport_certified = .true.
    state%cylindrical_backend_entries = 1
    state%legacy_backend_entries = 0
    state%chart_fallback_entries = 0
    state%real_field_amplitude_convention = &
        'real_field_amplitude_one_signed_n'
    state%conjugate_policy = 'conjugate_implicit'
    state%prefactor_convention = 'eq17_pi32_over_4_real_field'
    state%action_convention = GC_FULL_FOW_ACTION_CONVENTION
    state%phase_space_bound_method = GC_FULL_FOW_BOUND_METHOD
    state%perturbation_input_path = 'inputs/perturbation.dat'
    state%perturbation_provenance_certified = .true.
    state%quadrature_base_h0_order = 2
    state%quadrature_base_jk_order = 2
    state%quadrature_refined_h0_order = 4
    state%quadrature_refined_jk_order = 4
    state%quadrature_relative_tolerance = 1.0e-7_dp
    state%quadrature_absolute_tolerance = 1.0e-12_dp
    state%poloidal_harmonic_min = -8
    state%poloidal_harmonic_max = 8
    state%toroidal_harmonic = 3
    state%poloidal_harmonic_count = 17
    state%executed_harmonic_count = 17
    call format_gc_full_fow_frequency_convention(state%toroidal_harmonic, &
        expected_frequency)
    state%frequency_convention = expected_frequency
    state%quadrature_convergence_certified = .true.
    state%harmonic_batch_certified = .true.
    state%class_reconstruction_certified = .true.
    state%orbit_step_refinement_certified = .true.
    state%orbit_base_step = 2.0e-3_dp
    state%orbit_refined_step = 1.0e-3_dp

    statuses = -1
    messages = ''
    !$omp parallel do default(none) shared(base_path, state, statuses, messages) &
    !$omp private(i, surface_key, output_path)
    do i = 1, 20
        write (surface_key, '("s",I3.3)') i
        write (output_path, '(A,"/FULL_FOW_RUNTIME_METADATA.",A)') &
            trim(base_path), trim(surface_key)
        call emit_gc_full_fow_runtime_surface_record('phiI000', trim(surface_key), &
            real(i, dp)/25.0_dp, state, &
            'failure:0,success:1', &
            'success:1,no_return:0,radial_domain:0,wall_loss:0,error:0', &
            'not_hit:1,wall_loss:0,error:0', statuses(i), messages(i), output_path)
    end do
    !$omp end parallel do

    do i = 1, 20
        call require(statuses(i) == GC_FULL_FOW_DELIVERY_OK, &
            'parallel surface record did not publish')
        write (surface_key, '("s",I3.3)') i
        write (output_path, '(A,"/FULL_FOW_RUNTIME_METADATA.",A)') &
            trim(base_path), trim(surface_key)
        inquire (file=trim(output_path), exist=exists, iostat=io_status)
        call require(io_status == 0 .and. exists, &
            'parallel surface record is missing')
        call require(metadata_value(output_path, 'record_kind') == 'surface', &
            'surface record kind is wrong')
        call require(metadata_value(output_path, 'surface_key') == trim(surface_key), &
            'surface record key is wrong')
        call require(metadata_value(output_path, 'lane_kind') == 'torque', &
            'surface record lane is wrong')
        call require(metadata_value(output_path, 'prefactor_convention') == &
            'eq17_pi32_over_4_real_field', &
            'surface record lost the prefactor convention')
    end do

    write (surface_key, '("s",I3.3)') 1
    write (output_path, '(A,"/FULL_FOW_RUNTIME_METADATA.",A)') &
        trim(base_path), trim(surface_key)
    call emit_gc_full_fow_runtime_surface_record('phiI000', trim(surface_key), &
        1.0_dp/25.0_dp, state, 'success:1,failure:0', &
        'success:1,no_return:0,radial_domain:0,wall_loss:0,error:0', &
        'not_hit:1,wall_loss:0,error:0', status, message, output_path)
    call require(status == GC_FULL_FOW_DELIVERY_TARGET_EXISTS, &
        'immutable surface record was overwritten')

    state%nonlocal_transport_certified = .false.
    write (output_path, '(A,"/FULL_FOW_RUNTIME_METADATA.uncertified")') &
        trim(base_path)
    call emit_gc_full_fow_runtime_surface_record('phiI000', 'uncertified', &
        0.5_dp, state, 'success:1,failure:0', &
        'success:1,no_return:0,radial_domain:0,wall_loss:0,error:0', &
        'not_hit:1,wall_loss:0,error:0', status, message, output_path)
    call require(status == GC_FULL_FOW_DELIVERY_CERTIFICATION_REQUIRED, &
        'uncertified torque record was emitted')
    inquire (file=trim(output_path), exist=exists)
    call require(.not. exists, 'uncertified record left a file')

    write (output_path, '(A,"/FULL_FOW_RUNTIME_METADATA")') trim(base_path)
    call emit_gc_full_fow_runtime_surface_record('phiI000', 'generic', 0.5_dp, &
        state, 'success:1,failure:0', &
        'success:1,no_return:0,radial_domain:0,wall_loss:0,error:0', &
        'not_hit:1,wall_loss:0,error:0', status, message, output_path)
    call require(status == GC_FULL_FOW_DELIVERY_INVALID_INPUT, &
        'generic shared metadata path was accepted')

    print '(A)', 'test_gc_full_fow_runtime_delivery OK'

contains

    subroutine write_wall(path)
        character(len=*), intent(in) :: path
        integer :: unit

        open (newunit=unit, file=trim(path), status='replace', action='write', &
            iostat=io_status)
        call require(io_status == 0, 'could not create delivery test wall')
        write (unit, '(A)', iostat=io_status) '1.0 0.0'
        call require(io_status == 0, 'could not write delivery test wall')
        write (unit, '(A)', iostat=io_status) '1.1 0.0'
        call require(io_status == 0, 'could not write delivery test wall')
        write (unit, '(A)', iostat=io_status) '1.0 0.1'
        call require(io_status == 0, 'could not write delivery test wall')
        close (unit, iostat=io_status)
        call require(io_status == 0, 'could not close delivery test wall')
    end subroutine write_wall

    subroutine setenv_value(name, value, local_status)
        character(len=*), intent(in) :: name, value
        integer, intent(out) :: local_status
        character(kind=c_char), allocatable :: name_c(:), value_c(:)
        integer(c_int) :: code
        integer :: j, name_length, value_length

        name_length = len_trim(name)
        value_length = len_trim(value)
        allocate (name_c(name_length + 1), value_c(value_length + 1))
        do j = 1, name_length
            name_c(j) = name(j:j)
        end do
        do j = 1, value_length
            value_c(j) = value(j:j)
        end do
        name_c(name_length + 1) = c_null_char
        value_c(value_length + 1) = c_null_char
        code = c_setenv(name_c, value_c, 1_c_int)
        local_status = int(code)
        deallocate (name_c, value_c)
    end subroutine setenv_value

    character(len=4096) function metadata_value(path, key) result(value)
        character(len=*), intent(in) :: path, key
        character(len=4096) :: line, line_key
        integer :: unit, equal_position, local_status

        value = ''
        open (newunit=unit, file=trim(path), status='old', action='read', &
            iostat=local_status)
        if (local_status /= 0) return
        do
            read (unit, '(A)', iostat=local_status) line
            if (local_status /= 0) exit
            equal_position = index(line, '=')
            if (equal_position <= 1) cycle
            line_key = trim(line(:equal_position - 1))
            if (trim(line_key) == trim(key)) then
                value = trim(line(equal_position + 1:))
                exit
            end if
        end do
        close (unit)
    end function metadata_value

    subroutine require(condition, description)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: description

        if (.not. condition) then
            write (0, '(A)') trim(description)
            error stop 1
        end if
    end subroutine require

end program test_gc_full_fow_runtime_delivery
