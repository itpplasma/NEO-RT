module neort_main

    implicit none

    character(1024) :: runname

contains

    subroutine main
        use do_magfie_mod, only: s
        use driftorbit, only: frequency_model
        use neort_gc_full_fow_runtime_delivery, only: &
            GC_FULL_FOW_DELIVERY_OK, &
            emit_gc_full_fow_runtime_surface_record
        use neort_gc_frequency_provider, only: &
            get_gc_frequency_runtime_metadata, gc_frequency_runtime_metadata_t
        use neort_gc_wall_context, only: configured_wall_file, configured_wall_units
        use neort, only: check_magfie, write_magfie_data_to_files, write_transport_data_to_files
        use neort_datatypes, only: magfie_data_t
        use neort_lib
        use util, only: files_exist

        character(len=*), parameter :: boozer_file = "in_file"
        character(len=*), parameter :: boozer_pert_file = "in_file_pert"
        character(len=*), parameter :: plasma_file = "plasma.in"
        character(len=*), parameter :: profile_file = "profile.in"

        type(magfie_data_t) :: magfie_data
        type(transport_data_t) :: transport_data
        type(gc_frequency_runtime_metadata_t) :: runtime_metadata
        character(len=64) :: phase, lane
        character(len=256) :: delivery_message
        integer :: delivery_status, phase_length, lane_length, environment_status

        call get_command_argument(1, runname)

        call neort_init(trim(runname)//".in", boozer_file, boozer_pert_file)

        if (files_exist(plasma_file, profile_file)) then
            call neort_prepare_splines(plasma_file, profile_file)
            call neort_compute_at_s(s, transport_data)
        else
            call neort_compute_no_splines(transport_data)
        end if

        call check_magfie(magfie_data)  ! diagnostics

        ! The frequency diagnostic owns its one runtime record.  The main
        ! executable only delivers torque-lane records, one per surface, and
        ! the production provider has not yet exposed a nonlocal transport
        ! certificate.  Keeping this false is intentional: a model-2 orbit
        ! run must never be promoted to certified torque by the launcher.
        if (frequency_model == 2) then
            phase = ''
            lane = ''
            call get_environment_variable('NEORT_FULL_FOW_PHASE', value=phase, &
                length=phase_length, status=environment_status)
            if (environment_status == 0) then
                if (phase_length > 0 .and. phase_length <= len(phase)) then
                    phase = phase(:phase_length)
                end if
            end if
            call get_environment_variable('NEORT_FULL_FOW_LANE_KIND', value=lane, &
                length=lane_length, status=environment_status)
            if (environment_status == 0) then
                if (lane_length > 0 .and. lane_length <= len(lane)) then
                    lane = lane(:lane_length)
                end if
            end if
            if (trim(lane) == 'torque') then
                call get_gc_frequency_runtime_metadata(runtime_metadata)
                call deliver_torque_record(trim(phase), trim(runname), s, &
                    runtime_metadata, delivery_status, delivery_message)
                if (delivery_status /= GC_FULL_FOW_DELIVERY_OK) then
                    write (0, '(A,A)') 'full-FOW torque record was not emitted: ', &
                        trim(delivery_message)
                    error stop 'uncertified full-FOW torque lane'
                end if
            end if
        end if

        call write_magfie_data_to_files(magfie_data, trim(runname))
        call write_transport_data_to_files(transport_data, trim(runname))
    end subroutine main

    subroutine deliver_torque_record(phase, surface_key, surface, metadata, &
            status, message)
        use neort_gc_frequency_provider, only: gc_frequency_runtime_metadata_t
        use neort_gc_full_fow_runtime_metadata, only: &
            gc_full_fow_runtime_backend_state_t
        use neort_gc_full_fow_runtime_delivery, only: &
            emit_gc_full_fow_runtime_surface_record
        use neort_gc_wall_context, only: configured_wall_file, configured_wall_units
        use iso_fortran_env, only: dp => real64
        character(len=*), intent(in) :: phase, surface_key
        real(dp), intent(in) :: surface
        type(gc_frequency_runtime_metadata_t), intent(in) :: metadata
        integer, intent(out) :: status
        character(len=*), intent(out) :: message
        type(gc_full_fow_runtime_backend_state_t) :: state

        state = gc_full_fow_runtime_backend_state_t()
        state%backend = metadata%backend
        state%coordinates = metadata%coordinates
        state%model = 2
        state%frequency_model = 2
        state%wall_actual_path = configured_wall_file
        state%wall_units = configured_wall_units
        state%wall_sha256 = metadata%wall_hash
        state%runtime_execution_complete = metadata%cylindrical_entry_count > 0
        state%orbit_backend_certified = trim(metadata%backend) == &
            'cylindrical_full_fow' .and. metadata%cylindrical_entry_count > 0
        state%wall_certified = metadata%wall_certified
        state%canonical_measure_certified = metadata%canonical_measure_certified
        state%component_identity_certified = metadata%component_identity_certified
        ! No production nonlocal transport provider is connected to the
        ! transport executable yet.  Do not infer this from orbit success.
        state%nonlocal_transport_certified = .false.
        state%cylindrical_backend_entries = metadata%cylindrical_entry_count
        state%legacy_backend_entries = metadata%legacy_entry_count
        state%chart_fallback_entries = 0

        ! Empty coverage is deliberate until the provider supplies its own
        ! per-surface counters.  The delivery routine rejects it before any
        ! file can be published, even if a future caller accidentally sets a
        ! certification bit without wiring the counters.
        call emit_gc_full_fow_runtime_surface_record(phase, surface_key, surface, &
            state, '', '', '', status, message)
    end subroutine deliver_torque_record

end module neort_main
