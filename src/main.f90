module neort_main

    implicit none

    character(1024) :: runname

contains

    subroutine main
        use do_magfie_mod, only: do_magfie, s
        use driftorbit, only: bfac, efac, frequency_model, &
            FREQUENCY_MODEL_GC_FULL, m0, mph
        use neort_gc_full_fow_runtime_delivery, only: &
            GC_FULL_FOW_DELIVERY_OK
        use neort_gc_wall_context, only: configured_wall_file, configured_wall_units
        use neort, only: check_magfie, configure_model2_transport, &
            get_model2_transport_execution, write_magfie_data_to_files, &
            write_transport_data_to_files
        use neort_gc_eqdsk_nonlocal_transport, only: &
            GC_EQDSK_NONLOCAL_SUCCESS, &
            gc_eqdsk_nonlocal_factory_t, gc_eqdsk_nonlocal_options_t, &
            gc_eqdsk_nonlocal_species_t, &
            initialize_gc_eqdsk_nonlocal_transport
        use neort_gc_model2_transport_dispatch, only: &
            GC_MODEL2_DISPATCH_SUCCESS, emit_gc_model2_runtime_record, &
            gc_model2_backend_evidence_t, gc_model2_transport_execution_t, &
            gc_model2_transport_options_t
        use neort_gc_nonlocal_transport_types, only: &
            gc_nonlocal_transport_options_t
        use neort_datatypes, only: magfie_data_t
        use neort_lib
        use neort_profiles, only: vth
        use neort_gc_nonlocal_transport_types, only: &
            gc_nonlocal_transport_options_t
        use util, only: c, files_exist, mi, qi
        use iso_fortran_env, only: dp => real64

        character(len=*), parameter :: boozer_file = "in_file"
        character(len=*), parameter :: boozer_pert_file = "in_file_pert"
        character(len=*), parameter :: plasma_file = "plasma.in"
        character(len=*), parameter :: profile_file = "profile.in"

        type(magfie_data_t) :: magfie_data
        type(transport_data_t) :: transport_data
        type(gc_eqdsk_nonlocal_factory_t), pointer :: model2_factory => null()
        type(gc_model2_transport_execution_t) :: model2_execution
        character(len=64) :: phase, lane
        character(len=256) :: delivery_message
        integer :: delivery_status, phase_length, lane_length, environment_status
        integer :: model2_status

        call get_command_argument(1, runname)

        call neort_init(trim(runname)//".in", boozer_file, boozer_pert_file)

        if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
            call require_model2_inputs(boozer_file, boozer_pert_file, plasma_file, &
                profile_file, configured_wall_file)
            call neort_prepare_splines(plasma_file, profile_file)
            ! This initializes the common per-surface plasma/profile state.  The
            ! model-2 transport branch below does not enter the local eta loop.
            call neort_setup_at_s(s)
            call configure_model2_from_executable_inputs(model2_factory, &
                model2_status)
            if (model2_status /= GC_MODEL2_DISPATCH_SUCCESS) then
                write (0, '(A,I0)') 'model-2 factory configuration failed: ', &
                    model2_status
                error stop 'direct model-2 factory unavailable'
            end if
            call neort_compute_at_s(s, transport_data)
            call get_model2_transport_execution(model2_execution, model2_status)
            if (model2_status /= GC_MODEL2_DISPATCH_SUCCESS) then
                write (0, '(A,I0)') 'model-2 transport execution failed: ', &
                    model2_status
                error stop 'direct model-2 transport uncertified'
            end if
        else if (files_exist(plasma_file, profile_file)) then
            call neort_prepare_splines(plasma_file, profile_file)
            call neort_compute_at_s(s, transport_data)
        else
            call neort_compute_no_splines(transport_data)
        end if

        call check_magfie(magfie_data)  ! diagnostics

        if (frequency_model == FREQUENCY_MODEL_GC_FULL) then
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
                call emit_gc_model2_runtime_record(trim(phase), trim(runname), s, &
                    model2_execution, delivery_status, delivery_message)
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

    subroutine require_model2_inputs(eqdsk_file, perturbation_file, plasma_file, &
            profile_file, wall_file)
        character(len=*), intent(in) :: eqdsk_file, perturbation_file, plasma_file
        character(len=*), intent(in) :: profile_file, wall_file
        logical :: exists

        inquire (file=trim(eqdsk_file), exist=exists)
        if (.not. exists) error stop 'model-2 EQDSK input unavailable'
        inquire (file=trim(perturbation_file), exist=exists)
        if (.not. exists) error stop 'model-2 perturbation input unavailable'
        inquire (file=trim(plasma_file), exist=exists)
        if (.not. exists) error stop 'model-2 plasma input unavailable'
        inquire (file=trim(profile_file), exist=exists)
        if (.not. exists) error stop 'model-2 profile input unavailable'
        inquire (file=trim(wall_file), exist=exists)
        if (.not. exists) error stop 'model-2 wall input unavailable'
    end subroutine require_model2_inputs

    subroutine configure_model2_from_executable_inputs(factory, status)
        use do_magfie_mod, only: q, s
        use neort, only: configure_model2_transport, harmonic_bounds, mth_max_abs
        use neort_gc_nonlocal_transport_types, only: &
            gc_nonlocal_transport_options_t
        use neort_gc_model2_transport_dispatch, only: &
            GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE, &
            gc_model2_backend_evidence_t, gc_model2_transport_options_t
        use neort_gc_eqdsk_nonlocal_transport, only: &
            GC_EQDSK_NONLOCAL_SUCCESS, &
            gc_eqdsk_nonlocal_factory_t, gc_eqdsk_nonlocal_options_t, &
            gc_eqdsk_nonlocal_species_t, &
            initialize_gc_eqdsk_nonlocal_transport
        use neort_gc_wall_context, only: configured_wall_file, configured_wall_units
        use driftorbit, only: bfac, efac, m0, mph
        use neort_profiles, only: vth
        use util, only: c, mi, qi
        use iso_fortran_env, only: dp => real64

        type(gc_eqdsk_nonlocal_factory_t), pointer, intent(out) :: factory
        integer, intent(out) :: status

        type(gc_eqdsk_nonlocal_species_t) :: species
        type(gc_eqdsk_nonlocal_options_t) :: factory_options
        type(gc_model2_transport_options_t) :: dispatch_options
        type(gc_model2_backend_evidence_t) :: backend
        type(gc_nonlocal_transport_options_t) :: transport_options
        integer :: harmonic_index, mth_min, mth_max

        status = GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE
        allocate(factory)
        if (vth <= 0.0_dp .or. mi <= 0.0_dp .or. abs(qi) <= tiny(qi)) return

        species%name = 'configured_species'
        species%mass_g = mi
        species%charge_esu = qi
        species%reference_energy_erg = 0.5_dp*mi*vth**2
        species%reference_velocity_cm_s = vth

        factory_options = gc_eqdsk_nonlocal_options_t()
        factory_options%field_scale = bfac
        factory_options%profile_electric_factor = efac
        factory_options%profile_bfactor = bfac
        factory_options%reference_surface = s
        factory_options%residence_target_s_tor = s
        dispatch_options = gc_model2_transport_options_t()
        transport_options = gc_nonlocal_transport_options_t()
        transport_options%h0_order = 16
        transport_options%jk_order = 16
        transport_options%max_h0_nodes = 256
        transport_options%max_jperp_nodes = 256
        transport_options%max_total_nodes = 65536
        call propagate_model2_integral_options(transport_options, dispatch_options)
        dispatch_options%integral%resonance_options%force_count = 3
        dispatch_options%force_layout%d11_slot = 1
        dispatch_options%force_layout%d12_slot = 2
        dispatch_options%force_layout%torque_slot = 3
        call harmonic_bounds(mph, q, mth_max_abs, mth_min, mth_max)
        if (mth_max < mth_min) return
        allocate(dispatch_options%poloidal_harmonics(mth_min:mth_max))
        do harmonic_index = mth_min, mth_max
            dispatch_options%poloidal_harmonics(harmonic_index) = harmonic_index
        end do
        dispatch_options%toroidal_harmonic = mph

        call initialize_gc_eqdsk_nonlocal_transport('in_file', &
            trim(configured_wall_file), trim(configured_wall_units), &
            'in_file_pert', 'plasma.in', 'profile.in', species, mph, factory, &
            status, factory_options, transport_options)
        if (status /= GC_EQDSK_NONLOCAL_SUCCESS) return

        backend = gc_model2_backend_evidence_t()
        backend%factory_available = factory%initialized
        backend%field_certified = factory%field_ready
        backend%profile_certified = factory%profile_ready
        backend%perturbation_certified = factory%perturbation_ready
        backend%perturbation_amplitude_convention = &
            factory%perturbation_amplitude_convention
        backend%perturbation_input_path = factory%perturbation_path
        backend%perturbation_provenance_certified = &
            factory%perturbation_provenance_certified
        backend%phase_space_bound_method = &
            factory%options%phase_space_bound_method
        backend%wall_certified = factory%wall_ready
        backend%wall_actual_path = factory%wall_path
        backend%wall_units = factory%wall_units
        backend%wall_sha256 = factory%wall_hash
        ! Topology, canonical measure, and component identity are certified
        ! after the physical node executes in neort's dispatcher boundary.
        backend%topology_certified = .false.
        backend%canonical_measure_certified = .false.
        backend%component_identity_certified = .false.
        call configure_model2_transport(factory, m0, mph, backend, dispatch_options, &
            status)
    end subroutine configure_model2_from_executable_inputs

    subroutine propagate_model2_integral_options(transport_options, dispatch_options)
        use neort_gc_nonlocal_transport_types, only: &
            gc_nonlocal_transport_options_t
        use neort_gc_model2_transport_dispatch, only: &
            gc_model2_transport_options_t

        type(gc_nonlocal_transport_options_t), intent(in) :: transport_options
        type(gc_model2_transport_options_t), intent(inout) :: dispatch_options

        dispatch_options%integral = transport_options
    end subroutine propagate_model2_integral_options

end module neort_main
