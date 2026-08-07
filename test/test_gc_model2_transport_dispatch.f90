module gc_model2_transport_dispatch_test_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_CALLBACK_FAILURE, GC_NONLOCAL_SAMPLE_VALID, &
        gc_nonlocal_component_t, gc_nonlocal_orbit_sample_t
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES, GC_NONLOCAL_SUCCESS, &
        gc_nonlocal_transport_provider_t, gc_nonlocal_transport_quadrature_t, &
        gc_nonlocal_transport_reference_t

    implicit none
    private

    type, public, extends(gc_nonlocal_transport_provider_t) :: &
            mock_model2_provider_t
        logical :: fail_orbit = .false.
        integer :: callback_count = 0
    contains
        procedure :: get_section_reference => mock_get_section_reference
        procedure :: get_quadrature => mock_get_quadrature
        procedure :: get_components => mock_get_components
        procedure :: evaluate_orbit => mock_evaluate_orbit
        procedure :: evaluate_profiles => mock_evaluate_profiles
        procedure :: evaluate_outer_measure_factor => &
            mock_evaluate_outer_factor
    end type mock_model2_provider_t

contains

    subroutine mock_get_section_reference(provider, reference, status)
        class(mock_model2_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(out) :: reference
        integer, intent(out) :: status

        provider%callback_count = provider%callback_count + 1
        reference = gc_nonlocal_transport_reference_t()
        reference%section_id = 1
        reference%section_coordinate = 'x'
        reference%section_units = 'cm'
        reference%section_position = [1.0_dp, 0.0_dp, 0.0_dp]
        reference%section_flux = 0.5_dp
        reference%p_zeta_orientation = 1
        reference%frequency_semantics = GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES
        reference%hamiltonian_includes_phi = .true.
        reference%fixed = .true.
        status = GC_NONLOCAL_SUCCESS
    end subroutine mock_get_section_reference

    subroutine mock_get_quadrature(provider, h0_order, jk_order, quadrature, status)
        class(mock_model2_provider_t), intent(inout) :: provider
        integer, intent(in) :: h0_order, jk_order
        type(gc_nonlocal_transport_quadrature_t), intent(out) :: quadrature
        integer, intent(out) :: status

        provider%callback_count = provider%callback_count + 1
        quadrature = gc_nonlocal_transport_quadrature_t()
        allocate(quadrature%h0(h0_order*jk_order), quadrature%j_k(h0_order*jk_order), &
            quadrature%weight(h0_order*jk_order), &
            quadrature%j_k_upper_bound(h0_order*jk_order))
        quadrature%h0 = 1.0_dp
        quadrature%j_k = 0.5_dp
        quadrature%weight = 1.0_dp/real(h0_order*jk_order, dp)
        quadrature%j_k_upper_bound = 1.0_dp
        quadrature%h0_order = h0_order
        quadrature%jk_order = jk_order
        quadrature%n_nodes = h0_order*jk_order
        quadrature%h0_min = 0.0_dp
        quadrature%h0_scale = 1.0_dp
        quadrature%paired_domain = .true.
        quadrature%domain_certified = .true.
        quadrature%converged = .true.
        status = GC_NONLOCAL_SUCCESS
    end subroutine mock_get_quadrature

    subroutine mock_get_components(provider, reference, h0, jperp, components, &
            status)
        class(mock_model2_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp
        type(gc_nonlocal_component_t), allocatable, intent(out) :: components(:)
        integer, intent(out) :: status

        associate (unused_reference => reference, unused_h0 => h0, &
                unused_jperp => jperp)
        end associate
        provider%callback_count = provider%callback_count + 1
        allocate(components(1))
        components(1) = gc_nonlocal_component_t(component_id=1, sigma=1, &
            x_min=0.1_dp, x_max=0.9_dp)
        status = GC_NONLOCAL_SUCCESS
    end subroutine mock_get_components

    subroutine mock_evaluate_orbit(provider, reference, h0, jperp, x, sigma, &
            component_id, sample, status)
        class(mock_model2_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
        integer, intent(out) :: status

        associate (unused_reference => reference, unused_h0 => h0, &
                unused_jperp => jperp)
        end associate
        provider%callback_count = provider%callback_count + 1
        sample = gc_nonlocal_orbit_sample_t()
        if (provider%fail_orbit) then
            status = GC_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        sample%status = GC_NONLOCAL_SAMPLE_VALID
        sample%component_id = component_id
        sample%sigma = sigma
        sample%psi_star = 2.0_dp*x
        sample%dpsi_star_dx = 2.0_dp
        sample%tau_b = 1.0_dp
        sample%omega_b = x
        sample%omega_phi = -0.5_dp
        sample%domega_b_dx = 1.0_dp
        sample%domega_phi_dx = 0.0_dp
        sample%h_m = cmplx(1.0_dp, 0.0_dp, kind=dp)
        sample%derivatives_available = .true.
        status = GC_NONLOCAL_SUCCESS
    end subroutine mock_evaluate_orbit

    subroutine mock_evaluate_profiles(provider, reference, h0, jperp, x, sigma, &
            component_id, sample, force_count, force_values, status)
        class(mock_model2_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id, force_count
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        real(dp), intent(out) :: force_values(force_count)
        integer, intent(out) :: status

        associate (unused_reference => reference, unused_h0 => h0, &
                unused_jperp => jperp, unused_x => x, unused_sigma => sigma, &
                unused_component_id => component_id, unused_sample => sample)
        end associate
        provider%callback_count = provider%callback_count + 1
        force_values = 0.0_dp
        if (force_count /= 3) then
            status = GC_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        force_values = [2.0_dp, -3.0_dp, 5.0_dp]
        status = GC_NONLOCAL_SUCCESS
    end subroutine mock_evaluate_profiles

    subroutine mock_evaluate_outer_factor(provider, reference, h0, jperp, x, &
            sigma, component_id, sample, outer_factor, status)
        class(mock_model2_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        real(dp), intent(out) :: outer_factor
        integer, intent(out) :: status

        associate (unused_reference => reference, unused_h0 => h0, &
                unused_jperp => jperp, unused_x => x, unused_sigma => sigma, &
                unused_component_id => component_id, unused_sample => sample)
        end associate
        provider%callback_count = provider%callback_count + 1
        outer_factor = 1.0_dp
        status = GC_NONLOCAL_SUCCESS
    end subroutine mock_evaluate_outer_factor

end module gc_model2_transport_dispatch_test_support

program test_gc_model2_transport_dispatch
    use, intrinsic :: iso_c_binding, only: c_char, c_int, c_null_char
    use, intrinsic :: iso_fortran_env, only: dp => real64, int64
    use gc_model2_transport_dispatch_test_support, only: &
        mock_model2_provider_t
    use neort_main, only: propagate_model2_integral_options
    use neort_gc_model2_transport_dispatch, only: &
        GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE, &
        GC_MODEL2_DISPATCH_NOT_CERTIFIED, GC_MODEL2_DISPATCH_SUCCESS, &
        gc_model2_backend_evidence_t, gc_model2_force_layout_t, &
        gc_model2_transport_execution_t, gc_model2_transport_options_t, &
        emit_gc_model2_runtime_record, execute_gc_model2_transport, &
        finalize_gc_model2_transport_execution, gc_model2_dispatch_required, &
        gc_model2_observed_evidence_t
    use neort_transport, only: legacy_eta_transport_selected
    use neort_gc_nonlocal_transport_types, only: &
        gc_nonlocal_transport_options_t
    use neort_gc_full_fow_runtime_metadata, only: &
        GC_FULL_FOW_BOUND_METHOD, GC_FULL_FOW_REAL_FIELD_AMPLITUDE_CONVENTION
    use driftorbit, only: FREQUENCY_MODEL_GC_FULL, FREQUENCY_MODEL_GC_THIN, &
        FREQUENCY_MODEL_LEGACY
    implicit none

    interface
        function c_setenv(name, value, overwrite) bind(C, name='setenv') result(code)
            import c_char, c_int
            character(kind=c_char), intent(in) :: name(*), value(*)
            integer(c_int), value :: overwrite
            integer(c_int) :: code
        end function c_setenv
    end interface

    type(mock_model2_provider_t) :: provider
    type(gc_model2_backend_evidence_t) :: backend
    type(gc_model2_transport_options_t) :: options
    type(gc_model2_transport_execution_t) :: execution
    type(gc_model2_observed_evidence_t) :: observed
    character(len=4096) :: base_path, wall_path, output_path
    character(len=256) :: message
    integer(int64) :: clock_count
    integer :: status, io_status
    logical :: exists

    call require(.not. gc_model2_dispatch_required(0), &
        'model 0 was routed to model-2 dispatch')
    call require(.not. gc_model2_dispatch_required(1), &
        'model 1 was routed to model-2 dispatch')
    call require(gc_model2_dispatch_required(2), &
        'model 2 was not selected for full-FOW dispatch')
    call require(legacy_eta_transport_selected(FREQUENCY_MODEL_LEGACY), &
        'model 0 did not retain the legacy eta transport loop')
    call require(.not. legacy_eta_transport_selected(FREQUENCY_MODEL_GC_THIN), &
        'model 1 still selected the removed real-space thin transport loop')
    call require(.not. legacy_eta_transport_selected(FREQUENCY_MODEL_GC_FULL), &
        'model 2 still selected the legacy eta transport loop')

    call test_nondefault_options_propagation()

    call configure_options(options)
    call configure_backend(backend, base_path, wall_path)
    call execute_gc_model2_transport(provider, 1, 1, backend, options, execution, &
        status)
    call require(status == GC_MODEL2_DISPATCH_SUCCESS, &
        'mock transport aggregation did not complete')
    call require(.not. execution%certified, &
        'transport was certified before provider evidence was finalized')
    call require(abs(execution%d11 - 4.0_dp) < 1.0e-12_dp, &
        'D11 force slot was altered')
    call require(abs(execution%d12 + 6.0_dp) < 1.0e-12_dp, &
        'D12 force slot was altered')
    call require(abs(execution%torque - 10.0_dp) < 1.0e-12_dp, &
        'torque force slot was altered')
    call require(execution%integral%weighted_nodes == 4 .and. &
        execution%integral%certified_nodes == 4 .and. &
        execution%integral%unresolved_nodes == 0, &
        'dispatch did not expose actual node counters')
    call require(execution%force_slots_accepted, &
        'signed force slots were not accepted')
    call require(.not. execution%runtime%nonlocal_transport_certified, &
        'runtime state was certified before provider evidence was finalized')
    call require(execution%runtime%cylindrical_backend_entries == 4, &
        'runtime state did not report all attempted physical nodes')
    call configure_observed_evidence(observed)
    call finalize_gc_model2_transport_execution(execution, backend, observed, status)
    call require(status == GC_MODEL2_DISPATCH_SUCCESS .and. execution%certified, &
        'provider-observed evidence did not certify the aggregate')
    call require(execution%runtime%nonlocal_transport_certified, &
        'runtime state lost certified nonlocal transport')
    call require(index(execution%invariant_status_coverage, 'success:1') > 0, &
        'invariant coverage omitted the observed successful node')
    call require(trim(execution%return_status_coverage) == &
        'success:1,no_return:0,radial_domain:0,wall_loss:0,error:0', &
        'return coverage was not provider-observed')
    call require(trim(execution%wall_status_coverage) == &
        'not_hit:1,wall_loss:0,error:0', &
        'wall coverage was not conserved')

    call setenv_value('NEORT_FULL_FOW_LANE_KIND', 'torque', status)
    call require(status == 0, 'could not set torque delivery lane')
    call setenv_value('NEORT_FULL_FOW_WALL_RELATIVE_PATH', 'wall.dat', status)
    call require(status == 0, 'could not set wall delivery path')
    write (output_path, '(A,"/FULL_FOW_RUNTIME_METADATA.s050")') trim(base_path)
    call emit_gc_model2_runtime_record('phiI000', 's050', 0.5_dp, execution, &
        status, message, output_path)
    call require(status == 0, 'certified runtime record was not emitted')
    call require(metadata_value(output_path, 'nonlocal_transport_certified') == &
        'true', 'runtime record lost nonlocal certification')
    call require(metadata_value(output_path, 'cylindrical_backend_entries') == &
        '4', 'runtime record lost actual backend count')
    inquire (file=trim(output_path), exist=exists, iostat=io_status)
    call require(io_status == 0 .and. exists, 'runtime record is missing')

    call configure_backend(backend, base_path, wall_path)
    backend%canonical_measure_certified = .false.
    call execute_gc_model2_transport(provider, 1, 1, backend, options, execution, &
        status)
    call require(status == GC_MODEL2_DISPATCH_SUCCESS, &
        'static backend check rejected pre-execution configuration')
    call configure_observed_evidence(observed)
    call finalize_gc_model2_transport_execution(execution, backend, observed, status)
    call require(status /= GC_MODEL2_DISPATCH_SUCCESS .and. &
        .not. execution%certified, &
        'missing canonical certificate was promoted by finalization')

    call configure_backend(backend, base_path, wall_path)
    call execute_gc_model2_transport(provider, 1, 1, backend, options, execution, &
        status)
    call require(status == GC_MODEL2_DISPATCH_SUCCESS, &
        'second aggregate execution did not complete')
    observed = gc_model2_observed_evidence_t()
    call finalize_gc_model2_transport_execution(execution, backend, observed, status)
    call require(status /= GC_MODEL2_DISPATCH_SUCCESS .and. &
        .not. execution%certified, 'empty evidence was accepted')

    provider = mock_model2_provider_t()
    backend%factory_available = .false.
    call execute_gc_model2_transport(provider, 1, 1, backend, options, execution, &
        status)
    call require(status == GC_MODEL2_DISPATCH_FACTORY_UNAVAILABLE, &
        'unavailable physical factory was not rejected')
    call require(.not. execution%certified .and. provider%callback_count == 0, &
        'factory-unavailable path entered the provider')

    call configure_backend(backend, base_path, wall_path)
    provider%fail_orbit = .true.
    write (output_path, '(A,"/FULL_FOW_RUNTIME_METADATA.failed")') trim(base_path)
    call execute_gc_model2_transport(provider, 1, 1, backend, options, execution, &
        status)
    call require(status /= GC_MODEL2_DISPATCH_SUCCESS .and. &
        .not. execution%certified, 'failed physical run was certified')
    call emit_gc_model2_runtime_record('phiI000', 'failed', 0.5_dp, execution, &
        status, message, output_path)
    call require(status == GC_MODEL2_DISPATCH_NOT_CERTIFIED, &
        'failed physical run reached runtime delivery')
    inquire (file=trim(output_path), exist=exists)
    call require(.not. exists, 'failed physical run left a runtime record')

    print '(A)', 'test_gc_model2_transport_dispatch OK'

contains

    subroutine test_nondefault_options_propagation()
        type(gc_nonlocal_transport_options_t) :: transport_options
        type(gc_model2_transport_options_t) :: dispatch_options

        transport_options = gc_nonlocal_transport_options_t()
        transport_options%h0_order = 7
        transport_options%jk_order = 11
        transport_options%max_h0_nodes = 23
        transport_options%max_jperp_nodes = 29
        transport_options%max_total_nodes = 667
        transport_options%quadrature_relative_tolerance = 3.0e-8_dp
        transport_options%quadrature_absolute_tolerance = 4.0e-13_dp
        dispatch_options = gc_model2_transport_options_t()
        call propagate_model2_integral_options(transport_options, dispatch_options)
        call require(dispatch_options%integral%h0_order == 7 .and. &
            dispatch_options%integral%jk_order == 11 .and. &
            dispatch_options%integral%max_total_nodes == 667, &
            'non-default quadrature orders/capacity did not propagate')
        call require(abs(dispatch_options%integral%quadrature_relative_tolerance - &
            3.0e-8_dp) < 1.0e-20_dp .and. &
            abs(dispatch_options%integral%quadrature_absolute_tolerance - &
            4.0e-13_dp) < 1.0e-25_dp, &
            'non-default quadrature tolerances did not propagate')
    end subroutine test_nondefault_options_propagation

    subroutine configure_options(local_options)
        type(gc_model2_transport_options_t), intent(out) :: local_options

        local_options = gc_model2_transport_options_t()
        local_options%force_layout = gc_model2_force_layout_t(1, 2, 3)
        local_options%integral%max_h0_nodes = 2
        local_options%integral%max_jperp_nodes = 2
        local_options%integral%max_total_nodes = 4
        local_options%integral%h0_order = 2
        local_options%integral%jk_order = 2
        local_options%integral%require_converged = .false.
        local_options%integral%resonance_options%scan_intervals = 16
        local_options%integral%resonance_options%max_root_iterations = 64
        local_options%integral%resonance_options%max_roots = 4
        local_options%integral%resonance_options%force_count = 3
        local_options%poloidal_harmonics = [1]
        local_options%toroidal_harmonic = 1
    end subroutine configure_options

    subroutine configure_observed_evidence(local_observed)
        type(gc_model2_observed_evidence_t), intent(out) :: local_observed

        local_observed = gc_model2_observed_evidence_t()
        local_observed%physical_return_attempts = 1
        local_observed%invariant_successes = 1
        local_observed%return_successes = 1
        local_observed%wall_not_hit = 1
        local_observed%topology_certification_attempts = 1
        local_observed%topology_certification_successes = 1
        local_observed%return_accounting_complete = .true.
    end subroutine configure_observed_evidence

    subroutine configure_backend(local_backend, local_base, local_wall)
        type(gc_model2_backend_evidence_t), intent(out) :: local_backend
        character(len=*), intent(out) :: local_base, local_wall
        integer(int64) :: local_clock
        character(len=1024) :: perturbation_path
        integer :: local_status, perturbation_unit, wall_unit

        call system_clock(count=local_clock)
        write (local_base, '(A,I0)') '/var/tmp/ert/model2-dispatch-test-', &
            local_clock
        call execute_command_line('mkdir -p ' // trim(local_base), &
            exitstat=local_status)
        call require(local_status == 0, 'could not create dispatch test directory')
        local_wall = trim(local_base)//'/wall.dat'
        open (newunit=wall_unit, file=trim(local_wall), status='replace', &
            action='write', iostat=local_status)
        call require(local_status == 0, 'could not create dispatch test wall')
        write (wall_unit, '(A)', iostat=local_status) 'wall'
        close (wall_unit, iostat=local_status)
        call require(local_status == 0, 'could not write dispatch test wall')
        perturbation_path = trim(local_base)//'/perturbation.dat'
        open (newunit=perturbation_unit, file=trim(perturbation_path), &
            status='replace', action='write', iostat=local_status)
        call require(local_status == 0, 'could not create dispatch test perturbation')
        write (perturbation_unit, '(A)', iostat=local_status) 'perturbation'
        close (perturbation_unit, iostat=local_status)
        call require(local_status == 0, 'could not write dispatch test perturbation')

        local_backend = gc_model2_backend_evidence_t()
        local_backend%factory_available = .true.
        local_backend%field_certified = .true.
        local_backend%profile_certified = .true.
        local_backend%perturbation_certified = .true.
        local_backend%perturbation_amplitude_convention = &
            GC_FULL_FOW_REAL_FIELD_AMPLITUDE_CONVENTION
        local_backend%perturbation_input_path = trim(perturbation_path)
        local_backend%perturbation_provenance_certified = .true.
        local_backend%phase_space_bound_method = GC_FULL_FOW_BOUND_METHOD
        local_backend%orbit_step_refinement_certified = .true.
        local_backend%orbit_base_step = 0.1_dp
        local_backend%orbit_refined_step = 0.05_dp
        local_backend%topology_certified = .true.
        local_backend%wall_certified = .true.
        local_backend%canonical_measure_certified = .true.
        local_backend%component_identity_certified = .true.
        local_backend%wall_actual_path = trim(local_wall)
        local_backend%wall_units = 'm'
        local_backend%wall_sha256 = repeat('a', 64)
    end subroutine configure_backend

    subroutine setenv_value(name, value, local_status)
        character(len=*), intent(in) :: name, value
        integer, intent(out) :: local_status
        character(kind=c_char), allocatable :: name_c(:), value_c(:)
        integer(c_int) :: code
        integer :: j, name_length, value_length

        name_length = len_trim(name)
        value_length = len_trim(value)
        allocate(name_c(name_length + 1), value_c(value_length + 1))
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
        deallocate(name_c, value_c)
    end subroutine setenv_value

    character(len=4096) function metadata_value(path, key) result(value)
        character(len=*), intent(in) :: path, key
        character(len=4096) :: line
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
            if (trim(line(:equal_position - 1)) == trim(key)) then
                value = trim(line(equal_position + 1:))
                exit
            end if
        end do
        close (unit)
    end function metadata_value

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop trim(message)
    end subroutine require

end program test_gc_model2_transport_dispatch
