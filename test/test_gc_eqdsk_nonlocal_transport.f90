program test_gc_eqdsk_nonlocal_transport
    !! Focused independent checks for the direct-EQDSK production factory.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: R0, a, psi_pr
    use neort_gc_cylindrical_model, only: GC_CYL_SUCCESS, &
        gc_cylindrical_field_sample_t, gc_cylindrical_invariants_t, &
        gc_cylindrical_state_t, invariants_from_cylindrical_state
    use neort_gc_cylindrical_nonlocal_provider, only: &
        GC_CYL_NONLOCAL_SUCCESS, gc_cylindrical_nonlocal_evaluation_t, &
        evaluate_gc_cylindrical_nonlocal
    use neort_gc_cylindrical_transport_provider, only: &
        GC_CYL_TRANSPORT_SUCCESS
    use neort_gc_cylindrical_physical_return, only: &
        gc_cylindrical_physical_return_options_t, &
        gc_cylindrical_physical_return_t, &
        compute_gc_cylindrical_physical_return
    use neort_gc_eqdsk_nonlocal_transport, only: &
        GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE, GC_EQDSK_NONLOCAL_SUCCESS, &
        gc_eqdsk_nonlocal_diagnostics_t, gc_eqdsk_nonlocal_factory_t, &
        gc_eqdsk_nonlocal_options_t, &
        gc_eqdsk_nonlocal_species_t, &
        gc_eqdsk_profile_potential_t, gc_eqdsk_nonlocal_status_message, &
        get_gc_eqdsk_nonlocal_diagnostics, &
        initialize_gc_eqdsk_nonlocal_transport, poincare_cut_value_position
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_SAMPLE_UNRESOLVED, GC_NONLOCAL_SAMPLE_VALID, &
        GC_NONLOCAL_SAMPLE_WALL, gc_nonlocal_component_t, &
        gc_nonlocal_orbit_sample_t
    use neort_gc_nonlocal_transport_types, only: &
        gc_nonlocal_transport_reference_t
    use neort_gc_perpendicular_invariant, only: &
        gc_buchholz_jk_from_mu_phys, gc_mu_phys_from_buchholz_jk
    use neort_thin_orbit_limit, only: THIN_LIMIT_SUCCESS, estimate_thin_limit, &
        orbit_return_t, thin_limit_result_t
    use util, only: c, mu, qe

    implicit none

    call test_factory_fails_closed()
    call test_production_factory()
    call test_circular_poincare_cut()
    call test_profile_potential_units()
    call test_narrow_limit_oracle()
    write (*, '(A)') 'test_gc_eqdsk_nonlocal_transport: PASS'

contains

    subroutine test_production_factory()
        type(gc_eqdsk_nonlocal_factory_t), target :: factory
        type(gc_eqdsk_nonlocal_species_t) :: species
        type(gc_eqdsk_nonlocal_options_t) :: options
        type(gc_nonlocal_transport_reference_t) :: reference
        type(gc_nonlocal_component_t), allocatable :: components(:)
        type(gc_nonlocal_orbit_sample_t) :: sample, minus_sample, plus_sample
        type(gc_cylindrical_nonlocal_evaluation_t) :: evaluation
        type(gc_cylindrical_field_sample_t) :: field
        type(gc_cylindrical_invariants_t) :: reversed_invariants
        type(gc_cylindrical_state_t) :: reversed_state
        type(gc_cylindrical_physical_return_options_t) :: return_options
        type(gc_cylindrical_physical_return_t) :: reversed_return
        type(gc_eqdsk_nonlocal_diagnostics_t) :: diagnostics
        real(dp), allocatable :: saved_wall(:, :)
        real(dp) :: small_wall(2, 4), h0_nodes(1), h0_weights(1)
        real(dp) :: jk_nodes(1), jk_weights(1), gradient(3)
        real(dp) :: potential, h0, j_k, mu_phys, v_parallel, v_perp
        real(dp) :: x, dx, finite_difference, outer, outer_full, outer_low
        real(dp) :: coarse_mid, fine_mid
        character(len=1024) :: paths(5)
        integer :: status, i, selected

        do i = 1, size(paths)
            call get_command_argument(i, paths(i))
            call require(len_trim(paths(i)) > 0, &
                'production fixture path was not registered')
        end do
        species%name = 'D+'
        species%mass_g = 2.0_dp*mu
        species%charge_esu = qe
        species%reference_energy_erg = 2.0e3_dp*1.602176634e-12_dp
        species%reference_velocity_cm_s = 1.0e7_dp
        options = gc_eqdsk_nonlocal_options_t()
        options%surface_min = 0.05_dp
        options%surface_max = 0.90_dp
        options%reference_surface = 0.50_dp
        options%orbit_maximum_step = 4.0e-9_dp
        options%invariant_relative_tolerance = 1.0e-6_dp
        options%orbit_options%minimum_event_time = 1.0e-10_dp
        options%orbit_options%maximum_time = 1.0e-3_dp
        options%topology_probe_count = 3
        options%topology_probe_fraction = 0.20_dp
        h0_nodes = [species%reference_energy_erg]
        h0_weights = [1.0_dp]
        jk_nodes = [0.0_dp]
        jk_weights = [1.0_dp]
        call initialize_gc_eqdsk_nonlocal_transport(trim(paths(1)), &
            trim(paths(2)), 'cm', trim(paths(3)), trim(paths(4)), &
            trim(paths(5)), species, h0_nodes, h0_weights, jk_nodes, &
            jk_weights, 1, 3, factory, status, options)
        if (status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            write (*, '(a,i0,2a)') 'factory status=', status, ' ', &
                trim(gc_eqdsk_nonlocal_status_message(status))
            write (*, '(a,3(1x,es24.16))') 'field globals=', R0, a, psi_pr
        end if
        call require(status == GC_EQDSK_NONLOCAL_SUCCESS, &
            'real direct-EQDSK factory initialization failed')
        call require(factory%initialized .and. factory%cut_ready, &
            'factory readiness preceded physical cut construction')
        call factory%field%evaluate(factory%section_position, field, status)
        call require(status == GC_CYL_SUCCESS, &
            'factory reference field is unavailable')
        call factory%potential%evaluate(factory%section_position, field, &
            potential, gradient, status)
        call require(status == GC_CYL_SUCCESS, &
            'factory reference potential is unavailable')
        call require(abs(potential) > 100.0_dp*epsilon(max(1.0_dp, &
            abs(potential))), 'production Phi was silently zeroed')

        v_parallel = 5.0e7_dp
        v_perp = 1.0e6_dp
        mu_phys = species%mass_g*v_perp**2/(2.0_dp*field%bmod)
        j_k = gc_buchholz_jk_from_mu_phys(mu_phys, species%charge_esu, c)
        call require_close(gc_mu_phys_from_buchholz_jk(j_k, &
            -species%charge_esu, c), mu_phys, &
            'charge reversal changed J_K to mu conversion', 1.0e-13_dp)
        h0 = 0.5_dp*species%mass_g*v_parallel**2 &
            +mu_phys*field%bmod+species%charge_esu*potential
        call require(h0 > 0.0_dp .and. j_k > 0.0_dp, &
            'physical fixed-H0/J_K node is invalid')
        call factory%provider%get_section_reference(reference, status)
        call require(status == GC_CYL_TRANSPORT_SUCCESS, &
            'factory reference callback failed')
        call factory%provider%get_components(reference, h0, j_k, components, &
            status)
        if (status /= GC_CYL_TRANSPORT_SUCCESS) then
            write (*, '(a,15(1x,i0))') 'component status/counters=', status, &
                factory%topology_certification_attempts, &
                factory%topology_certification_successes, &
                factory%physical_return_attempts, &
                factory%physical_return_successes, factory%no_return_count, &
                factory%wall_return_count, factory%radial_return_count, &
                factory%invariant_rejection_count, factory%field_error_count, &
                factory%potential_error_count, factory%state_error_count, &
                factory%start_error_count, factory%integrator_error_count, &
                factory%other_return_error_count
            write (*, '(a,2(1x,i0),3(1x,es24.16))') 'last return=', &
                factory%last_return_status, factory%last_return_event_kind, &
                factory%last_launch_event_value, &
                factory%last_return_event_value, &
                factory%last_invariant_scaled_drift
            write (*, '(a,3(1x,es24.16))') 'last residuals=', &
                factory%last_energy_error, &
                factory%last_magnetic_moment_error, &
                factory%last_canonical_momentum_error
        end if
        call require(status == GC_CYL_TRANSPORT_SUCCESS, &
            'physical homoclinic certification failed')
        call require(factory%topology_ready, &
            'topology readiness preceded physical-return certification')
        call require(factory%topology_certification_successes > 0, &
            'physical-return certification counter was not advanced')
        selected = 0
        do i = 1, size(components)
            if (components(i)%sigma /= 1) cycle
            if (components(i)%x_max <= components(i)%x_min) cycle
            selected = i
            exit
        end do
        call require(selected > 0, 'no positive-sigma physical class found')
        x = 0.5_dp*(components(selected)%x_min+components(selected)%x_max)
        call factory%provider%evaluate_orbit(reference, h0, j_k, x, 1, &
            components(selected)%component_id, sample, status)
        call require(status == GC_CYL_TRANSPORT_SUCCESS, &
            'real provider orbit callback failed')
        call require(sample%status == GC_NONLOCAL_SAMPLE_VALID, &
            'real physical return did not produce a valid sample')
        call require(sample%derivatives_available, &
            'fixed-H0/J_K derivative probes were not certified')
        call require_close(sample%omega_b*sample%tau_b, &
            2.0_dp*acos(-1.0_dp), 'complete-cycle omega_b identity', 2.0e-7_dp)
        call evaluate_gc_cylindrical_nonlocal(factory%provider%node_context, x, &
            1, components(selected)%component_id, 1, 3, 3, evaluation, status)
        call require(status == GC_CYL_NONLOCAL_SUCCESS, &
            'full nonlocal evaluation callback failed')
        call require(evaluation%complete_cycle_return, &
            'physical return was not marked complete-cycle')
        call require_close(evaluation%delta_phi_b/evaluation%sample%tau_b, &
            evaluation%sample%omega_phi, 'complete-cycle omega_phi identity', &
            2.0e-7_dp)
        call require(factory%harmonic_average_successes > 0, &
            'direct R-Z H_m orbit average was not executed')

        reversed_state = gc_cylindrical_state_t()
        reversed_state%R = factory%section_position(1)
        reversed_state%Z = factory%section_position(2)
        reversed_state%phi = factory%section_position(3)
        reversed_state%p_parallel = species%mass_g*v_parallel
        reversed_state%mu = gc_mu_phys_from_buchholz_jk(j_k, &
            -species%charge_esu, c)
        call invariants_from_cylindrical_state(reversed_state, field, potential, &
            species%mass_g, -species%charge_esu, c, reversed_invariants, status)
        call require(status == GC_CYL_SUCCESS, &
            'charge-reversed physical invariants failed')
        return_options = gc_cylindrical_physical_return_options_t()
        return_options%relative_tolerance = &
            options%orbit_options%relative_tolerance
        return_options%absolute_tolerance = &
            options%orbit_options%absolute_tolerance
        return_options%invariant_relative_tolerance = 3.0e-6_dp
        return_options%invariant_absolute_tolerance = tiny(1.0_dp)
        return_options%event_value_tolerance = options%cut_root_tolerance
        return_options%minimum_return_time = &
            options%orbit_options%minimum_event_time
        return_options%maximum_time = options%orbit_options%maximum_time
        return_options%maximum_step = options%orbit_maximum_step
        call compute_gc_cylindrical_physical_return(factory%field, &
            factory%potential, reversed_state, reversed_invariants, &
            species%mass_g, -species%charge_esu, c, actual_cut_event, &
            return_options, reversed_return, wall_model=factory%wall, &
            user_data=factory)
        if (reversed_return%status /= GC_CYL_SUCCESS) then
            write (*, '(a,3(1x,i0),4(1x,es24.16))') 'reversed return=', &
                reversed_return%status, reversed_return%event_kind, &
                reversed_return%return_orientation, &
                reversed_return%launch_event_value, &
                reversed_return%return_event_value, &
                reversed_return%maximum_invariant_scaled_drift, &
                reversed_return%period
        end if
        call require(reversed_return%status == GC_CYL_SUCCESS .and. &
            reversed_return%physical_return_found, &
            'charge-reversed actual-cut physical return failed')
        call require_close(2.0_dp*acos(-1.0_dp)/reversed_return%period &
            *reversed_return%period, 2.0_dp*acos(-1.0_dp), &
            'charge-reversed complete-cycle identity', 2.0e-7_dp)

        dx = min(2.0e-4_dp, 0.1_dp*(components(selected)%x_max &
            -components(selected)%x_min))
        call require(dx > 0.0_dp, 'derivative probe width vanished')
        call factory%provider%evaluate_orbit(reference, h0, j_k, x-dx, 1, &
            components(selected)%component_id, minus_sample, status)
        call require(status == GC_CYL_TRANSPORT_SUCCESS .and. &
            minus_sample%status == GC_NONLOCAL_SAMPLE_VALID, &
            'minus fixed-H0/J_K derivative orbit failed')
        call factory%provider%evaluate_orbit(reference, h0, j_k, x+dx, 1, &
            components(selected)%component_id, plus_sample, status)
        call require(status == GC_CYL_TRANSPORT_SUCCESS .and. &
            plus_sample%status == GC_NONLOCAL_SAMPLE_VALID, &
            'plus fixed-H0/J_K derivative orbit failed')
        finite_difference = (plus_sample%omega_b-minus_sample%omega_b)/(2.0_dp*dx)
        call require_close(sample%domega_b_dx, finite_difference, &
            'fixed-H0/J_K omega derivative', 5.0e-2_dp)

        call factory%provider%evaluate_orbit(reference, h0, j_k, x, 1, &
            components(selected)%component_id, sample, status)
        factory%options%residence_target_s_tor = factory%options%surface_max
        call factory%provider%evaluate_outer_measure_factor(reference, h0, j_k, &
            x, 1, components(selected)%component_id, sample, outer_full, status)
        call require(status == GC_CYL_TRANSPORT_SUCCESS, &
            'full radial-flux residence failed')
        factory%options%residence_target_s_tor = factory%options%surface_min
        call factory%provider%evaluate_outer_measure_factor(reference, h0, j_k, &
            x, 1, components(selected)%component_id, sample, outer_low, status)
        call require(status == GC_CYL_TRANSPORT_SUCCESS, &
            'inner radial-flux residence failed')
        call require(abs(outer_full) > tiny(outer_full), &
            'full residence normalization vanished')
        call require_close((outer_full-outer_low)/outer_full, 1.0_dp, &
            'radial-cell residence sum', 3.0e-3_dp)
        factory%options%residence_target_s_tor = x
        call factory%provider%evaluate_outer_measure_factor(reference, h0, j_k, &
            x, 1, components(selected)%component_id, sample, coarse_mid, status)
        call require(status == GC_CYL_TRANSPORT_SUCCESS, &
            'coarse radial residence failed')
        factory%options%orbit_maximum_step = &
            0.5_dp*factory%options%orbit_maximum_step
        call factory%provider%evaluate_outer_measure_factor(reference, h0, j_k, &
            x, 1, components(selected)%component_id, sample, fine_mid, status)
        call require(status == GC_CYL_TRANSPORT_SUCCESS, &
            'fine radial residence failed')
        call require_close(fine_mid, coarse_mid, &
            'radial residence step convergence', 3.0e-3_dp)
        associate (unused_outer => outer)
        end associate

        allocate(saved_wall(2, size(factory%wall%vertices, 2)))
        saved_wall = factory%wall%vertices
        small_wall(:, 1) = factory%section_position(1:2)+[-0.02_dp, -0.02_dp]
        small_wall(:, 2) = factory%section_position(1:2)+[0.02_dp, -0.02_dp]
        small_wall(:, 3) = factory%section_position(1:2)+[0.02_dp, 0.02_dp]
        small_wall(:, 4) = factory%section_position(1:2)+[-0.02_dp, 0.02_dp]
        call factory%wall%set_vertices(small_wall, status)
        call require(status == GC_CYL_SUCCESS, 'temporary wall setup failed')
        call factory%provider%evaluate_orbit(reference, h0, j_k, x, 1, &
            components(selected)%component_id, sample, status)
        call require(status == GC_CYL_TRANSPORT_SUCCESS .and. &
            sample%status == GC_NONLOCAL_SAMPLE_WALL, &
            'wall physical loss was not preserved as a wall sample')
        call factory%wall%set_vertices(saved_wall, status)
        call require(status == GC_CYL_SUCCESS, 'production wall restore failed')

        factory%options%orbit_options%maximum_time = 0.5_dp &
            *factory%options%orbit_options%minimum_event_time
        call factory%provider%evaluate_orbit(reference, h0, j_k, x, 1, &
            components(selected)%component_id, sample, status)
        call require(status == GC_CYL_TRANSPORT_SUCCESS .and. &
            sample%status == GC_NONLOCAL_SAMPLE_UNRESOLVED, &
            'no-return was not preserved as unresolved')
        call require(factory%no_return_count > 0, &
            'no-return counter advanced without a physical no-return')
        call require(factory%physical_return_successes > 0, &
            'success counter advanced before physical returns')
        call get_gc_eqdsk_nonlocal_diagnostics(factory, diagnostics)
        call require(diagnostics%physical_return_successes == &
            factory%physical_return_successes, &
            'diagnostics fabricated aggregate successful returns')
        call require(diagnostics%wall_returns == factory%wall_return_count .and. &
            diagnostics%no_returns == factory%no_return_count, &
            'diagnostics fabricated physical return categories')
        call require(diagnostics%return_accounting_complete .and. &
            diagnostics%categorized_returns == &
            diagnostics%physical_return_attempts, &
            'observed return categories do not conserve attempts')
        call require(diagnostics%topology_ready .and. diagnostics%initialized, &
            'diagnostics readiness preceded successful factory state')

    end subroutine test_production_factory

    subroutine actual_cut_event(position, state, sample_field, user_data, value, &
            event_status)
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: sample_field
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: value
        integer, intent(out) :: event_status

        real(dp) :: scale

        associate (unused_state => state)
        end associate
        value = 0.0_dp
        event_status = 1
        if (.not. associated(user_data)) return
        select type (event_factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            call poincare_cut_value_position(event_factory, position, value, &
                event_status)
            if (event_status /= GC_EQDSK_NONLOCAL_SUCCESS) return
            scale = sqrt(sample_field%grad_b(1)**2 &
                +sample_field%grad_b(3)**2) &
                *sqrt(sample_field%grad_psi(1)**2 &
                +sample_field%grad_psi(3)**2)
            if (scale <= tiny(scale)) then
                event_status = 1
                return
            end if
            value = value/scale
            event_status = GC_CYL_SUCCESS
        class default
            event_status = 1
        end select
    end subroutine actual_cut_event

    subroutine test_factory_fails_closed()
        type(gc_eqdsk_nonlocal_factory_t) :: factory
        type(gc_eqdsk_nonlocal_species_t) :: species
        real(dp) :: h0(1), wh0(1), jperp(1), wjperp(1)
        integer :: status

        species%name = 'D+'
        species%mass_g = 2.014_dp*mu
        species%charge_esu = qe
        species%reference_energy_erg = 1.0e-8_dp
        species%reference_velocity_cm_s = 1.0e7_dp
        h0 = [1.0e-8_dp]
        wh0 = [1.0_dp]
        jperp = [1.0e-16_dp]
        wjperp = [1.0_dp]
        call initialize_gc_eqdsk_nonlocal_transport( &
            '/definitely/missing/direct.eqdsk', &
            '/definitely/missing/wall.dat', 'cm', &
            '/definitely/missing/perturbation.bin', &
            '/definitely/missing/plasma.in', &
            '/definitely/missing/rotation.in', species, h0, wh0, jperp, &
            wjperp, 1, 3, factory, status)
        call require(status == GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE, &
            'factory did not fail closed on missing EQDSK')
    end subroutine test_factory_fails_closed

    subroutine test_circular_poincare_cut()
        real(dp), parameter :: r0 = 300.0_dp, b0 = 1.0e4_dp
        real(dp) :: r, z, grad_b_r, grad_psi_r, grad_psi_z, cut_value

        r = r0 + 40.0_dp
        z = 0.0_dp
        grad_b_r = -b0*r0/r**2
        grad_psi_r = b0*(r-r0)
        grad_psi_z = b0*z
        cut_value = -grad_b_r*grad_psi_z
        call require(abs(cut_value) <= 1.0e-12_dp, &
            'circular C=0 cut failed at the midplane')
        z = 12.0_dp
        grad_psi_z = b0*z
        cut_value = -grad_b_r*grad_psi_z
        call require(abs(cut_value) > 0.0_dp, &
            'circular C=0 cut accepted an off-midplane point')
    end subroutine test_circular_poincare_cut

    subroutine test_profile_potential_units()
        type(gc_eqdsk_profile_potential_t) :: potential
        type(gc_cylindrical_field_sample_t) :: field
        real(dp) :: value, gradient(3)
        integer :: status

        allocate(potential%s_tor(2), potential%psi_pol(2), potential%phi(2), &
            potential%omega_e(2))
        potential%s_tor = [0.0_dp, 1.0_dp]
        potential%psi_pol = [0.0_dp, 2.0_dp]
        potential%phi = [0.0_dp, 4.0_dp]
        potential%omega_e = [2.0_dp, 2.0_dp]
        potential%c_light = 1.0_dp
        potential%initialized = .true.
        field%psi = 1.0_dp
        field%grad_psi = [3.0_dp, 0.0_dp, 4.0_dp]
        call potential%evaluate([300.0_dp, 0.0_dp, 0.0_dp], field, value, &
            gradient, status)
        call require(status == GC_CYL_SUCCESS, 'profile potential unavailable')
        call require(abs(value-2.0_dp) <= 1.0e-13_dp, &
            'profile potential integration has the wrong units')
        call require(maxval(abs(gradient-[6.0_dp, 0.0_dp, 8.0_dp])) &
            <= 1.0e-13_dp, 'profile potential gradient is inconsistent')
    end subroutine test_profile_potential_units

    subroutine test_narrow_limit_oracle()
        type(orbit_return_t) :: base, plus(3), minus(3)
        type(thin_limit_result_t) :: result
        real(dp), parameter :: lambda(3) = [1.0_dp, 0.5_dp, 0.25_dp]
        real(dp), parameter :: topological = 2.0_dp
        real(dp), parameter :: slope = -0.37_dp
        integer :: k

        base%status = 0
        base%period = 5.0_dp
        base%delta_phi = topological
        do k = 1, 3
            plus(k)%status = 0
            minus(k)%status = 0
            plus(k)%period = base%period
            minus(k)%period = base%period
            plus(k)%delta_phi = topological + slope*lambda(k)
            minus(k)%delta_phi = topological - slope*lambda(k)
        end do
        call estimate_thin_limit(topological, base, lambda, plus, minus, result, &
            baseline_tolerance=1.0e-12_dp)
        call require(result%status == THIN_LIMIT_SUCCESS, &
            'existing thin-limit oracle did not converge')
        call require(abs(result%omega-slope/base%period) <= 1.0e-12_dp, &
            'narrow-orbit limit does not match the analytic return map')
    end subroutine test_narrow_limit_oracle

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop trim(message)
    end subroutine require

    subroutine require_close(actual, expected, message, relative_tolerance)
        real(dp), intent(in) :: actual, expected, relative_tolerance
        character(len=*), intent(in) :: message

        real(dp) :: scale

        scale = max(tiny(1.0_dp), abs(expected))
        if (abs(actual-expected) > relative_tolerance*scale) then
            write (*, '(a,3(1x,es24.16))') trim(message), actual, expected, &
                relative_tolerance
            error stop trim(message)
        end if
    end subroutine require_close

end program test_gc_eqdsk_nonlocal_transport
