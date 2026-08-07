program test_gc_cylindrical_nonlocal_provider
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_SAMPLE_UNRESOLVED, GC_NONLOCAL_SAMPLE_VALID, &
        GC_NONLOCAL_SAMPLE_WALL, gc_nonlocal_component_t, &
        gc_nonlocal_orbit_sample_t, gc_nonlocal_resonance_options_t, &
        gc_nonlocal_resonance_result_t
    use neort_gc_cylindrical_nonlocal_provider, only: &
        GC_CYL_NONLOCAL_COMPONENT_IDENTITY, &
        GC_CYL_NONLOCAL_FREQUENCY_REPRESENTATION, &
        GC_CYL_NONLOCAL_INTEGRAL_UNAVAILABLE, &
        GC_CYL_NONLOCAL_NORMALIZATION_UNAVAILABLE, &
        GC_CYL_NONLOCAL_ORBIT_UNRESOLVED, GC_CYL_NONLOCAL_ORBIT_VALID, &
        GC_CYL_NONLOCAL_ORBIT_WALL, GC_CYL_NONLOCAL_REFERENCE_MISMATCH, &
        GC_CYL_NONLOCAL_SUCCESS, GC_CYL_NONLOCAL_WALL_CLEAR, &
        GC_CYL_NONLOCAL_WALL_HIT, gc_cylindrical_nonlocal_context_t, &
        gc_cylindrical_nonlocal_evaluation_t, &
        gc_cylindrical_nonlocal_orbit_t, &
        gc_cylindrical_nonlocal_resonance_residuals, &
        gc_cylindrical_nonlocal_resonance_weights, &
        initialize_gc_cylindrical_nonlocal_provider, &
        enumerate_gc_cylindrical_nonlocal_components, &
        evaluate_gc_cylindrical_nonlocal, &
        integrate_gc_cylindrical_nonlocal_resonance
    implicit none

    real(dp), parameter :: H0_REFERENCE = 17.25_dp
    real(dp), parameter :: JPERP_REFERENCE = 2.75_dp
    real(dp), parameter :: PARTICLE_CHARGE = 2.0_dp
    real(dp), parameter :: C_LIGHT = 10.0_dp
    real(dp), parameter :: SECTION_REFERENCE(3) = [1.7_dp, 0.2_dp, 0.4_dp]
    real(dp), parameter :: X_ROOT = 0.2_dp
    real(dp), parameter :: TWO_PI = 6.28318530717958647692528676656_dp

    type :: manufactured_data_t
        real(dp) :: expected_h0 = H0_REFERENCE
        real(dp) :: expected_jperp = JPERP_REFERENCE
        real(dp) :: expected_reference(3) = SECTION_REFERENCE
        integer :: invariant_calls = 0
        logical :: invariant_arguments_exact = .true.
        logical :: bad_identity = .false.
        logical :: bad_reference = .false.
        logical :: wall_mode = .false.
        logical :: unresolved_mode = .false.
    end type manufactured_data_t

    type(gc_cylindrical_nonlocal_context_t) :: context
    type(gc_cylindrical_nonlocal_context_t) :: conversion_context
    type(gc_cylindrical_nonlocal_evaluation_t) :: evaluation
    type(gc_cylindrical_nonlocal_evaluation_t) :: evaluation_plus
    type(gc_cylindrical_nonlocal_evaluation_t) :: evaluation_minus
    type(gc_nonlocal_component_t), allocatable :: components(:)
    type(gc_nonlocal_resonance_options_t) :: options
    type(gc_nonlocal_resonance_result_t) :: result
    type(manufactured_data_t), target :: data
    integer :: status
    real(dp) :: phase_residual, phase_derivative
    real(dp) :: frequency_residual, frequency_derivative
    real(dp) :: phase_weight, frequency_weight
    integer :: frequency_prefactor
    real(dp) :: step, central_derivative

    call initialize_provider(context, data, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'full manufactured provider initialization failed')
    call enumerate_gc_cylindrical_nonlocal_components(context, components, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'disconnected component enumeration failed')
    call require(size(components) == 3, &
        'the three manufactured Poincare components were not retained')
    call require(components(1)%sigma == 1 .and. components(2)%sigma == -1, &
        'sigma identity was not retained')
    call require(abs(components(1)%x_min - components(2)%x_min) <= 0.0_dp .and. &
        abs(components(1)%x_max - components(2)%x_max) <= 0.0_dp, &
        'opposite-sigma branches were incorrectly treated as overlapping')
    call require(components(3)%sigma == 1 .and. components(3)%x_min > 1.0_dp, &
        'disconnected same-sigma component was not retained')

    call evaluate_gc_cylindrical_nonlocal(context, X_ROOT, 1, 11, 2, 3, 2, &
        evaluation, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'valid manufactured cylindrical orbit was rejected')
    call require(evaluation%sample%status == GC_NONLOCAL_SAMPLE_VALID, &
        'valid orbit did not produce a valid nonlocal sample')
    call require(evaluation%winding == 2, 'complete-cycle winding was lost')
    call require(evaluation%section_return_crossings == 2, &
        'section return multiplicity was lost')
    call require(evaluation%complete_cycle_return, &
        'frequency was not certified as a complete-cycle return frequency')
    call require(evaluation%frequency_representation == &
        GC_CYL_NONLOCAL_FREQUENCY_REPRESENTATION, &
        'frequency representation was not declared')
    call require(evaluation%electric_potential_included, &
        'electric-potential inclusion was not locked in orbit metadata')
    call require(evaluation%p_phi_mapping_certified, &
        'one-to-one P_phi section mapping was not certified')
    call require(evaluation%section%locked, 'section reference was not locked')
    call require(trim(evaluation%section%reference_id) == 'manufactured-section', &
        'section reference identity was not preserved')
    call require(evaluation%canonical_normalization_certified, &
        'canonical normalization was not certified')
    call require(trim(evaluation%psi_star_units) == 'psi_star=(c/q)*p_phi', &
        'canonical normalization units were not recorded')
    call require_close(evaluation%sample%psi_star, 5.0_dp*evaluation%p_phi, &
        'psi_star was not normalized from p_phi exactly once', 1.0e-13_dp)
    call require_close(evaluation%sample%dpsi_star_dx, &
        5.0_dp*evaluation%dp_phi_dx, &
        'dpsi_star/dx was not normalized from dp_phi/dx exactly once', 1.0e-13_dp)
    call require_close(evaluation%sample%omega_phi, &
        evaluation%delta_phi_b/evaluation%sample%tau_b, &
        'electric frequency was shifted independently of the full orbit', 1.0e-12_dp)

    step = 1.0e-5_dp
    call evaluate_gc_cylindrical_nonlocal(context, X_ROOT + step, 1, 11, 2, 3, 2, &
        evaluation_plus, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'positive derivative sample failed')
    call evaluate_gc_cylindrical_nonlocal(context, X_ROOT - step, 1, 11, 2, 3, 2, &
        evaluation_minus, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'negative derivative sample failed')
    central_derivative = (evaluation_plus%p_phi - evaluation_minus%p_phi) &
        /(2.0_dp*step)
    call require_close(central_derivative, evaluation%dp_phi_dx, &
        'P_phi derivative is not section-symmetric', 1.0e-8_dp)
    central_derivative = (evaluation_plus%sample%psi_star &
        -evaluation_minus%sample%psi_star)/(2.0_dp*step)
    call require_close(central_derivative, evaluation%sample%dpsi_star_dx, &
        'psi_star derivative is not section-symmetric', 1.0e-8_dp)
    central_derivative = (evaluation_plus%sample%omega_b &
        -evaluation_minus%sample%omega_b)/(2.0_dp*step)
    call require_close(central_derivative, evaluation%sample%domega_b_dx, &
        'complete-cycle omega_b derivative is not section-symmetric', 1.0e-8_dp)
    central_derivative = (evaluation_plus%delta_phi_b &
        -evaluation_minus%delta_phi_b)/(2.0_dp*step)
    call require_close(central_derivative, evaluation%ddelta_phi_b_dx, &
        'phase-advance derivative is not section-symmetric', 1.0e-8_dp)
    call require(data%invariant_arguments_exact, &
        'a callback did not receive fixed H0 and Jperp unchanged')
    call require(data%invariant_calls >= 4, &
        'fixed invariant callback contract was not exercised')

    call gc_cylindrical_nonlocal_resonance_residuals(evaluation, 2, 3, &
        phase_residual, phase_derivative, frequency_residual, frequency_derivative, &
        status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'phase/frequency residual conversion failed')
    call require_close(phase_residual, 0.0_dp, 'manufactured phase root', &
        1.0e-12_dp)
    call require_close(frequency_residual, 0.0_dp, 'manufactured frequency root', &
        1.0e-12_dp)
    call require_close(frequency_derivative, &
        3.0_dp/evaluation%sample%tau_b*phase_derivative, &
        'phase/frequency derivative equivalence', 1.0e-12_dp)
    call gc_cylindrical_nonlocal_resonance_weights(evaluation, 2, 3, &
        phase_weight, frequency_weight, frequency_prefactor, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'phase/frequency weight conversion failed')
    call require(frequency_prefactor == 9, &
        'frequency representation did not expose n squared')
    call require_close(phase_weight, 9.0_dp*frequency_weight, &
        'phase and frequency delta weights disagree', 1.0e-11_dp)
    call require_close(frequency_weight, &
        evaluation%sample%tau_b/abs(frequency_derivative), &
        'frequency inner weight contains more than one tau_b', 1.0e-12_dp)

    data%wall_mode = .true.
    call evaluate_gc_cylindrical_nonlocal(context, 0.75_dp, -1, 22, 2, 3, 2, &
        evaluation, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'wall status was turned into a callback failure')
    call require(evaluation%sample%status == GC_NONLOCAL_SAMPLE_WALL, &
        'wall identity was not preserved')
    call require(evaluation%wall_status == GC_CYL_NONLOCAL_WALL_HIT, &
        'wall hit status was not preserved')
    data%wall_mode = .false.

    data%unresolved_mode = .true.
    call evaluate_gc_cylindrical_nonlocal(context, 2.5_dp, 1, 33, 2, 3, 2, &
        evaluation, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'unresolved topology was turned into a callback failure')
    call require(evaluation%sample%status == GC_NONLOCAL_SAMPLE_UNRESOLVED, &
        'unresolved topology was silently accepted')
    data%unresolved_mode = .false.

    data%bad_identity = .true.
    call evaluate_gc_cylindrical_nonlocal(context, X_ROOT, 1, 11, 2, 3, 2, &
        evaluation, status)
    call require(status == GC_CYL_NONLOCAL_COMPONENT_IDENTITY, &
        'component identity mismatch was accepted')
    data%bad_identity = .false.

    data%bad_reference = .true.
    call evaluate_gc_cylindrical_nonlocal(context, X_ROOT, 1, 11, 2, 3, 2, &
        evaluation, status)
    call require(status == GC_CYL_NONLOCAL_REFERENCE_MISMATCH, &
        'section reference mismatch was accepted')
    data%bad_reference = .false.

    call initialize_conversion_provider(conversion_context, data, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'certified conversion provider initialization failed')
    call enumerate_gc_cylindrical_nonlocal_components(conversion_context, &
        components, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'conversion-provider component enumeration failed')
    call evaluate_gc_cylindrical_nonlocal(conversion_context, X_ROOT, 1, 11, 2, 3, &
        2, evaluation, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'certified conversion callback was not accepted')
    call require(trim(evaluation%psi_star_units) == 'manufactured_psi_star', &
        'certified conversion units were not retained')
    call require_close(evaluation%sample%psi_star, 7.0_dp*evaluation%p_phi, &
        'certified conversion callback was not applied exactly once', 1.0e-13_dp)

    call initialize_missing_normalization(context, data, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'missing-normalization provider initialization failed')
    call enumerate_gc_cylindrical_nonlocal_components(context, components, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'missing-normalization enumeration failed')
    call evaluate_gc_cylindrical_nonlocal(context, X_ROOT, 1, 11, 2, 3, 2, &
        evaluation, status)
    call require(status == GC_CYL_NONLOCAL_NORMALIZATION_UNAVAILABLE, &
        'raw p_phi was accepted without a certified normalization')

    call initialize_provider(context, data, status, omit_harmonic=.true.)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'missing-harmonic provider initialization failed')
    call enumerate_gc_cylindrical_nonlocal_components(context, components, status)
    call require(status == GC_CYL_NONLOCAL_SUCCESS, &
        'missing-harmonic enumeration failed')
    call evaluate_gc_cylindrical_nonlocal(context, X_ROOT, 1, 11, 2, 3, 2, &
        evaluation, status)
    call require(status /= GC_CYL_NONLOCAL_SUCCESS, &
        'missing phase-averaged H_m was silently fabricated')

    call initialize_provider(context, data, status)
    options = gc_nonlocal_resonance_options_t()
    call integrate_gc_cylindrical_nonlocal_resonance(context, 2, 3, options, &
        result, status)
    call require(status == GC_CYL_NONLOCAL_INTEGRAL_UNAVAILABLE, &
        'bounded provider integration seam was accidentally enabled')
    call require(.not. result%certified, &
        'uncorrected inner resonance representation was certified')

    write (*, '(a)') 'test_gc_cylindrical_nonlocal_provider OK'

contains

    subroutine initialize_provider(local_context, local_data, local_status, &
            omit_harmonic)
        type(gc_cylindrical_nonlocal_context_t), intent(out) :: local_context
        type(manufactured_data_t), target, intent(inout) :: local_data
        integer, intent(out) :: local_status
        logical, intent(in), optional :: omit_harmonic
        logical :: skip_harmonic

        skip_harmonic = .false.
        if (present(omit_harmonic)) skip_harmonic = omit_harmonic
        if (skip_harmonic) then
            call initialize_gc_cylindrical_nonlocal_provider(H0_REFERENCE, &
                JPERP_REFERENCE, local_context, local_status, &
                particle_charge=PARTICLE_CHARGE, c_light=C_LIGHT, &
                component_provider=manufactured_components, &
                orbit_provider=manufactured_orbit, &
                force_provider=manufactured_force, section_coordinate='R_c', &
                section_reference=SECTION_REFERENCE, &
                section_reference_id='manufactured-section', user_data=local_data)
        else
            call initialize_gc_cylindrical_nonlocal_provider(H0_REFERENCE, &
                JPERP_REFERENCE, local_context, local_status, &
                particle_charge=PARTICLE_CHARGE, c_light=C_LIGHT, &
                component_provider=manufactured_components, &
                orbit_provider=manufactured_orbit, &
                harmonic_provider=manufactured_harmonic, &
                force_provider=manufactured_force, section_coordinate='R_c', &
                section_reference=SECTION_REFERENCE, &
                section_reference_id='manufactured-section', user_data=local_data)
        end if
    end subroutine initialize_provider

    subroutine initialize_conversion_provider(local_context, local_data, local_status)
        type(gc_cylindrical_nonlocal_context_t), intent(out) :: local_context
        type(manufactured_data_t), target, intent(inout) :: local_data
        integer, intent(out) :: local_status

        call initialize_gc_cylindrical_nonlocal_provider(H0_REFERENCE, &
            JPERP_REFERENCE, local_context, local_status, &
            component_provider=manufactured_components, &
            orbit_provider=manufactured_orbit, &
            harmonic_provider=manufactured_harmonic, &
            force_provider=manufactured_force, &
            canonical_conversion_provider=manufactured_conversion, &
            section_coordinate='R_c', section_reference=SECTION_REFERENCE, &
            section_reference_id='manufactured-section', user_data=local_data)
    end subroutine initialize_conversion_provider

    subroutine initialize_missing_normalization(local_context, local_data, local_status)
        type(gc_cylindrical_nonlocal_context_t), intent(out) :: local_context
        type(manufactured_data_t), target, intent(inout) :: local_data
        integer, intent(out) :: local_status

        call initialize_gc_cylindrical_nonlocal_provider(H0_REFERENCE, &
            JPERP_REFERENCE, local_context, local_status, &
            component_provider=manufactured_components, &
            orbit_provider=manufactured_orbit, &
            harmonic_provider=manufactured_harmonic, &
            force_provider=manufactured_force, &
            section_coordinate='R_c', section_reference=SECTION_REFERENCE, &
            section_reference_id='manufactured-section', user_data=local_data)
    end subroutine initialize_missing_normalization

    subroutine manufactured_components(h0, jperp, user_data, components, status)
        real(dp), intent(in) :: h0, jperp
        class(*), pointer, intent(inout) :: user_data
        type(gc_nonlocal_component_t), allocatable, intent(out) :: components(:)
        integer, intent(out) :: status

        type(manufactured_data_t), pointer :: state

        status = 1
        if (.not. associated(user_data)) return
        select type (user_data)
            type is (manufactured_data_t)
            state => user_data
            call record_invariants(state, h0, jperp)
            allocate(components(3))
            components(1) = gc_nonlocal_component_t(11, 1, -1.0_dp, 1.0_dp)
            components(2) = gc_nonlocal_component_t(22, -1, -1.0_dp, 1.0_dp)
            components(3) = gc_nonlocal_component_t(33, 1, 2.0_dp, 3.0_dp)
            status = 0
        class default
            nullify(state)
        end select
    end subroutine manufactured_components

    subroutine manufactured_orbit(h0, jperp, x, sigma, component_id, user_data, &
            orbit, status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        class(*), pointer, intent(inout) :: user_data
        type(gc_cylindrical_nonlocal_orbit_t), intent(out) :: orbit
        integer, intent(out) :: status

        type(manufactured_data_t), pointer :: state

        orbit = gc_cylindrical_nonlocal_orbit_t()
        status = 1
        if (.not. associated(user_data)) return
        select type (user_data)
            type is (manufactured_data_t)
            state => user_data
            call record_invariants(state, h0, jperp)
            call fill_orbit(state, x, sigma, component_id, orbit)
            status = 0
        class default
            nullify(state)
        end select
    end subroutine manufactured_orbit

    subroutine manufactured_harmonic(h0, jperp, x, sigma, component_id, harmonic_m, &
            harmonic_n, orbit, user_data, h_m, status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id, harmonic_m, harmonic_n
        type(gc_cylindrical_nonlocal_orbit_t), intent(in) :: orbit
        class(*), pointer, intent(inout) :: user_data
        complex(dp), intent(out) :: h_m
        integer, intent(out) :: status

        type(manufactured_data_t), pointer :: state

        h_m = cmplx(0.0_dp, 0.0_dp, kind=dp)
        status = 1
        if (.not. associated(user_data)) return
        select type (user_data)
            type is (manufactured_data_t)
            state => user_data
            call record_invariants(state, h0, jperp)
            if (orbit%component_id /= component_id) return
            if (orbit%sigma /= sigma) return
            h_m = cmplx(1.25_dp + 0.05_dp*x, -0.75_dp + &
                0.1_dp*real(sigma + harmonic_m - harmonic_n, dp), kind=dp)
            status = 0
        class default
            nullify(state)
        end select
    end subroutine manufactured_harmonic

    subroutine manufactured_force(h0, jperp, x, sigma, component_id, orbit, &
            user_data, force, status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_cylindrical_nonlocal_orbit_t), intent(in) :: orbit
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: force(:)
        integer, intent(out) :: status

        type(manufactured_data_t), pointer :: state

        force = 0.0_dp
        status = 1
        if (.not. associated(user_data)) return
        select type (user_data)
            type is (manufactured_data_t)
            state => user_data
            call record_invariants(state, h0, jperp)
            if (orbit%component_id /= component_id) return
            if (orbit%sigma /= sigma) return
            if (size(force) >= 1) force(1) = 2.0_dp*x + real(sigma, dp)
            if (size(force) >= 2) force(2) = -0.5_dp*x
            status = 0
        class default
            nullify(state)
        end select
    end subroutine manufactured_force

    subroutine manufactured_conversion(p_phi, dp_phi_dx, user_data, psi_star, &
            dpsi_star_dx, units, certified, status)
        real(dp), intent(in) :: p_phi, dp_phi_dx
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: psi_star, dpsi_star_dx
        character(len=*), intent(out) :: units
        logical, intent(out) :: certified
        integer, intent(out) :: status

        associate (unused_user_data => user_data)
        end associate
        psi_star = 7.0_dp*p_phi
        dpsi_star_dx = 7.0_dp*dp_phi_dx
        units = 'manufactured_psi_star'
        certified = .true.
        status = 0
    end subroutine manufactured_conversion

    subroutine fill_orbit(state, x, sigma, component_id, orbit)
        type(manufactured_data_t), intent(in) :: state
        real(dp), intent(in) :: x
        integer, intent(in) :: sigma, component_id
        type(gc_cylindrical_nonlocal_orbit_t), intent(out) :: orbit

        real(dp) :: orientation, tau, dtau, delta, ddelta

        orbit = gc_cylindrical_nonlocal_orbit_t()
        orientation = real(sigma, dp)
        tau = 4.0_dp + 0.1_dp*(x - X_ROOT)
        dtau = 0.1_dp
        delta = -TWO_PI*2.0_dp/3.0_dp + 0.3_dp*(x - X_ROOT) &
            +0.02_dp*(x - X_ROOT)**2
        ddelta = 0.3_dp + 0.04_dp*(x - X_ROOT)
        orbit%status = GC_CYL_NONLOCAL_ORBIT_VALID
        orbit%component_id = component_id
        orbit%sigma = sigma
        orbit%winding = 2*sigma
        orbit%section_return_crossings = 2
        orbit%winding_available = .true.
        orbit%section_return_available = .true.
        orbit%complete_cycle_return = .true.
        orbit%p_phi_mapping_certified = .true.
        orbit%mapping_orientation = sigma
        orbit%electric_potential_included = .true.
        orbit%p_phi = 10.0_dp + 0.2_dp*real(component_id, dp) &
            +orientation*(0.7_dp*x + 0.1_dp*x**2)
        orbit%dp_phi_dx = orientation*(0.7_dp + 0.2_dp*x)
        orbit%tau_b = tau
        orbit%dtau_b_dx = dtau
        orbit%delta_phi_b = delta
        orbit%ddelta_phi_b_dx = ddelta
        orbit%phase_advance_available = .true.
        orbit%period_derivative_available = .true.
        orbit%omega_b = TWO_PI/tau
        orbit%omega_phi = delta/tau
        orbit%domega_b_dx = -TWO_PI*dtau/tau**2
        orbit%domega_phi_dx = (ddelta*tau - delta*dtau)/tau**2
        orbit%derivatives_available = .true.
        orbit%wall_checked = .true.
        orbit%wall_status = GC_CYL_NONLOCAL_WALL_CLEAR
        orbit%section%coordinate = 'R_c'
        orbit%section%reference = state%expected_reference
        orbit%section%reference_id = 'manufactured-section'
        orbit%section%locked = .true.
        if (state%bad_identity) orbit%component_id = component_id + 100
        if (state%bad_reference) orbit%section%reference(1) = &
            orbit%section%reference(1) + 0.01_dp
        if (state%wall_mode) then
            if (component_id == 22 .and. x > 0.5_dp) then
                orbit%status = GC_CYL_NONLOCAL_ORBIT_WALL
                orbit%wall_status = GC_CYL_NONLOCAL_WALL_HIT
            end if
        end if
        if (state%unresolved_mode) then
            if (component_id == 33) orbit%status = GC_CYL_NONLOCAL_ORBIT_UNRESOLVED
        end if
    end subroutine fill_orbit

    subroutine record_invariants(state, h0, jperp)
        type(manufactured_data_t), intent(inout) :: state
        real(dp), intent(in) :: h0, jperp

        state%invariant_calls = state%invariant_calls + 1
        if (abs(h0 - state%expected_h0) > 0.0_dp) then
            state%invariant_arguments_exact = .false.
        end if
        if (abs(jperp - state%expected_jperp) > 0.0_dp) then
            state%invariant_arguments_exact = .false.
        end if
    end subroutine record_invariants

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) then
            write (*, '(a)') trim(message)
            error stop 1
        end if
    end subroutine require

    subroutine require_close(value, reference, message, tolerance)
        real(dp), intent(in) :: value, reference, tolerance
        character(len=*), intent(in) :: message

        if (abs(value - reference) > tolerance*max(1.0_dp, &
            max(abs(value), abs(reference)))) then
            write (*, '(a,2(1x,es24.16))') trim(message), value, reference
            error stop 1
        end if
    end subroutine require_close

end program test_gc_cylindrical_nonlocal_provider
