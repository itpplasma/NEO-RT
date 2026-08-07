module neort_gc_cylindrical_orbit
    !! Fixed-invariant real-space return maps and temporal perturbation phases.
    !! No flux-chart coordinate, Boozer angle, or radial-domain sentinel is
    !! used by this module.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use fortnum_ode, only: ODE_EVENT_FALLING
    use fortnum_ode_vode, only: vode_init, vode_integrate_to, vode_state_t
    use fortnum_status, only: FORTNUM_OK, fortnum_status_t
    use neort_gc_cylindrical_dynamics, only: gc_cylindrical_rhs
    use neort_gc_cylindrical_model, only: GC_CYL_FIELD_ERROR, &
        GC_CYL_INTEGRATOR_ERROR, GC_CYL_INVALID_INPUT, &
        GC_CYL_NO_RETURN_MODEL => GC_CYL_NO_RETURN, &
        GC_CYL_EQUILIBRIUM_DOMAIN, GC_CYL_INVARIANT_ERROR, &
        GC_CYL_PERTURBATION_ERROR, &
        GC_CYL_POTENTIAL_ERROR, &
        GC_CYL_SINGULAR_BSTAR, &
        GC_CYL_SECTION_PHI, GC_CYL_START_ERROR, GC_CYL_STATE_ERROR, &
        GC_CYL_SUCCESS, GC_CYL_WALL_ERROR, &
        GC_CYL_WALL_LOSS_MODEL => GC_CYL_WALL_LOSS, &
        gc_cylindrical_field_sample_t, &
        gc_cylindrical_field_t, gc_cylindrical_invariants_t, &
        gc_cylindrical_potential_t, gc_cylindrical_section_t, &
        gc_cylindrical_state_t, gc_cylindrical_wall_t, &
        gc_cylindrical_section_rate, gc_cylindrical_section_value, &
        gc_cylindrical_invariant_residuals, &
        invariants_from_cylindrical_state, &
        state_from_cylindrical_invariants

    implicit none
    private

    real(dp), parameter :: PI = 3.1415926535897932384626433832795_dp

    ! Kept as local public aliases so callers can classify orbit results
    ! without importing the implementation/model layer separately.
    integer, parameter, public :: GC_CYL_WALL_LOSS = GC_CYL_WALL_LOSS_MODEL
    integer, parameter, public :: GC_CYL_NO_RETURN = GC_CYL_NO_RETURN_MODEL

    type, public :: gc_cylindrical_orbit_options_t
        real(dp) :: relative_tolerance = 3.0e-10_dp
        real(dp) :: absolute_tolerance = 1.0e-11_dp
        real(dp) :: invariant_relative_tolerance = 1.0e-8_dp
        real(dp) :: event_time_tolerance = 1.0e-10_dp
        real(dp) :: minimum_event_time = 1.0e-9_dp
        real(dp) :: maximum_time = 1.0e3_dp
    end type gc_cylindrical_orbit_options_t

    type, public :: gc_cylindrical_return_t
        integer :: status = GC_CYL_INTEGRATOR_ERROR
        real(dp) :: period = 0.0_dp
        real(dp) :: delta_phi = 0.0_dp
        real(dp) :: omega_return = 0.0_dp
        real(dp) :: omega_b = 0.0_dp
        real(dp) :: omega_phi = 0.0_dp
        real(dp) :: omega_prec = 0.0_dp
        real(dp) :: precession = 0.0_dp
        real(dp) :: topology_delta_phi = 0.0_dp
        real(dp) :: fieldline_q = 0.0_dp
        integer :: frequency_winding = 0
        real(dp) :: energy_error = 0.0_dp
        real(dp) :: magnetic_moment_error = 0.0_dp
        real(dp) :: canonical_momentum_error = 0.0_dp
        type(gc_cylindrical_state_t) :: state_at_return
        logical :: physical_loss = .false.
        logical :: numerical_failure = .true.
        logical :: invariant_rejected = .false.
    end type gc_cylindrical_return_t

    type, public :: gc_cylindrical_phase_result_t
        integer :: status = GC_CYL_INTEGRATOR_ERROR
        real(dp) :: period = 0.0_dp
        real(dp) :: delta_phi = 0.0_dp
        real(dp) :: omega_b = 0.0_dp
        real(dp) :: omega_phi = 0.0_dp
        real(dp) :: omega_prec = 0.0_dp
        complex(dp) :: perturbation_average = (0.0_dp, 0.0_dp)
        real(dp) :: inverse_b_average = 0.0_dp
        real(dp) :: b_average = 0.0_dp
        real(dp) :: energy_error = 0.0_dp
        real(dp) :: magnetic_moment_error = 0.0_dp
        real(dp) :: canonical_momentum_error = 0.0_dp
        type(gc_cylindrical_state_t) :: state_at_return
        logical :: physical_loss = .false.
        logical :: numerical_failure = .true.
        logical :: invariant_rejected = .false.
    end type gc_cylindrical_phase_result_t

    abstract interface
        subroutine gc_cylindrical_perturbation_i(position, field, amplitude, status)
            import :: dp, gc_cylindrical_field_sample_t
            real(dp), intent(in) :: position(3)
            type(gc_cylindrical_field_sample_t), intent(in) :: field
            complex(dp), intent(out) :: amplitude
            integer, intent(out) :: status
        end subroutine gc_cylindrical_perturbation_i
    end interface

    public :: compute_gc_cylindrical_return
    public :: compute_gc_cylindrical_phase_average
    public :: gc_cylindrical_perturbation_i
    public :: gc_cylindrical_status_is_physical_loss
    public :: gc_cylindrical_status_is_numerical

contains

    pure logical function gc_cylindrical_status_is_physical_loss(status)
        integer, intent(in) :: status

        gc_cylindrical_status_is_physical_loss = status == GC_CYL_WALL_LOSS
    end function gc_cylindrical_status_is_physical_loss

    pure logical function gc_cylindrical_status_is_numerical(status)
        integer, intent(in) :: status

        gc_cylindrical_status_is_numerical = status == GC_CYL_EQUILIBRIUM_DOMAIN
        if (.not. gc_cylindrical_status_is_numerical) then
            gc_cylindrical_status_is_numerical = status == GC_CYL_FIELD_ERROR
        end if
        if (.not. gc_cylindrical_status_is_numerical) then
            gc_cylindrical_status_is_numerical = status == GC_CYL_POTENTIAL_ERROR
        end if
        if (.not. gc_cylindrical_status_is_numerical) then
            gc_cylindrical_status_is_numerical = status == GC_CYL_STATE_ERROR
        end if
        if (.not. gc_cylindrical_status_is_numerical) then
            gc_cylindrical_status_is_numerical = status == GC_CYL_INTEGRATOR_ERROR
        end if
        if (.not. gc_cylindrical_status_is_numerical) then
            gc_cylindrical_status_is_numerical = status == GC_CYL_NO_RETURN_MODEL
        end if
        if (.not. gc_cylindrical_status_is_numerical) then
            gc_cylindrical_status_is_numerical = status == GC_CYL_WALL_ERROR
        end if
        if (.not. gc_cylindrical_status_is_numerical) then
            gc_cylindrical_status_is_numerical = status == GC_CYL_PERTURBATION_ERROR
        end if
        if (.not. gc_cylindrical_status_is_numerical) then
            gc_cylindrical_status_is_numerical = status == GC_CYL_START_ERROR
        end if
        if (.not. gc_cylindrical_status_is_numerical) then
            gc_cylindrical_status_is_numerical = status == GC_CYL_SINGULAR_BSTAR
        end if
        if (.not. gc_cylindrical_status_is_numerical) then
            gc_cylindrical_status_is_numerical = status == GC_CYL_INVARIANT_ERROR
        end if
    end function gc_cylindrical_status_is_numerical

    subroutine compute_gc_cylindrical_return(field_model, potential_model, &
            invariants, reference_position, parallel_sign, mass, charge, c_light, &
            section, options, result, wall_model)
        class(gc_cylindrical_field_t), intent(in) :: field_model
        class(gc_cylindrical_potential_t), intent(in) :: potential_model
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: reference_position(3), mass, charge, c_light
        integer, intent(in) :: parallel_sign
        type(gc_cylindrical_section_t), intent(in) :: section
        type(gc_cylindrical_orbit_options_t), intent(in) :: options
        type(gc_cylindrical_return_t), intent(out) :: result
        class(gc_cylindrical_wall_t), intent(in), optional :: wall_model

        type(vode_state_t) :: integrator
        type(fortnum_status_t) :: integration_status
        real(dp) :: initial_state(5)
        real(dp), allocatable :: final_state(:)
        real(dp) :: atol(5), target_time, root_time, event_tolerance
        integer :: status, rhs_status, event_direction, event_index
        logical :: found

        result = gc_cylindrical_return_t()
        call initialize_orbit_state(field_model, potential_model, invariants, &
            reference_position, parallel_sign, mass, charge, c_light, &
            initial_state, status)
        if (status /= GC_CYL_SUCCESS) then
            result%status = status
            call classify_result(result)
            return
        end if
        if (present(wall_model)) then
            call configure_event(section, initial_state, mass, charge, c_light, &
                field_model, potential_model, event_direction, rhs_status, &
                wall_model)
        else
            call configure_event(section, initial_state, mass, charge, c_light, &
                field_model, potential_model, event_direction, rhs_status)
        end if
        if (rhs_status /= GC_CYL_SUCCESS) then
            result%status = rhs_status
            call classify_result(result)
            return
        end if
        target_time = options%maximum_time
        if (target_time <= 0.0_dp) then
            result%status = GC_CYL_INVALID_INPUT
            call classify_result(result)
            return
        end if
        call set_component_tolerances(field_model, potential_model, initial_state, &
            invariants, mass, charge, options, atol, rhs_status)
        if (rhs_status /= GC_CYL_SUCCESS) then
            result%status = rhs_status
            call classify_result(result)
            return
        end if
        event_tolerance = max(options%event_time_tolerance, &
            100.0_dp*epsilon(target_time)*max(1.0_dp, target_time))
        call vode_init(integrator, 5, 0.0_dp, initial_state)
        root_time = 0.0_dp
        call vode_integrate_to(orbit_rhs, integrator, target_time, &
            options%relative_tolerance, atol, final_state, integration_status, &
            event=return_event, event_dir=event_direction, &
            event_tol=event_tolerance, t_root=root_time, root_found=found, &
            event_index=event_index, event2=wall_event, &
            event_dir2=ODE_EVENT_FALLING)
        if (rhs_status /= GC_CYL_SUCCESS) then
            result%status = rhs_status
            call classify_result(result)
            return
        end if
        if (integration_status%code /= FORTNUM_OK) then
            result%status = GC_CYL_INTEGRATOR_ERROR
            call classify_result(result)
            return
        end if
        if (event_index == 2 .or. rhs_status == GC_CYL_WALL_LOSS_MODEL) then
            result%status = GC_CYL_WALL_LOSS_MODEL
            call classify_result(result)
            return
        end if
        if (.not. found .or. root_time <= 0.0_dp) then
            result%status = GC_CYL_NO_RETURN
            call classify_result(result)
            return
        end if
        call complete_return(field_model, potential_model, invariants, &
            initial_state(3), final_state, root_time, mass, charge, c_light, &
            section, options, result)

    contains

        subroutine orbit_rhs(time, state, derivative, context)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            real(dp), intent(out) :: derivative(:)
            class(*), intent(in), optional :: context

            associate (unused_time => time, unused_context => context)
            end associate
            if (present(wall_model)) then
                call evaluate_rhs(field_model, potential_model, state, mass, &
                    charge, c_light, derivative, rhs_status, wall_model)
            else
                call evaluate_rhs(field_model, potential_model, state, mass, &
                    charge, c_light, derivative, rhs_status)
            end if
        end subroutine orbit_rhs

        function return_event(time, state, context) result(value)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            class(*), intent(in), optional :: context
            real(dp) :: value
            type(gc_cylindrical_state_t) :: local_state

            associate (unused_time => time, unused_context => context)
            end associate
            local_state = state_from_array(state)
            value = gc_cylindrical_section_value(section, local_state)
        end function return_event

        function wall_event(time, state, context) result(value)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            class(*), intent(in), optional :: context
            real(dp) :: value, distance
            integer :: wall_status

            associate (unused_time => time, unused_context => context)
            end associate
            value = huge(1.0_dp)
            if (.not. present(wall_model)) return
            call wall_model%evaluate(state(1:3), distance, wall_status)
            if (wall_status /= GC_CYL_SUCCESS) then
                rhs_status = GC_CYL_WALL_ERROR
                return
            end if
            value = distance
        end function wall_event

    end subroutine compute_gc_cylindrical_return

    subroutine compute_gc_cylindrical_phase_average(field_model, potential_model, &
            invariants, reference_position, parallel_sign, mass, charge, c_light, &
            omega_b, omega_phi, bounce_harmonic, toroidal_harmonic, perturbation, &
            options, section, result, wall_model)
        class(gc_cylindrical_field_t), intent(in) :: field_model
        class(gc_cylindrical_potential_t), intent(in) :: potential_model
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: reference_position(3), mass, charge, c_light
        real(dp), intent(in) :: omega_b, omega_phi
        integer, intent(in) :: parallel_sign, bounce_harmonic, toroidal_harmonic
        procedure(gc_cylindrical_perturbation_i) :: perturbation
        type(gc_cylindrical_orbit_options_t), intent(in) :: options
        type(gc_cylindrical_section_t), intent(in) :: section
        type(gc_cylindrical_phase_result_t), intent(out) :: result
        class(gc_cylindrical_wall_t), intent(in), optional :: wall_model

        type(vode_state_t) :: integrator
        type(fortnum_status_t) :: integration_status
        real(dp) :: initial_state(9)
        real(dp), allocatable :: final_state(:)
        real(dp) :: atol(9), target_time, root_time, event_tolerance
        integer :: status, rhs_status, event_direction, event_index
        logical :: found

        result = gc_cylindrical_phase_result_t()
        call initialize_orbit_state(field_model, potential_model, invariants, &
            reference_position, parallel_sign, mass, charge, c_light, &
            initial_state(1:5), status)
        if (status /= GC_CYL_SUCCESS) then
            result%status = status
            call classify_phase_result(result)
            return
        end if
        initial_state(6:9) = 0.0_dp
        if (present(wall_model)) then
            call configure_event(section, initial_state(1:5), mass, charge, c_light, &
                field_model, potential_model, event_direction, rhs_status, &
                wall_model)
        else
            call configure_event(section, initial_state(1:5), mass, charge, c_light, &
                field_model, potential_model, event_direction, rhs_status)
        end if
        if (rhs_status /= GC_CYL_SUCCESS) then
            result%status = rhs_status
            call classify_phase_result(result)
            return
        end if
        target_time = options%maximum_time
        if (target_time <= 0.0_dp) then
            result%status = GC_CYL_INVALID_INPUT
            call classify_phase_result(result)
            return
        end if
        call set_component_tolerances(field_model, potential_model, initial_state(1:5), &
            invariants, mass, charge, options, atol(1:5), rhs_status)
        if (rhs_status /= GC_CYL_SUCCESS) then
            result%status = rhs_status
            call classify_phase_result(result)
            return
        end if
        atol(6:7) = options%absolute_tolerance
        atol(8:9) = options%absolute_tolerance
        event_tolerance = max(options%event_time_tolerance, &
            100.0_dp*epsilon(target_time)*max(1.0_dp, target_time))
        call vode_init(integrator, 9, 0.0_dp, initial_state)
        root_time = 0.0_dp
        call vode_integrate_to(phase_rhs, integrator, target_time, &
            options%relative_tolerance, atol, final_state, integration_status, &
            event=phase_return_event, event_dir=event_direction, &
            event_tol=event_tolerance, t_root=root_time, root_found=found, &
            event_index=event_index, event2=wall_event, &
            event_dir2=ODE_EVENT_FALLING)
        if (rhs_status /= GC_CYL_SUCCESS) then
            result%status = rhs_status
            call classify_phase_result(result)
            return
        end if
        if (integration_status%code /= FORTNUM_OK) then
            result%status = GC_CYL_INTEGRATOR_ERROR
            call classify_phase_result(result)
            return
        end if
        if (event_index == 2 .or. rhs_status == GC_CYL_WALL_LOSS_MODEL) then
            result%status = GC_CYL_WALL_LOSS_MODEL
            call classify_phase_result(result)
            return
        end if
        if (.not. found .or. root_time <= 0.0_dp) then
            result%status = GC_CYL_NO_RETURN
            call classify_phase_result(result)
            return
        end if
        result%period = root_time
        result%delta_phi = final_state(3) - initial_state(3)
        result%omega_phi = result%delta_phi/result%period
        result%omega_b = omega_b
        result%omega_prec = result%omega_phi
        if (section%frequency_winding /= 0) then
            result%omega_prec = result%omega_phi &
                -2.0_dp*PI*real(section%frequency_winding, dp) &
                *section%fieldline_q/result%period
        end if
        result%perturbation_average = cmplx(final_state(6), final_state(7), dp) &
            /result%period
        result%inverse_b_average = final_state(8)/result%period
        result%b_average = final_state(9)/result%period
        result%state_at_return = state_from_array(final_state(1:5))
        call validate_invariant_return(field_model, potential_model, invariants, &
            final_state(1:5), mass, charge, c_light, options, &
            result%energy_error, result%magnetic_moment_error, &
            result%canonical_momentum_error, result%invariant_rejected, status)
        if (status /= GC_CYL_SUCCESS) then
            result%status = status
        else
            result%status = GC_CYL_SUCCESS
        end if
        call classify_phase_result(result)

    contains

        subroutine phase_rhs(time, state, derivative, context)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            real(dp), intent(out) :: derivative(:)
            class(*), intent(in), optional :: context

            type(gc_cylindrical_state_t) :: local_state
            type(gc_cylindrical_field_sample_t) :: field
            complex(dp) :: amplitude, phase
            real(dp) :: gradient(3), potential, phase_argument
            integer :: local_status, field_status, potential_status, perturbation_status

            associate (unused_context => context)
            end associate
            derivative = 0.0_dp
            if (present(wall_model)) then
                call evaluate_rhs(field_model, potential_model, state(1:5), mass, &
                    charge, c_light, derivative(1:5), rhs_status, wall_model)
            else
                call evaluate_rhs(field_model, potential_model, state(1:5), mass, &
                    charge, c_light, derivative(1:5), rhs_status)
            end if
            if (rhs_status /= GC_CYL_SUCCESS) return
            local_state = state_from_array(state(1:5))
            call field_model%evaluate(state(1:3), field, field_status)
            if (field_status /= GC_CYL_SUCCESS) then
                rhs_status = map_field_status(field_status)
                return
            end if
            call potential_model%evaluate(state(1:3), field, potential, gradient, &
                potential_status)
            if (potential_status /= GC_CYL_SUCCESS) then
                rhs_status = GC_CYL_POTENTIAL_ERROR
                return
            end if
            call perturbation(state(1:3), field, amplitude, perturbation_status)
            if (perturbation_status /= GC_CYL_SUCCESS) then
                rhs_status = GC_CYL_PERTURBATION_ERROR
                return
            end if
            phase_argument = real(toroidal_harmonic, dp) &
                *(local_state%phi - section%point(3)) &
                -(real(bounce_harmonic, dp)*omega_b &
                +real(toroidal_harmonic, dp)*omega_phi)*time
            phase = cmplx(cos(phase_argument), sin(phase_argument), dp)
            derivative(6) = real(amplitude*phase, dp)
            derivative(7) = aimag(amplitude*phase)
            derivative(8) = 1.0_dp/field%bmod
            derivative(9) = field%bmod
        end subroutine phase_rhs

        function phase_return_event(time, state, context) result(value)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            class(*), intent(in), optional :: context
            real(dp) :: value
            type(gc_cylindrical_state_t) :: local_state

            associate (unused_time => time, unused_context => context)
            end associate
            local_state = state_from_array(state(1:5))
            value = gc_cylindrical_section_value(section, local_state)
        end function phase_return_event

        function wall_event(time, state, context) result(value)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            class(*), intent(in), optional :: context
            real(dp) :: value, distance
            integer :: wall_status

            associate (unused_time => time, unused_context => context)
            end associate
            value = huge(1.0_dp)
            if (.not. present(wall_model)) return
            call wall_model%evaluate(state(1:3), distance, wall_status)
            if (wall_status /= GC_CYL_SUCCESS) then
                rhs_status = GC_CYL_WALL_ERROR
                return
            end if
            value = distance
        end function wall_event

    end subroutine compute_gc_cylindrical_phase_average

    subroutine initialize_orbit_state(field_model, potential_model, invariants, &
            reference_position, parallel_sign, mass, charge, c_light, state_array, &
            status)
        class(gc_cylindrical_field_t), intent(in) :: field_model
        class(gc_cylindrical_potential_t), intent(in) :: potential_model
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: reference_position(3), mass, charge, c_light
        integer, intent(in) :: parallel_sign
        real(dp), intent(out) :: state_array(5)
        integer, intent(out) :: status

        type(gc_cylindrical_field_sample_t) :: field
        type(gc_cylindrical_state_t) :: state
        real(dp) :: potential, gradient(3)
        integer :: field_status, potential_status

        status = GC_CYL_INVALID_INPUT
        state_array = 0.0_dp
        call field_model%evaluate(reference_position, field, field_status)
        if (field_status /= GC_CYL_SUCCESS) then
            status = map_field_status(field_status)
            return
        end if
        call potential_model%evaluate(reference_position, field, potential, gradient, &
            potential_status)
        if (potential_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_POTENTIAL_ERROR
            return
        end if
        call state_from_cylindrical_invariants(reference_position, field, potential, &
            invariants, parallel_sign, mass, charge, c_light, state, status)
        if (status /= GC_CYL_SUCCESS) return
        state_array = state_to_array(state)
    end subroutine initialize_orbit_state

    subroutine evaluate_rhs(field_model, potential_model, state_array, mass, charge, &
            c_light, derivative, status, wall_model)
        class(gc_cylindrical_field_t), intent(in) :: field_model
        class(gc_cylindrical_potential_t), intent(in) :: potential_model
        real(dp), intent(in) :: state_array(5), mass, charge, c_light
        real(dp), intent(out) :: derivative(5)
        integer, intent(out) :: status
        class(gc_cylindrical_wall_t), intent(in), optional :: wall_model

        type(gc_cylindrical_field_sample_t) :: field
        type(gc_cylindrical_state_t) :: state
        real(dp) :: position(3), potential, gradient(3), distance
        integer :: field_status, potential_status, wall_status, dynamics_status

        derivative = 0.0_dp
        state = state_from_array(state_array)
        position = state_array(1:3)
        if (present(wall_model)) then
            call wall_model%evaluate(position, distance, wall_status)
            if (wall_status /= GC_CYL_SUCCESS) then
                status = GC_CYL_WALL_ERROR
                return
            end if
            if (distance <= 0.0_dp) then
                status = GC_CYL_WALL_LOSS
                return
            end if
        end if
        call field_model%evaluate(position, field, field_status)
        if (field_status /= GC_CYL_SUCCESS) then
            status = map_field_status(field_status)
            return
        end if
        call potential_model%evaluate(position, field, potential, gradient, &
            potential_status)
        if (potential_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_POTENTIAL_ERROR
            return
        end if
        call gc_cylindrical_rhs(field, gradient, mass, charge, c_light, state, &
            derivative, dynamics_status)
        status = dynamics_status
    end subroutine evaluate_rhs

    subroutine set_component_tolerances(field_model, potential_model, state_array, &
            invariants, mass, charge, options, atol, status)
        class(gc_cylindrical_field_t), intent(in) :: field_model
        class(gc_cylindrical_potential_t), intent(in) :: potential_model
        real(dp), intent(in) :: state_array(5)
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: mass, charge
        type(gc_cylindrical_orbit_options_t), intent(in) :: options
        real(dp), intent(out) :: atol(5)
        integer, intent(out) :: status

        type(gc_cylindrical_field_sample_t) :: field
        type(gc_cylindrical_state_t) :: state
        real(dp) :: potential, gradient(3), coordinate_scale
        real(dp) :: momentum_scale, moment_scale, base
        integer :: field_status, potential_status

        atol = 0.0_dp
        status = GC_CYL_INVALID_INPUT
        if (options%relative_tolerance <= 0.0_dp .or. &
                options%absolute_tolerance <= 0.0_dp .or. &
                options%invariant_relative_tolerance <= 0.0_dp) return
        if (mass <= 0.0_dp .or. charge == 0.0_dp) return
        state = state_from_array(state_array)
        call field_model%evaluate(state_array(1:3), field, field_status)
        if (field_status /= GC_CYL_SUCCESS) then
            status = map_field_status(field_status)
            return
        end if
        call potential_model%evaluate(state_array(1:3), field, potential, gradient, &
            potential_status)
        if (potential_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_POTENTIAL_ERROR
            return
        end if

        base = options%absolute_tolerance
        ! VODE's absolute tolerances are physical component scales.  This is
        ! deliberately not a single CGS tolerance and does not compare p^2
        ! against an arbitrary max(1, ...) floor.
        coordinate_scale = max(abs(state%R), abs(state%Z), tiny(state%R))
        momentum_scale = sqrt(2.0_dp*mass*max(abs(invariants%energy), &
            tiny(invariants%energy)))
        momentum_scale = max(momentum_scale, abs(state%p_parallel), &
            tiny(state%p_parallel))
        moment_scale = max(abs(invariants%magnetic_moment), &
            abs(invariants%energy)/max(field%bmod, tiny(field%bmod)), &
            tiny(invariants%magnetic_moment))
        atol(1:2) = base*coordinate_scale
        atol(3) = base
        atol(4) = base*momentum_scale
        atol(5) = base*moment_scale
        if (.not. all(ieee_is_finite(atol))) then
            status = GC_CYL_INVALID_INPUT
            return
        end if
        status = GC_CYL_SUCCESS
    end subroutine set_component_tolerances

    subroutine validate_invariant_return(field_model, potential_model, invariants, &
            state_array, mass, charge, c_light, options, energy_error, &
            magnetic_moment_error, canonical_momentum_error, rejected, status)
        class(gc_cylindrical_field_t), intent(in) :: field_model
        class(gc_cylindrical_potential_t), intent(in) :: potential_model
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: state_array(5), mass, charge, c_light
        type(gc_cylindrical_orbit_options_t), intent(in) :: options
        real(dp), intent(out) :: energy_error, magnetic_moment_error
        real(dp), intent(out) :: canonical_momentum_error
        logical, intent(out) :: rejected
        integer, intent(out) :: status

        type(gc_cylindrical_field_sample_t) :: field
        type(gc_cylindrical_state_t) :: state
        real(dp) :: potential, gradient(3), tolerance
        integer :: field_status, potential_status, residual_status

        energy_error = 0.0_dp
        magnetic_moment_error = 0.0_dp
        canonical_momentum_error = 0.0_dp
        rejected = .false.
        status = GC_CYL_INVALID_INPUT
        if (options%invariant_relative_tolerance <= 0.0_dp) return
        state = state_from_array(state_array)
        call field_model%evaluate(state_array(1:3), field, field_status)
        if (field_status /= GC_CYL_SUCCESS) then
            status = map_field_status(field_status)
            return
        end if
        call potential_model%evaluate(state_array(1:3), field, potential, gradient, &
            potential_status)
        if (potential_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_POTENTIAL_ERROR
            return
        end if
        call gc_cylindrical_invariant_residuals(state, field, potential, mass, charge, &
            c_light, invariants, energy_error, magnetic_moment_error, &
            canonical_momentum_error, residual_status)
        if (residual_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_STATE_ERROR
            return
        end if
        tolerance = options%invariant_relative_tolerance
        rejected = abs(energy_error) > tolerance*max(abs(invariants%energy), &
            tiny(invariants%energy)) .or. &
            abs(magnetic_moment_error) > tolerance*max(&
            abs(invariants%magnetic_moment), tiny(invariants%magnetic_moment)) .or. &
            abs(canonical_momentum_error) > tolerance*max(&
            abs(invariants%canonical_toroidal_momentum), &
            tiny(invariants%canonical_toroidal_momentum))
        if (rejected) then
            status = GC_CYL_INVARIANT_ERROR
            return
        end if
        status = GC_CYL_SUCCESS
    end subroutine validate_invariant_return

    subroutine configure_event(section, initial_array, mass, charge, c_light, &
            field_model, potential_model, event_direction, status, wall_model)
        type(gc_cylindrical_section_t), intent(in) :: section
        real(dp), intent(in) :: initial_array(5), mass, charge, c_light
        class(gc_cylindrical_field_t), intent(in) :: field_model
        class(gc_cylindrical_potential_t), intent(in) :: potential_model
        integer, intent(out) :: event_direction, status
        class(gc_cylindrical_wall_t), intent(in), optional :: wall_model

        real(dp) :: derivative(5), event_rate
        if (present(wall_model)) then
            call evaluate_rhs(field_model, potential_model, initial_array, mass, &
                charge, c_light, derivative, status, wall_model)
        else
            call evaluate_rhs(field_model, potential_model, initial_array, mass, &
                charge, c_light, derivative, status)
        end if
        if (status /= GC_CYL_SUCCESS) return
        if (section%return_crossings < 1) then
            status = GC_CYL_INVALID_INPUT
            return
        end if
        event_rate = gc_cylindrical_section_rate(section, derivative)
        event_direction = section%direction
        if (event_direction == 0) then
            if (section%kind == GC_CYL_SECTION_PHI) then
                event_direction = sign_integer(section%winding)
                if (event_direction == 0) event_direction = sign_as_integer(event_rate)
            else
                ! The default section orientation follows the launch.  This
                ! requests the next same-oriented crossing of a full return,
                ! rather than the opposite crossing or a synthetic initial
                ! side that can manufacture a root at t=0.
                event_direction = sign_as_integer(event_rate)
            end if
        end if
        if (event_direction == 0) then
            status = GC_CYL_START_ERROR
            return
        end if
    end subroutine configure_event

    integer function map_field_status(field_status)
        integer, intent(in) :: field_status

        if (field_status == GC_CYL_EQUILIBRIUM_DOMAIN) then
            map_field_status = GC_CYL_EQUILIBRIUM_DOMAIN
        else
            map_field_status = GC_CYL_FIELD_ERROR
        end if
    end function map_field_status

    integer function sign_as_integer(value)
        real(dp), intent(in) :: value

        sign_as_integer = 0
        if (value > 0.0_dp) sign_as_integer = 1
        if (value < 0.0_dp) sign_as_integer = -1
    end function sign_as_integer

    integer function sign_integer(value)
        integer, intent(in) :: value

        sign_integer = 0
        if (value > 0) sign_integer = 1
        if (value < 0) sign_integer = -1
    end function sign_integer

    subroutine complete_return(field_model, potential_model, invariants, &
            initial_phi, final_array, period, mass, charge, c_light, section, &
            options, result)
        class(gc_cylindrical_field_t), intent(in) :: field_model
        class(gc_cylindrical_potential_t), intent(in) :: potential_model
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: initial_phi, final_array(:), period
        real(dp), intent(in) :: mass, charge, c_light
        type(gc_cylindrical_section_t), intent(in) :: section
        type(gc_cylindrical_orbit_options_t), intent(in) :: options
        type(gc_cylindrical_return_t), intent(inout) :: result

        type(gc_cylindrical_state_t) :: state
        integer :: validation_status

        state = state_from_array(final_array(1:5))
        call validate_invariant_return(field_model, potential_model, invariants, &
            final_array(1:5), mass, charge, c_light, options, &
            result%energy_error, result%magnetic_moment_error, &
            result%canonical_momentum_error, result%invariant_rejected, &
            validation_status)
        if (validation_status /= GC_CYL_SUCCESS) then
            result%status = validation_status
            call classify_result(result)
            return
        end if
        result%period = period
        result%delta_phi = final_array(3) - initial_phi
        result%omega_return = 2.0_dp*PI/period
        result%frequency_winding = section%frequency_winding
        result%fieldline_q = section%fieldline_q
        result%topology_delta_phi = 2.0_dp*PI*real(section%winding, dp) &
            *section%fieldline_q
        result%omega_b = 2.0_dp*PI*real(section%frequency_winding, dp)/period
        result%omega_phi = result%delta_phi/period
        result%omega_prec = result%omega_phi &
            -2.0_dp*PI*real(section%frequency_winding, dp) &
            *section%fieldline_q/period
        result%precession = result%omega_prec
        result%state_at_return = state
        result%status = GC_CYL_SUCCESS
        call classify_result(result)
    end subroutine complete_return

    subroutine classify_result(result)
        type(gc_cylindrical_return_t), intent(inout) :: result

        result%physical_loss = gc_cylindrical_status_is_physical_loss(result%status)
        result%numerical_failure = gc_cylindrical_status_is_numerical(result%status)
    end subroutine classify_result

    subroutine classify_phase_result(result)
        type(gc_cylindrical_phase_result_t), intent(inout) :: result

        result%physical_loss = gc_cylindrical_status_is_physical_loss(result%status)
        result%numerical_failure = gc_cylindrical_status_is_numerical(result%status)
    end subroutine classify_phase_result

    pure function state_to_array(state) result(array)
        type(gc_cylindrical_state_t), intent(in) :: state
        real(dp) :: array(5)

        array = [state%R, state%Z, state%phi, state%p_parallel, state%mu]
    end function state_to_array

    pure function state_from_array(array) result(state)
        real(dp), intent(in) :: array(5)
        type(gc_cylindrical_state_t) :: state

        state%R = array(1)
        state%Z = array(2)
        state%phi = array(3)
        state%p_parallel = array(4)
        state%mu = array(5)
    end function state_from_array

end module neort_gc_cylindrical_orbit
