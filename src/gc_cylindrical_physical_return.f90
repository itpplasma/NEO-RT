module neort_gc_cylindrical_physical_return
    !! Physical same-oriented return map in cylindrical real space.
    !!
    !! The caller supplies the scalar event that defines the physical
    !! Poincare cut.  The event is evaluated on the field sample at the
    !! current (R,Z,phi) state; no plane, Boozer angle, or flux coordinate is
    !! introduced here.  VODE supplies the RHS integration and dense
    !! Nordsieck/Illinois event interpolation.  Its internal step is bounded
    !! by maximum_step so an accepted step cannot silently contain the two
    !! opposite-oriented crossings of a periodic cut.
    !!
    !! A launch on C=0 is disarmed by integrating to minimum_return_time and
    !! restarting the same physical state with a fresh VODE history.  The
    !! restart is a numerical continuation only: it does not alter the RHS or
    !! the state.  The first subsequent crossing with the launch orientation
    !! is the returned complete cycle.  A wall and a radial-domain margin are
    !! combined as a signed positive-inside domain event because fortnum VODE
    !! exposes at most two simultaneous event functions.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode, only: ODE_EVENT_FALLING, ODE_EVENT_RISING
    use fortnum_ode_vode, only: vode_init, vode_integrate_to, vode_state_t
    use fortnum_status, only: FORTNUM_OK, fortnum_status_t
    use neort_gc_cylindrical_dynamics, only: gc_cylindrical_rhs
    use neort_gc_cylindrical_model, only: &
        GC_CYL_EQUILIBRIUM_DOMAIN, GC_CYL_FIELD_ERROR, &
        GC_CYL_INTEGRATOR_ERROR, GC_CYL_INVARIANT_ERROR, GC_CYL_INVALID_INPUT, &
        GC_CYL_NO_RETURN, GC_CYL_POTENTIAL_ERROR, GC_CYL_START_ERROR, &
        GC_CYL_STATE_ERROR, GC_CYL_SUCCESS, GC_CYL_WALL_ERROR, &
        GC_CYL_WALL_LOSS, gc_cylindrical_field_sample_t, &
        gc_cylindrical_field_t, gc_cylindrical_invariants_t, &
        gc_cylindrical_potential_t, gc_cylindrical_state_t, &
        gc_cylindrical_wall_t, gc_cylindrical_invariant_residuals

    implicit none
    private

    integer, parameter, public :: GC_CYL_PHYSICAL_EVENT_NONE = 0
    integer, parameter, public :: GC_CYL_PHYSICAL_EVENT_RETURN = 1
    integer, parameter, public :: GC_CYL_PHYSICAL_EVENT_WALL = 2
    integer, parameter, public :: GC_CYL_PHYSICAL_EVENT_RADIAL_DOMAIN = 3

    ! These are deliberately outside the model status range.  A caller can
    ! distinguish a callback contract failure from a field/RHS failure while
    ! still treating both as fail-closed numerical outcomes.
    integer, parameter, public :: GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR = 32
    integer, parameter, public :: GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE = 33

    type, public :: gc_cylindrical_physical_return_options_t
        real(dp) :: relative_tolerance = 3.0e-10_dp
        real(dp) :: absolute_tolerance(5) = 1.0e-11_dp
        real(dp) :: invariant_relative_tolerance = 1.0e-8_dp
        real(dp) :: invariant_absolute_tolerance(3) = 1.0e-11_dp
        real(dp) :: event_time_tolerance = 1.0e-10_dp
        real(dp) :: event_value_tolerance = 1.0e-8_dp
        real(dp) :: minimum_return_time = 1.0e-9_dp
        real(dp) :: maximum_time = 1.0e3_dp
        ! Zero selects maximum_time/1000.  Production callers should set this
        ! from the shortest physical cut period they intend to resolve.
        real(dp) :: maximum_step = 0.0_dp
        real(dp) :: wall_distance_scale = 1.0_dp
        real(dp) :: radial_distance_scale = 1.0_dp
        integer :: return_orientation = 0
        integer :: maximum_steps = 500000
    end type gc_cylindrical_physical_return_options_t

    type, public :: gc_cylindrical_physical_return_t
        integer :: status = GC_CYL_INTEGRATOR_ERROR
        integer :: event_kind = GC_CYL_PHYSICAL_EVENT_NONE
        integer :: return_orientation = 0
        real(dp) :: period = 0.0_dp
        real(dp) :: delta_phi = 0.0_dp
        real(dp) :: event_time = 0.0_dp
        real(dp) :: launch_event_value = 0.0_dp
        real(dp) :: return_event_value = 0.0_dp
        real(dp) :: energy_error = 0.0_dp
        real(dp) :: magnetic_moment_error = 0.0_dp
        real(dp) :: canonical_momentum_error = 0.0_dp
        real(dp) :: maximum_invariant_scaled_drift = 0.0_dp
        integer :: accepted_steps = 0
        integer :: rhs_evaluations = 0
        type(gc_cylindrical_state_t) :: state_at_event
        logical :: physical_return_found = .false.
        logical :: wall_hit = .false.
        logical :: radial_domain_hit = .false.
        logical :: invariant_rejected = .false.
        logical :: numerical_failure = .true.
    end type gc_cylindrical_physical_return_t

    abstract interface
        subroutine gc_cylindrical_physical_event_i(position, state, field, &
                user_data, value, status)
            import :: dp, gc_cylindrical_field_sample_t, &
                gc_cylindrical_state_t
            real(dp), intent(in) :: position(3)
            type(gc_cylindrical_state_t), intent(in) :: state
            type(gc_cylindrical_field_sample_t), intent(in) :: field
            class(*), pointer, intent(inout) :: user_data
            real(dp), intent(out) :: value
            integer, intent(out) :: status
        end subroutine gc_cylindrical_physical_event_i
    end interface

    abstract interface
        subroutine gc_cylindrical_radial_domain_i(position, state, field, &
                user_data, margin, status)
            import :: dp, gc_cylindrical_field_sample_t, &
                gc_cylindrical_state_t
            real(dp), intent(in) :: position(3)
            type(gc_cylindrical_state_t), intent(in) :: state
            type(gc_cylindrical_field_sample_t), intent(in) :: field
            class(*), pointer, intent(inout) :: user_data
            real(dp), intent(out) :: margin
            integer, intent(out) :: status
        end subroutine gc_cylindrical_radial_domain_i
    end interface

    public :: gc_cylindrical_physical_event_i
    public :: gc_cylindrical_radial_domain_i
    public :: compute_gc_cylindrical_physical_return

contains

    subroutine compute_gc_cylindrical_physical_return(field_model, &
            potential_model, initial_state, invariants, mass, charge, c_light, &
            return_event, options, result, wall_model, radial_domain, user_data)
        class(gc_cylindrical_field_t), intent(in) :: field_model
        class(gc_cylindrical_potential_t), intent(in) :: potential_model
        type(gc_cylindrical_state_t), intent(in) :: initial_state
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: mass, charge, c_light
        procedure(gc_cylindrical_physical_event_i) :: return_event
        type(gc_cylindrical_physical_return_options_t), intent(in) :: options
        type(gc_cylindrical_physical_return_t), intent(out) :: result
        class(gc_cylindrical_wall_t), intent(in), optional :: wall_model
        procedure(gc_cylindrical_radial_domain_i), optional :: radial_domain
        class(*), target, intent(inout), optional :: user_data

        type(vode_state_t) :: integrator
        type(fortnum_status_t) :: integration_status
        type(gc_cylindrical_state_t) :: state_at_start
        real(dp) :: y_initial(5), y_start(5)
        real(dp), allocatable :: y_final(:)
        real(dp) :: derivative(5), target_time, pre_time
        real(dp) :: root_time, event_tolerance, maximum_step
        real(dp) :: event_value, probe_value, probe_time, slope
        real(dp) :: domain_value
        integer :: callback_status, event_index, pre_steps, pre_nfev
        integer :: event_orientation, domain_kind
        logical :: found, have_domain, have_wall, have_radial, valid
        class(*), pointer :: callback_data
        procedure(gc_cylindrical_radial_domain_i), pointer :: radial_domain_proc

        result = gc_cylindrical_physical_return_t()
        callback_status = GC_CYL_SUCCESS
        pre_steps = 0
        pre_nfev = 0
        nullify(callback_data)
        nullify(radial_domain_proc)
        if (present(user_data)) callback_data => user_data
        have_wall = .false.
        if (present(wall_model)) have_wall = .true.
        have_radial = .false.
        if (present(radial_domain)) then
            have_radial = .true.
            radial_domain_proc => radial_domain
        end if

        y_initial = state_to_array(initial_state)
        if (.not. valid_inputs(y_initial, invariants, mass, charge, c_light, &
                options)) then
            result%status = GC_CYL_INVALID_INPUT
            return
        end if
        if (have_radial) then
            have_domain = .true.
        else
            have_domain = .false.
        end if
        if (have_wall) have_domain = .true.
        if (options%maximum_step > 0.0_dp) then
            maximum_step = options%maximum_step
        else
            maximum_step = options%maximum_time/1000.0_dp
        end if
        if (maximum_step <= 0.0_dp) then
            result%status = GC_CYL_INVALID_INPUT
            return
        end if
        if (options%minimum_return_time < 0.0_dp) then
            result%status = GC_CYL_INVALID_INPUT
            return
        end if
        if (options%minimum_return_time >= options%maximum_time) then
            result%status = GC_CYL_NO_RETURN
            result%numerical_failure = .true.
            return
        end if

        call evaluate_state_data(y_initial, state_at_start, derivative, valid)
        if (.not. valid) then
            result%status = callback_status
            if (result%status == GC_CYL_SUCCESS) then
                result%status = GC_CYL_FIELD_ERROR
            end if
            return
        end if
        call evaluate_invariant_drift(state_at_start, derivative, valid)
        if (.not. valid) then
            result%status = callback_status
            if (result%status == GC_CYL_SUCCESS) then
                result%status = GC_CYL_STATE_ERROR
            end if
            return
        end if

        call evaluate_physical_value(y_initial, event_value, valid)
        if (.not. valid) then
            result%status = callback_status
            if (result%status == GC_CYL_SUCCESS) then
                result%status = GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR
            end if
            return
        end if
        result%launch_event_value = event_value
        if (abs(event_value) > options%event_value_tolerance) then
            result%status = GC_CYL_START_ERROR
            return
        end if

        event_orientation = options%return_orientation
        if (event_orientation == 0) then
            probe_time = choose_probe_time(options, maximum_step)
            if (probe_time <= 0.0_dp) then
                result%status = GC_CYL_START_ERROR
                return
            end if
            call evaluate_physical_value(y_initial + probe_time*derivative, &
                probe_value, valid)
            if (.not. valid) then
                result%status = callback_status
                if (result%status == GC_CYL_SUCCESS) then
                    result%status = GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR
                end if
                return
            end if
            slope = probe_value - event_value
            if (abs(slope) <= max(options%event_value_tolerance, &
                    100.0_dp*epsilon(max(1.0_dp, abs(probe_value))))) then
                result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
                return
            end if
            event_orientation = sign_integer(slope)
        end if
        if (abs(event_orientation) /= 1) then
            result%status = GC_CYL_START_ERROR
            return
        end if
        result%return_orientation = event_orientation

        if (have_domain) then
            call evaluate_domain(y_initial, domain_value, domain_kind, valid)
            if (.not. valid) then
                result%status = callback_status
                if (result%status == GC_CYL_SUCCESS) then
                    result%status = GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR
                end if
                return
            end if
            if (domain_value <= 0.0_dp) then
                call finish_domain_event(y_initial, 0.0_dp, domain_kind)
                return
            end if
        end if

        target_time = options%maximum_time
        pre_time = min(options%minimum_return_time, target_time)
        y_start = y_initial

        ! Integrate through the launch-disarm interval.  Physical cut events
        ! are deliberately absent in this call, while domain loss remains
        ! terminal and is still refined by VODE.
        call vode_init(integrator, 5, 0.0_dp, y_initial)
        integrator%hmax = maximum_step
        integrator%max_steps = options%maximum_steps
        if (pre_time > 0.0_dp) then
            root_time = 0.0_dp
            found = .false.
            event_index = 0
            if (have_domain) then
                call vode_integrate_to(physical_rhs, integrator, pre_time, &
                    options%relative_tolerance, options%absolute_tolerance, &
                    y_final, integration_status, event=domain_event, &
                    event_dir=ODE_EVENT_FALLING, &
                    event_tol=make_event_tolerance(pre_time, options), &
                    t_root=root_time, root_found=found, &
                    event_index=event_index)
            else
                call vode_integrate_to(physical_rhs, integrator, pre_time, &
                    options%relative_tolerance, options%absolute_tolerance, &
                    y_final, integration_status)
            end if
            pre_steps = integrator%nsteps
            pre_nfev = integrator%nfev
            if (callback_status /= GC_CYL_SUCCESS) then
                result%status = callback_status
                return
            end if
            if (integration_status%code /= FORTNUM_OK) then
                result%status = GC_CYL_INTEGRATOR_ERROR
                return
            end if
            if (have_domain) then
                if (found) then
                    call finish_domain_event(y_final, root_time, domain_kind)
                    return
                end if
            end if
            y_start = y_final
        end if

        ! Restarting at the exact disarmed state avoids using the VODE event
        ! detector's t=0 value as a physical return.  The RHS itself is
        ! unchanged, and the period is measured from the original launch.
        call vode_init(integrator, 5, pre_time, y_start)
        integrator%hmax = maximum_step
        integrator%max_steps = options%maximum_steps
        root_time = pre_time
        found = .false.
        event_index = 0
        event_tolerance = make_event_tolerance(target_time, options)
        if (have_domain) then
            call vode_integrate_to(physical_rhs, integrator, target_time, &
                options%relative_tolerance, options%absolute_tolerance, y_final, &
                integration_status, event=physical_event, &
                event_dir=map_event_direction(event_orientation), &
                event_tol=event_tolerance, t_root=root_time, root_found=found, &
                event_index=event_index, event2=domain_event, &
                event_dir2=ODE_EVENT_FALLING)
        else
            call vode_integrate_to(physical_rhs, integrator, target_time, &
                options%relative_tolerance, options%absolute_tolerance, y_final, &
                integration_status, event=physical_event, &
                event_dir=map_event_direction(event_orientation), &
                event_tol=event_tolerance, t_root=root_time, root_found=found, &
                event_index=event_index)
        end if
        result%accepted_steps = pre_steps+integrator%nsteps
        result%rhs_evaluations = pre_nfev+integrator%nfev
        if (callback_status /= GC_CYL_SUCCESS) then
            result%status = callback_status
            return
        end if
        if (integration_status%code /= FORTNUM_OK) then
            result%status = GC_CYL_INTEGRATOR_ERROR
            return
        end if
        if (.not. found) then
            result%status = GC_CYL_NO_RETURN
            result%numerical_failure = .true.
            return
        end if
        if (event_index == 2) then
            call evaluate_domain(y_final, domain_value, domain_kind, valid)
            if (.not. valid) then
                result%status = callback_status
                if (result%status == GC_CYL_SUCCESS) then
                    result%status = GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR
                end if
                return
            end if
            call finish_domain_event(y_final, root_time, domain_kind)
            return
        end if
        if (event_index /= 1) then
            result%status = GC_CYL_INTEGRATOR_ERROR
            return
        end if
        call finish_physical_event(y_final, root_time, valid)
        if (.not. valid) return

    contains

        subroutine physical_rhs(time, state_array, derivative_array, context)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state_array(:)
            real(dp), intent(out) :: derivative_array(:)
            class(*), intent(in), optional :: context

            type(gc_cylindrical_state_t) :: local_state
            type(gc_cylindrical_field_sample_t) :: local_field
            real(dp) :: potential, gradient(3)
            integer :: field_status, potential_status, dynamics_status
            logical :: local_valid

            associate (unused_time => time, unused_context => context)
            end associate
            derivative_array = 0.0_dp
            if (size(state_array) /= 5) then
                call note_callback_failure(GC_CYL_STATE_ERROR)
                return
            end if
            local_state = state_from_array(state_array(1:5))
            call field_model%evaluate(state_array(1:3), local_field, field_status)
            if (field_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(map_field_status(field_status))
                return
            end if
            call potential_model%evaluate(state_array(1:3), local_field, potential, &
                gradient, potential_status)
            if (potential_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(GC_CYL_POTENTIAL_ERROR)
                return
            end if
            call gc_cylindrical_rhs(local_field, gradient, mass, charge, c_light, &
                local_state, derivative_array(1:5), dynamics_status)
            if (dynamics_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(dynamics_status)
                return
            end if
            call update_drift(local_state, local_field, potential, local_valid)
            if (.not. local_valid) return
        end subroutine physical_rhs

        function physical_event(time, state_array, context) result(value)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state_array(:)
            class(*), intent(in), optional :: context
            real(dp) :: value
            logical :: local_valid

            associate (unused_time => time, unused_context => context)
            end associate
            value = huge(1.0_dp)
            if (size(state_array) /= 5) then
                call note_callback_failure(GC_CYL_STATE_ERROR)
                return
            end if
            call evaluate_physical_value(state_array(1:5), value, local_valid)
            if (.not. local_valid) value = huge(1.0_dp)
        end function physical_event

        function domain_event(time, state_array, context) result(value)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state_array(:)
            class(*), intent(in), optional :: context
            real(dp) :: value
            logical :: local_valid

            associate (unused_time => time, unused_context => context)
            end associate
            value = huge(1.0_dp)
            if (size(state_array) /= 5) then
                call note_callback_failure(GC_CYL_STATE_ERROR)
                return
            end if
            call evaluate_domain(state_array(1:5), value, domain_kind, local_valid)
            if (.not. local_valid) value = huge(1.0_dp)
        end function domain_event

        subroutine evaluate_state_data(state_array, state, derivative_array, ok)
            real(dp), intent(in) :: state_array(5)
            type(gc_cylindrical_state_t), intent(out) :: state
            real(dp), intent(out) :: derivative_array(5)
            logical, intent(out) :: ok

            type(gc_cylindrical_field_sample_t) :: local_field
            real(dp) :: potential, gradient(3)
            integer :: field_status, potential_status, dynamics_status

            state = state_from_array(state_array)
            derivative_array = 0.0_dp
            ok = .false.
            call field_model%evaluate(state_array(1:3), local_field, field_status)
            if (field_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(map_field_status(field_status))
                return
            end if
            call potential_model%evaluate(state_array(1:3), local_field, potential, &
                gradient, potential_status)
            if (potential_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(GC_CYL_POTENTIAL_ERROR)
                return
            end if
            call gc_cylindrical_rhs(local_field, gradient, mass, charge, c_light, &
                state, derivative_array, dynamics_status)
            if (dynamics_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(dynamics_status)
                return
            end if
            ok = .true.
        end subroutine evaluate_state_data

        subroutine evaluate_physical_value(state_array, value, ok)
            real(dp), intent(in) :: state_array(5)
            real(dp), intent(out) :: value
            logical, intent(out) :: ok

            type(gc_cylindrical_field_sample_t) :: local_field
            type(gc_cylindrical_state_t) :: local_state
            integer :: field_status, event_status

            value = huge(1.0_dp)
            ok = .false.
            local_state = state_from_array(state_array)
            call field_model%evaluate(state_array(1:3), local_field, field_status)
            if (field_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(map_field_status(field_status))
                return
            end if
            call return_event(state_array(1:3), local_state, local_field, &
                callback_data, value, event_status)
            if (event_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR)
                return
            end if
            if (.not. ieee_is_finite(value)) then
                call note_callback_failure(GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR)
                return
            end if
            ok = .true.
        end subroutine evaluate_physical_value

        subroutine evaluate_domain(state_array, value, kind, ok)
            real(dp), intent(in) :: state_array(5)
            real(dp), intent(out) :: value
            integer, intent(out) :: kind
            logical, intent(out) :: ok

            type(gc_cylindrical_field_sample_t) :: local_field
            type(gc_cylindrical_state_t) :: local_state
            real(dp) :: wall_value, radial_value, wall_scaled, radial_scaled
            integer :: wall_status, radial_status, field_status
            logical :: need_field

            value = huge(1.0_dp)
            kind = GC_CYL_PHYSICAL_EVENT_NONE
            ok = .false.
            wall_value = huge(1.0_dp)
            radial_value = huge(1.0_dp)
            need_field = .false.
            if (have_radial) need_field = .true.
            local_state = state_from_array(state_array)
            if (have_wall) then
                call wall_model%evaluate(state_array(1:3), wall_value, wall_status)
                if (wall_status /= GC_CYL_SUCCESS) then
                    call note_callback_failure(GC_CYL_WALL_ERROR)
                    return
                end if
            end if
            if (need_field) then
                call field_model%evaluate(state_array(1:3), local_field, field_status)
                if (field_status /= GC_CYL_SUCCESS) then
                    call note_callback_failure(map_field_status(field_status))
                    return
                end if
                call radial_domain_proc(state_array(1:3), local_state, local_field, &
                    callback_data, radial_value, radial_status)
                if (radial_status /= GC_CYL_SUCCESS) then
                    call note_callback_failure(GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR)
                    return
                end if
            end if
            if (.not. ieee_is_finite(wall_value)) then
                call note_callback_failure(GC_CYL_WALL_ERROR)
                return
            end if
            if (.not. ieee_is_finite(radial_value)) then
                call note_callback_failure(GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR)
                return
            end if
            wall_scaled = huge(1.0_dp)
            radial_scaled = huge(1.0_dp)
            if (have_wall) wall_scaled = &
                wall_value/options%wall_distance_scale
            if (have_radial) radial_scaled = &
                radial_value/options%radial_distance_scale
            if (wall_scaled <= radial_scaled) then
                value = wall_scaled
                kind = GC_CYL_PHYSICAL_EVENT_WALL
            else
                value = radial_scaled
                kind = GC_CYL_PHYSICAL_EVENT_RADIAL_DOMAIN
            end if
            ok = .true.
        end subroutine evaluate_domain

        subroutine update_drift(state, field, potential, ok)
            type(gc_cylindrical_state_t), intent(in) :: state
            type(gc_cylindrical_field_sample_t), intent(in) :: field
            real(dp), intent(in) :: potential
            logical, intent(out) :: ok

            real(dp) :: energy_error, moment_error, canonical_error
            real(dp) :: scaled_energy, scaled_moment, scaled_canonical
            integer :: residual_status

            ok = .false.
            call gc_cylindrical_invariant_residuals(state, field, potential, mass, &
                charge, c_light, invariants, energy_error, moment_error, &
                canonical_error, residual_status)
            if (residual_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(GC_CYL_STATE_ERROR)
                return
            end if
            scaled_energy = scaled_residual(energy_error, invariants%energy, 1, &
                options)
            scaled_moment = scaled_residual(moment_error, &
                invariants%magnetic_moment, 2, options)
            scaled_canonical = scaled_residual(canonical_error, &
                invariants%canonical_toroidal_momentum, 3, options)
            result%maximum_invariant_scaled_drift = max(&
                result%maximum_invariant_scaled_drift, scaled_energy, &
                scaled_moment, scaled_canonical)
            ok = .true.
        end subroutine update_drift

        subroutine finish_physical_event(state_array, event_time_value, ok)
            real(dp), intent(in) :: state_array(:), event_time_value
            logical, intent(out) :: ok

            type(gc_cylindrical_field_sample_t) :: local_field
            type(gc_cylindrical_state_t) :: local_state
            real(dp) :: potential, gradient(3), event_value_final
            real(dp) :: energy_error, moment_error, canonical_error
            integer :: field_status, potential_status, residual_status

            ok = .false.
            local_state = state_from_array(state_array(1:5))
            call field_model%evaluate(state_array(1:3), local_field, field_status)
            if (field_status /= GC_CYL_SUCCESS) then
                result%status = map_field_status(field_status)
                return
            end if
            call potential_model%evaluate(state_array(1:3), local_field, potential, &
                gradient, potential_status)
            if (potential_status /= GC_CYL_SUCCESS) then
                result%status = GC_CYL_POTENTIAL_ERROR
                return
            end if
            call gc_cylindrical_invariant_residuals(local_state, local_field, &
                potential, mass, charge, c_light, invariants, energy_error, &
                moment_error, canonical_error, residual_status)
            if (residual_status /= GC_CYL_SUCCESS) then
                result%status = GC_CYL_STATE_ERROR
                return
            end if
            result%energy_error = energy_error
            result%magnetic_moment_error = moment_error
            result%canonical_momentum_error = canonical_error
            if (scaled_residual(energy_error, invariants%energy, 1, options) > &
                    1.0_dp) then
                result%invariant_rejected = .true.
            end if
            if (scaled_residual(moment_error, invariants%magnetic_moment, 2, &
                    options) > 1.0_dp) then
                result%invariant_rejected = .true.
            end if
            if (scaled_residual(canonical_error, &
                    invariants%canonical_toroidal_momentum, 3, options) > &
                    1.0_dp) then
                result%invariant_rejected = .true.
            end if
            if (result%maximum_invariant_scaled_drift > 1.0_dp) then
                result%invariant_rejected = .true.
            end if
            if (result%invariant_rejected) then
                result%status = GC_CYL_INVARIANT_ERROR
                result%event_kind = GC_CYL_PHYSICAL_EVENT_RETURN
                result%event_time = event_time_value
                result%state_at_event = local_state
                return
            end if
            call evaluate_physical_value(state_array(1:5), event_value_final, ok)
            if (.not. ok) then
                result%status = callback_status
                if (result%status == GC_CYL_SUCCESS) then
                    result%status = GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR
                end if
                return
            end if
            result%status = GC_CYL_SUCCESS
            result%event_kind = GC_CYL_PHYSICAL_EVENT_RETURN
            result%period = event_time_value
            result%delta_phi = state_array(3) - y_initial(3)
            result%event_time = event_time_value
            result%return_event_value = event_value_final
            result%state_at_event = local_state
            result%physical_return_found = .true.
            result%numerical_failure = .false.
            ok = .true.
        end subroutine finish_physical_event

        subroutine finish_domain_event(state_array, event_time_value, kind)
            real(dp), intent(in) :: state_array(:), event_time_value
            integer, intent(in) :: kind

            result%event_kind = kind
            result%event_time = event_time_value
            result%state_at_event = state_from_array(state_array(1:5))
            if (kind == GC_CYL_PHYSICAL_EVENT_WALL) then
                result%status = GC_CYL_WALL_LOSS
                result%wall_hit = .true.
                result%numerical_failure = .false.
            else
                result%status = GC_CYL_EQUILIBRIUM_DOMAIN
                result%radial_domain_hit = .true.
                result%numerical_failure = .true.
            end if
        end subroutine finish_domain_event

        subroutine evaluate_invariant_drift(state, derivative_array, ok)
            type(gc_cylindrical_state_t), intent(in) :: state
            real(dp), intent(in) :: derivative_array(5)
            logical, intent(out) :: ok

            type(gc_cylindrical_field_sample_t) :: local_field
            real(dp) :: state_array(5)
            real(dp) :: potential, gradient(3)
            integer :: field_status, potential_status

            associate (unused_derivative => derivative_array)
            end associate
            ok = .false.
            state_array = state_to_array(state)
            call field_model%evaluate(state_array(1:3), local_field, field_status)
            if (field_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(map_field_status(field_status))
                return
            end if
            call potential_model%evaluate(state_array(1:3), local_field, &
                potential, gradient, potential_status)
            if (potential_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(GC_CYL_POTENTIAL_ERROR)
                return
            end if
            call update_drift(state, local_field, potential, ok)
        end subroutine evaluate_invariant_drift

        subroutine note_callback_failure(new_status)
            integer, intent(in) :: new_status

            if (callback_status == GC_CYL_SUCCESS) callback_status = new_status
        end subroutine note_callback_failure

    end subroutine compute_gc_cylindrical_physical_return

    logical function valid_inputs(state_array, invariants, mass, charge, c_light, &
            options) result(valid)
        real(dp), intent(in) :: state_array(5)
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: mass, charge, c_light
        type(gc_cylindrical_physical_return_options_t), intent(in) :: options

        valid = .false.
        if (.not. all(ieee_is_finite(state_array))) return
        if (.not. all(ieee_is_finite([invariants%energy, &
                invariants%magnetic_moment, &
                invariants%canonical_toroidal_momentum]))) return
        if (.not. all(ieee_is_finite(options%absolute_tolerance))) return
        if (.not. all(ieee_is_finite(options%invariant_absolute_tolerance))) return
        if (.not. all(ieee_is_finite([options%relative_tolerance, &
                options%invariant_relative_tolerance, &
                options%event_time_tolerance, options%event_value_tolerance, &
                options%minimum_return_time, options%maximum_time, &
                options%maximum_step, options%wall_distance_scale, &
                options%radial_distance_scale]))) return
        if (state_array(1) <= 0.0_dp) return
        if (state_array(5) < 0.0_dp) return
        if (mass <= 0.0_dp) return
        if (charge == 0.0_dp) return
        if (c_light <= 0.0_dp) return
        if (options%relative_tolerance <= 0.0_dp) return
        if (any(options%absolute_tolerance <= 0.0_dp)) return
        if (options%invariant_relative_tolerance <= 0.0_dp) return
        if (any(options%invariant_absolute_tolerance <= 0.0_dp)) return
        if (options%event_time_tolerance <= 0.0_dp) return
        if (options%event_value_tolerance < 0.0_dp) return
        if (options%maximum_time <= 0.0_dp) return
        if (options%maximum_step < 0.0_dp) return
        if (options%wall_distance_scale <= 0.0_dp) return
        if (options%radial_distance_scale <= 0.0_dp) return
        if (options%maximum_steps < 1) return
        if (abs(options%return_orientation) > 1) return
        valid = .true.
    end function valid_inputs

    real(dp) function choose_probe_time(options, maximum_step) result(value)
        type(gc_cylindrical_physical_return_options_t), intent(in) :: options
        real(dp), intent(in) :: maximum_step

        value = max(100.0_dp*options%event_time_tolerance, &
            1.0e-6_dp*options%maximum_time)
        value = min(value, 0.1_dp*options%maximum_time)
        value = min(value, 0.1_dp*maximum_step)
    end function choose_probe_time

    real(dp) function make_event_tolerance(time_scale, options) result(value)
        real(dp), intent(in) :: time_scale
        type(gc_cylindrical_physical_return_options_t), intent(in) :: options

        value = max(options%event_time_tolerance, &
            100.0_dp*epsilon(max(1.0_dp, time_scale))*max(1.0_dp, time_scale))
    end function make_event_tolerance

    integer function map_event_direction(orientation) result(value)
        integer, intent(in) :: orientation

        if (orientation > 0) then
            value = ODE_EVENT_RISING
        else
            value = ODE_EVENT_FALLING
        end if
    end function map_event_direction

    integer function sign_integer(value) result(sign_value)
        real(dp), intent(in) :: value

        sign_value = 0
        if (value > 0.0_dp) sign_value = 1
        if (value < 0.0_dp) sign_value = -1
    end function sign_integer

    integer function map_field_status(field_status) result(status)
        integer, intent(in) :: field_status

        if (field_status == GC_CYL_EQUILIBRIUM_DOMAIN) then
            status = GC_CYL_EQUILIBRIUM_DOMAIN
        else
            status = GC_CYL_FIELD_ERROR
        end if
    end function map_field_status

    real(dp) function scaled_residual(error_value, invariant_value, index, &
            options) &
            result(value)
        real(dp), intent(in) :: error_value, invariant_value
        integer, intent(in) :: index
        type(gc_cylindrical_physical_return_options_t), intent(in) :: options
        real(dp) :: scale

        ! This function is only called with indices 1..3 by this module.
        scale = options%invariant_absolute_tolerance(index) &
            +options%invariant_relative_tolerance*abs(invariant_value)
        value = abs(error_value)/max(scale, tiny(scale))
    end function scaled_residual

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

end module neort_gc_cylindrical_physical_return
