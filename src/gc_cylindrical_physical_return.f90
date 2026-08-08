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
    use neort_gc_callback_context, only: gc_callback_context_t
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

    ! These statuses describe the independent multiplicity-certification
    ! layer, not the numerical ODE result%status.
    integer, parameter, public :: &
        GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_UNKNOWN = 0
    integer, parameter, public :: &
        GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_CERTIFIED = 1
    integer, parameter, public :: &
        GC_CYL_PHYSICAL_RETURN_CERTIFICATE_UNAVAILABLE = 2
    integer, parameter, public :: &
        GC_CYL_PHYSICAL_RETURN_CERTIFICATE_INVALID = 3

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
        logical :: require_opposite_intersection = .false.
        logical :: require_transverse_intersection = .false.
        ! Retained for source compatibility with callers that have not yet
        ! been migrated to a theorem provider.  It is intentionally ignored:
        ! a caller flag cannot certify global cut multiplicity.
        logical :: complete_atlas_multiplicity_certified = .false.
        real(dp) :: cut_rate_tolerance = 1.0e-12_dp
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
        integer :: intersection_count = 0
        integer :: intersection_orientations(2) = 0
        real(dp) :: intersection_times(2) = 0.0_dp
        real(dp) :: intersection_rates(2) = 0.0_dp
        logical :: intersection_multiplicity_certified = .false.
        integer :: multiplicity_status = &
            GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_UNKNOWN
        type(gc_cylindrical_state_t) :: state_at_event
        logical :: physical_return_found = .false.
        logical :: wall_hit = .false.
        logical :: radial_domain_hit = .false.
        logical :: invariant_rejected = .false.
        logical :: numerical_failure = .true.
    end type gc_cylindrical_physical_return_t

    ! The numerical crossing list is evidence only.  This typed value is
    ! accepted as multiplicity proof only when an independent provider fills
    ! it with a nonzero identity and an exactly-two theorem result.
    type, public :: gc_cylindrical_physical_return_certificate_t
        integer :: certificate_id = 0
        integer :: crossing_count = 0
        logical :: exactly_two_proved = .false.
    end type gc_cylindrical_physical_return_certificate_t

    abstract interface
        subroutine gc_cylindrical_physical_event_i(position, state, field, &
                user_data, value, status)
            import :: dp, gc_cylindrical_field_sample_t, &
                gc_cylindrical_state_t, gc_callback_context_t
            real(dp), intent(in) :: position(3)
            type(gc_cylindrical_state_t), intent(in) :: state
            type(gc_cylindrical_field_sample_t), intent(in) :: field
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            real(dp), intent(out) :: value
            integer, intent(out) :: status
        end subroutine gc_cylindrical_physical_event_i
    end interface

    abstract interface
        subroutine gc_cylindrical_physical_return_multiplicity_provider( &
                evidence, certificate, status)
            import :: gc_cylindrical_physical_return_t
            import :: gc_cylindrical_physical_return_certificate_t
            type(gc_cylindrical_physical_return_t), intent(in) :: evidence
            type(gc_cylindrical_physical_return_certificate_t), intent(out) :: &
                certificate
            integer, intent(out) :: status
        end subroutine gc_cylindrical_physical_return_multiplicity_provider
    end interface

    abstract interface
        subroutine gc_cylindrical_physical_event_rate_i(position, state, field, &
                user_data, rate, status)
            import :: dp, gc_cylindrical_field_sample_t, &
                gc_cylindrical_state_t, gc_callback_context_t
            real(dp), intent(in) :: position(3)
            type(gc_cylindrical_state_t), intent(in) :: state
            type(gc_cylindrical_field_sample_t), intent(in) :: field
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            real(dp), intent(out) :: rate
            integer, intent(out) :: status
        end subroutine gc_cylindrical_physical_event_rate_i
    end interface

    abstract interface
        subroutine gc_cylindrical_radial_domain_i(position, state, field, &
                user_data, margin, status)
            import :: dp, gc_cylindrical_field_sample_t, &
                gc_cylindrical_state_t, gc_callback_context_t
            real(dp), intent(in) :: position(3)
            type(gc_cylindrical_state_t), intent(in) :: state
            type(gc_cylindrical_field_sample_t), intent(in) :: field
            class(gc_callback_context_t), pointer, intent(inout) :: user_data
            real(dp), intent(out) :: margin
            integer, intent(out) :: status
        end subroutine gc_cylindrical_radial_domain_i
    end interface

    public :: gc_cylindrical_physical_event_i
    public :: gc_cylindrical_physical_event_rate_i
    public :: gc_cylindrical_radial_domain_i
    public :: gc_cylindrical_physical_return_multiplicity_provider
    public :: compute_gc_cylindrical_physical_return
    public :: certify_gc_cylindrical_physical_return
    public :: attach_gc_cylindrical_physical_return_certificate

contains

    subroutine compute_gc_cylindrical_physical_return(field_model, &
            potential_model, initial_state, invariants, mass, charge, c_light, &
            return_event, options, result, wall_model, radial_domain, user_data, &
            return_event_rate)
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
        ! Borrowed callback context: the caller-owned target must remain alive
        ! for this call and every callback invocation it performs.
        class(gc_callback_context_t), target, intent(inout), optional :: user_data
        procedure(gc_cylindrical_physical_event_rate_i), optional :: &
            return_event_rate

        type(vode_state_t) :: integrator
        type(fortnum_status_t) :: integration_status
        type(gc_cylindrical_state_t) :: state_at_start
        real(dp) :: y_initial(5), y_start(5)
        real(dp), allocatable :: y_final(:)
        real(dp), allocatable :: y_stage_start(:)
        real(dp) :: derivative(5), target_time, pre_time
        real(dp) :: root_time, maximum_step, cut_rate
        real(dp) :: opposite_root_time, disarm_time, disarmed_time
        real(dp) :: event_value, probe_value, probe_time, slope
        real(dp) :: disarmed_event_value
        real(dp) :: domain_value
        integer :: callback_status, event_index, pre_steps, pre_nfev
        integer :: stage_steps, stage_nfev
        integer :: disarm_steps, disarm_nfev, disarm_event_index
        integer :: event_orientation, domain_kind, local_status
        logical :: found, have_domain, have_wall, have_radial, valid
        logical :: disarm_found, disarmed_event_valid
        class(gc_callback_context_t), pointer :: callback_data
        procedure(gc_cylindrical_physical_event_i), pointer :: return_event_proc
        procedure(gc_cylindrical_radial_domain_i), pointer :: radial_domain_proc
        procedure(gc_cylindrical_physical_event_rate_i), pointer :: &
            return_event_rate_proc

        result = gc_cylindrical_physical_return_t()
        callback_status = GC_CYL_SUCCESS
        pre_steps = 0
        pre_nfev = 0
        nullify(callback_data)
        nullify(return_event_proc)
        nullify(radial_domain_proc)
        nullify(return_event_rate_proc)
        if (present(user_data)) callback_data => user_data
        return_event_proc => return_event
        if (present(return_event_rate)) return_event_rate_proc => return_event
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
        ! Always disarm the launch root before the first event restart.  The
        ! configured minimum return time is retained, but a zero value cannot
        ! disable this root-separation contract.
        pre_time = min(max(options%minimum_return_time, &
            choose_probe_time(options, maximum_step)), target_time)
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
            if (pre_time >= target_time) then
                result%status = GC_CYL_NO_RETURN
                result%numerical_failure = .true.
                return
            end if
        end if

        ! A launch at C=0 must be outside the event root before the first
        ! directed search.  This is the launch counterpart of the explicit
        ! post-opposite disarm below; it prevents VODE from returning t=0 as
        ! either directed crossing.
        call evaluate_physical_value(y_start, event_value, valid)
        if (.not. valid) then
            result%status = GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR
            result%numerical_failure = .true.
            return
        end if
        if (abs(event_value) <= options%event_value_tolerance) then
            result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
            result%numerical_failure = .true.
            return
        end if
        if (sign_integer(event_value) /= event_orientation) then
            result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
            result%numerical_failure = .true.
            return
        end if

        ! A Buchholz full-bounce return has two physical cut crossings: the
        ! opposite-oriented crossing first, followed by the same-oriented
        ! crossing that closes the period.  A single same-oriented event is
        ! not sufficient evidence and must not be repaired by a scalar factor.
        if (options%require_opposite_intersection) then
            call integrate_to_cut(y_start, pre_time, target_time, &
                -event_orientation, y_final, root_time, event_index, found, &
                stage_steps, stage_nfev, local_status)
            if (local_status /= GC_CYL_SUCCESS) then
                result%status = local_status
                result%numerical_failure = .true.
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
                result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
                return
            end if
            call evaluate_physical_rate(y_final, cut_rate, valid)
            if (.not. valid .or. abs(cut_rate) <= options%cut_rate_tolerance) then
                result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
                return
            end if
            pre_steps = pre_steps+stage_steps
            pre_nfev = pre_nfev+stage_nfev
            result%intersection_count = 1
            result%intersection_orientations(1) = -event_orientation
            result%intersection_times(1) = root_time
            result%intersection_rates(1) = cut_rate
            opposite_root_time = root_time

            ! The next event search must not restart on the exact previous
            ! root.  Advance the unchanged physical state with cut detection
            ! disabled, retaining wall/domain events, then restart from the
            ! disarmed state.  The saved opposite_root_time remains the true
            ! event time in the evidence.
            disarm_time = min(choose_probe_time(options, maximum_step), &
                0.25_dp*(target_time-opposite_root_time))
            if (disarm_time <= options%event_time_tolerance) then
                result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
                return
            end if
            y_stage_start = y_final
            call advance_without_cut(y_stage_start, opposite_root_time, &
                disarm_time, y_final, disarmed_time, disarm_found, &
                disarm_event_index, disarm_steps, disarm_nfev, local_status)
            if (local_status /= GC_CYL_SUCCESS) then
                result%status = local_status
                result%numerical_failure = .true.
                return
            end if
            pre_steps = pre_steps+disarm_steps
            pre_nfev = pre_nfev+disarm_nfev
            if (disarm_found) then
                call evaluate_domain(y_final, domain_value, domain_kind, valid)
                if (.not. valid) then
                    result%status = callback_status
                    if (result%status == GC_CYL_SUCCESS) then
                        result%status = GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR
                    end if
                    return
                end if
                call finish_domain_event(y_final, disarmed_time, domain_kind)
                return
            end if

            ! The disarm interval is part of the event contract, not merely a
            ! delay.  Verify that the unchanged state has left the root on the
            ! side implied by the first crossing before enabling the cut event
            ! again.  This makes a t=start rediscovery an explicit failure.
            call evaluate_physical_value(y_final, disarmed_event_value, &
                disarmed_event_valid)
            if (.not. disarmed_event_valid) then
                result%status = GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR
                result%numerical_failure = .true.
                return
            end if
            if (abs(disarmed_event_value) <= options%event_value_tolerance) then
                result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
                result%numerical_failure = .true.
                return
            end if
            if (sign_integer(disarmed_event_value) /= -event_orientation) then
                result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
                result%numerical_failure = .true.
                return
            end if

            y_stage_start = y_final
            call integrate_to_cut(y_stage_start, disarmed_time, target_time, &
                event_orientation, y_final, root_time, event_index, found, &
                stage_steps, stage_nfev, local_status)
            if (local_status /= GC_CYL_SUCCESS) then
                result%status = local_status
                result%numerical_failure = .true.
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
                result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
                return
            end if
            call evaluate_physical_rate(y_final, cut_rate, valid)
            if (.not. valid .or. abs(cut_rate) <= options%cut_rate_tolerance) then
                result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
                return
            end if
            pre_steps = pre_steps+stage_steps
            pre_nfev = pre_nfev+stage_nfev
            result%intersection_count = 2
            result%intersection_orientations(2) = event_orientation
            result%intersection_times(2) = root_time
            result%intersection_rates(2) = cut_rate
            ! Two events are numerical evidence only.  Global exactly-two
            ! multiplicity belongs to a separate theorem/provider and is
            ! attached by certify_gc_cylindrical_physical_return.
            result%intersection_multiplicity_certified = .false.
            result%multiplicity_status = &
                GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_UNKNOWN
        else
            call integrate_to_cut(y_start, pre_time, target_time, &
                event_orientation, y_final, root_time, event_index, found, &
                stage_steps, stage_nfev, local_status)
            if (local_status /= GC_CYL_SUCCESS) then
                result%status = local_status
                result%numerical_failure = .true.
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
                result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
                return
            end if
            if (options%require_transverse_intersection) then
                call evaluate_physical_rate(y_final, cut_rate, valid)
                if (.not. valid .or. &
                        abs(cut_rate) <= options%cut_rate_tolerance) then
                    result%status = GC_CYL_PHYSICAL_EVENT_NONTRANSVERSE
                    return
                end if
            end if
            pre_steps = pre_steps+stage_steps
            pre_nfev = pre_nfev+stage_nfev
        end if
        result%accepted_steps = pre_steps
        result%rhs_evaluations = pre_nfev
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
            call gc_cylindrical_rhs(local_field, potential, gradient, mass, &
                charge, c_light, local_state, derivative_array(1:5), &
                dynamics_status)
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
            call gc_cylindrical_rhs(local_field, potential, gradient, mass, &
                charge, c_light, state, derivative_array, dynamics_status)
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
            call return_event_proc(state_array(1:3), local_state, local_field, &
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

        subroutine evaluate_physical_rate(state_array, rate, ok)
            real(dp), intent(in) :: state_array(5)
            real(dp), intent(out) :: rate
            logical, intent(out) :: ok

            type(gc_cylindrical_state_t) :: local_state
            type(gc_cylindrical_field_sample_t) :: local_field
            integer :: field_status, rate_status

            rate = 0.0_dp
            ok = .false.
            if (.not. associated(return_event_rate_proc)) then
                if (options%require_transverse_intersection) then
                    call note_callback_failure(GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR)
                    return
                end if
                ok = .true.
                return
            end if
            local_state = state_from_array(state_array)
            call field_model%evaluate(state_array(1:3), local_field, field_status)
            if (field_status /= GC_CYL_SUCCESS) then
                call note_callback_failure(map_field_status(field_status))
                return
            end if
            call return_event_rate_proc(state_array(1:3), local_state, local_field, &
                callback_data, rate, rate_status)
            if (rate_status /= GC_CYL_SUCCESS .or. .not. ieee_is_finite(rate)) then
                call note_callback_failure(GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR)
                return
            end if
            ok = .true.
        end subroutine evaluate_physical_rate

        subroutine advance_without_cut(state_in, start_time, delta_time, state_out, &
                end_time, domain_hit, domain_index, accepted_steps_out, &
                rhs_evaluations_out, stage_status)
            real(dp), intent(in) :: state_in(:), start_time, delta_time
            real(dp), allocatable, intent(out) :: state_out(:)
            real(dp), intent(out) :: end_time
            logical, intent(out) :: domain_hit
            integer, intent(out) :: domain_index
            integer, intent(out) :: accepted_steps_out, rhs_evaluations_out
            integer, intent(out) :: stage_status

            type(vode_state_t) :: disarm_integrator
            type(fortnum_status_t) :: disarm_integration_status
            real(dp) :: stop_time, root_time
            logical :: found_domain

            stage_status = GC_CYL_INTEGRATOR_ERROR
            end_time = start_time
            domain_hit = .false.
            domain_index = 0
            accepted_steps_out = 0
            rhs_evaluations_out = 0
            if (size(state_in) /= 5 .or. delta_time <= 0.0_dp) then
                stage_status = GC_CYL_INVALID_INPUT
                return
            end if
            stop_time = start_time+delta_time
            if (stop_time >= target_time) then
                stage_status = GC_CYL_NO_RETURN
                return
            end if
            call vode_init(disarm_integrator, 5, start_time, state_in)
            disarm_integrator%hmax = maximum_step
            disarm_integrator%max_steps = options%maximum_steps
            root_time = stop_time
            found_domain = .false.
            domain_index = 0
            if (have_domain) then
                call vode_integrate_to(physical_rhs, disarm_integrator, stop_time, &
                    options%relative_tolerance, options%absolute_tolerance, &
                    state_out, disarm_integration_status, event=domain_event, &
                    event_dir=ODE_EVENT_FALLING, &
                    event_tol=make_event_tolerance(stop_time, options), &
                    t_root=root_time, root_found=found_domain, &
                    event_index=domain_index)
            else
                call vode_integrate_to(physical_rhs, disarm_integrator, stop_time, &
                    options%relative_tolerance, options%absolute_tolerance, &
                    state_out, disarm_integration_status)
            end if
            accepted_steps_out = disarm_integrator%nsteps
            rhs_evaluations_out = disarm_integrator%nfev
            if (callback_status /= GC_CYL_SUCCESS) then
                stage_status = callback_status
                return
            end if
            if (disarm_integration_status%code /= FORTNUM_OK) then
                stage_status = GC_CYL_INTEGRATOR_ERROR
                return
            end if
            if (.not. allocated(state_out)) then
                stage_status = GC_CYL_STATE_ERROR
                return
            end if
            if (size(state_out) /= 5) then
                stage_status = GC_CYL_STATE_ERROR
                return
            end if
            if (.not. all(ieee_is_finite(state_out))) then
                stage_status = GC_CYL_STATE_ERROR
                return
            end if
            end_time = root_time
            domain_hit = found_domain
            stage_status = GC_CYL_SUCCESS
        end subroutine advance_without_cut

        subroutine integrate_to_cut(state_in, start_time, stop_time, orientation, &
                state_out, root_time_out, event_index_out, root_found_out, &
                accepted_steps_out, rhs_evaluations_out, stage_status)
            real(dp), intent(in) :: state_in(:), start_time, stop_time
            integer, intent(in) :: orientation
            real(dp), allocatable, intent(out) :: state_out(:)
            real(dp), intent(out) :: root_time_out
            integer, intent(out) :: event_index_out
            logical, intent(out) :: root_found_out
            integer, intent(out) :: accepted_steps_out, rhs_evaluations_out
            integer, intent(out) :: stage_status

            type(vode_state_t) :: stage_integrator
            type(fortnum_status_t) :: stage_integration_status
            real(dp) :: stage_event_tolerance

            stage_status = GC_CYL_INTEGRATOR_ERROR
            root_time_out = start_time
            event_index_out = 0
            root_found_out = .false.
            accepted_steps_out = 0
            rhs_evaluations_out = 0
            if (size(state_in) /= 5) then
                stage_status = GC_CYL_STATE_ERROR
                return
            end if
            if (stop_time <= start_time .or. abs(orientation) /= 1) then
                stage_status = GC_CYL_INVALID_INPUT
                return
            end if
            call vode_init(stage_integrator, 5, start_time, state_in)
            stage_integrator%hmax = maximum_step
            stage_integrator%max_steps = options%maximum_steps
            stage_event_tolerance = make_event_tolerance(stop_time, options)
            if (have_domain) then
                call vode_integrate_to(physical_rhs, stage_integrator, stop_time, &
                    options%relative_tolerance, options%absolute_tolerance, &
                    state_out, stage_integration_status, event=physical_event, &
                    event_dir=map_event_direction(orientation), &
                    event_tol=stage_event_tolerance, t_root=root_time_out, &
                    root_found=root_found_out, event_index=event_index_out, &
                    event2=domain_event, event_dir2=ODE_EVENT_FALLING)
            else
                call vode_integrate_to(physical_rhs, stage_integrator, stop_time, &
                    options%relative_tolerance, options%absolute_tolerance, &
                    state_out, stage_integration_status, event=physical_event, &
                    event_dir=map_event_direction(orientation), &
                    event_tol=stage_event_tolerance, t_root=root_time_out, &
                    root_found=root_found_out, event_index=event_index_out)
            end if
            accepted_steps_out = stage_integrator%nsteps
            rhs_evaluations_out = stage_integrator%nfev
            if (callback_status /= GC_CYL_SUCCESS) then
                stage_status = callback_status
                return
            end if
            if (stage_integration_status%code /= FORTNUM_OK) then
                stage_status = GC_CYL_INTEGRATOR_ERROR
                return
            end if
            if (.not. root_found_out) then
                stage_status = GC_CYL_NO_RETURN
                return
            end if
            if (.not. allocated(state_out)) then
                stage_status = GC_CYL_STATE_ERROR
                return
            end if
            if (size(state_out) /= 5) then
                stage_status = GC_CYL_STATE_ERROR
                return
            end if
            if (.not. all(ieee_is_finite(state_out))) then
                stage_status = GC_CYL_STATE_ERROR
                return
            end if
            stage_status = GC_CYL_SUCCESS
        end subroutine integrate_to_cut

    end subroutine compute_gc_cylindrical_physical_return

    subroutine certify_gc_cylindrical_physical_return(result, provider, status, &
            message)
        type(gc_cylindrical_physical_return_t), intent(inout) :: result
        procedure(gc_cylindrical_physical_return_multiplicity_provider) :: provider
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        type(gc_cylindrical_physical_return_certificate_t) :: certificate
        integer :: provider_status

        result%intersection_multiplicity_certified = .false.
        result%multiplicity_status = GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_UNKNOWN
        certificate = gc_cylindrical_physical_return_certificate_t()
        status = GC_CYL_PHYSICAL_RETURN_CERTIFICATE_INVALID
        message = 'numerical return has no valid multiplicity certificate'
        if (result%status /= GC_CYL_SUCCESS .or. &
                .not. result%physical_return_found) return
        call provider(result, certificate, provider_status)
        if (provider_status == GC_CYL_PHYSICAL_RETURN_CERTIFICATE_UNAVAILABLE) then
            status = provider_status
            message = 'exactly-two theorem provider is unavailable'
            return
        end if
        if (provider_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_PHYSICAL_RETURN_CERTIFICATE_INVALID
            message = 'exactly-two theorem provider failed'
            return
        end if
        call attach_gc_cylindrical_physical_return_certificate(result, certificate, &
            status, message)
    end subroutine certify_gc_cylindrical_physical_return

    subroutine attach_gc_cylindrical_physical_return_certificate(result, &
            certificate, status, message)
        type(gc_cylindrical_physical_return_t), intent(inout) :: result
        type(gc_cylindrical_physical_return_certificate_t), intent(in) :: certificate
        integer, intent(out) :: status
        character(len=*), intent(out) :: message

        status = GC_CYL_PHYSICAL_RETURN_CERTIFICATE_INVALID
        message = 'exactly-two certificate does not match return evidence'
        result%intersection_multiplicity_certified = .false.
        result%multiplicity_status = GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_UNKNOWN
        if (result%status /= GC_CYL_SUCCESS .or. &
                .not. result%physical_return_found) return
        if (result%event_kind /= GC_CYL_PHYSICAL_EVENT_RETURN .or. &
                result%wall_hit .or. result%radial_domain_hit .or. &
                result%numerical_failure .or. result%return_orientation == 0) return
        if (result%intersection_count /= 2) return
        if (certificate%certificate_id <= 0 .or. &
                certificate%crossing_count /= 2 .or. &
                .not. certificate%exactly_two_proved) return
        if (abs(result%intersection_orientations(1)) /= 1 .or. &
                abs(result%intersection_orientations(2)) /= 1 .or. &
                result%intersection_orientations(1) == &
                result%intersection_orientations(2) .or. &
                result%intersection_orientations(2) /= result%return_orientation .or. &
                result%intersection_times(2) <= result%intersection_times(1) .or. &
                result%intersection_times(1) <= 0.0_dp .or. &
                result%period < result%intersection_times(2) .or. &
                .not. all(ieee_is_finite(result%intersection_times)) .or. &
                .not. all(ieee_is_finite(result%intersection_rates)) .or. &
                any(abs(result%intersection_rates) <= tiny(1.0_dp))) return

        result%intersection_multiplicity_certified = .true.
        result%multiplicity_status = GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_CERTIFIED
        status = GC_CYL_SUCCESS
        message = 'independent exactly-two theorem certificate attached'
    end subroutine attach_gc_cylindrical_physical_return_certificate

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
                options%radial_distance_scale, options%cut_rate_tolerance]))) return
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
        if (options%cut_rate_tolerance < 0.0_dp) return
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
