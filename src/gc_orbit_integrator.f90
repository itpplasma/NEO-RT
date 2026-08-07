module neort_gc_orbit_integrator
    !! Full axisymmetric guiding-center return map and strict thin-orbit limit.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode, only: ODE_EVENT_RISING, ODE_EVENT_FALLING
    use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use neort_gc_dynamics, only: GC_SUCCESS, gc_field_sample_t, gc_rhs
    use neort_gc_models, only: GC_MODEL_SUCCESS, gc_field_t, gc_potential_t, &
        gc_invariants_t, state_from_invariants, canonical_flux_from_state
    use neort_thin_orbit_limit, only: orbit_return_t, thin_limit_result_t, &
        estimate_thin_limit, THIN_LIMIT_SUCCESS, THIN_LIMIT_RETURN_ERROR, &
        THIN_LIMIT_BASELINE_ERROR, THIN_LIMIT_CONVERGENCE_ERROR
    use util, only: pi

    implicit none
    private

    integer, parameter, public :: GC_ORBIT_SUCCESS = 0
    integer, parameter, public :: GC_ORBIT_FIELD_ERROR = 1
    integer, parameter, public :: GC_ORBIT_STATE_ERROR = 2
    integer, parameter, public :: GC_ORBIT_START_ROOT_ERROR = 3
    integer, parameter, public :: GC_ORBIT_INTEGRATOR_ERROR = 4
    integer, parameter, public :: GC_ORBIT_NO_RETURN = 5
    integer, parameter, public :: GC_ORBIT_PERTURBATION_ERROR = 6

    integer, parameter, public :: GC_ORBIT_TRAPPED = 1
    integer, parameter, public :: GC_ORBIT_PASSING = 2

    type, public :: gc_orbit_options_t
        real(dp) :: relative_tolerance = 3.0e-11_dp
        real(dp) :: absolute_tolerance = 1.0e-12_dp
        !! Root-location tolerance as a fraction of the estimated return
        !! period.  The ODE solver expects an absolute tolerance in its
        !! independent variable, which here is distance v_ref*t.
        real(dp) :: event_relative_tolerance = 1.0e-11_dp
        real(dp) :: invariant_root_tolerance = 3.0e-12_dp
        real(dp) :: baseline_relative_tolerance = 3.0e-8_dp
        real(dp) :: radial_min = 1.0e-10_dp
        real(dp) :: radial_max = 1.0_dp - 1.0e-10_dp
        real(dp) :: max_periods = 3.0_dp
        real(dp) :: lambda(3) = [1.0_dp, 0.5_dp, 0.25_dp]
        !! A strict limit is accepted only when the centered ladder displays
        !! its expected second-order behavior and the two Richardson values
        !! agree.  Failed ladders are retried at successively smaller width.
        !! Near-passing real-space returns can require a smaller centered
        !! ladder than the trapped branch. Keep the strict order/tolerance
        !! gate; this only gives the adaptive limit more rungs to resolve the
        !! same derivative instead of rejecting a valid return map.
        integer :: max_limit_refinements = 12
        real(dp) :: minimum_observed_order = 1.25_dp
        real(dp) :: limit_relative_tolerance = 5.0e-3_dp
        real(dp) :: limit_absolute_tolerance = 1.0e-2_dp
        !! The derivative at zero orbit width is evaluated from the linear
        !! variational return map.  The Richardson map remains available as
        !! an independent oracle by setting this switch false.
        logical :: use_variational_limit = .false.
        logical :: use_variational_fallback = .true.
        real(dp) :: variational_state_relative_step = 1.0e-6_dp
        real(dp) :: variational_parameter_step = 1.0e-5_dp
        !! In a numerically generated real-space chart, the nominal fixed-s
        !! surface need not coincide pointwise with the local psi(R,Z) surface.
        !! Use the lambda=0 return as 2*pi*W*q(psi*) in that case.
        logical :: topology_from_zero_width_return = .false.
    end type gc_orbit_options_t

    type, public :: gc_orbit_average_t
        integer :: status = GC_ORBIT_INTEGRATOR_ERROR
        real(dp) :: period = 0.0_dp
        complex(dp) :: perturbation_average = (0.0_dp, 0.0_dp)
        real(dp) :: inverse_b_average = 0.0_dp
        real(dp) :: b_average = 0.0_dp
    end type gc_orbit_average_t

    abstract interface
        subroutine gc_orbit_perturbation_i(position, bmod, amplitude, status)
            import dp
            real(dp), intent(in) :: position(3), bmod
            complex(dp), intent(out) :: amplitude
            integer, intent(out) :: status
        end subroutine gc_orbit_perturbation_i
    end interface

    public :: compute_return_map, compute_thin_precession
    public :: compute_gc_orbit_average
    public :: gc_orbit_perturbation_i

contains

    subroutine compute_thin_precession(field_model, potential_model, &
            invariants, reference_position, parallel_sign, rho0, &
            reference_velocity, q_reference, orbit_class, winding, &
            period_estimate, options, result, base_return, plus_return, &
            minus_return)
        class(gc_field_t), intent(in) :: field_model
        class(gc_potential_t), intent(in) :: potential_model
        type(gc_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: reference_position(3)
        integer, intent(in) :: parallel_sign
        real(dp), intent(in) :: rho0, reference_velocity, q_reference
        integer, intent(in) :: orbit_class, winding
        real(dp), intent(in) :: period_estimate
        type(gc_orbit_options_t), intent(in) :: options
        type(thin_limit_result_t), intent(out) :: result
        type(orbit_return_t), intent(out), optional :: base_return
        type(orbit_return_t), intent(out), optional :: plus_return(3), minus_return(3)

        type(orbit_return_t) :: base, plus(3), minus(3)
        real(dp) :: topological_delta, baseline_tolerance, lambda(3)
        real(dp) :: convergence_tolerance
        integer :: attempt, k
        logical :: converged

        if (options%use_variational_limit) then
            call compute_variational_precession(field_model, potential_model, &
                invariants, reference_position, parallel_sign, rho0, &
                reference_velocity, q_reference, orbit_class, winding, &
                period_estimate, options, result, base)
            if (result%status == THIN_LIMIT_SUCCESS) then
                if (present(base_return)) base_return = base
                if (present(plus_return)) plus_return = orbit_return_t()
                if (present(minus_return)) minus_return = orbit_return_t()
                return
            end if
        end if

        call compute_return_map(field_model, potential_model, invariants, &
            reference_position, parallel_sign, rho0, 0.0_dp, &
            reference_velocity, orbit_class, winding, period_estimate, &
            options, base)
        topological_delta = 2.0_dp*pi*real(winding, dp)*q_reference
        if (options%topology_from_zero_width_return &
            .and. orbit_class == GC_ORBIT_PASSING) then
            topological_delta = base%delta_phi
        end if
        baseline_tolerance = options%baseline_relative_tolerance &
            *max(1.0_dp, abs(topological_delta))
        result = thin_limit_result_t()
        converged = .false.
        do attempt = 0, options%max_limit_refinements
            lambda = options%lambda/2.0_dp**attempt
            do k = 1, 3
                call compute_return_map(field_model, potential_model, invariants, &
                    reference_position, parallel_sign, rho0, lambda(k), &
                    reference_velocity, orbit_class, winding, period_estimate, &
                    options, plus(k))
                call compute_return_map(field_model, potential_model, invariants, &
                    reference_position, parallel_sign, rho0, -lambda(k), &
                    reference_velocity, orbit_class, winding, period_estimate, &
                    options, minus(k))
            end do
            call estimate_thin_limit(topological_delta, base, lambda, plus, &
                minus, result, baseline_tolerance)
            result%refinement_count = attempt
            if (result%status /= THIN_LIMIT_SUCCESS) cycle
            convergence_tolerance = options%limit_absolute_tolerance &
                + options%limit_relative_tolerance*abs(result%omega)
            converged = result%error_estimate <= convergence_tolerance &
                .and. result%observed_order >= options%minimum_observed_order
            if (converged) exit
        end do
        if (result%status == THIN_LIMIT_SUCCESS .and. .not. converged) then
            result%status = THIN_LIMIT_CONVERGENCE_ERROR
        end if
        if (options%use_variational_fallback &
            .and. result%status /= THIN_LIMIT_SUCCESS) then
            block
                type(thin_limit_result_t) :: finite_width_result

                finite_width_result = result
                call compute_variational_precession(field_model, potential_model, &
                    invariants, reference_position, parallel_sign, rho0, &
                    reference_velocity, q_reference, orbit_class, winding, &
                    period_estimate, options, result, base)
                if (result%status /= THIN_LIMIT_SUCCESS) then
                    result = finite_width_result
                else
                    if (present(base_return)) base_return = base
                    if (present(plus_return)) plus_return = orbit_return_t()
                    if (present(minus_return)) minus_return = orbit_return_t()
                    return
                end if
            end block
        end if
        if (present(base_return)) base_return = base
        if (present(plus_return)) plus_return = plus
        if (present(minus_return)) minus_return = minus
    end subroutine compute_thin_precession

    subroutine compute_variational_precession(field_model, potential_model, &
            invariants, reference_position, parallel_sign, rho0, &
            reference_velocity, q_reference, orbit_class, winding, &
            period_estimate, options, result, base)
        !! Evaluate d[(Delta phi - topology)/period]/d(lambda) at lambda=0.
        !! The augmented ODE is the exact first variation of the same
        !! coordinate-general guiding-center RHS used by compute_return_map.
        !! Finite differences are used only for the local Jacobian and the
        !! fixed-invariant initial tangent; no orbit-width cancellation occurs.
        class(gc_field_t), intent(in) :: field_model
        class(gc_potential_t), intent(in) :: potential_model
        type(gc_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: reference_position(3)
        integer, intent(in) :: parallel_sign
        real(dp), intent(in) :: rho0, reference_velocity, q_reference
        integer, intent(in) :: orbit_class, winding
        real(dp), intent(in) :: period_estimate
        type(gc_orbit_options_t), intent(in) :: options
        type(thin_limit_result_t), intent(out) :: result
        type(orbit_return_t), intent(out) :: base

        type(vode_state_t) :: integrator
        type(fortnum_status_t) :: integration_status
        real(dp) :: initial_state(5), initial_plus(5), initial_minus(5)
        real(dp) :: initial_augmented(10), derivative(10), atol(10)
        real(dp), allocatable :: final_state(:)
        real(dp) :: final_rhs(5), return_time, maximum_time
        real(dp) :: event_time_tolerance, topology, time_tangent
        real(dp) :: delta_phi_tangent, state_step, parameter_step
        real(dp) :: convergence_tolerance
        integer :: status, k, event_direction, rhs_status
        logical :: found

        result = thin_limit_result_t()
        base = orbit_return_t()
        base%orbit_class = orbit_class
        base%winding = winding
        if (reference_velocity <= 0.0_dp .or. period_estimate <= 0.0_dp &
            .or. rho0 == 0.0_dp .or. parallel_sign == 0) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if
        if (orbit_class /= GC_ORBIT_TRAPPED &
            .and. orbit_class /= GC_ORBIT_PASSING) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if
        if ((orbit_class == GC_ORBIT_TRAPPED .and. winding /= 0) &
            .or. (orbit_class == GC_ORBIT_PASSING .and. abs(winding) /= 1)) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if

        call initialize_fixed_invariants(field_model, potential_model, &
            invariants, reference_position, parallel_sign, rho0, 0.0_dp, &
            options, initial_state, status)
        if (status /= GC_ORBIT_SUCCESS) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if
        parameter_step = options%variational_parameter_step
        if (parameter_step <= 0.0_dp) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if
        call initialize_fixed_invariants(field_model, potential_model, &
            invariants, reference_position, parallel_sign, rho0, parameter_step, &
            options, initial_plus, status)
        if (status /= GC_ORBIT_SUCCESS) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if
        call initialize_fixed_invariants(field_model, potential_model, &
            invariants, reference_position, parallel_sign, rho0, -parameter_step, &
            options, initial_minus, status)
        if (status /= GC_ORBIT_SUCCESS) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if
        initial_augmented(1:5) = initial_state
        initial_augmented(6:10) = (initial_plus - initial_minus) &
            /(2.0_dp*parameter_step)

        rhs_status = GC_ORBIT_SUCCESS
        call variational_rhs(0.0_dp, initial_augmented, derivative)
        if (rhs_status /= GC_ORBIT_SUCCESS) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if
        if (orbit_class == GC_ORBIT_TRAPPED) then
            if (abs(derivative(3)) <= tiny(derivative(3))) then
                result%status = THIN_LIMIT_RETURN_ERROR
                return
            end if
            event_direction = merge(ODE_EVENT_RISING, ODE_EVENT_FALLING, &
                derivative(3) > 0.0_dp)
        else
            event_direction = merge(ODE_EVENT_RISING, ODE_EVENT_FALLING, &
                winding > 0)
        end if

        maximum_time = options%max_periods*period_estimate*reference_velocity
        event_time_tolerance = options%event_relative_tolerance &
            *period_estimate*reference_velocity
        event_time_tolerance = max(event_time_tolerance, &
            100.0_dp*epsilon(maximum_time)*max(1.0_dp, maximum_time))
        atol = options%absolute_tolerance
        call vode_init(integrator, 10, 0.0_dp, initial_augmented)
        call vode_integrate_to(variational_rhs, integrator, maximum_time, &
            options%relative_tolerance, atol, final_state, integration_status, &
            event=return_event, event_dir=event_direction, &
            event_tol=event_time_tolerance, t_root=return_time, &
            root_found=found)
        if (rhs_status /= GC_ORBIT_SUCCESS &
            .or. integration_status%code /= FORTNUM_OK .or. .not. found) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if

        call evaluate_rhs(final_state(1:5), 0.0_dp, final_rhs, status)
        if (status /= GC_ORBIT_SUCCESS .or. abs(final_rhs(3)) <= tiny(final_rhs(3))) then
            result%status = THIN_LIMIT_RETURN_ERROR
            return
        end if
        time_tangent = -final_state(8)/final_rhs(3)
        delta_phi_tangent = final_state(7) + final_rhs(2)*time_tangent
        base%period = return_time/reference_velocity
        base%delta_phi = final_state(2) - initial_state(2)
        base%status = GC_ORBIT_SUCCESS
        topology = 2.0_dp*pi*real(winding, dp)*q_reference
        if (options%topology_from_zero_width_return &
            .and. orbit_class == GC_ORBIT_PASSING) then
            topology = base%delta_phi
        end if
        result%baseline_residual = base%delta_phi - topology
        result%omega = reference_velocity*delta_phi_tangent/return_time
        result%centered = result%omega
        result%richardson_coarse = result%omega
        result%richardson_fine = result%omega
        result%error_estimate = 0.0_dp
        result%observed_order = 2.0_dp
        result%lambda_used = [parameter_step, parameter_step/2.0_dp, &
            parameter_step/4.0_dp]
        result%refinement_count = 0
        convergence_tolerance = options%baseline_relative_tolerance &
            *max(1.0_dp, abs(topology))
        if (abs(result%baseline_residual) > convergence_tolerance) then
            result%status = THIN_LIMIT_BASELINE_ERROR
        else
            result%status = THIN_LIMIT_SUCCESS
        end if

    contains

        subroutine variational_rhs(time, state, derivative, context)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            real(dp), intent(out) :: derivative(:)
            class(*), intent(in), optional :: context

            real(dp) :: f0(5), fplus(5), fminus(5), jacobian(5, 5)
            real(dp) :: perturbed(5), step
            integer :: local_status, j

            associate (unused_time => time, unused_context => context)
            end associate
            derivative = 0.0_dp
            call evaluate_rhs(state(1:5), 0.0_dp, f0, local_status)
            if (local_status /= GC_ORBIT_SUCCESS) then
                rhs_status = local_status
                return
            end if
            call evaluate_rhs(state(1:5), parameter_step, fplus, local_status)
            if (local_status /= GC_ORBIT_SUCCESS) then
                rhs_status = local_status
                return
            end if
            call evaluate_rhs(state(1:5), -parameter_step, fminus, local_status)
            if (local_status /= GC_ORBIT_SUCCESS) then
                rhs_status = local_status
                return
            end if
            derivative(1:5) = f0
            derivative(6:10) = (fplus - fminus)/(2.0_dp*parameter_step)
            do j = 1, 5
                perturbed = state(1:5)
                step = options%variational_state_relative_step &
                    *max(1.0_dp, abs(perturbed(j)))
                if (step <= 0.0_dp) then
                    rhs_status = GC_ORBIT_STATE_ERROR
                    return
                end if
                perturbed(j) = perturbed(j) + step
                call evaluate_rhs(perturbed, 0.0_dp, fplus, local_status)
                if (local_status /= GC_ORBIT_SUCCESS) then
                    rhs_status = local_status
                    return
                end if
                perturbed(j) = state(j) - step
                call evaluate_rhs(perturbed, 0.0_dp, fminus, local_status)
                if (local_status /= GC_ORBIT_SUCCESS) then
                    rhs_status = local_status
                    return
                end if
                jacobian(:, j) = (fplus - fminus)/(2.0_dp*step)
            end do
            derivative(6:10) = derivative(6:10) &
                + matmul(jacobian, state(6:10))
        end subroutine variational_rhs

        subroutine evaluate_rhs(state, orbit_width_scale, derivative, local_status)
            real(dp), intent(in) :: state(5), orbit_width_scale
            real(dp), intent(out) :: derivative(5)
            integer, intent(out) :: local_status
            type(gc_field_sample_t) :: sample
            real(dp) :: potential, gradient(3), xdot(3), pdot, xidot
            integer :: field_status, potential_status, dynamics_status

            derivative = 0.0_dp
            call field_model%evaluate(state(1:3), sample, field_status)
            if (field_status /= GC_MODEL_SUCCESS) then
                local_status = GC_ORBIT_FIELD_ERROR
                return
            end if
            call potential_model%evaluate(state(1:3), sample, potential, &
                gradient, potential_status)
            if (potential_status /= GC_MODEL_SUCCESS) then
                local_status = GC_ORBIT_FIELD_ERROR
                return
            end if
            call gc_rhs(sample, gradient, rho0, orbit_width_scale, state(4), &
                state(5), xdot, pdot, xidot, dynamics_status)
            if (dynamics_status /= GC_SUCCESS) then
                local_status = GC_ORBIT_STATE_ERROR
                return
            end if
            derivative(1:3) = xdot
            derivative(4) = pdot
            derivative(5) = xidot
            local_status = GC_ORBIT_SUCCESS
        end subroutine evaluate_rhs

        function return_event(time, state, context) result(value)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            class(*), intent(in), optional :: context
            real(dp) :: value

            associate (unused_time => time, unused_context => context)
            end associate
            if (orbit_class == GC_ORBIT_PASSING) then
                value = state(3) - reference_position(3) &
                    - 2.0_dp*pi*real(winding, dp)
            else
                value = state(3) - reference_position(3)
            end if
        end function return_event

    end subroutine compute_variational_precession

    subroutine compute_return_map(field_model, potential_model, invariants, &
            reference_position, parallel_sign, rho0, orbit_width_scale, &
            reference_velocity, orbit_class, winding, period_estimate, &
            options, result)
        class(gc_field_t), intent(in) :: field_model
        class(gc_potential_t), intent(in) :: potential_model
        type(gc_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: reference_position(3)
        integer, intent(in) :: parallel_sign
        real(dp), intent(in) :: rho0, orbit_width_scale, reference_velocity
        integer, intent(in) :: orbit_class, winding
        real(dp), intent(in) :: period_estimate
        type(gc_orbit_options_t), intent(in) :: options
        type(orbit_return_t), intent(out) :: result

        type(vode_state_t) :: integrator
        type(fortnum_status_t) :: integration_status
        real(dp) :: initial_state(5), initial_derivative(5), atol(5)
        real(dp), allocatable :: final_state(:)
        real(dp) :: maximum_time, return_time, event_time_tolerance
        integer :: rhs_status, event_direction, start_status
        logical :: found

        result = orbit_return_t()
        result%orbit_class = orbit_class
        result%winding = winding
        if (reference_velocity <= 0.0_dp .or. period_estimate <= 0.0_dp &
            .or. rho0 == 0.0_dp .or. parallel_sign == 0) then
            result%status = GC_ORBIT_STATE_ERROR
            return
        end if
        if (orbit_class /= GC_ORBIT_TRAPPED &
            .and. orbit_class /= GC_ORBIT_PASSING) then
            result%status = GC_ORBIT_STATE_ERROR
            return
        end if
        if (orbit_class == GC_ORBIT_TRAPPED .and. winding /= 0) then
            result%status = GC_ORBIT_STATE_ERROR
            return
        end if
        if (orbit_class == GC_ORBIT_PASSING .and. abs(winding) /= 1) then
            result%status = GC_ORBIT_STATE_ERROR
            return
        end if

        call initialize_fixed_invariants(field_model, potential_model, &
            invariants, reference_position, parallel_sign, rho0, &
            orbit_width_scale, options, initial_state, start_status)
        if (start_status /= GC_ORBIT_SUCCESS) then
            result%status = start_status
            return
        end if

        rhs_status = GC_ORBIT_SUCCESS
        call orbit_rhs(0.0_dp, initial_state, initial_derivative)
        if (rhs_status /= GC_ORBIT_SUCCESS) then
            result%status = rhs_status
            return
        end if
        if (orbit_class == GC_ORBIT_TRAPPED) then
            if (abs(initial_derivative(3)) <= tiny(initial_derivative(3))) then
                result%status = GC_ORBIT_STATE_ERROR
                return
            end if
            event_direction = merge(ODE_EVENT_RISING, ODE_EVENT_FALLING, &
                initial_derivative(3) > 0.0_dp)
        else
            event_direction = merge(ODE_EVENT_RISING, ODE_EVENT_FALLING, &
                winding > 0)
        end if

        atol = options%absolute_tolerance
        maximum_time = options%max_periods*period_estimate*reference_velocity
        event_time_tolerance = options%event_relative_tolerance &
            *period_estimate*reference_velocity
        event_time_tolerance = max(event_time_tolerance, &
            100.0_dp*epsilon(maximum_time)*max(1.0_dp, maximum_time))
        call vode_init(integrator, 5, 0.0_dp, initial_state)
        call vode_integrate_to(orbit_rhs, integrator, maximum_time, &
            options%relative_tolerance, atol, final_state, integration_status, &
            event=return_event, event_dir=event_direction, &
            event_tol=event_time_tolerance, t_root=return_time, &
            root_found=found)
        if (rhs_status /= GC_ORBIT_SUCCESS) then
            result%status = rhs_status
            return
        end if
        if (integration_status%code /= FORTNUM_OK) then
            result%status = GC_ORBIT_INTEGRATOR_ERROR
            return
        end if
        if (.not. found) then
            result%status = GC_ORBIT_NO_RETURN
            return
        end if

        result%period = return_time/reference_velocity
        result%delta_phi = final_state(2) - initial_state(2)
        result%status = GC_ORBIT_SUCCESS

    contains

        subroutine orbit_rhs(time, state, derivative, context)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            real(dp), intent(out) :: derivative(:)
            class(*), intent(in), optional :: context

            type(gc_field_sample_t) :: sample
            real(dp) :: potential, grad_potential(3), xdot(3), pdot, xidot
            integer :: field_status, potential_status, dynamics_status

            associate (unused_time => time, unused_context => context)
            end associate
            derivative = 0.0_dp
            call field_model%evaluate(state(1:3), sample, field_status)
            if (field_status /= GC_MODEL_SUCCESS) then
                rhs_status = GC_ORBIT_FIELD_ERROR
                return
            end if
            call potential_model%evaluate(state(1:3), sample, potential, &
                grad_potential, potential_status)
            if (potential_status /= GC_MODEL_SUCCESS) then
                rhs_status = GC_ORBIT_FIELD_ERROR
                return
            end if
            call gc_rhs(sample, grad_potential, rho0, orbit_width_scale, &
                state(4), state(5), xdot, pdot, xidot, dynamics_status)
            if (dynamics_status /= GC_SUCCESS) then
                rhs_status = GC_ORBIT_STATE_ERROR
                return
            end if
            derivative(1:3) = xdot
            derivative(4) = pdot
            derivative(5) = xidot
        end subroutine orbit_rhs

        function return_event(time, state, context) result(value)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            class(*), intent(in), optional :: context
            real(dp) :: value

            associate (unused_time => time, unused_context => context)
            end associate
            if (orbit_class == GC_ORBIT_PASSING) then
                value = state(3) - reference_position(3) &
                    - 2.0_dp*pi*real(winding, dp)
            else
                value = state(3) - reference_position(3)
            end if
        end function return_event

    end subroutine compute_return_map

    subroutine compute_gc_orbit_average(field_model, potential_model, invariants, &
            reference_position, parallel_sign, rho0, reference_velocity, eta, &
            orbit_class, winding, period_estimate, omega_b, omega_phi, &
            q_fieldline, mth, mph, &
            perturbation, options, result)
        !! Average a native real-space perturbation along the zero-width
        !! guiding-center trajectory.  The trajectory is the same coordinate-
        !! general RHS used by the return map; only its first-order orbit-width
        !! terms are set to zero.  The drift/precession frequency is supplied by
        !! the caller, so this routine never introduces a Boozer reduction or a
        !! second magnetic/electric-drift model.
        class(gc_field_t), intent(in) :: field_model
        class(gc_potential_t), intent(in) :: potential_model
        type(gc_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: reference_position(3)
        integer, intent(in) :: parallel_sign
        real(dp), intent(in) :: rho0, reference_velocity, eta
        integer, intent(in) :: orbit_class, winding, mth, mph
        real(dp), intent(in) :: period_estimate, omega_b, omega_phi, q_fieldline
        procedure(gc_orbit_perturbation_i) :: perturbation
        type(gc_orbit_options_t), intent(in) :: options
        type(gc_orbit_average_t), intent(out) :: result

        type(vode_state_t) :: integrator
        type(fortnum_status_t) :: integration_status
        real(dp) :: initial_state(5), initial_derivative(5)
        real(dp) :: initial_augmented(9), atol(9)
        real(dp), allocatable :: final_state(:)
        real(dp) :: maximum_time, return_time, event_time_tolerance
        real(dp) :: delta_omega_phi
        integer :: start_status, rhs_status, event_direction
        logical :: found

        result = gc_orbit_average_t()
        if (reference_velocity <= 0.0_dp .or. period_estimate <= 0.0_dp &
            .or. rho0 == 0.0_dp .or. parallel_sign == 0) return
        if (orbit_class /= GC_ORBIT_TRAPPED &
            .and. orbit_class /= GC_ORBIT_PASSING) return
        if ((orbit_class == GC_ORBIT_TRAPPED .and. winding /= 0) &
            .or. (orbit_class == GC_ORBIT_PASSING &
            .and. abs(winding) /= 1)) return
        if (orbit_class == GC_ORBIT_PASSING) then
            delta_omega_phi = omega_phi - q_fieldline*omega_b
        else
            delta_omega_phi = omega_phi
        end if

        call initialize_fixed_invariants(field_model, potential_model, &
            invariants, reference_position, parallel_sign, rho0, 0.0_dp, &
            options, initial_state, start_status)
        if (start_status /= GC_ORBIT_SUCCESS) then
            result%status = start_status
            return
        end if

        rhs_status = GC_ORBIT_SUCCESS
        call evaluate_state_rhs(initial_state, initial_derivative, rhs_status)
        if (rhs_status /= GC_ORBIT_SUCCESS) then
            result%status = rhs_status
            return
        end if
        if (orbit_class == GC_ORBIT_TRAPPED) then
            if (abs(initial_derivative(3)) <= tiny(initial_derivative(3))) then
                result%status = GC_ORBIT_STATE_ERROR
                return
            end if
            event_direction = merge(ODE_EVENT_RISING, ODE_EVENT_FALLING, &
                initial_derivative(3) > 0.0_dp)
        else
            event_direction = merge(ODE_EVENT_RISING, ODE_EVENT_FALLING, &
                winding > 0)
        end if

        initial_augmented = 0.0_dp
        initial_augmented(1:5) = initial_state
        maximum_time = options%max_periods*period_estimate*reference_velocity
        event_time_tolerance = options%event_relative_tolerance &
            *period_estimate*reference_velocity
        event_time_tolerance = max(event_time_tolerance, &
            100.0_dp*epsilon(maximum_time)*max(1.0_dp, maximum_time))
        atol = options%absolute_tolerance
        call vode_init(integrator, 9, 0.0_dp, initial_augmented)
        call vode_integrate_to(average_rhs, integrator, maximum_time, &
            options%relative_tolerance, atol, final_state, integration_status, &
            event=return_event, event_dir=event_direction, &
            event_tol=event_time_tolerance, t_root=return_time, &
            root_found=found)
        if (rhs_status /= GC_ORBIT_SUCCESS) then
            result%status = rhs_status
            return
        end if
        if (integration_status%code /= FORTNUM_OK) then
            result%status = GC_ORBIT_INTEGRATOR_ERROR
            return
        end if
        if (.not. found .or. return_time <= 0.0_dp) then
            result%status = GC_ORBIT_NO_RETURN
            return
        end if

        result%period = return_time/reference_velocity
        result%perturbation_average = cmplx(final_state(6), final_state(7), dp) &
            /return_time
        result%inverse_b_average = final_state(8)/return_time
        result%b_average = final_state(9)/return_time
        result%status = GC_ORBIT_SUCCESS

    contains

        subroutine evaluate_state_rhs(state, derivative, local_status)
            real(dp), intent(in) :: state(5)
            real(dp), intent(out) :: derivative(5)
            integer, intent(out) :: local_status

            type(gc_field_sample_t) :: sample
            real(dp) :: potential, grad_potential(3), xdot(3), pdot, xidot
            integer :: field_status, potential_status, dynamics_status

            derivative = 0.0_dp
            call field_model%evaluate(state(1:3), sample, field_status)
            if (field_status /= GC_MODEL_SUCCESS) then
                local_status = GC_ORBIT_FIELD_ERROR
                return
            end if
            call potential_model%evaluate(state(1:3), sample, potential, &
                grad_potential, potential_status)
            if (potential_status /= GC_MODEL_SUCCESS) then
                local_status = GC_ORBIT_FIELD_ERROR
                return
            end if
            call gc_rhs(sample, grad_potential, rho0, 0.0_dp, state(4), &
                state(5), xdot, pdot, xidot, dynamics_status)
            if (dynamics_status /= GC_SUCCESS) then
                local_status = GC_ORBIT_STATE_ERROR
                return
            end if
            derivative(1:3) = xdot
            derivative(4) = pdot
            derivative(5) = xidot
            local_status = GC_ORBIT_SUCCESS
        end subroutine evaluate_state_rhs

        subroutine average_rhs(time, state, derivative, context)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            real(dp), intent(out) :: derivative(:)
            class(*), intent(in), optional :: context

            type(gc_field_sample_t) :: sample
            real(dp) :: base_derivative(5), phase_argument
            complex(dp) :: amplitude, phase, hamiltonian
            integer :: local_status, field_status, perturbation_status

            associate (unused_context => context)
            end associate
            derivative = 0.0_dp
            call evaluate_state_rhs(state(1:5), base_derivative, local_status)
            if (local_status /= GC_ORBIT_SUCCESS) then
                rhs_status = local_status
                return
            end if
            call field_model%evaluate(state(1:3), sample, field_status)
            if (field_status /= GC_MODEL_SUCCESS) then
                rhs_status = GC_ORBIT_FIELD_ERROR
                return
            end if
            call perturbation(state(1:3), sample%bmod, amplitude, &
                perturbation_status)
            if (perturbation_status /= 0) then
                rhs_status = GC_ORBIT_PERTURBATION_ERROR
                return
            end if
            if (orbit_class == GC_ORBIT_PASSING) then
                phase_argument = real(mph, dp) &
                    *(state(2) + delta_omega_phi*time/reference_velocity) &
                    -(real(mth, dp) + real(mph, dp)*q_fieldline) &
                    *omega_b*time/reference_velocity
            else
                phase_argument = real(mph, dp) &
                    *(state(2) + delta_omega_phi*time/reference_velocity) &
                    -real(mth, dp)*omega_b*time/reference_velocity
            end if
            phase = cmplx(cos(phase_argument), sin(phase_argument), dp)
            hamiltonian = (2.0_dp - eta*sample%bmod)*amplitude*phase
            derivative(1:5) = base_derivative
            derivative(6) = real(hamiltonian)
            derivative(7) = aimag(hamiltonian)
            derivative(8) = 1.0_dp/sample%bmod
            derivative(9) = sample%bmod
        end subroutine average_rhs

        function return_event(time, state, context) result(value)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            class(*), intent(in), optional :: context
            real(dp) :: value

            associate (unused_time => time, unused_context => context)
            end associate
            if (orbit_class == GC_ORBIT_PASSING) then
                value = state(3) - reference_position(3) &
                    -2.0_dp*pi*real(winding, dp)
            else
                value = state(3) - reference_position(3)
            end if
        end function return_event

    end subroutine compute_gc_orbit_average

    subroutine initialize_fixed_invariants(field_model, potential_model, &
            invariants, reference_position, parallel_sign, rho0, &
            orbit_width_scale, options, state, status)
        class(gc_field_t), intent(in) :: field_model
        class(gc_potential_t), intent(in) :: potential_model
        type(gc_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: reference_position(3)
        integer, intent(in) :: parallel_sign
        real(dp), intent(in) :: rho0, orbit_width_scale
        type(gc_orbit_options_t), intent(in) :: options
        real(dp), intent(out) :: state(5)
        integer, intent(out) :: status

        type(gc_field_sample_t) :: reference_sample, sample
        real(dp) :: potential, gradient(3), p, xi, residual
        real(dp) :: radial, left, right, fleft, fright, derivative, step
        real(dp) :: finite_step, residual_scale
        integer :: field_status, potential_status, state_status, iteration

        state = 0.0_dp
        status = GC_ORBIT_START_ROOT_ERROR
        call field_model%evaluate(reference_position, reference_sample, field_status)
        if (field_status /= GC_MODEL_SUCCESS) then
            status = GC_ORBIT_FIELD_ERROR
            return
        end if
        call potential_model%evaluate(reference_position, reference_sample, &
            potential, gradient, potential_status)
        if (potential_status /= GC_MODEL_SUCCESS) then
            status = GC_ORBIT_FIELD_ERROR
            return
        end if
        call state_from_invariants(reference_sample, potential, invariants, &
            parallel_sign, p, xi, state_status)
        if (state_status /= GC_MODEL_SUCCESS) then
            status = GC_ORBIT_STATE_ERROR
            return
        end if

        radial = reference_position(1)
        if (orbit_width_scale /= 0.0_dp) then
            if (abs(reference_sample%grad_psi(1)) <= &
                tiny(reference_sample%grad_psi(1))) return
            radial = radial - orbit_width_scale*rho0*p*xi &
                *reference_sample%hcov(2)/reference_sample%grad_psi(1)
        end if
        finite_step = 1.0e-6_dp*(options%radial_max - options%radial_min)
        residual_scale = max(1.0_dp, abs(invariants%canonical_flux))

        do iteration = 1, 30
            call residual_at(radial, residual, sample, potential, p, xi, state_status)
            if (state_status /= GC_ORBIT_SUCCESS) then
                status = state_status
                return
            end if
            if (abs(residual) <= options%invariant_root_tolerance*residual_scale) then
                state(1:3) = [radial, reference_position(2), reference_position(3)]
                state(4:5) = [p, xi]
                status = GC_ORBIT_SUCCESS
                return
            end if

            left = max(options%radial_min, radial - finite_step)
            right = min(options%radial_max, radial + finite_step)
            if (right <= left) return
            call residual_value(left, fleft, state_status)
            if (state_status /= GC_ORBIT_SUCCESS) return
            call residual_value(right, fright, state_status)
            if (state_status /= GC_ORBIT_SUCCESS) return
            derivative = (fright - fleft)/(right - left)
            if (abs(derivative) <= tiny(derivative)) return
            step = residual/derivative
            step = sign(min(abs(step), 0.2_dp*(options%radial_max &
                - options%radial_min)), step)
            radial = radial - step
            if (radial <= options%radial_min .or. radial >= options%radial_max) return
        end do

    contains

        subroutine residual_value(radial_value, value, local_status)
            real(dp), intent(in) :: radial_value
            real(dp), intent(out) :: value
            integer, intent(out) :: local_status
            type(gc_field_sample_t) :: unused_sample
            real(dp) :: unused_potential, unused_p, unused_xi

            call residual_at(radial_value, value, unused_sample, &
                unused_potential, unused_p, unused_xi, local_status)
        end subroutine residual_value

        subroutine residual_at(radial_value, value, local_sample, &
                local_potential, local_p, local_xi, local_status)
            real(dp), intent(in) :: radial_value
            real(dp), intent(out) :: value, local_potential, local_p, local_xi
            type(gc_field_sample_t), intent(out) :: local_sample
            integer, intent(out) :: local_status
            real(dp) :: position(3), local_gradient(3)
            integer :: local_field_status, local_potential_status

            position = [radial_value, reference_position(2), reference_position(3)]
            call field_model%evaluate(position, local_sample, local_field_status)
            if (local_field_status /= GC_MODEL_SUCCESS) then
                local_status = GC_ORBIT_FIELD_ERROR
                return
            end if
            call potential_model%evaluate(position, local_sample, &
                local_potential, local_gradient, local_potential_status)
            if (local_potential_status /= GC_MODEL_SUCCESS) then
                local_status = GC_ORBIT_FIELD_ERROR
                return
            end if
            call state_from_invariants(local_sample, local_potential, invariants, &
                parallel_sign, local_p, local_xi, local_status)
            if (local_status /= GC_MODEL_SUCCESS) then
                local_status = GC_ORBIT_STATE_ERROR
                return
            end if
            value = canonical_flux_from_state(local_sample, rho0, &
                orbit_width_scale, local_p, local_xi) &
                - invariants%canonical_flux
            local_status = GC_ORBIT_SUCCESS
        end subroutine residual_at

    end subroutine initialize_fixed_invariants

end module neort_gc_orbit_integrator
