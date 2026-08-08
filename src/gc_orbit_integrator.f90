module neort_gc_orbit_integrator
    !! Full axisymmetric guiding-center return map and strict thin-orbit limit.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode, only: ODE_EVENT_RISING, ODE_EVENT_FALLING
    use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use neort_gc_dynamics, only: GC_SUCCESS, gc_field_sample_t, gc_rhs
    use neort_gc_models, only: GC_MODEL_SUCCESS, gc_field_t, gc_potential_t, &
        gc_invariants_t, state_from_invariants, canonical_flux_from_state
    use neort_full_fow_cycle_average_symbolic, only: &
        evaluate_neort_full_fow_cycle_average
    use neort_full_fow_resonance_scalar_symbolic, only: &
        evaluate_neort_full_fow_resonance_scalar
    use neort_thin_orbit_limit, only: orbit_return_t
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
    !! This status is reserved for a separately certified open-orbit result.
    !! A max-periods timeout is GC_ORBIT_NO_RETURN and is never unconfined.
    integer, parameter, public :: GC_ORBIT_UNCONFINED = 7
    !! This status may be returned only by a wall-polygon classifier.
    integer, parameter, public :: GC_ORBIT_WALL_LOSS = 8
    !! This status denotes the computational equilibrium-domain guard.  It is
    !! unresolved unless a physical wall/domain contract certifies it.
    integer, parameter, public :: GC_ORBIT_RADIAL_DOMAIN = 9

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

    public :: compute_return_map
    public :: compute_zero_width_passing_cycle
    public :: compute_gc_full_orbit_average
    public :: normalized_full_hamiltonian_factor
    public :: gc_orbit_perturbation_i
    public :: gc_orbit_status_is_physical_loss
    public :: gc_orbit_status_is_unconfined
    public :: gc_orbit_status_is_wall_loss
    public :: gc_orbit_status_is_radial_domain

contains

    pure logical function gc_orbit_status_is_unconfined(status)
        integer, intent(in) :: status

        gc_orbit_status_is_unconfined = status == GC_ORBIT_UNCONFINED
    end function gc_orbit_status_is_unconfined

    pure logical function gc_orbit_status_is_wall_loss(status)
        integer, intent(in) :: status

        gc_orbit_status_is_wall_loss = status == GC_ORBIT_WALL_LOSS
    end function gc_orbit_status_is_wall_loss

    pure logical function gc_orbit_status_is_radial_domain(status)
        integer, intent(in) :: status

        gc_orbit_status_is_radial_domain = status == GC_ORBIT_RADIAL_DOMAIN
    end function gc_orbit_status_is_radial_domain

    pure logical function gc_orbit_status_is_physical_loss(status)
        integer, intent(in) :: status

        gc_orbit_status_is_physical_loss = gc_orbit_status_is_wall_loss(status)
    end function gc_orbit_status_is_physical_loss

    pure integer function radial_domain_status(radial, options)
        real(dp), intent(in) :: radial
        type(gc_orbit_options_t), intent(in) :: options

        radial_domain_status = GC_ORBIT_SUCCESS
        if (radial <= options%radial_min .or. radial >= options%radial_max) then
            radial_domain_status = GC_ORBIT_RADIAL_DOMAIN
        end if
    end function radial_domain_status

    pure function normalized_full_hamiltonian_factor(invariants, potential, &
            bmod, reference_kinetic_energy) result(value)
        !! POTATO normalization: [2(H-Phi)-J_perp B]/p_ref^2.
        type(gc_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: potential, bmod, reference_kinetic_energy
        real(dp) :: value

        value = (2.0_dp*(invariants%energy - potential) &
            -invariants%magnetic_moment*bmod)/reference_kinetic_energy
    end function normalized_full_hamiltonian_factor

    subroutine compute_zero_width_passing_cycle(field_model, potential_model, &
            invariants, reference_position, parallel_sign, reference_velocity, &
            winding, result, convergence_error)
        !! Independent lambda=0 full-cycle oracle for a direct flux-surface
        !! section.  The theta quadrature is deliberately separate from the
        !! VODE event return map used by the finite-width ladder.
        !! If requested, convergence_error reports the final maximum relative
        !! change of the period and toroidal increment under panel doubling.
        class(gc_field_t), intent(in) :: field_model
        class(gc_potential_t), intent(in) :: potential_model
        type(gc_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: reference_position(3), reference_velocity
        integer, intent(in) :: parallel_sign, winding
        type(orbit_return_t), intent(out) :: result
        real(dp), intent(out), optional :: convergence_error

        integer, parameter :: initial_panels = 64
        integer, parameter :: maximum_refinement = 9
        real(dp), parameter :: quadrature_tolerance = 2.0e-12_dp
        real(dp) :: previous_period_u, previous_delta_phi
        real(dp) :: period_u, delta_phi
        real(dp) :: nested_error
        integer :: refinement, n_panels, local_status
        logical :: converged

        result = orbit_return_t()
        result%orbit_class = GC_ORBIT_PASSING
        result%winding = winding
        if (present(convergence_error)) convergence_error = huge(1.0_dp)
        if (reference_velocity <= 0.0_dp .or. parallel_sign == 0 &
            .or. abs(winding) /= 1) then
            result%status = GC_ORBIT_STATE_ERROR
            return
        end if

        previous_period_u = 0.0_dp
        previous_delta_phi = 0.0_dp
        converged = .false.
        do refinement = 0, maximum_refinement
            n_panels = initial_panels*2**refinement
            call integrate_cycle(n_panels, period_u, delta_phi, local_status)
            if (local_status /= GC_ORBIT_SUCCESS) then
                result%status = local_status
                return
            end if
            if (refinement > 0) then
                nested_error = max( &
                    abs(period_u - previous_period_u) &
                    /max(1.0_dp, abs(period_u)), &
                    abs(delta_phi - previous_delta_phi) &
                    /max(1.0_dp, abs(delta_phi)))
                if (present(convergence_error)) convergence_error = nested_error
                converged = nested_error <= quadrature_tolerance
                if (converged) exit
            end if
            previous_period_u = period_u
            previous_delta_phi = delta_phi
        end do
        if (.not. converged) then
            result%status = GC_ORBIT_INTEGRATOR_ERROR
            return
        end if

        result%period = period_u/reference_velocity
        result%delta_phi = delta_phi
        result%status = GC_ORBIT_SUCCESS

    contains

        subroutine integrate_cycle(n, period_value, delta_phi_value, status)
            integer, intent(in) :: n
            real(dp), intent(out) :: period_value, delta_phi_value
            integer, intent(out) :: status

            type(gc_field_sample_t) :: sample
            real(dp) :: dtheta, theta, p, xi, potential, gradient(3)
            real(dp) :: time_integrand, phi_integrand, time_sum, phi_sum
            integer :: j, field_status, potential_status, state_status
            integer :: weight

            period_value = 0.0_dp
            delta_phi_value = 0.0_dp
            status = GC_ORBIT_SUCCESS
            dtheta = 2.0_dp*pi*real(winding, dp)/real(n, dp)
            time_sum = 0.0_dp
            phi_sum = 0.0_dp
            do j = 0, n
                theta = reference_position(3) + real(j, dp)*dtheta
                call field_model%evaluate([reference_position(1), &
                    reference_position(2), theta], sample, field_status)
                if (field_status /= GC_MODEL_SUCCESS) then
                    status = GC_ORBIT_FIELD_ERROR
                    return
                end if
                ! The zero-width expression is the tangent pullback on the
                ! fixed-s section and uses only h^phi/h^theta.  A smooth
                ! cubic EQDSK interpolation leaves a small h^s residual on
                ! that section; rejecting it at machine epsilon incorrectly
                ! turns a valid tangent quadrature into a state failure.
                call potential_model%evaluate([reference_position(1), &
                    reference_position(2), theta], sample, potential, gradient, &
                    potential_status)
                if (potential_status /= GC_MODEL_SUCCESS) then
                    status = GC_ORBIT_FIELD_ERROR
                    return
                end if
                call state_from_invariants(sample, potential, invariants, &
                    parallel_sign, p, xi, state_status)
                if (state_status /= GC_MODEL_SUCCESS &
                    .or. abs(xi*sample%hcon(3)) <= tiny(xi)) then
                    status = GC_ORBIT_STATE_ERROR
                    return
                end if
                time_integrand = 1.0_dp/(p*xi*sample%hcon(3))
                phi_integrand = sample%hcon(2)/sample%hcon(3)
                weight = merge(4, 2, mod(j, 2) == 1)
                if (j == 0 .or. j == n) weight = 1
                time_sum = time_sum + real(weight, dp)*time_integrand
                phi_sum = phi_sum + real(weight, dp)*phi_integrand
            end do
            period_value = abs(dtheta*time_sum/3.0_dp)
            delta_phi_value = dtheta*phi_sum/3.0_dp
            if (period_value <= 0.0_dp) status = GC_ORBIT_STATE_ERROR
        end subroutine integrate_cycle

    end subroutine compute_zero_width_passing_cycle

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
            rhs_status = radial_domain_status(state(1), options)
            if (rhs_status /= GC_ORBIT_SUCCESS) return
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

    subroutine compute_gc_full_orbit_average(field_model, potential_model, invariants, &
            reference_position, parallel_sign, rho0, reference_velocity, eta, &
            orbit_class, winding, period_estimate, omega_b, omega_phi, mth, mph, perturbation, &
            options, result)
        !! Average a perturbation and field moments on the physical finite-width
        !! return map. The Fourier phase uses the generated temporal
        !! bounce-harmonic resonance scalar. This routine deliberately has no
        !! q reduction; a failed return remains a failed orbit.
        class(gc_field_t), intent(in) :: field_model
        class(gc_potential_t), intent(in) :: potential_model
        type(gc_invariants_t), intent(in) :: invariants
        real(dp), intent(in) :: reference_position(3)
        integer, intent(in) :: parallel_sign
        real(dp), intent(in) :: rho0, reference_velocity, eta, period_estimate
        real(dp), intent(in) :: omega_b, omega_phi
        integer, intent(in) :: orbit_class, winding, mth, mph
        procedure(gc_orbit_perturbation_i) :: perturbation
        type(gc_orbit_options_t), intent(in) :: options
        type(gc_orbit_average_t), intent(out) :: result

        type(vode_state_t) :: integrator
        type(fortnum_status_t) :: integration_status
        type(gc_field_sample_t) :: reference_sample
        real(dp) :: initial_state(5), initial_derivative(5), initial_augmented(9)
        real(dp), allocatable :: final_state(:)
        real(dp) :: atol(9), maximum_time, return_time, event_time_tolerance
        real(dp) :: reference_potential, reference_gradient(3), reference_kinetic_energy
        real(dp) :: average_real, average_imag, residence_average, field_average
        real(dp) :: phase_rate, unused_phase_rate_derivative
        integer :: start_status, rhs_status, event_direction, reference_status
        logical :: found

        result = gc_orbit_average_t()
        if (reference_velocity <= 0.0_dp .or. period_estimate <= 0.0_dp &
            .or. rho0 == 0.0_dp .or. parallel_sign == 0) return
        if (orbit_class /= GC_ORBIT_TRAPPED &
            .and. orbit_class /= GC_ORBIT_PASSING) return
        if (orbit_class == GC_ORBIT_TRAPPED .and. winding /= 0) return
        if (orbit_class == GC_ORBIT_PASSING .and. abs(winding) /= 1) return

        call field_model%evaluate(reference_position, reference_sample, reference_status)
        if (reference_status /= GC_MODEL_SUCCESS) return
        call potential_model%evaluate(reference_position, reference_sample, &
            reference_potential, reference_gradient, reference_status)
        if (reference_status /= GC_MODEL_SUCCESS) return
        reference_kinetic_energy = invariants%energy - reference_potential
        if (reference_kinetic_energy <= 0.0_dp) return

        call initialize_fixed_invariants(field_model, potential_model, &
            invariants, reference_position, parallel_sign, rho0, 1.0_dp, &
            options, initial_state, start_status)
        if (start_status /= GC_ORBIT_SUCCESS) then
            result%status = start_status
            return
        end if
        call evaluate_neort_full_fow_resonance_scalar(real(mth, dp), real(mph, dp), &
            omega_b, omega_phi, 0.0_dp, 0.0_dp, phase_rate, &
            unused_phase_rate_derivative)
        if (.not. ieee_is_finite(phase_rate)) then
            result%status = GC_ORBIT_STATE_ERROR
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
        call vode_integrate_to(full_average_rhs, integrator, maximum_time, &
            options%relative_tolerance, atol, final_state, integration_status, &
            event=full_return_event, event_dir=event_direction, &
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
        if (.not. found .or. .not. ieee_is_finite(return_time) .or. &
                return_time <= 0.0_dp) then
            result%status = GC_ORBIT_NO_RETURN
            return
        end if

        result%period = return_time/reference_velocity
        call evaluate_neort_full_fow_cycle_average(final_state(6), final_state(7), &
            final_state(8), final_state(9), return_time, average_real, average_imag, &
            residence_average, field_average)
        result%perturbation_average = cmplx(average_real, average_imag, dp)
        result%inverse_b_average = residence_average
        result%b_average = field_average
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
            local_status = radial_domain_status(state(1), options)
            if (local_status /= GC_ORBIT_SUCCESS) return
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
            call gc_rhs(sample, grad_potential, rho0, 1.0_dp, state(4), &
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

        subroutine full_average_rhs(time, state, derivative, context)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state(:)
            real(dp), intent(out) :: derivative(:)
            class(*), intent(in), optional :: context

            type(gc_field_sample_t) :: sample
            real(dp) :: base_derivative(5), phase_argument, local_potential
            real(dp) :: local_gradient(3), hamiltonian_factor
            complex(dp) :: amplitude, phase, hamiltonian
            integer :: local_status, field_status, potential_status, perturbation_status

            associate (unused_time => time, unused_context => context)
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
            if (perturbation_status /= GC_ORBIT_SUCCESS) then
                rhs_status = GC_ORBIT_PERTURBATION_ERROR
                return
            end if
            call potential_model%evaluate(state(1:3), sample, local_potential, &
                local_gradient, potential_status)
            if (potential_status /= GC_MODEL_SUCCESS) then
                rhs_status = GC_ORBIT_FIELD_ERROR
                return
            end if
            ! Bounce harmonics are temporal Fourier coefficients, not spatial
            ! poloidal Fourier modes. This is POTATO's canonical phase and is
            ! valid for trapped, passing, shaped, and finite-width orbits.
            phase_argument = real(mph, dp) &
                *(state(2) - reference_position(2)) &
                -phase_rate*time/reference_velocity
            phase = cmplx(cos(phase_argument), sin(phase_argument), dp)
            hamiltonian_factor = normalized_full_hamiltonian_factor(invariants, &
                local_potential, sample%bmod, reference_kinetic_energy)
            hamiltonian = hamiltonian_factor*amplitude*phase
            derivative(1:5) = base_derivative
            derivative(6) = real(hamiltonian)
            derivative(7) = aimag(hamiltonian)
            derivative(8) = 1.0_dp/sample%bmod
            derivative(9) = sample%bmod
        end subroutine full_average_rhs

        function full_return_event(time, state, context) result(value)
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
        end function full_return_event

    end subroutine compute_gc_full_orbit_average

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
        ! The finite-width first guess can lie outside the adapter's
        ! computational chart.  Clamp the numerical launch guess and solve
        ! for a represented fixed-invariant state; this is not a loss claim.
        radial = min(options%radial_max, max(options%radial_min, radial))
        finite_step = 1.0e-5_dp*(options%radial_max - options%radial_min)
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
            radial = min(options%radial_max, max(options%radial_min, radial - step))
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
