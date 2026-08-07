module neort_gc_frequency_provider
    !! Real-space canonical-frequency provider for the direct GEQDSK backend.
    !! The transport layer sees only strict first-order frequencies; finite-
    !! orbit trajectories are used as a centered numerical limiting device.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_dynamics, only: gc_field_sample_t
    use neort_gc_eqdsk_adapter, only: eqdsk_gc_field_t
    use neort_gc_models, only: GC_MODEL_SUCCESS, gc_invariants_t, &
        gc_zero_potential_t, gc_linear_flux_potential_t, &
        invariants_from_state, make_linear_flux_potential
    use neort_gc_orbit_integrator, only: GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, &
        GC_ORBIT_SUCCESS, &
        gc_orbit_options_t, gc_orbit_average_t, gc_orbit_perturbation_i, &
        compute_return_map, compute_thin_precession, compute_gc_orbit_average, &
        compute_gc_full_orbit_average
    use neort_thin_orbit_limit, only: THIN_LIMIT_SUCCESS, orbit_return_t, &
        thin_limit_result_t
    use util, only: pi, c

    implicit none
    private

    integer, parameter, public :: GC_FREQUENCY_SUCCESS = 0
    integer, parameter, public :: GC_FREQUENCY_INVALID_INPUT = 1
    integer, parameter, public :: GC_FREQUENCY_FIELD_ERROR = 2
    integer, parameter, public :: GC_FREQUENCY_LIMIT_ERROR = 3
    integer, parameter, public :: GC_FREQUENCY_ORBIT_ERROR = 4

    type, public :: gc_frequency_context_t
        type(eqdsk_gc_field_t) :: field
        type(gc_zero_potential_t) :: zero_potential
        type(gc_linear_flux_potential_t) :: electric_potential
        type(gc_orbit_options_t) :: orbit_options
        type(gc_field_sample_t) :: reference_sample
        real(dp) :: reference_position(3) = 0.0_dp
        real(dp) :: reference_velocity = 0.0_dp
        real(dp) :: rho0 = 0.0_dp
        real(dp) :: q_fieldline = 0.0_dp
        integer :: htheta_sign = 0
        logical :: initialized = .false.
    end type gc_frequency_context_t

    ! Multiplier of the physical guiding-centre width used only by the
    ! finite-width model.  One is production physics; zero is the exact
    ! real-space thin-orbit limiting oracle.
    real(dp), public :: gc_full_orbit_width_scale = 1.0_dp

    type, public :: gc_frequency_result_t
        real(dp) :: omega_b = 0.0_dp
        real(dp) :: omega_magnetic = 0.0_dp
        real(dp) :: omega_electric = 0.0_dp
        real(dp) :: q_fieldline = 0.0_dp
        real(dp) :: magnetic_error = 0.0_dp
        real(dp) :: electric_error = 0.0_dp
        real(dp) :: magnetic_order = 0.0_dp
        real(dp) :: total_order = 0.0_dp
        real(dp) :: baseline_residual = 0.0_dp
        real(dp) :: lambda_used(4) = 0.0_dp
        integer :: maximum_refinements = 0
        integer :: magnetic_limit_status = 0
        integer :: total_limit_status = 0
    end type gc_frequency_result_t

    type, public :: gc_full_orbit_frequency_result_t
        !! Native finite-width canonical frequencies.  Unlike the thin-limit
        !! result, Omega_phi includes field-line transit, magnetic drift, and
        !! electric drift exactly once through the full return map.
        real(dp) :: omega_b = 0.0_dp
        real(dp) :: omega_phi = 0.0_dp
        real(dp) :: period = 0.0_dp
        real(dp) :: delta_phi = 0.0_dp
        integer :: orbit_status = 0
    end type gc_full_orbit_frequency_result_t

    public :: initialize_gc_frequency_context, evaluate_gc_frequency
    public :: evaluate_gc_full_orbit_frequency
    public :: evaluate_gc_phase_average
    public :: evaluate_gc_full_orbit_phase_average

contains

    subroutine evaluate_gc_full_orbit_frequency(context, eta, &
            parallel_direction, orbit_class, period_estimate, result, status, velocity)
        !! Evaluate one physical-width guiding-center return at fixed
        !! (H, mu, P_phi).  There is deliberately no pitch spline or thin-
        !! orbit velocity scaling here: callers must retain non-return status.
        type(gc_frequency_context_t), intent(in) :: context
        real(dp), intent(in) :: eta, period_estimate
        integer, intent(in) :: parallel_direction, orbit_class
        type(gc_full_orbit_frequency_result_t), intent(out) :: result
        integer, intent(out) :: status
        real(dp), intent(in), optional :: velocity

        type(gc_invariants_t) :: invariants
        type(orbit_return_t) :: orbit_return
        real(dp) :: xi_squared, potential, grad_potential(3), speed_ratio
        integer :: invariant_status, parallel_sign, winding, potential_status

        result = gc_full_orbit_frequency_result_t()
        status = GC_FREQUENCY_INVALID_INPUT
        if (.not. context%initialized .or. eta <= 0.0_dp &
            .or. period_estimate <= 0.0_dp) return
        if (abs(parallel_direction) /= 1) return
        if (orbit_class /= GC_ORBIT_TRAPPED &
            .and. orbit_class /= GC_ORBIT_PASSING) return

        speed_ratio = 1.0_dp
        if (present(velocity)) speed_ratio = velocity/context%reference_velocity
        if (speed_ratio <= 0.0_dp) return

        xi_squared = 1.0_dp - eta*context%reference_sample%bmod
        if (xi_squared <= 0.0_dp) return
        parallel_sign = parallel_direction*context%htheta_sign
        winding = merge(parallel_direction, 0, orbit_class == GC_ORBIT_PASSING)
        call context%electric_potential%evaluate(context%reference_position, &
            context%reference_sample, potential, grad_potential, &
            potential_status)
        if (potential_status /= GC_MODEL_SUCCESS) then
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if
        call invariants_from_state(context%reference_sample, potential, &
            context%rho0, gc_full_orbit_width_scale, speed_ratio, &
            real(parallel_sign, dp)*sqrt(xi_squared), invariants, &
            invariant_status)
        if (invariant_status /= GC_MODEL_SUCCESS) return

        call compute_return_map(context%field, context%electric_potential, &
            invariants, context%reference_position, parallel_sign, &
            context%rho0, gc_full_orbit_width_scale, context%reference_velocity, orbit_class, &
            winding, period_estimate, context%orbit_options, orbit_return)
        result%orbit_status = orbit_return%status
        if (orbit_return%status /= GC_ORBIT_SUCCESS) then
            status = GC_FREQUENCY_ORBIT_ERROR
            return
        end if
        result%period = orbit_return%period
        result%delta_phi = orbit_return%delta_phi
        result%omega_b = 2.0_dp*pi/orbit_return%period
        if (orbit_class == GC_ORBIT_PASSING) then
            result%omega_b = real(parallel_direction, dp)*result%omega_b
        end if
        result%omega_phi = orbit_return%delta_phi/orbit_return%period
        status = GC_FREQUENCY_SUCCESS
    end subroutine evaluate_gc_full_orbit_frequency

    subroutine initialize_gc_frequency_context(surface, reference_theta, &
            field_scale, omega_e, &
            particle_mass, particle_charge, reference_velocity, context, status)
        real(dp), intent(in) :: surface, reference_theta, field_scale, omega_e
        real(dp), intent(in) :: particle_mass, particle_charge, reference_velocity
        type(gc_frequency_context_t), intent(out) :: context
        integer, intent(out) :: status

        integer, parameter :: ntheta = 2048
        type(gc_field_sample_t) :: sample
        real(dp) :: theta, dtheta
        integer :: k, field_status

        context%field%field_scale = 1.0_dp
        context%electric_potential = gc_linear_flux_potential_t()
        context%orbit_options = gc_orbit_options_t()
        context%reference_sample = gc_field_sample_t()
        context%reference_position = 0.0_dp
        context%reference_velocity = 0.0_dp
        context%rho0 = 0.0_dp
        context%q_fieldline = 0.0_dp
        context%htheta_sign = 0
        context%initialized = .false.
        status = GC_FREQUENCY_INVALID_INPUT
        if (surface <= 0.0_dp .or. surface >= 1.0_dp) return
        if (field_scale <= 0.0_dp .or. particle_mass <= 0.0_dp &
            .or. particle_charge == 0.0_dp &
            .or. reference_velocity <= 0.0_dp) return

        context%field%field_scale = field_scale
        ! Use one fixed physical Poincare section.  In shaped equilibria the
        ! field-strength minimum is generally not at geometric theta=0.
        context%reference_position = [surface, 0.0_dp, reference_theta]
        context%reference_velocity = reference_velocity
        context%rho0 = particle_mass*c*reference_velocity/particle_charge
        call context%field%evaluate(context%reference_position, &
            context%reference_sample, field_status)
        if (field_status /= GC_MODEL_SUCCESS) then
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if
        if (abs(context%reference_sample%hcon(3)) &
            <= tiny(context%reference_sample%hcon(3))) then
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if
        context%htheta_sign = merge(1, -1, &
            context%reference_sample%hcon(3) > 0.0_dp)

        ! q(s*) is defined by the same zero-width field-line return used in
        ! the limiting expression, not by a separately interpolated profile.
        dtheta = 2.0_dp*pi/real(ntheta, dp)
        context%q_fieldline = 0.0_dp
        do k = 1, ntheta
            theta = -pi + (real(k, dp) - 0.5_dp)*dtheta
            call context%field%evaluate([surface, 0.0_dp, theta], sample, &
                field_status)
            if (field_status /= GC_MODEL_SUCCESS &
                .or. abs(sample%hcon(3)) <= tiny(sample%hcon(3))) then
                status = GC_FREQUENCY_FIELD_ERROR
                return
            end if
            context%q_fieldline = context%q_fieldline &
                + sample%hcon(2)/sample%hcon(3)*dtheta/(2.0_dp*pi)
        end do
        if (abs(context%q_fieldline) <= tiny(context%q_fieldline)) then
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if

        context%electric_potential = make_linear_flux_potential(omega_e, &
            particle_charge, particle_mass, reference_velocity, &
            context%reference_sample%psi)
        context%orbit_options = gc_orbit_options_t()
        context%orbit_options%radial_min = 1.0e-10_dp
        context%orbit_options%radial_max = 1.0_dp - 1.0e-10_dp
        context%orbit_options%topology_from_zero_width_return = .true.
        context%initialized = .true.
        status = GC_FREQUENCY_SUCCESS
    end subroutine initialize_gc_frequency_context

    subroutine evaluate_gc_frequency(context, eta, parallel_direction, &
            orbit_class, period_estimate, result, status)
        type(gc_frequency_context_t), intent(in) :: context
        real(dp), intent(in) :: eta, period_estimate
        integer, intent(in) :: parallel_direction, orbit_class
        type(gc_frequency_result_t), intent(out) :: result
        integer, intent(out) :: status

        type(gc_invariants_t) :: invariants
        type(thin_limit_result_t) :: magnetic, total
        type(orbit_return_t) :: base
        real(dp) :: xi_squared
        integer :: invariant_status, parallel_sign, winding

        result = gc_frequency_result_t()
        status = GC_FREQUENCY_INVALID_INPUT
        if (.not. context%initialized .or. eta <= 0.0_dp &
            .or. period_estimate <= 0.0_dp) return
        if (abs(parallel_direction) /= 1) return
        if (orbit_class /= GC_ORBIT_TRAPPED &
            .and. orbit_class /= GC_ORBIT_PASSING) return

        xi_squared = 1.0_dp - eta*context%reference_sample%bmod
        if (xi_squared <= 0.0_dp) return
        parallel_sign = parallel_direction*context%htheta_sign
        winding = merge(parallel_direction, 0, orbit_class == GC_ORBIT_PASSING)
        call invariants_from_state(context%reference_sample, 0.0_dp, &
            context%rho0, 0.0_dp, 1.0_dp, &
            real(parallel_sign, dp)*sqrt(xi_squared), invariants, &
            invariant_status)
        if (invariant_status /= GC_MODEL_SUCCESS) return

        call compute_thin_precession(context%field, context%zero_potential, &
            invariants, context%reference_position, parallel_sign, &
            context%rho0, context%reference_velocity, context%q_fieldline, &
            orbit_class, winding, period_estimate, context%orbit_options, &
            magnetic, base)
        call compute_thin_precession(context%field, context%electric_potential, &
            invariants, context%reference_position, parallel_sign, &
            context%rho0, context%reference_velocity, context%q_fieldline, &
            orbit_class, winding, period_estimate, context%orbit_options, total)

        result%magnetic_limit_status = magnetic%status
        result%total_limit_status = total%status
        result%q_fieldline = context%q_fieldline
        if (orbit_class == GC_ORBIT_PASSING) then
            result%q_fieldline = base%delta_phi &
                /(2.0_dp*pi*real(winding, dp))
        end if
        result%baseline_residual = magnetic%baseline_residual
        result%magnetic_error = magnetic%error_estimate
        if (magnetic%status == THIN_LIMIT_SUCCESS &
            .and. total%status == THIN_LIMIT_SUCCESS) then
            result%electric_error = magnetic%error_estimate &
                + total%error_estimate
        else
            ! Keep failed-limit diagnostics finite; adding two huge error
            ! sentinels can overflow before the caller applies its boundary
            ! separatrix policy.
            result%electric_error = huge(1.0_dp)
        end if
        result%magnetic_order = magnetic%observed_order
        result%total_order = total%observed_order
        result%lambda_used(1:2) = [minval(magnetic%lambda_used), &
            maxval(magnetic%lambda_used)]
        result%lambda_used(3:4) = [minval(total%lambda_used), &
            maxval(total%lambda_used)]
        result%maximum_refinements = max(magnetic%refinement_count, &
            total%refinement_count)
        if (magnetic%status /= THIN_LIMIT_SUCCESS &
            .or. total%status /= THIN_LIMIT_SUCCESS) then
            status = GC_FREQUENCY_LIMIT_ERROR
            return
        end if

        result%omega_b = 2.0_dp*pi/base%period
        if (orbit_class == GC_ORBIT_PASSING) then
            result%omega_b = real(parallel_direction, dp)*result%omega_b
        end if
        result%omega_magnetic = magnetic%omega
        result%omega_electric = total%omega - magnetic%omega
        status = GC_FREQUENCY_SUCCESS
    end subroutine evaluate_gc_frequency

    subroutine evaluate_gc_phase_average(context, velocity, eta, &
            parallel_direction, orbit_class, period_estimate, omega_b, omega_phi, &
            q_fieldline, mth, mph, &
            perturbation, result, status)
        !! Evaluate the perturbation Hamiltonian on the direct real-space
        !! zero-width orbit.  The caller supplies the native perturbation
        !! evaluator; this layer supplies the same field, potential, invariant,
        !! and sign conventions as the direct frequency provider.
        type(gc_frequency_context_t), intent(in) :: context
        real(dp), intent(in) :: velocity, eta, period_estimate, omega_b
        real(dp), intent(in) :: omega_phi, q_fieldline
        integer, intent(in) :: parallel_direction, orbit_class, mth, mph
        procedure(gc_orbit_perturbation_i) :: perturbation
        type(gc_orbit_average_t), intent(out) :: result
        integer, intent(out) :: status

        type(gc_invariants_t) :: invariants
        real(dp) :: speed_ratio, xi_squared, rho0
        integer :: invariant_status, parallel_sign, winding

        result = gc_orbit_average_t()
        status = GC_FREQUENCY_INVALID_INPUT
        if (.not. context%initialized .or. velocity <= 0.0_dp &
            .or. eta <= 0.0_dp .or. period_estimate <= 0.0_dp) return
        if (abs(parallel_direction) /= 1) return
        if (orbit_class /= GC_ORBIT_TRAPPED &
            .and. orbit_class /= GC_ORBIT_PASSING) return

        speed_ratio = velocity/context%reference_velocity
        rho0 = context%rho0*speed_ratio
        xi_squared = 1.0_dp - eta*context%reference_sample%bmod
        if (xi_squared <= 0.0_dp) return
        parallel_sign = parallel_direction*context%htheta_sign
        winding = merge(parallel_direction, 0, orbit_class == GC_ORBIT_PASSING)
        call invariants_from_state(context%reference_sample, 0.0_dp, rho0, &
            0.0_dp, speed_ratio, real(parallel_sign, dp)*sqrt(xi_squared), &
            invariants, invariant_status)
        if (invariant_status /= GC_MODEL_SUCCESS) return

        call compute_gc_orbit_average(context%field, context%electric_potential, &
            invariants, context%reference_position, parallel_sign, rho0, &
            context%reference_velocity, eta, orbit_class, winding, &
            period_estimate, omega_b, omega_phi, q_fieldline, mth, mph, perturbation, &
            context%orbit_options, result)
        status = result%status
    end subroutine evaluate_gc_phase_average

    subroutine evaluate_gc_full_orbit_phase_average(context, eta, &
            parallel_direction, orbit_class, period_estimate, omega_b, omega_phi, mth, mph, &
            perturbation, result, status, velocity)
        !! Full-width counterpart of evaluate_gc_phase_average.  It shares the
        !! fixed-H, mu, P_phi construction and physical return-map convention
        !! of evaluate_gc_full_orbit_frequency, while retaining the actual
        !! canonical Fourier phase mth*theta + mph*phi.
        type(gc_frequency_context_t), intent(in) :: context
        real(dp), intent(in) :: eta, period_estimate, omega_b, omega_phi
        integer, intent(in) :: parallel_direction, orbit_class, mth, mph
        procedure(gc_orbit_perturbation_i) :: perturbation
        type(gc_orbit_average_t), intent(out) :: result
        integer, intent(out) :: status
        real(dp), intent(in), optional :: velocity

        type(gc_invariants_t) :: invariants
        real(dp) :: xi_squared, potential, grad_potential(3), speed_ratio
        integer :: invariant_status, parallel_sign, winding, potential_status

        result = gc_orbit_average_t()
        status = GC_FREQUENCY_INVALID_INPUT
        if (.not. context%initialized .or. eta <= 0.0_dp &
            .or. period_estimate <= 0.0_dp) return
        if (abs(parallel_direction) /= 1) return
        if (orbit_class /= GC_ORBIT_TRAPPED &
            .and. orbit_class /= GC_ORBIT_PASSING) return

        speed_ratio = 1.0_dp
        if (present(velocity)) speed_ratio = velocity/context%reference_velocity
        if (speed_ratio <= 0.0_dp) return

        xi_squared = 1.0_dp - eta*context%reference_sample%bmod
        if (xi_squared <= 0.0_dp) return
        parallel_sign = parallel_direction*context%htheta_sign
        winding = merge(parallel_direction, 0, orbit_class == GC_ORBIT_PASSING)
        call context%electric_potential%evaluate(context%reference_position, &
            context%reference_sample, potential, grad_potential, &
            potential_status)
        if (potential_status /= GC_MODEL_SUCCESS) then
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if
        call invariants_from_state(context%reference_sample, potential, &
            context%rho0, gc_full_orbit_width_scale, speed_ratio, &
            real(parallel_sign, dp)*sqrt(xi_squared), invariants, &
            invariant_status)
        if (invariant_status /= GC_MODEL_SUCCESS) then
            status = GC_FREQUENCY_INVALID_INPUT
            return
        end if

        call compute_gc_full_orbit_average(context%field, &
            context%electric_potential, invariants, context%reference_position, &
            parallel_sign, context%rho0, gc_full_orbit_width_scale, &
            context%reference_velocity, eta, &
            orbit_class, winding, period_estimate, omega_b, omega_phi, mth, mph, perturbation, &
            context%orbit_options, result)
        status = result%status
    end subroutine evaluate_gc_full_orbit_phase_average

end module neort_gc_frequency_provider
