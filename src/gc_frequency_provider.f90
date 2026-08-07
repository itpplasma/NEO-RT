module neort_gc_frequency_provider
    !! Real-space canonical-frequency provider for the direct GEQDSK backend.
    !! The transport layer sees only strict first-order frequencies; finite-
    !! orbit trajectories are used as a centered numerical limiting device.
    use, intrinsic :: iso_fortran_env, only: dp => real64
#ifdef NEO_RT_USE_STANDALONE
    use do_magfie_mod, only: inp_swi
    use neort_gc_cylindrical_backend, only: &
        GC_CYL_BACKEND_SUCCESS, gc_cylindrical_backend_context_t, &
        GC_CYL_BACKEND_WALL_ERROR, &
        gc_cylindrical_backend_result_t, initialize_gc_cylindrical_backend, &
        evaluate_gc_cylindrical_backend, evaluate_gc_cylindrical_phase_average
    use neort_gc_cylindrical_model, only: GC_CYL_EQUILIBRIUM_DOMAIN, &
        GC_CYL_FIELD_ERROR, GC_CYL_INTEGRATOR_ERROR, GC_CYL_INVARIANT_ERROR, &
        GC_CYL_NO_RETURN, GC_CYL_PERTURBATION_ERROR, GC_CYL_POTENTIAL_ERROR, &
        GC_CYL_SINGULAR_BSTAR, GC_CYL_START_ERROR, GC_CYL_STATE_ERROR, &
        GC_CYL_WALL_ERROR, GC_CYL_WALL_LOSS
#endif
    use neort_gc_dynamics, only: gc_field_sample_t
    use neort_gc_eqdsk_adapter, only: eqdsk_gc_field_t
    use neort_gc_models, only: GC_MODEL_SUCCESS, gc_invariants_t, &
        gc_zero_potential_t, gc_linear_flux_potential_t, &
        invariants_from_state, make_linear_flux_potential
    use neort_gc_orbit_integrator, only: GC_ORBIT_FIELD_ERROR, &
        GC_ORBIT_INTEGRATOR_ERROR, GC_ORBIT_NO_RETURN, &
        GC_ORBIT_PASSING, GC_ORBIT_PERTURBATION_ERROR, GC_ORBIT_RADIAL_DOMAIN, &
        GC_ORBIT_START_ROOT_ERROR, GC_ORBIT_STATE_ERROR, GC_ORBIT_SUCCESS, &
        GC_ORBIT_TRAPPED, GC_ORBIT_WALL_LOSS, &
        gc_orbit_options_t, gc_orbit_average_t, gc_orbit_perturbation_i, &
        compute_return_map, compute_thin_precession, compute_gc_orbit_average, &
        compute_gc_full_orbit_average, compute_zero_width_passing_cycle
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
    integer, parameter, public :: GC_FREQUENCY_WALL_ERROR = 5
    integer, parameter, public :: GC_FREQUENCY_BACKEND_LEGACY = 1
    integer, parameter, public :: GC_FREQUENCY_BACKEND_CYLINDRICAL = 2

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
        integer :: frequency_model = 0
#ifdef NEO_RT_USE_STANDALONE
        type(gc_cylindrical_backend_context_t) :: cylindrical_backend
#endif
        integer :: htheta_sign = 0
        logical :: initialized = .false.
    end type gc_frequency_context_t

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
        real(dp) :: omega_prec = 0.0_dp
        real(dp) :: energy_error = 0.0_dp
        real(dp) :: magnetic_moment_error = 0.0_dp
        real(dp) :: canonical_momentum_error = 0.0_dp
        integer :: backend_id = GC_FREQUENCY_BACKEND_LEGACY
        integer :: orbit_status = GC_ORBIT_SUCCESS
    end type gc_full_orbit_frequency_result_t

    type, public :: gc_frequency_runtime_metadata_t
        character(len=32) :: backend = 'uninitialized'
        character(len=16) :: coordinates = ''
        logical :: wall_certified = .false.
        character(len=64) :: wall_hash = ''
        character(len=16) :: wall_units = ''
        character(len=16) :: wall_backend_units = ''
        logical :: canonical_measure_certified = .false.
        logical :: component_identity_certified = .false.
        logical :: nonlocal_transport_required = .false.
        integer :: cylindrical_entry_count = 0
        integer :: legacy_entry_count = 0
    end type gc_frequency_runtime_metadata_t

    integer, save :: runtime_cylindrical_entry_count = 0
    integer, save :: runtime_legacy_entry_count = 0
    logical, save :: runtime_cylindrical_initialized = .false.
    logical, save :: runtime_wall_certified = .false.
    character(len=64), save :: runtime_wall_hash = ''
    character(len=16), save :: runtime_wall_units = ''
    character(len=16), save :: runtime_wall_backend_units = ''

    public :: initialize_gc_frequency_context, evaluate_gc_frequency
    public :: evaluate_gc_full_orbit_frequency
    public :: evaluate_gc_phase_average
    public :: evaluate_gc_full_orbit_phase_average
    public :: reset_gc_frequency_runtime_metadata
    public :: get_gc_frequency_runtime_metadata
#ifdef NEO_RT_USE_STANDALONE
    public :: note_gc_cylindrical_initialization
#endif

contains

    subroutine reset_gc_frequency_runtime_metadata()
        runtime_cylindrical_entry_count = 0
        runtime_legacy_entry_count = 0
        runtime_cylindrical_initialized = .false.
        runtime_wall_certified = .false.
        runtime_wall_hash = ''
        runtime_wall_units = ''
        runtime_wall_backend_units = ''
    end subroutine reset_gc_frequency_runtime_metadata

    subroutine get_gc_frequency_runtime_metadata(metadata)
        type(gc_frequency_runtime_metadata_t), intent(out) :: metadata

        metadata = gc_frequency_runtime_metadata_t()
        if (runtime_cylindrical_initialized .or. runtime_cylindrical_entry_count > 0) then
            metadata%backend = 'cylindrical_full_fow'
            metadata%coordinates = 'R,Z,phi'
            metadata%wall_certified = runtime_wall_certified
            metadata%wall_hash = runtime_wall_hash
            metadata%wall_units = runtime_wall_units
            metadata%wall_backend_units = runtime_wall_backend_units
            metadata%nonlocal_transport_required = .true.
        else if (runtime_legacy_entry_count > 0) then
            metadata%backend = 'legacy'
        end if
        metadata%cylindrical_entry_count = runtime_cylindrical_entry_count
        metadata%legacy_entry_count = runtime_legacy_entry_count
    end subroutine get_gc_frequency_runtime_metadata

#ifdef NEO_RT_USE_STANDALONE
    subroutine note_gc_cylindrical_initialization(context)
        type(gc_cylindrical_backend_context_t), intent(in) :: context

        runtime_cylindrical_initialized = context%initialized
        runtime_wall_certified = context%initialized .and. context%wall_initialized
        runtime_wall_hash = context%wall_hash
        runtime_wall_units = context%wall_units
        runtime_wall_backend_units = context%wall_backend_units
    end subroutine note_gc_cylindrical_initialization
#endif

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
#ifdef NEO_RT_USE_STANDALONE
        type(gc_cylindrical_backend_result_t) :: cylindrical_result
        integer :: cylindrical_status
#endif

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

#ifdef NEO_RT_USE_STANDALONE
        if (context%frequency_model == 2) then
            if (.not. context%cylindrical_backend%initialized) then
                result%backend_id = GC_FREQUENCY_BACKEND_CYLINDRICAL
                result%orbit_status = GC_ORBIT_FIELD_ERROR
                status = GC_FREQUENCY_FIELD_ERROR
                return
            end if
            runtime_cylindrical_entry_count = runtime_cylindrical_entry_count + 1
            call evaluate_gc_cylindrical_backend(context%cylindrical_backend, eta, &
                parallel_direction, orbit_class, period_estimate, cylindrical_result, &
                cylindrical_status, velocity)
            result%backend_id = GC_FREQUENCY_BACKEND_CYLINDRICAL
            result%orbit_status = map_cylindrical_orbit_status( &
                cylindrical_result%orbit_status)
            if (cylindrical_status /= GC_CYL_BACKEND_SUCCESS) then
                status = GC_FREQUENCY_ORBIT_ERROR
                return
            end if
            result%period = cylindrical_result%period
            result%delta_phi = cylindrical_result%delta_phi
            result%omega_b = cylindrical_result%omega_b
            result%omega_phi = cylindrical_result%omega_phi
            result%omega_prec = cylindrical_result%omega_prec
            result%energy_error = cylindrical_result%energy_error
            result%magnetic_moment_error = cylindrical_result%magnetic_moment_error
            result%canonical_momentum_error = &
                cylindrical_result%canonical_momentum_error
            status = GC_FREQUENCY_SUCCESS
            return
        end if
#else
        if (context%frequency_model == 2) then
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if
#endif

        xi_squared = 1.0_dp - eta*context%reference_sample%bmod
        if (xi_squared < 0.0_dp) return
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
            context%rho0, 1.0_dp, speed_ratio, &
            real(parallel_sign, dp)*sqrt(xi_squared), invariants, &
            invariant_status)
        if (invariant_status /= GC_MODEL_SUCCESS) return

        runtime_legacy_entry_count = runtime_legacy_entry_count + 1
        call compute_return_map(context%field, context%electric_potential, &
            invariants, context%reference_position, parallel_sign, &
            context%rho0, 1.0_dp, context%reference_velocity, orbit_class, &
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
        ! The precession part is defined against the same full-cycle field-line
        ! winding as the return map: omega_prec = omega_phi - q*omega_b.
        result%omega_prec = result%omega_phi &
            - context%q_fieldline*result%omega_b
        status = GC_FREQUENCY_SUCCESS
    end subroutine evaluate_gc_full_orbit_frequency

    subroutine initialize_gc_frequency_context(surface, reference_theta, &
            field_scale, omega_e, &
            particle_mass, particle_charge, reference_velocity, context, status, &
            selected_frequency_model, wall_file, wall_units)
        real(dp), intent(in) :: surface, reference_theta, field_scale, omega_e
        real(dp), intent(in) :: particle_mass, particle_charge, reference_velocity
        type(gc_frequency_context_t), intent(out) :: context
        integer, intent(out) :: status
        integer, intent(in), optional :: selected_frequency_model
        character(len=*), intent(in), optional :: wall_file, wall_units

        type(gc_invariants_t) :: q_invariants
        type(orbit_return_t) :: q_return
        integer :: field_status, q_invariant_status
#ifdef NEO_RT_USE_STANDALONE
        integer :: cylindrical_status
#endif

        context%field%field_scale = 1.0_dp
        context%electric_potential = gc_linear_flux_potential_t()
        context%orbit_options = gc_orbit_options_t()
        context%reference_sample = gc_field_sample_t()
        context%reference_position = 0.0_dp
        context%reference_velocity = 0.0_dp
        context%rho0 = 0.0_dp
        context%q_fieldline = 0.0_dp
        context%frequency_model = 0
        if (present(selected_frequency_model)) then
            context%frequency_model = selected_frequency_model
        end if
#ifdef NEO_RT_USE_STANDALONE
        context%cylindrical_backend = gc_cylindrical_backend_context_t()
#endif
        context%htheta_sign = 0
        context%initialized = .false.
        status = GC_FREQUENCY_INVALID_INPUT
        if (surface <= 0.0_dp .or. surface >= 1.0_dp) return
        if (field_scale <= 0.0_dp .or. particle_mass <= 0.0_dp &
            .or. particle_charge == 0.0_dp &
            .or. reference_velocity <= 0.0_dp) return

#ifdef NEO_RT_USE_STANDALONE
        if (context%frequency_model == 2 .and. inp_swi /= 11) then
            ! Model 2 is a direct cylindrical backend contract, not an
            ! inference from the input switch.  Refuse a chart mismatch
            ! instead of allowing a legacy fallback.
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if
#else
        if (context%frequency_model == 2) then
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if
#endif

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

        ! q(s*) is the full-cycle field-line winding, Delta_phi/(2*pi),
        ! evaluated by the independent lambda=0 expression used by the
        ! passing thin-limit base.  This keeps finite-width precession and
        ! the literature definitions omega_b=2*pi/tau and
        ! omega_phi=Delta_phi/tau on one Poincare section.
        call invariants_from_state(context%reference_sample, 0.0_dp, 0.0_dp, &
            0.0_dp, 1.0_dp, 1.0_dp, q_invariants, q_invariant_status)
        if (q_invariant_status /= GC_MODEL_SUCCESS) then
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if
        call compute_zero_width_passing_cycle(context%field, context%zero_potential, &
            q_invariants, context%reference_position, context%htheta_sign, &
            context%reference_velocity, 1, q_return)
        if (q_return%status /= GC_ORBIT_SUCCESS) then
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if
        context%q_fieldline = q_return%delta_phi/(2.0_dp*pi)
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
#ifdef NEO_RT_USE_STANDALONE
        if (context%frequency_model == 2) then
            if (.not. present(wall_file) .or. .not. present(wall_units)) then
                status = GC_FREQUENCY_WALL_ERROR
                return
            end if
            call initialize_gc_cylindrical_backend(surface, reference_theta, &
                field_scale, omega_e, particle_mass, particle_charge, &
                reference_velocity, c, context%htheta_sign, &
                wall_file, wall_units, context%cylindrical_backend, &
                cylindrical_status)
            if (cylindrical_status /= GC_CYL_BACKEND_SUCCESS) then
                if (cylindrical_status == GC_CYL_BACKEND_WALL_ERROR) then
                    status = GC_FREQUENCY_WALL_ERROR
                else
                    status = GC_FREQUENCY_FIELD_ERROR
                end if
                return
            end if
            call note_gc_cylindrical_initialization(context%cylindrical_backend)
        end if
#endif
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
        !! of evaluate_gc_full_orbit_frequency. The bounce harmonic is the
        !! temporal canonical phase, with no Boozer-angle reduction.
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

#ifdef NEO_RT_USE_STANDALONE
        if (context%frequency_model == 2) then
            if (.not. context%cylindrical_backend%initialized) then
                status = GC_FREQUENCY_FIELD_ERROR
                return
            end if
            runtime_cylindrical_entry_count = runtime_cylindrical_entry_count + 1
            call evaluate_gc_cylindrical_phase_average( &
                context%cylindrical_backend, eta, parallel_direction, orbit_class, &
                period_estimate, omega_b, omega_phi, mth, mph, perturbation, &
                result, status, velocity)
            return
        end if
#else
        if (context%frequency_model == 2) then
            status = GC_FREQUENCY_FIELD_ERROR
            return
        end if
#endif

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
            context%rho0, 1.0_dp, speed_ratio, &
            real(parallel_sign, dp)*sqrt(xi_squared), invariants, &
            invariant_status)
        if (invariant_status /= GC_MODEL_SUCCESS) then
            status = GC_FREQUENCY_INVALID_INPUT
            return
        end if

        call compute_gc_full_orbit_average(context%field, &
            context%electric_potential, invariants, context%reference_position, &
            parallel_sign, context%rho0, context%reference_velocity, eta, &
            orbit_class, winding, period_estimate, omega_b, omega_phi, mth, mph, perturbation, &
            context%orbit_options, result)
        status = result%status
    end subroutine evaluate_gc_full_orbit_phase_average

#ifdef NEO_RT_USE_STANDALONE
    integer function map_cylindrical_orbit_status(cylindrical_status)
        integer, intent(in) :: cylindrical_status

        select case (cylindrical_status)
        case (0)
            map_cylindrical_orbit_status = GC_ORBIT_SUCCESS
        case (GC_CYL_EQUILIBRIUM_DOMAIN)
            map_cylindrical_orbit_status = GC_ORBIT_RADIAL_DOMAIN
        case (GC_CYL_FIELD_ERROR, GC_CYL_POTENTIAL_ERROR)
            map_cylindrical_orbit_status = GC_ORBIT_FIELD_ERROR
        case (GC_CYL_STATE_ERROR, GC_CYL_START_ERROR)
            map_cylindrical_orbit_status = GC_ORBIT_STATE_ERROR
        case (GC_CYL_NO_RETURN)
            map_cylindrical_orbit_status = GC_ORBIT_NO_RETURN
        case (GC_CYL_PERTURBATION_ERROR)
            map_cylindrical_orbit_status = GC_ORBIT_PERTURBATION_ERROR
        case (GC_CYL_INVARIANT_ERROR)
            map_cylindrical_orbit_status = GC_ORBIT_STATE_ERROR
        case (GC_CYL_WALL_LOSS)
            map_cylindrical_orbit_status = GC_ORBIT_WALL_LOSS
        case default
            map_cylindrical_orbit_status = GC_ORBIT_INTEGRATOR_ERROR
        end select
    end function map_cylindrical_orbit_status
#endif

end module neort_gc_frequency_provider
