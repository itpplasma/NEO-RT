module neort_gc_cylindrical_backend
    !! Production entry for the direct-GEQDSK cylindrical return map.
    !!
    !! The transport-facing context still owns the legacy chart provider for
    !! thin model 0.  This module is the finite-width model-2 seam: it maps a
    !! launch surface only once to physical (R,Z,phi), constructs fixed
    !! physical (H,mu,P_phi), and delegates all evolution to the cylindrical
    !! Hamiltonian backend.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: GC_CYL_EQUILIBRIUM_DOMAIN, &
        GC_CYL_FIELD_ERROR, GC_CYL_INTEGRATOR_ERROR, &
        GC_CYL_NO_RETURN, GC_CYL_PERTURBATION_ERROR, GC_CYL_POTENTIAL_ERROR, &
        GC_CYL_START_ERROR, GC_CYL_STATE_ERROR, GC_CYL_SUCCESS, &
        GC_CYL_WALL_ERROR, GC_CYL_WALL_LOSS, &
        gc_cylindrical_field_sample_t, gc_cylindrical_invariants_t, &
        gc_cylindrical_linear_flux_potential_t, gc_cylindrical_section_t, &
        gc_cylindrical_polygon_wall_t, gc_cylindrical_state_t, &
        invariants_from_cylindrical_state
    use neort_wall_io, only: WALL_IO_OK, load_wall_polygon, wall_polygon_t
    use neort_gc_cylindrical_orbit, only: compute_gc_cylindrical_phase_average, &
        compute_gc_cylindrical_return, gc_cylindrical_orbit_options_t, &
        gc_cylindrical_phase_result_t, gc_cylindrical_return_t
    use neort_gc_orbit_integrator, only: GC_ORBIT_FIELD_ERROR, &
        GC_ORBIT_INTEGRATOR_ERROR, &
        GC_ORBIT_NO_RETURN, GC_ORBIT_PASSING, GC_ORBIT_PERTURBATION_ERROR, &
        GC_ORBIT_RADIAL_DOMAIN, GC_ORBIT_STATE_ERROR, GC_ORBIT_SUCCESS, &
        GC_ORBIT_TRAPPED, GC_ORBIT_WALL_LOSS, gc_orbit_average_t, &
        gc_orbit_perturbation_i
    use geoflux_coordinates, only: cyl_to_geoflux
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        configure_eqdsk_cylindrical_field, eqdsk_cylindrical_field_t, &
        map_eqdsk_flux_position

    implicit none
    private

    real(dp), parameter, public :: GC_CYL_BACKEND_C = 2.99792458e10_dp
    integer, parameter, public :: GC_CYL_BACKEND_SUCCESS = 0
    integer, parameter, public :: GC_CYL_BACKEND_INVALID_INPUT = 1
    integer, parameter, public :: GC_CYL_BACKEND_FIELD_ERROR = 2
    integer, parameter, public :: GC_CYL_BACKEND_ORBIT_ERROR = 3
    integer, parameter, public :: GC_CYL_BACKEND_WALL_ERROR = 4
    integer, parameter, public :: GC_CYL_BACKEND_MEASURE_UNAVAILABLE = 5

    type, public :: gc_cylindrical_backend_context_t
        type(eqdsk_cylindrical_field_t) :: field
        type(gc_cylindrical_linear_flux_potential_t) :: potential
        type(gc_cylindrical_orbit_options_t) :: orbit_options
        type(gc_cylindrical_section_t) :: section
        real(dp) :: reference_position(3) = 0.0_dp
        real(dp) :: reference_velocity = 0.0_dp
        real(dp) :: particle_mass = 0.0_dp
        real(dp) :: particle_charge = 0.0_dp
        real(dp) :: c_light = GC_CYL_BACKEND_C
        real(dp) :: surface = 0.0_dp
        real(dp) :: reference_theta = 0.0_dp
        real(dp) :: q_fieldline = 0.0_dp
        integer :: htheta_sign = 1
        type(gc_cylindrical_polygon_wall_t) :: wall
        character(len=1024) :: wall_path = ''
        character(len=16) :: wall_units = ''
        character(len=16) :: wall_backend_units = 'cm'
        character(len=64) :: wall_hash = ''
        integer :: wall_io_status = WALL_IO_OK
        character(len=128) :: wall_error_message = ''
        logical :: wall_initialized = .false.
        logical :: initialized = .false.
    end type gc_cylindrical_backend_context_t

    type, public :: gc_cylindrical_backend_result_t
        integer :: status = GC_CYL_BACKEND_INVALID_INPUT
        integer :: orbit_status = 0
        real(dp) :: period = 0.0_dp
        real(dp) :: delta_phi = 0.0_dp
        real(dp) :: omega_b = 0.0_dp
        real(dp) :: omega_phi = 0.0_dp
        real(dp) :: omega_prec = 0.0_dp
        real(dp) :: energy_error = 0.0_dp
        real(dp) :: magnetic_moment_error = 0.0_dp
        real(dp) :: canonical_momentum_error = 0.0_dp
        logical :: used_cylindrical_backend = .false.
    end type gc_cylindrical_backend_result_t

    type, public :: gc_cylindrical_nonlocal_measure_result_t
        integer :: status = GC_CYL_BACKEND_INVALID_INPUT
        integer :: component_id = 0
        integer :: sigma = 0
        logical :: available = .false.
        logical :: canonical_measure_certified = .false.
        logical :: component_identity_certified = .false.
    end type gc_cylindrical_nonlocal_measure_result_t

    public :: initialize_gc_cylindrical_backend
    public :: evaluate_gc_cylindrical_backend
    public :: evaluate_gc_cylindrical_phase_average
    public :: evaluate_gc_cylindrical_nonlocal_measure

contains

    subroutine initialize_gc_cylindrical_backend(surface, reference_theta, &
            field_scale, omega_e, particle_mass, particle_charge, &
            reference_velocity, c_light, htheta_sign, wall_path, wall_units, &
            context, status)
        real(dp), intent(in) :: surface, reference_theta, field_scale, omega_e
        real(dp), intent(in) :: particle_mass, particle_charge, reference_velocity
        real(dp), intent(in) :: c_light
        integer, intent(in) :: htheta_sign
        character(len=*), intent(in) :: wall_path, wall_units
        type(gc_cylindrical_backend_context_t), intent(out) :: context
        integer, intent(out) :: status

        type(gc_cylindrical_field_sample_t) :: sample
        type(wall_polygon_t) :: input_wall
        character(len=128) :: wall_message
        integer :: field_status, map_status, wall_status

        context = gc_cylindrical_backend_context_t()
        status = GC_CYL_BACKEND_INVALID_INPUT
        if (surface <= 0.0_dp .or. surface >= 1.0_dp) return
        if (field_scale <= 0.0_dp .or. particle_mass <= 0.0_dp .or. &
                particle_charge == 0.0_dp .or. reference_velocity <= 0.0_dp .or. &
                c_light <= 0.0_dp .or. abs(htheta_sign) /= 1) return
        if (len_trim(wall_path) == 0 .or. len_trim(wall_units) == 0) then
            status = GC_CYL_BACKEND_WALL_ERROR
            return
        end if

        call load_wall_polygon(trim(wall_path), input_wall, wall_status, &
            wall_message, trim(wall_units))
        if (wall_status /= WALL_IO_OK) then
            status = GC_CYL_BACKEND_WALL_ERROR
            context%wall_io_status = wall_status
            context%wall_error_message = wall_message
            context%wall_path = wall_path
            context%wall_units = wall_units
            return
        end if

        call configure_eqdsk_cylindrical_field(field_scale, context%field, &
            status)
        if (status /= GC_CYL_SUCCESS) then
            status = GC_CYL_BACKEND_FIELD_ERROR
            return
        end if
        call map_eqdsk_flux_position(surface, reference_theta, 0.0_dp, &
            context%reference_position, map_status)
        if (map_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_BACKEND_FIELD_ERROR
            return
        end if
        call context%field%evaluate(context%reference_position, sample, field_status)
        if (field_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_BACKEND_FIELD_ERROR
            return
        end if

        ! The electrostatic potential is Phi=(omega_e/c)*psi in the
        ! cylindrical backend.  The factor two sometimes used in the thin
        ! normalized drift formula is already accounted for by that model's
        ! normalization; carrying it here doubles the physical E x B drift.
        context%potential%coefficient = omega_e/c_light
        context%potential%psi_reference = sample%psi
        context%orbit_options = gc_cylindrical_orbit_options_t()
        context%reference_velocity = reference_velocity
        context%particle_mass = particle_mass
        context%particle_charge = particle_charge
        context%c_light = c_light
        context%htheta_sign = htheta_sign
        call context%wall%set_vertices(input_wall%vertices, wall_status)
        if (wall_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_BACKEND_WALL_ERROR
            context%wall_io_status = wall_status
            context%wall_error_message = 'validated wall could not enter polygon object'
            return
        end if
        context%wall_path = wall_path
        context%wall_units = input_wall%input_units
        context%wall_backend_units = 'cm'
        context%wall_hash = input_wall%hash
        context%wall_io_status = WALL_IO_OK
        context%wall_error_message = ''
        context%wall_initialized = .true.
        context%surface = surface
        context%reference_theta = reference_theta
        call make_physical_section(context, status)
        if (status /= GC_CYL_BACKEND_SUCCESS) return
        call compute_direct_fieldline_q(context, status)
        if (status /= GC_CYL_BACKEND_SUCCESS) return
        context%section%fieldline_q = context%q_fieldline
        context%initialized = .true.
        status = GC_CYL_BACKEND_SUCCESS
    end subroutine initialize_gc_cylindrical_backend

    subroutine evaluate_gc_cylindrical_backend(context, eta, parallel_direction, &
            orbit_class, period_estimate, result, status, velocity)
        type(gc_cylindrical_backend_context_t), intent(in) :: context
        real(dp), intent(in) :: eta, period_estimate
        integer, intent(in) :: parallel_direction, orbit_class
        type(gc_cylindrical_backend_result_t), intent(out) :: result
        integer, intent(out) :: status
        real(dp), intent(in), optional :: velocity

        type(gc_cylindrical_field_sample_t) :: field
        type(gc_cylindrical_state_t) :: state
        type(gc_cylindrical_invariants_t) :: invariants
        type(gc_cylindrical_section_t) :: section
        type(gc_cylindrical_return_t) :: orbit
        real(dp) :: potential, gradient(3), speed, xi_squared
        integer :: field_status, potential_status, invariant_status
        integer :: parallel_sign, winding

        result = gc_cylindrical_backend_result_t()
        status = GC_CYL_BACKEND_INVALID_INPUT
        if (.not. context%initialized .or. eta <= 0.0_dp .or. &
                period_estimate <= 0.0_dp) return
        if (abs(parallel_direction) /= 1) return
        if (orbit_class /= GC_ORBIT_TRAPPED .and. &
                orbit_class /= GC_ORBIT_PASSING) return

        speed = context%reference_velocity
        if (present(velocity)) speed = velocity
        if (speed <= 0.0_dp) return
        call context%field%evaluate(context%reference_position, field, field_status)
        if (field_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_BACKEND_FIELD_ERROR
            return
        end if
        call context%potential%evaluate(context%reference_position, field, potential, &
            gradient, potential_status)
        if (potential_status /= GC_CYL_SUCCESS) then
            status = GC_CYL_BACKEND_FIELD_ERROR
            return
        end if

        xi_squared = 1.0_dp - eta*field%bmod
        if (xi_squared < 0.0_dp) return
        state = gc_cylindrical_state_t()
        state%R = context%reference_position(1)
        state%Z = context%reference_position(2)
        state%phi = context%reference_position(3)
        parallel_sign = parallel_direction*context%htheta_sign
        state%p_parallel = context%particle_mass*speed*real(parallel_sign, dp) &
            *sqrt(max(0.0_dp, xi_squared))
        state%mu = 0.5_dp*context%particle_mass*speed**2*eta
        call invariants_from_cylindrical_state(state, field, potential, &
            context%particle_mass, context%particle_charge, context%c_light, &
            invariants, invariant_status)
        if (invariant_status /= GC_CYL_SUCCESS) return

        section = context%section
        if (orbit_class == GC_ORBIT_PASSING) then
            winding = parallel_direction
            section%winding = winding
            section%frequency_winding = winding
            section%fieldline_q = context%q_fieldline
        else
            section%winding = 0
            section%frequency_winding = 1
            section%fieldline_q = 0.0_dp
        end if
        call compute_gc_cylindrical_return(context%field, context%potential, invariants, &
            context%reference_position, parallel_sign, context%particle_mass, &
            context%particle_charge, context%c_light, section, &
            orbit_options_with_period(context, period_estimate), orbit, context%wall)
        result%used_cylindrical_backend = .true.
        result%orbit_status = orbit%status
        if (orbit%status /= GC_CYL_SUCCESS) then
            result%status = GC_CYL_BACKEND_ORBIT_ERROR
            return
        end if
        result%period = orbit%period
        result%delta_phi = orbit%delta_phi
        result%omega_b = orbit%omega_b
        result%omega_phi = orbit%omega_phi
        result%omega_prec = orbit%omega_prec
        result%energy_error = orbit%energy_error
        result%magnetic_moment_error = orbit%magnetic_moment_error
        result%canonical_momentum_error = orbit%canonical_momentum_error
        result%status = GC_CYL_BACKEND_SUCCESS
        status = GC_CYL_BACKEND_SUCCESS
    end subroutine evaluate_gc_cylindrical_backend

    subroutine evaluate_gc_cylindrical_phase_average(context, eta, &
            parallel_direction, orbit_class, period_estimate, omega_b, omega_phi, &
            bounce_harmonic, toroidal_harmonic, perturbation, result, status, velocity)
        type(gc_cylindrical_backend_context_t), intent(in) :: context
        real(dp), intent(in) :: eta, period_estimate, omega_b, omega_phi
        integer, intent(in) :: parallel_direction, orbit_class
        integer, intent(in) :: bounce_harmonic, toroidal_harmonic
        procedure(gc_orbit_perturbation_i) :: perturbation
        type(gc_orbit_average_t), intent(out) :: result
        integer, intent(out) :: status
        real(dp), intent(in), optional :: velocity

        type(gc_cylindrical_field_sample_t) :: field
        type(gc_cylindrical_invariants_t) :: invariants
        type(gc_cylindrical_phase_result_t) :: phase
        type(gc_cylindrical_section_t) :: section
        real(dp) :: potential, gradient(3), speed, xi_squared
        integer :: field_status, potential_status, invariant_status
        integer :: parallel_sign, winding

        result = gc_orbit_average_t()
        status = GC_ORBIT_STATE_ERROR
        if (.not. context%initialized .or. eta <= 0.0_dp .or. &
                period_estimate <= 0.0_dp) return
        if (abs(parallel_direction) /= 1) return
        if (orbit_class /= GC_ORBIT_TRAPPED .and. &
                orbit_class /= GC_ORBIT_PASSING) return
        speed = context%reference_velocity
        if (present(velocity)) speed = velocity
        if (speed <= 0.0_dp) return

        call context%field%evaluate(context%reference_position, field, field_status)
        if (field_status /= GC_CYL_SUCCESS) then
            status = GC_ORBIT_FIELD_ERROR
            return
        end if
        call context%potential%evaluate(context%reference_position, field, potential, &
            gradient, potential_status)
        if (potential_status /= GC_CYL_SUCCESS) then
            status = GC_ORBIT_FIELD_ERROR
            return
        end if
        xi_squared = 1.0_dp - eta*field%bmod
        if (xi_squared < 0.0_dp) return

        parallel_sign = parallel_direction*context%htheta_sign
        call make_launch_invariants(context, field, potential, speed, &
            xi_squared, parallel_sign, eta, invariants, invariant_status)
        if (invariant_status /= GC_CYL_SUCCESS) then
            status = GC_ORBIT_STATE_ERROR
            return
        end if
        section = context%section
        if (orbit_class == GC_ORBIT_PASSING) then
            winding = parallel_direction
            section%winding = winding
            section%frequency_winding = winding
            section%fieldline_q = context%q_fieldline
        else
            section%winding = 0
            section%frequency_winding = 1
            section%fieldline_q = 0.0_dp
        end if
        call compute_gc_cylindrical_phase_average(context%field, context%potential, &
            invariants, context%reference_position, parallel_sign, &
            context%particle_mass, context%particle_charge, context%c_light, &
            omega_b, omega_phi, bounce_harmonic, toroidal_harmonic, &
            cylindrical_perturbation, orbit_options_with_period(context, period_estimate), &
            section, phase, context%wall)
        result%status = map_cylindrical_orbit_status(phase%status)
        status = result%status
        if (phase%status /= GC_CYL_SUCCESS) return
        result%period = phase%period
        result%perturbation_average = phase%perturbation_average
        result%inverse_b_average = phase%inverse_b_average
        result%b_average = phase%b_average

    contains

        subroutine cylindrical_perturbation(position, cylindrical_field, amplitude, &
                perturbation_status)
            real(dp), intent(in) :: position(3)
            type(gc_cylindrical_field_sample_t), intent(in) :: cylindrical_field
            complex(dp), intent(out) :: amplitude
            integer, intent(out) :: perturbation_status

            real(dp) :: geoflux_position(3), chart_position(3)

            ! The orbit remains entirely in cylindrical space.  This is only
            ! the adapter boundary for the existing perturbation callback,
            ! which is still keyed by the direct GEQDSK chart.
            call cyl_to_geoflux([position(1), position(3), position(2)], &
                geoflux_position)
            chart_position = [geoflux_position(1), geoflux_position(3), &
                geoflux_position(2)]
            call perturbation(chart_position, cylindrical_field%bmod, amplitude, &
                perturbation_status)
        end subroutine cylindrical_perturbation
    end subroutine evaluate_gc_cylindrical_phase_average

    subroutine evaluate_gc_cylindrical_nonlocal_measure(context, eta, velocity, &
            parallel_direction, orbit_class, result, status)
        !! Architectural placeholder for the nonlocal invariant measure.
        !!
        !! A local eta callback cannot be certified here: Eq. 14--17 of
        !! Buchholz et al. require the fixed-(H,J_perp) Poincare classes and
        !! the change of variable from the physical section coordinate to
        !! P_phi.  An energy-allowed interval with a varying psi_star is not
        !! that transform.  Keep the API explicit and fail closed until the
        !! Eq. 17 nonlocal transport integral is implemented.
        type(gc_cylindrical_backend_context_t), intent(in) :: context
        real(dp), intent(in) :: eta, velocity
        integer, intent(in) :: parallel_direction, orbit_class
        type(gc_cylindrical_nonlocal_measure_result_t), intent(out) :: result
        integer, intent(out) :: status

        associate (unused_context => context, unused_eta => eta, &
            unused_velocity => velocity, unused_parallel_direction => parallel_direction, &
            unused_orbit_class => orbit_class)
        end associate
        result = gc_cylindrical_nonlocal_measure_result_t()
        result%status = GC_CYL_BACKEND_MEASURE_UNAVAILABLE
        status = GC_CYL_BACKEND_MEASURE_UNAVAILABLE
        return

    end subroutine evaluate_gc_cylindrical_nonlocal_measure

    subroutine make_launch_invariants(context, field, potential, speed, xi_squared, &
            parallel_sign, eta, invariants, status)
        type(gc_cylindrical_backend_context_t), intent(in) :: context
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(in) :: potential, speed, xi_squared
        integer, intent(in) :: parallel_sign
        real(dp), intent(in) :: eta
        type(gc_cylindrical_invariants_t), intent(out) :: invariants
        integer, intent(out) :: status

        type(gc_cylindrical_state_t) :: state

        state = gc_cylindrical_state_t()
        state%R = context%reference_position(1)
        state%Z = context%reference_position(2)
        state%phi = context%reference_position(3)
        state%p_parallel = context%particle_mass*speed*real(parallel_sign, dp) &
            *sqrt(max(0.0_dp, xi_squared))
        state%mu = 0.5_dp*context%particle_mass*speed**2*eta
        call invariants_from_cylindrical_state(state, field, potential, &
            context%particle_mass, context%particle_charge, context%c_light, &
            invariants, status)
    end subroutine make_launch_invariants

    function orbit_options_with_period(context, period_estimate) result(options)
        type(gc_cylindrical_backend_context_t), intent(in) :: context
        type(gc_cylindrical_orbit_options_t) :: options
        real(dp), intent(in) :: period_estimate

        options = context%orbit_options
        options%maximum_time = max(20.0_dp*period_estimate, &
            10.0_dp*options%event_time_tolerance)
        options%event_time_tolerance = min(options%event_time_tolerance, &
            1.0e-3_dp*period_estimate)
    end function orbit_options_with_period

    subroutine make_physical_section(context, status)
        type(gc_cylindrical_backend_context_t), intent(inout) :: context
        integer, intent(out) :: status

        real(dp) :: poloidal_field(2), norm
        type(gc_cylindrical_field_sample_t) :: sample
        integer :: field_status

        context%section = gc_cylindrical_section_t()
        context%section%point = context%reference_position
        status = GC_CYL_BACKEND_FIELD_ERROR
        call context%field%evaluate(context%reference_position, sample, field_status)
        if (field_status /= GC_CYL_SUCCESS) return
        ! The section must be transverse to the poloidal trajectory.  The
        ! flux-surface tangent is parallel to B_p and would define a section
        ! that the orbit never crosses; use B_p itself as the section normal.
        poloidal_field = [sample%b(1), sample%b(3)]
        norm = sqrt(dot_product(poloidal_field, poloidal_field))
        if (norm <= tiny(norm)) return
        context%section%normal = poloidal_field/norm
        context%section%direction = 0
        status = GC_CYL_BACKEND_SUCCESS
    end subroutine make_physical_section

    subroutine compute_direct_fieldline_q(context, status)
        type(gc_cylindrical_backend_context_t), intent(inout) :: context
        integer, intent(out) :: status

        integer, parameter :: samples = 1024
        real(dp), parameter :: two_pi = 6.28318530717958647692528676656_dp
        real(dp) :: step, theta, previous_position(3), current_position(3)
        real(dp) :: next_position(3), tangent(2), bp(2), bp_squared
        real(dp) :: integrand, q_sum
        type(gc_cylindrical_field_sample_t) :: field
        integer :: k, map_status, field_status

        status = GC_CYL_BACKEND_FIELD_ERROR
        step = two_pi/real(samples, dp)
        q_sum = 0.0_dp
        do k = 1, samples
            theta = -0.5_dp*two_pi + (real(k, dp) - 0.5_dp)*step
            call map_eqdsk_flux_position(context%surface, theta - step*0.5_dp, &
                0.0_dp, previous_position, map_status)
            if (map_status /= GC_CYL_SUCCESS) return
            call map_eqdsk_flux_position(context%surface, theta, 0.0_dp, &
                current_position, map_status)
            if (map_status /= GC_CYL_SUCCESS) return
            call map_eqdsk_flux_position(context%surface, theta + step*0.5_dp, &
                0.0_dp, next_position, map_status)
            if (map_status /= GC_CYL_SUCCESS) return
            call context%field%evaluate(current_position, field, field_status)
            if (field_status /= GC_CYL_SUCCESS) return
            tangent = [(next_position(1) - previous_position(1))/step, &
                (next_position(2) - previous_position(2))/step]
            bp = [field%b(1), field%b(3)]
            bp_squared = dot_product(bp, bp)
            if (bp_squared <= tiny(bp_squared)) return
            integrand = dot_product(tangent, bp)*field%b(2) &
                /(current_position(1)*bp_squared)
            q_sum = q_sum + integrand*step
        end do
        context%q_fieldline = q_sum/two_pi
        status = GC_CYL_BACKEND_SUCCESS
    end subroutine compute_direct_fieldline_q

    integer function map_cylindrical_orbit_status(cylindrical_status)
        integer, intent(in) :: cylindrical_status

        select case (cylindrical_status)
        case (GC_CYL_SUCCESS)
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
        case (GC_CYL_WALL_LOSS)
            map_cylindrical_orbit_status = GC_ORBIT_WALL_LOSS
        case (GC_CYL_INTEGRATOR_ERROR, GC_CYL_WALL_ERROR)
            map_cylindrical_orbit_status = GC_ORBIT_INTEGRATOR_ERROR
        case default
            map_cylindrical_orbit_status = GC_ORBIT_INTEGRATOR_ERROR
        end select
    end function map_cylindrical_orbit_status

end module neort_gc_cylindrical_backend
