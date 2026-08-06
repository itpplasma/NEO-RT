program test_gc_orbit_limit
    !! Independent circular-tokamak oracle for the real-space return map.
    !! The field kernel is generated from a Cartesian chart by fortsym.  The
    !! expected trapped precession is the large-aspect-ratio elliptic-integral
    !! result, so agreement improves as epsilon=r/R0 tends to zero.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use manufactured_gc_field, only: manufactured_torus_field_t
    use neort_gc_dynamics, only: gc_field_sample_t
    use neort_gc_models, only: GC_MODEL_SUCCESS, gc_invariants_t, &
        gc_zero_potential_t, gc_linear_flux_potential_t, &
        invariants_from_state, make_linear_flux_potential
    use neort_gc_orbit_integrator, only: GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, &
        gc_orbit_options_t, compute_thin_precession
    use neort_thin_orbit_limit, only: THIN_LIMIT_SUCCESS, orbit_return_t, &
        thin_limit_result_t
    use util, only: pi, qe, mu, c
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: major_radius = 1000.0_dp
    real(dp), parameter :: field_axis = 2.0e6_dp
    real(dp), parameter :: q_axis = 2.0_dp
    real(dp), parameter :: target_shear = 0.5_dp
    real(dp), parameter :: reference_velocity = 1.0e8_dp
    real(dp), parameter :: particle_mass = 2.014_dp*mu
    real(dp), parameter :: particle_charge = qe
    real(dp), parameter :: rho0 = particle_mass*c*reference_velocity/particle_charge
    real(dp), parameter :: electric_frequency = 4.0e3_dp

    real(dp) :: magnetic_numeric(2), magnetic_analytic(2), relative_error(2)
    real(dp) :: electric_trapped, electric_passing, magnetic_reversed

    call trapped_case(0.02_dp, particle_charge, rho0, magnetic_numeric(1), &
        magnetic_analytic(1), electric_trapped)
    call trapped_case(0.01_dp, particle_charge, rho0, magnetic_numeric(2), &
        magnetic_analytic(2), electric_trapped)
    relative_error = abs(magnetic_numeric - magnetic_analytic) &
        /max(abs(magnetic_analytic), tiny(1.0_dp))
    write(*, '(a,2es14.5)') 'manufactured trapped magnetic numeric:  ', magnetic_numeric
    write(*, '(a,2es14.5)') 'manufactured trapped magnetic analytic: ', magnetic_analytic
    write(*, '(a,2es14.5)') 'manufactured trapped relative error:     ', relative_error
    if (relative_error(2) >= relative_error(1) .or. relative_error(2) > 6.0e-2_dp) then
        error stop 'thin-orbit magnetic precession missed its circular limit'
    end if
    if (abs(electric_trapped - electric_frequency) &
        > 3.0e-3_dp*abs(electric_frequency)) then
        write(*,*) 'trapped electric/reference:', electric_trapped, electric_frequency
        error stop 'trapped electric precession has wrong sign or magnitude'
    end if

    call passing_case(0.02_dp, particle_charge, rho0, electric_passing)
    if (abs(electric_passing - electric_frequency) &
        > 3.0e-3_dp*abs(electric_frequency)) then
        write(*,*) 'passing electric/reference:', electric_passing, electric_frequency
        error stop 'passing electric precession has wrong sign or magnitude'
    end if

    call trapped_case(0.02_dp, -particle_charge, -rho0, magnetic_reversed, &
        magnetic_analytic(1), electric_trapped)
    if (abs(magnetic_reversed + magnetic_numeric(1)) &
        > 3.0e-3_dp*abs(magnetic_numeric(1))) then
        write(*,*) 'charge reversal magnetic:', magnetic_numeric(1), magnetic_reversed
        error stop 'magnetic precession did not reverse with charge'
    end if
    if (abs(electric_trapped - electric_frequency) &
        > 3.0e-3_dp*abs(electric_frequency)) then
        error stop 'electric precession changed under charge reversal'
    end if

    call pass_test

contains

    subroutine trapped_case(epsilon, charge, signed_rho0, omega_magnetic, &
            omega_analytic, omega_electric)
        real(dp), intent(in) :: epsilon, charge, signed_rho0
        real(dp), intent(out) :: omega_magnetic, omega_analytic, omega_electric
        type(manufactured_torus_field_t) :: field
        type(gc_zero_potential_t) :: zero_potential
        type(gc_linear_flux_potential_t) :: total_potential
        type(gc_field_sample_t) :: sample_min, sample_max
        type(gc_invariants_t) :: invariants
        type(gc_orbit_options_t) :: options
        type(thin_limit_result_t) :: magnetic, total
        type(orbit_return_t) :: base
        real(dp) :: reference(3), radius, q_safety, magnetic_shear
        real(dp) :: eta, xi, k_squared, elliptic_k, elliptic_e
        real(dp) :: b_center, c_pitch, period_estimate, cyclotron
        integer :: status

        call configure_field(epsilon, field, radius)
        reference = [radius, 0.0_dp, 0.0_dp]
        call field%evaluate(reference, sample_min, status)
        if (status /= GC_MODEL_SUCCESS) error stop 'manufactured Bmin failed'
        call field%evaluate([radius, 0.0_dp, pi], sample_max, status)
        if (status /= GC_MODEL_SUCCESS) error stop 'manufactured Bmax failed'
        call field%radial_profiles(radius, q_safety, magnetic_shear)

        k_squared = 0.4_dp
        b_center = 0.5_dp*(sample_min%bmod + sample_max%bmod) &
            *(1.0_dp - epsilon**2)
        eta = 1.0_dp/(sample_min%bmod + 2.0_dp*epsilon*b_center*k_squared)
        xi = sqrt(1.0_dp - eta*sample_min%bmod)
        call invariants_from_state(sample_min, 0.0_dp, signed_rho0, 0.0_dp, &
            1.0_dp, xi, invariants, status)
        if (status /= GC_MODEL_SUCCESS) error stop 'trapped invariants failed'

        call complete_elliptic(k_squared, elliptic_k, elliptic_e)
        c_pitch = 2.0_dp*epsilon*eta*b_center
        period_estimate = 8.0_dp*q_safety*major_radius*elliptic_k &
            /(reference_velocity*sqrt(c_pitch))
        call configure_options(radius, options)
        call compute_thin_precession(field, zero_potential, invariants, &
            reference, 1, signed_rho0, reference_velocity, q_safety, &
            GC_ORBIT_TRAPPED, 0, period_estimate, options, magnetic, base)
        call require_limit('trapped magnetic', magnetic, base)

        total_potential = make_linear_flux_potential(electric_frequency, &
            charge, particle_mass, reference_velocity, sample_min%psi)
        call compute_thin_precession(field, total_potential, invariants, &
            reference, 1, signed_rho0, reference_velocity, q_safety, &
            GC_ORBIT_TRAPPED, 0, period_estimate, options, total)
        call require_limit('trapped total', total)

        cyclotron = reference_velocity*b_center/signed_rho0
        ! This chart has psi increasing outwards and alpha=phi-q*theta.
        ! With psi_star=psi+rho*p_parallel*h_phi, the circular-limit
        ! precession is the negative of the positive-frequency expression
        ! often quoted for the opposite canonical-flux convention.
        omega_analytic = -q_safety*reference_velocity**2 &
            /(cyclotron*radius*major_radius)*(elliptic_e/elliptic_k - 0.5_dp &
            + 2.0_dp*magnetic_shear*(elliptic_e/elliptic_k &
            + k_squared - 1.0_dp))
        omega_magnetic = magnetic%omega
        omega_electric = total%omega - magnetic%omega
    end subroutine trapped_case

    subroutine passing_case(epsilon, charge, signed_rho0, omega_electric)
        real(dp), intent(in) :: epsilon, charge, signed_rho0
        real(dp), intent(out) :: omega_electric
        type(manufactured_torus_field_t) :: field
        type(gc_zero_potential_t) :: zero_potential
        type(gc_linear_flux_potential_t) :: total_potential
        type(gc_field_sample_t) :: sample_min
        type(gc_invariants_t) :: invariants
        type(gc_orbit_options_t) :: options
        type(thin_limit_result_t) :: magnetic, total
        type(orbit_return_t) :: base
        real(dp) :: reference(3), radius, q_safety, magnetic_shear
        real(dp) :: xi, eta, period_estimate
        integer :: status

        call configure_field(epsilon, field, radius)
        reference = [radius, 0.0_dp, 0.0_dp]
        call field%evaluate(reference, sample_min, status)
        if (status /= GC_MODEL_SUCCESS) error stop 'passing sample failed'
        call field%radial_profiles(radius, q_safety, magnetic_shear)
        xi = 0.8_dp
        eta = (1.0_dp - xi**2)/sample_min%bmod
        call invariants_from_state(sample_min, 0.0_dp, signed_rho0, 0.0_dp, &
            1.0_dp, xi, invariants, status)
        if (status /= GC_MODEL_SUCCESS) error stop 'passing invariants failed'
        period_estimate = passing_period(field, radius, eta)
        call configure_options(radius, options)

        call compute_thin_precession(field, zero_potential, invariants, &
            reference, 1, signed_rho0, reference_velocity, q_safety, &
            GC_ORBIT_PASSING, 1, period_estimate, options, magnetic, base)
        call require_limit('passing magnetic', magnetic, base)
        if (abs(base%delta_phi - 2.0_dp*pi*q_safety) &
            > options%baseline_relative_tolerance*2.0_dp*pi*q_safety) then
            error stop 'zero-width passing return does not recover q'
        end if

        total_potential = make_linear_flux_potential(electric_frequency, &
            charge, particle_mass, reference_velocity, sample_min%psi)
        call compute_thin_precession(field, total_potential, invariants, &
            reference, 1, signed_rho0, reference_velocity, q_safety, &
            GC_ORBIT_PASSING, 1, period_estimate, options, total)
        call require_limit('passing total', total)
        omega_electric = total%omega - magnetic%omega
    end subroutine passing_case

    subroutine configure_field(epsilon, field, radius)
        real(dp), intent(in) :: epsilon
        type(manufactured_torus_field_t), intent(out) :: field
        real(dp), intent(out) :: radius
        real(dp) :: u, geometric_shear, shape

        radius = epsilon*major_radius
        u = epsilon**2
        geometric_shear = u/(1.0_dp - u)
        shape = -(target_shear - geometric_shear) &
            /(u*(2.0_dp + target_shear - geometric_shear))
        field%major_radius = major_radius
        field%toroidal_flux = field_axis*major_radius
        field%q_axis = q_axis
        field%shear_shape = shape
    end subroutine configure_field

    subroutine configure_options(radius, options)
        real(dp), intent(in) :: radius
        type(gc_orbit_options_t), intent(out) :: options

        options = gc_orbit_options_t()
        options%radial_min = 0.2_dp*radius
        options%radial_max = 1.8_dp*radius
        options%relative_tolerance = 1.0e-11_dp
        options%absolute_tolerance = 3.0e-12_dp
        options%event_relative_tolerance = 3.0e-12_dp
        options%baseline_relative_tolerance = 2.0e-8_dp
    end subroutine configure_options

    function passing_period(field, radius, eta) result(period)
        type(manufactured_torus_field_t), intent(in) :: field
        real(dp), intent(in) :: radius, eta
        real(dp) :: period
        integer, parameter :: n = 20000
        type(gc_field_sample_t) :: sample
        real(dp) :: theta, dtheta
        integer :: k, status

        dtheta = 2.0_dp*pi/real(n, dp)
        period = 0.0_dp
        do k = 1, n
            theta = (real(k, dp) - 0.5_dp)*dtheta
            call field%evaluate([radius, 0.0_dp, theta], sample, status)
            if (status /= GC_MODEL_SUCCESS) error stop 'passing quadrature field failed'
            period = period + dtheta/(reference_velocity*sample%hcon(3) &
                *sqrt(1.0_dp - eta*sample%bmod))
        end do
    end function passing_period

    subroutine complete_elliptic(parameter, elliptic_k, elliptic_e)
        real(dp), intent(in) :: parameter
        real(dp), intent(out) :: elliptic_k, elliptic_e
        integer, parameter :: n = 200000
        real(dp) :: angle, dangle, denominator
        integer :: k

        dangle = 0.5_dp*pi/real(n, dp)
        elliptic_k = 0.0_dp
        elliptic_e = 0.0_dp
        do k = 1, n
            angle = (real(k, dp) - 0.5_dp)*dangle
            denominator = sqrt(1.0_dp - parameter*sin(angle)**2)
            elliptic_k = elliptic_k + dangle/denominator
            elliptic_e = elliptic_e + dangle*denominator
        end do
    end subroutine complete_elliptic

    subroutine require_limit(label, result, base)
        character(*), intent(in) :: label
        type(thin_limit_result_t), intent(in) :: result
        type(orbit_return_t), intent(in), optional :: base

        if (result%status /= THIN_LIMIT_SUCCESS) then
            write(*,*) trim(label)//' thin-limit status:', result%status
            if (present(base)) write(*,*) 'base return:', base
            error stop 'manufactured thin-limit evaluation failed'
        end if
        if (result%observed_order < 1.5_dp) then
            write(*,*) trim(label)//' observed order:', result%observed_order
            write(*,*) 'centered ladder:', result%centered
            error stop 'manufactured thin-limit did not converge quadratically'
        end if
    end subroutine require_limit

end program test_gc_orbit_limit
