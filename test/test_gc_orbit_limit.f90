program test_gc_orbit_limit
    !! Retained orbit-integrator oracle for the direct/full model-2 entry points.
    !!
    !! This test deliberately exercises the independent zero-width passing-cycle
    !! expression and the shared return-map integrator.  The removed
    !! real-space-thin precession and phase-average APIs are not part of this
    !! contract.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use manufactured_gc_field, only: manufactured_torus_field_t
    use neort_gc_dynamics, only: gc_field_sample_t
    use neort_gc_models, only: GC_MODEL_SUCCESS, gc_invariants_t, &
        gc_zero_potential_t, invariants_from_state
    use neort_gc_orbit_integrator, only: GC_ORBIT_PASSING, GC_ORBIT_SUCCESS, &
        gc_orbit_options_t, compute_return_map, &
        compute_zero_width_passing_cycle
    use neort_thin_orbit_limit, only: orbit_return_t
    use util, only: pi, qe, mu, c
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: major_radius = 1000.0_dp
    real(dp), parameter :: field_axis = 2.0e6_dp
    real(dp), parameter :: q_axis = 2.0_dp
    real(dp), parameter :: reference_velocity = 1.0e8_dp
    real(dp), parameter :: particle_mass = 2.014_dp*mu
    real(dp), parameter :: particle_charge = qe
    real(dp), parameter :: rho0 = particle_mass*c*reference_velocity/particle_charge

    type(manufactured_torus_field_t) :: field
    type(gc_zero_potential_t) :: zero_potential
    type(gc_field_sample_t) :: sample
    type(gc_invariants_t) :: invariants
    type(gc_orbit_options_t) :: options
    type(orbit_return_t) :: zero_width, return_map
    real(dp) :: radius, q_safety, magnetic_shear, xi, eta, period_estimate
    real(dp) :: relative_period_error, relative_delta_phi_error
    integer :: status

    call configure_field(0.02_dp, field, radius)
    call field%evaluate([radius, 0.0_dp, 0.0_dp], sample, status)
    if (status /= GC_MODEL_SUCCESS) error stop 'manufactured passing sample failed'
    call field%radial_profiles(radius, q_safety, magnetic_shear)
    xi = 0.8_dp
    eta = (1.0_dp - xi**2)/sample%bmod
    call invariants_from_state(sample, 0.0_dp, rho0, 0.0_dp, 1.0_dp, xi, &
        invariants, status)
    if (status /= GC_MODEL_SUCCESS) error stop 'manufactured passing invariants failed'

    period_estimate = passing_period(field, radius, eta)
    call configure_options(radius, options)
    call compute_zero_width_passing_cycle(field, zero_potential, invariants, &
        [radius, 0.0_dp, 0.0_dp], 1, reference_velocity, 1, zero_width)
    if (zero_width%status /= GC_ORBIT_SUCCESS) then
        error stop 'zero-width passing cycle failed'
    end if

    call compute_return_map(field, zero_potential, invariants, &
        [radius, 0.0_dp, 0.0_dp], 1, rho0, 0.0_dp, reference_velocity, &
        GC_ORBIT_PASSING, 1, period_estimate, options, return_map)
    if (return_map%status /= GC_ORBIT_SUCCESS) then
        error stop 'zero-width return map failed'
    end if

    relative_period_error = abs(zero_width%period - period_estimate) &
        /max(abs(period_estimate), tiny(1.0_dp))
    relative_delta_phi_error = abs(zero_width%delta_phi &
        -2.0_dp*pi*q_safety)/max(1.0_dp, abs(2.0_dp*pi*q_safety))
    if (relative_period_error > 3.0e-8_dp) then
        error stop 'zero-width cycle disagrees with independent period quadrature'
    end if
    if (relative_delta_phi_error > 3.0e-8_dp) then
        error stop 'zero-width cycle has incorrect toroidal increment'
    end if
    if (abs(return_map%period - zero_width%period) &
            > 3.0e-8_dp*max(1.0_dp, abs(zero_width%period))) then
        error stop 'return map period disagrees with zero-width cycle'
    end if
    if (abs(return_map%delta_phi - zero_width%delta_phi) &
            > 3.0e-8_dp*max(1.0_dp, abs(zero_width%delta_phi))) then
        error stop 'return map increment disagrees with zero-width cycle'
    end if

    call pass_test

contains

    subroutine configure_field(epsilon, local_field, local_radius)
        real(dp), intent(in) :: epsilon
        type(manufactured_torus_field_t), intent(out) :: local_field
        real(dp), intent(out) :: local_radius

        local_radius = epsilon*major_radius
        local_field%major_radius = major_radius
        local_field%toroidal_flux = field_axis*major_radius
        local_field%q_axis = q_axis
        local_field%shear_shape = 0.0_dp
    end subroutine configure_field

    subroutine configure_options(local_radius, local_options)
        real(dp), intent(in) :: local_radius
        type(gc_orbit_options_t), intent(out) :: local_options

        local_options = gc_orbit_options_t()
        local_options%radial_min = 0.2_dp*local_radius
        local_options%radial_max = 1.8_dp*local_radius
        local_options%relative_tolerance = 1.0e-11_dp
        local_options%absolute_tolerance = 3.0e-12_dp
        local_options%event_relative_tolerance = 3.0e-12_dp
    end subroutine configure_options

    function passing_period(local_field, local_radius, local_eta) result(period)
        type(manufactured_torus_field_t), intent(in) :: local_field
        real(dp), intent(in) :: local_radius, local_eta
        real(dp) :: period

        integer, parameter :: n = 20000
        type(gc_field_sample_t) :: local_sample
        real(dp) :: theta, dtheta, time_sum
        integer :: k, local_status

        dtheta = 2.0_dp*pi/real(n, dp)
        time_sum = 0.0_dp
        do k = 1, n
            theta = (real(k, dp) - 0.5_dp)*dtheta
            call local_field%evaluate([local_radius, 0.0_dp, theta], &
                local_sample, local_status)
            if (local_status /= GC_MODEL_SUCCESS) then
                error stop 'passing period quadrature field failed'
            end if
            time_sum = time_sum + dtheta/(reference_velocity &
                *local_sample%hcon(3) &
                *sqrt(1.0_dp - local_eta*local_sample%bmod))
        end do
        period = time_sum
    end function passing_period

end program test_gc_orbit_limit
