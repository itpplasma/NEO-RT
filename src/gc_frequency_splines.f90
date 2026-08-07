module neort_gc_frequency_splines
    !! Surface adapter from the invariant-preserving GC frequency provider to
    !! NEO-RT's velocity/pitch full-orbit interface.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_frequency_provider, only: GC_FREQUENCY_SUCCESS, &
        gc_frequency_context_t, &
        initialize_gc_frequency_context, &
        evaluate_gc_full_orbit_frequency, &
        evaluate_gc_full_orbit_phase_average, gc_full_orbit_frequency_result_t
    use neort_gc_orbit_integrator, only: gc_orbit_average_t, &
        gc_orbit_perturbation_i

    implicit none
    private

    integer, parameter, public :: GC_SPLINE_SUCCESS = 0
    integer, parameter, public :: GC_SPLINE_PROVIDER_ERROR = 2
    integer, parameter, public :: GC_SPLINE_NOT_INITIALIZED = 3
    type(gc_frequency_context_t) :: surface_context
    logical :: surface_initialized = .false.
    !$omp threadprivate (surface_context, surface_initialized)

    public :: initialize_gc_spline_surface
    public :: evaluate_gc_full_orbit_frequency_surface
    public :: evaluate_gc_full_orbit_phase_average_surface

contains

    subroutine initialize_gc_spline_surface(surface, reference_theta, &
            field_scale, omega_e, particle_mass, particle_charge, velocity, &
            status, selected_frequency_model, wall_path, wall_units)
        real(dp), intent(in) :: surface, reference_theta, field_scale, omega_e
        real(dp), intent(in) :: particle_mass, particle_charge, velocity
        integer, intent(out) :: status
        integer, intent(in), optional :: selected_frequency_model
        character(len=*), intent(in), optional :: wall_path, wall_units

        integer :: provider_status

        surface_initialized = .false.
        call initialize_gc_frequency_context(surface, reference_theta, &
            field_scale, omega_e, particle_mass, particle_charge, velocity, &
            surface_context, provider_status, selected_frequency_model, &
            wall_path, wall_units)
        if (provider_status /= GC_FREQUENCY_SUCCESS) then
            status = GC_SPLINE_PROVIDER_ERROR
            return
        end if
        surface_initialized = .true.
        status = GC_SPLINE_SUCCESS
    end subroutine initialize_gc_spline_surface

    subroutine evaluate_gc_full_orbit_frequency_surface(velocity, eta, &
            parallel_direction, orbit_class, period_estimate, result, status)
        real(dp), intent(in) :: velocity, eta, period_estimate
        integer, intent(in) :: parallel_direction, orbit_class
        type(gc_full_orbit_frequency_result_t), intent(out) :: result
        integer, intent(out) :: status

        result = gc_full_orbit_frequency_result_t()
        status = GC_SPLINE_NOT_INITIALIZED
        if (.not. surface_initialized) return
        call evaluate_gc_full_orbit_frequency(surface_context, eta, &
            parallel_direction, orbit_class, period_estimate, result, status, &
            velocity)
    end subroutine evaluate_gc_full_orbit_frequency_surface

    subroutine evaluate_gc_full_orbit_phase_average_surface(velocity, eta, &
            parallel_direction, orbit_class, period_estimate, omega_b, omega_phi, mth, mph, &
            perturbation, result, status)
        real(dp), intent(in) :: velocity, eta, period_estimate, omega_b, omega_phi
        integer, intent(in) :: parallel_direction, orbit_class, mth, mph
        procedure(gc_orbit_perturbation_i) :: perturbation
        type(gc_orbit_average_t), intent(out) :: result
        integer, intent(out) :: status

        result = gc_orbit_average_t()
        status = GC_SPLINE_NOT_INITIALIZED
        if (.not. surface_initialized) return
        call evaluate_gc_full_orbit_phase_average(surface_context, eta, &
            parallel_direction, orbit_class, period_estimate, omega_b, omega_phi, mth, mph, &
            perturbation, result, status, velocity)
    end subroutine evaluate_gc_full_orbit_phase_average_surface

end module neort_gc_frequency_splines
