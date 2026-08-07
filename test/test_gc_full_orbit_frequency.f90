program test_gc_full_orbit_frequency
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_frequency_provider, only: GC_FREQUENCY_SUCCESS, &
        GC_FREQUENCY_BACKEND_LEGACY, &
        GC_FREQUENCY_ORBIT_ERROR, gc_frequency_context_t, &
        gc_full_orbit_frequency_result_t, initialize_gc_frequency_context, &
        evaluate_gc_full_orbit_frequency
    use neort_gc_orbit_integrator, only: GC_ORBIT_PASSING
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s, R0
    use neort_magfie, only: init_flux_surface_average
    use util, only: pi, qe, mu

    implicit none

    type(gc_frequency_context_t) :: context, narrow_context, tiny_context
    type(gc_full_orbit_frequency_result_t) :: full, reverse, narrow, tiny_orbit
    character(len=1024) :: eqdsk_file
    real(dp), parameter :: surface = 0.368_dp
    real(dp), parameter :: theta0 = 0.0_dp
    real(dp), parameter :: speed = 6.92e7_dp
    real(dp) :: eta, period_estimate
    real(dp) :: residual_full, residual_narrow, residual_tiny
    integer :: status

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    if (len_trim(eqdsk_file) == 0) error stop 'EQDSK_FILE is required'
    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))
    call set_s(surface)
    call init_magfie_at_s()
    call init_flux_surface_average(surface)

    call initialize_gc_frequency_context(surface, theta0, 1.0_dp, 0.0_dp, &
        2.014_dp*mu, qe, speed, context, status, &
        selected_frequency_model=0)
    call require(status == GC_FREQUENCY_SUCCESS, 'context initialization')
    eta = (1.0_dp - 0.8_dp**2)/context%reference_sample%bmod
    period_estimate = 12.0_dp*abs(context%q_fieldline)*R0/speed
    call evaluate_gc_full_orbit_frequency(context, eta, 1, &
        GC_ORBIT_PASSING, period_estimate, full, status)
    call require(status == GC_FREQUENCY_SUCCESS, 'physical-width return')
    call require(full%backend_id == GC_FREQUENCY_BACKEND_LEGACY, &
        'legacy regression did not stay on legacy backend')
    call require(full%omega_b > 0.0_dp, 'positive bounce frequency')
    call require(abs(full%omega_phi/full%omega_b &
        - full%delta_phi/(2.0_dp*pi)) < 2.0e-12_dp, &
        'canonical return-map identity')
    call evaluate_gc_full_orbit_frequency(context, eta, -1, &
        GC_ORBIT_PASSING, period_estimate, reverse, status)
    call require(status == GC_FREQUENCY_SUCCESS, 'reverse physical-width return')
    call require(reverse%omega_b < 0.0_dp, 'signed reverse bounce frequency')
    call require(abs(reverse%omega_phi/reverse%omega_b &
        - reverse%delta_phi/(-2.0_dp*pi)) < 2.0e-12_dp, &
        'reverse canonical return-map identity')

    ! Limiting expression for the physical launch family: the local surface
    ! and pitch stay fixed while rho0, and therefore the launch P_phi shift,
    ! shrink together.  Successive residual/rho0 values must converge.
    narrow_context = context
    narrow_context%rho0 = 0.1_dp*context%rho0
    call evaluate_gc_full_orbit_frequency(narrow_context, eta, 1, &
        GC_ORBIT_PASSING, period_estimate, narrow, status)
    call require(status == GC_FREQUENCY_SUCCESS, 'narrow-width return')
    tiny_context = context
    tiny_context%rho0 = 0.01_dp*context%rho0
    call evaluate_gc_full_orbit_frequency(tiny_context, eta, 1, &
        GC_ORBIT_PASSING, period_estimate, tiny_orbit, status)
    call require(status == GC_FREQUENCY_SUCCESS, 'tiny-width return')
    residual_full = full%omega_phi &
        - context%q_fieldline*full%omega_b
    residual_narrow = narrow%omega_phi &
        - narrow_context%q_fieldline*narrow%omega_b
    residual_tiny = tiny_orbit%omega_phi &
        - tiny_context%q_fieldline*tiny_orbit%omega_b
    call require(abs(residual_tiny/0.01_dp - residual_narrow/0.1_dp) &
        < abs(residual_narrow/0.1_dp - residual_full), &
        'finite-width limiting expression')

    ! A physical orbit that cannot return is data, not a zero-frequency row.
    narrow_context%orbit_options%max_periods = 1.0e-6_dp
    call evaluate_gc_full_orbit_frequency(narrow_context, eta, 1, &
        GC_ORBIT_PASSING, period_estimate, narrow, status)
    call require(status == GC_FREQUENCY_ORBIT_ERROR, 'non-return status')
    call require(narrow%omega_b == 0.0_dp .and. narrow%omega_phi == 0.0_dp, &
        'non-return carries no fabricated frequency')

contains

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label
        if (.not. condition) then
            write(*, '(A)') 'FAIL: '//label
            error stop 1
        end if
    end subroutine require

end program test_gc_full_orbit_frequency
