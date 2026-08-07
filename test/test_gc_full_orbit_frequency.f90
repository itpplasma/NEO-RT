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
    real(dp) :: scaled_full, scaled_narrow, scaled_tiny
    real(dp) :: extrapolated_coarse, extrapolated_fine
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
    call require(abs(full%omega_prec &
        - (full%omega_phi - context%q_fieldline*full%omega_b)) &
        < 2.0e-12_dp*max(1.0_dp, abs(full%omega_phi)), &
        'full-cycle precession identity')
    call evaluate_gc_full_orbit_frequency(context, eta, -1, &
        GC_ORBIT_PASSING, period_estimate, reverse, status)
    call require(status == GC_FREQUENCY_SUCCESS, 'reverse physical-width return')
    call require(reverse%omega_b < 0.0_dp, 'signed reverse bounce frequency')
    call require(abs(reverse%omega_phi/reverse%omega_b &
        - reverse%delta_phi/(-2.0_dp*pi)) < 2.0e-12_dp, &
        'reverse canonical return-map identity')

    ! Limiting expression for the physical launch family.  The local surface
    ! and pitch stay fixed while rho0, and therefore the launch P_phi shift,
    ! shrink together.  At fixed (H, mu, P_phi) for each orbit, smoothness of
    ! the return map gives F(a)=a*C0+O(a**2), where
    ! F=Omega_phi-q*Omega_b.  Richardson extrapolation of F(a)/a therefore
    ! supplies a one-sided independent limit oracle.
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
    residual_full = full%omega_prec
    residual_narrow = narrow%omega_prec
    residual_tiny = tiny_orbit%omega_prec
    scaled_full = residual_full
    scaled_narrow = residual_narrow/0.1_dp
    scaled_tiny = residual_tiny/0.01_dp
    extrapolated_coarse = (10.0_dp*scaled_narrow - scaled_full)/9.0_dp
    extrapolated_fine = (10.0_dp*scaled_tiny - scaled_narrow)/9.0_dp
    call require(abs(extrapolated_fine - scaled_tiny) &
        < abs(extrapolated_coarse - scaled_narrow), &
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
