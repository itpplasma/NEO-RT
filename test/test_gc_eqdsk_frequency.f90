program test_gc_eqdsk_frequency
    !! Regression for the provider's model-independent full-orbit path.
    !! The obsolete real-space-thin frequency contract is intentionally absent.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s, R0
    use neort_gc_frequency_provider, only: GC_FREQUENCY_SUCCESS, &
        gc_frequency_context_t, gc_full_orbit_frequency_result_t, &
        initialize_gc_frequency_context, evaluate_gc_full_orbit_frequency
    use neort_gc_orbit_integrator, only: GC_ORBIT_PASSING
    use neort_magfie, only: init_flux_surface_average
    use util, only: pi, qe, mu
    use util_for_test, only: pass_test

    implicit none

    type(gc_frequency_context_t) :: context
    type(gc_full_orbit_frequency_result_t) :: forward, reverse
    character(len=1024) :: eqdsk_file
    real(dp), parameter :: surface = 0.368_dp
    real(dp), parameter :: theta0 = 0.0_dp
    real(dp), parameter :: speed = 6.92e7_dp
    real(dp) :: eta, period_estimate
    integer :: status

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    if (len_trim(eqdsk_file) == 0) error stop 'EQDSK_FILE is required'
    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))
    call set_s(surface)
    call init_magfie_at_s()
    call init_flux_surface_average(surface)

    call initialize_gc_frequency_context(surface, theta0, 1.0_dp, 0.0_dp, &
        2.014_dp*mu, qe, speed, context, status, selected_frequency_model=0)
    call require(status == GC_FREQUENCY_SUCCESS, 'context initialization')
    eta = (1.0_dp - 0.8_dp**2)/context%reference_sample%bmod
    period_estimate = 12.0_dp*abs(context%q_fieldline)*R0/speed

    call evaluate_gc_full_orbit_frequency(context, eta, 1, GC_ORBIT_PASSING, &
        period_estimate, forward, status)
    call require(status == GC_FREQUENCY_SUCCESS, 'forward full-orbit return')
    call require(forward%period > 0.0_dp, 'positive return period')
    call require(abs(forward%omega_phi/forward%omega_b &
        - forward%delta_phi/(2.0_dp*pi)) < 2.0e-12_dp, &
        'forward canonical return-map identity')
    call require(abs(forward%omega_prec &
        - (forward%omega_phi - context%q_fieldline*forward%omega_b)) &
        < 2.0e-12_dp*max(1.0_dp, abs(forward%omega_phi)), &
        'forward full-cycle precession identity')

    call evaluate_gc_full_orbit_frequency(context, eta, -1, GC_ORBIT_PASSING, &
        period_estimate, reverse, status)
    call require(status == GC_FREQUENCY_SUCCESS, 'reverse full-orbit return')
    call require(reverse%period > 0.0_dp, 'reverse return period')
    call require(abs(reverse%omega_phi/reverse%omega_b &
        - reverse%delta_phi/(-2.0_dp*pi)) < 2.0e-12_dp, &
        'reverse canonical return-map identity')

    call pass_test

contains

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (.not. condition) then
            write(*, '(A)') 'FAIL: '//label
            error stop 1
        end if
    end subroutine require

end program test_gc_eqdsk_frequency
