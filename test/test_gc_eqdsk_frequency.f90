program test_gc_eqdsk_frequency
    !! End-to-end real-space frequency gate on an independently generated
    !! circular GEQDSK.  The old fixed-s bounce integrator is used only as an
    !! independent oracle for the lambda=0 return period.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s, R0
    use driftorbit, only: etatp, etadt
    use neort_gc_dynamics, only: gc_field_sample_t
    use neort_gc_frequency_provider, only: GC_FREQUENCY_SUCCESS, &
        gc_frequency_context_t, gc_frequency_result_t, &
        initialize_gc_frequency_context, evaluate_gc_frequency
    use neort_gc_models, only: GC_MODEL_SUCCESS
    use neort_gc_orbit_integrator, only: GC_ORBIT_TRAPPED, GC_ORBIT_PASSING
    use neort_magfie, only: init_flux_surface_average
    use neort_orbit, only: bounce_time, th0
    use util, only: pi, qe, mu
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: surface = 0.25_dp
    real(dp), parameter :: velocity = 1.0e8_dp
    real(dp), parameter :: mass = 2.014_dp*mu
    real(dp), parameter :: electric_frequency = 4.0e3_dp
    type(gc_frequency_context_t) :: positive, negative, near_separatrix
    type(gc_frequency_result_t) :: trapped, passing, reversed, near_result
    type(gc_field_sample_t) :: sample_min, sample_max
    character(len=1024) :: eqdsk_file
    real(dp) :: eta_trapped, eta_passing, epsilon
    real(dp) :: period_estimate, legacy_period, legacy_frequency
    real(dp) :: quadrature_period, quadrature_frequency
    integer :: status

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    if (len_trim(eqdsk_file) == 0) then
        error stop 'GEQDSK test fixture is not configured'
    end if
    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))
    call set_s(surface)
    call init_magfie_at_s()
    call init_flux_surface_average(surface)

    call initialize_gc_frequency_context(surface, th0, 1.0_dp, electric_frequency, &
        mass, qe, velocity, positive, status, selected_frequency_model=1)
    if (status /= GC_FREQUENCY_SUCCESS) error stop 'positive GC context failed'
    call positive%field%evaluate([surface, 0.0_dp, 0.0_dp], sample_min, status)
    if (status /= GC_MODEL_SUCCESS) error stop 'GC Bmin sample failed'
    call positive%field%evaluate([surface, 0.0_dp, pi], sample_max, status)
    if (status /= GC_MODEL_SUCCESS) error stop 'GC Bmax sample failed'
    if (abs(sample_min%bmod - 1.0_dp/etadt) &
        > 2.0e-5_dp*sample_min%bmod) then
        error stop 'GC and legacy Bmin disagree'
    end if
    if (abs(sample_max%bmod - 1.0_dp/etatp) &
        > 2.0e-5_dp*sample_max%bmod) then
        error stop 'GC and legacy Bmax disagree'
    end if

    eta_trapped = 1.0_dp/(sample_min%bmod &
        + 0.4_dp*(sample_max%bmod - sample_min%bmod))
    epsilon = (sample_max%bmod - sample_min%bmod) &
        /(sample_max%bmod + sample_min%bmod)
    period_estimate = 30.0_dp*abs(positive%q_fieldline)*R0/velocity &
        /sqrt(epsilon)
    call evaluate_gc_frequency(positive, eta_trapped, 1, &
        GC_ORBIT_TRAPPED, period_estimate, trapped, status)
    call require_frequency('trapped', trapped, status)
    legacy_period = bounce_time(velocity, eta_trapped)
    legacy_frequency = 2.0_dp*pi/legacy_period
    if (abs(trapped%omega_b/legacy_frequency - 1.0_dp) > 3.0e-4_dp) then
        write(*,*) 'trapped GC/legacy bounce:', trapped%omega_b, legacy_frequency
        error stop 'GC trapped zero-width period disagrees with independent solver'
    end if
    call require_electric('trapped', trapped%omega_electric)

    eta_passing = (1.0_dp - 0.8_dp**2)/sample_min%bmod
    period_estimate = 12.0_dp*abs(positive%q_fieldline)*R0/velocity
    call evaluate_gc_frequency(positive, eta_passing, 1, &
        GC_ORBIT_PASSING, period_estimate, passing, status)
    call require_frequency('passing', passing, status)
    legacy_period = bounce_time(velocity, eta_passing)
    legacy_frequency = 2.0_dp*pi/legacy_period
    quadrature_period = zero_width_passing_period(positive, eta_passing, velocity)
    quadrature_frequency = 2.0_dp*pi/quadrature_period
    ! Independent full-cycle expression: omega_b=2*pi/tau, with
    ! tau=(1/v)*integral[dtheta/(xi*h^theta)] on the passing section.
    if (abs(passing%omega_b/quadrature_frequency - 1.0_dp) > 2.0e-6_dp) then
        write(*,*) 'passing GC/analytic quadrature:', passing%omega_b, &
            quadrature_frequency
        error stop 'GC passing full-cycle limit disagrees with analytic expression'
    end if
    if (abs(passing%omega_b/legacy_frequency - 1.0_dp) > 3.0e-4_dp) then
        write(*,*) 'passing GC/legacy transit:', passing%omega_b, legacy_frequency
        error stop 'GC passing zero-width period disagrees with independent solver'
    end if
    call require_electric('passing', passing%omega_electric)

    call initialize_gc_frequency_context(surface, th0, 1.0_dp, electric_frequency, &
        mass, -qe, velocity, negative, status, selected_frequency_model=1)
    if (status /= GC_FREQUENCY_SUCCESS) error stop 'negative GC context failed'
    period_estimate = 30.0_dp*abs(negative%q_fieldline)*R0/velocity/sqrt(epsilon)
    call evaluate_gc_frequency(negative, eta_trapped, 1, &
        GC_ORBIT_TRAPPED, period_estimate, reversed, status)
    call require_frequency('charge-reversed', reversed, status)
    if (abs(reversed%omega_magnetic + trapped%omega_magnetic) &
        > 3.0e-3_dp*abs(trapped%omega_magnetic)) then
        write(*,*) 'magnetic charge reversal:', trapped%omega_magnetic, &
            reversed%omega_magnetic
        error stop 'GEQDSK magnetic precession did not reverse with charge'
    end if
    call require_electric('charge-reversed', reversed%omega_electric)

    ! This surface is close to the passing/trapped boundary in the ITER
    ! GEQDSK and exercises the adaptive small-lambda ladder that the
    ! production real-space sweep needs.  The assertion is deliberately on
    ! the provider contract; the circular legacy fixture is not a substitute
    ! for the real-space return-map result at this shaped surface.
    call set_s(0.368_dp)
    call init_magfie_at_s()
    call init_flux_surface_average(0.368_dp)
    call initialize_gc_frequency_context(0.368_dp, th0, 1.0_dp, &
        electric_frequency, mass, qe, velocity, near_separatrix, status, &
        selected_frequency_model=1)
    if (status /= GC_FREQUENCY_SUCCESS) then
        error stop 'near-separatrix GC context failed'
    end if
    call near_separatrix%field%evaluate([0.368_dp, 0.0_dp, pi], &
        sample_max, status)
    if (status /= GC_MODEL_SUCCESS) error stop 'near-separatrix Bmax failed'
    eta_passing = 0.98_dp/sample_max%bmod
    period_estimate = 12.0_dp*abs(near_separatrix%q_fieldline)*R0/velocity
    call evaluate_gc_frequency(near_separatrix, eta_passing, 1, &
        GC_ORBIT_PASSING, period_estimate, near_result, status)
    call require_frequency('near-separatrix passing', near_result, status)
    call require_electric('near-separatrix passing', near_result%omega_electric)

    write(*, '(a,3es14.5)') 'GEQDSK trapped omega_b, omega_B, omega_E: ', &
        trapped%omega_b, trapped%omega_magnetic, trapped%omega_electric
    write(*, '(a,3es14.5)') 'GEQDSK passing omega_b, omega_B, omega_E: ', &
        passing%omega_b, passing%omega_magnetic, passing%omega_electric
    call pass_test

contains

    subroutine require_frequency(label, value, local_status)
        character(*), intent(in) :: label
        type(gc_frequency_result_t), intent(in) :: value
        integer, intent(in) :: local_status

        if (local_status /= GC_FREQUENCY_SUCCESS) then
            write(*,*) trim(label)//' provider/statuses:', local_status, &
                value%magnetic_limit_status, value%total_limit_status, &
                value%baseline_residual
            error stop 'GEQDSK GC frequency evaluation failed'
        end if
        if (abs(value%baseline_residual) > 2.0e-7_dp) then
            error stop 'GEQDSK zero-width topological residual is nonzero'
        end if
    end subroutine require_frequency

    subroutine require_electric(label, value)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value

        if (abs(value/electric_frequency - 1.0_dp) > 5.0e-3_dp) then
            write(*,*) trim(label)//' electric/reference:', value, electric_frequency
            error stop 'GEQDSK electric precession has wrong sign or magnitude'
        end if
    end subroutine require_electric

    function zero_width_passing_period(context, eta, velocity) result(period)
        type(gc_frequency_context_t), intent(in) :: context
        real(dp), intent(in) :: eta, velocity
        real(dp) :: period
        integer, parameter :: ntheta = 32768
        type(gc_field_sample_t) :: sample
        real(dp) :: dtheta, theta, xi_squared
        integer :: k, local_status

        period = 0.0_dp
        dtheta = 2.0_dp*pi/real(ntheta, dp)
        do k = 1, ntheta
            theta = context%reference_position(3) &
                + (real(k, dp) - 0.5_dp)*dtheta
            call context%field%evaluate([context%reference_position(1), &
                context%reference_position(2), theta], sample, local_status)
            if (local_status /= GC_MODEL_SUCCESS) error stop 'quadrature field failed'
            xi_squared = 1.0_dp - eta*sample%bmod
            if (xi_squared <= 0.0_dp) error stop 'quadrature entered trapped region'
            period = period + dtheta/(velocity*sqrt(xi_squared) &
                *abs(sample%hcon(3)))
        end do
    end function zero_width_passing_period

end program test_gc_eqdsk_frequency
