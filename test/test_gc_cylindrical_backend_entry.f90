program test_gc_cylindrical_backend_entry
    !! Behavioral provider gate: model 2 must enter the cylindrical backend
    !! for both return frequencies and phase averages, while model 0 must
    !! leave it completely uninitialized.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s, R0, a
    use driftorbit, only: FREQUENCY_MODEL_LEGACY, FREQUENCY_MODEL_GC_FULL
    use neort_gc_frequency_provider, only: GC_FREQUENCY_SUCCESS, &
        GC_FREQUENCY_BACKEND_LEGACY, GC_FREQUENCY_BACKEND_CYLINDRICAL, &
        gc_frequency_context_t, gc_full_orbit_frequency_result_t, &
        initialize_gc_frequency_context, evaluate_gc_full_orbit_frequency, &
        evaluate_gc_full_orbit_phase_average, gc_frequency_runtime_metadata_t, &
        reset_gc_frequency_runtime_metadata, get_gc_frequency_runtime_metadata
    use neort_gc_cylindrical_backend, only: &
        GC_CYL_BACKEND_MEASURE_UNAVAILABLE, &
        gc_cylindrical_nonlocal_measure_result_t, &
        evaluate_gc_cylindrical_nonlocal_measure
    use neort_gc_orbit_integrator, only: GC_ORBIT_PASSING, GC_ORBIT_WALL_LOSS, &
        gc_orbit_average_t, gc_orbit_perturbation_i
    use neort_magfie, only: init_flux_surface_average
    use neort_orbit, only: th0
    use util, only: mu, pi, qe
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: surface = 0.368_dp
    real(dp), parameter :: speed = 6.92e7_dp
    real(dp), parameter :: mass = 2.014_dp*mu
    type(gc_frequency_context_t) :: legacy_context, full_context
    type(gc_full_orbit_frequency_result_t) :: full, reverse
    type(gc_orbit_average_t) :: average
    type(gc_frequency_runtime_metadata_t) :: metadata
    type(gc_cylindrical_nonlocal_measure_result_t) :: measure
    character(len=1024) :: eqdsk_file
    character(len=*), parameter :: large_wall_file = 'backend_entry_wall_large.dat'
    character(len=*), parameter :: narrow_wall_file = 'backend_entry_wall_narrow.dat'
    real(dp) :: eta, period_estimate
    integer :: status

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    if (len_trim(eqdsk_file) == 0) error stop 'EQDSK_FILE is required'
    call reset_gc_frequency_runtime_metadata()

    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))
    call set_s(surface)
    call init_magfie_at_s()
    call init_flux_surface_average(surface)
    call write_rectangle(large_wall_file, R0 - 2.0_dp*a, R0 + 2.0_dp*a, &
        -2.0_dp*a, 2.0_dp*a)
    call write_rectangle(narrow_wall_file, R0 - 0.1_dp*a, R0 + 0.1_dp*a, &
        -0.1_dp*a, 0.1_dp*a)

    call initialize_gc_frequency_context(surface, th0, 1.0_dp, 0.0_dp, &
        mass, qe, speed, legacy_context, status, &
        selected_frequency_model=FREQUENCY_MODEL_LEGACY)
    call require(status == GC_FREQUENCY_SUCCESS, 'model 0 initialization')
    call require(.not. legacy_context%cylindrical_backend%initialized, &
        'model 0 entered cylindrical backend')

    call initialize_gc_frequency_context(surface, th0, 1.0_dp, 0.0_dp, &
        mass, qe, speed, full_context, status, &
        selected_frequency_model=FREQUENCY_MODEL_GC_FULL, &
        wall_file=large_wall_file, wall_units='cm')
    call require(status == GC_FREQUENCY_SUCCESS, 'model 2 initialization')
    call require(full_context%cylindrical_backend%initialized, &
        'model 2 did not initialize cylindrical backend')
    call require(full_context%frequency_model == FREQUENCY_MODEL_GC_FULL, &
        'model 2 dispatch flag was not retained')
    call require(abs(full_context%cylindrical_backend%htheta_sign) == 1, &
        'physical htheta sign was not initialized')

    eta = (1.0_dp - 0.8_dp**2)/full_context%reference_sample%bmod
    period_estimate = 12.0_dp*abs(full_context%q_fieldline)*R0/speed
    call evaluate_gc_full_orbit_frequency(full_context, eta, 1, GC_ORBIT_PASSING, &
        period_estimate, full, status)
    call require(status == GC_FREQUENCY_SUCCESS, &
        'model 2 frequency provider failed')
    call require(full%backend_id == GC_FREQUENCY_BACKEND_CYLINDRICAL, &
        'model 2 frequency used legacy backend')
    call require(full%period > 0.0_dp .and. full%omega_b > 0.0_dp, &
        'model 2 frequency was not populated')
    call require(abs(full%energy_error) < 1.0e-7_dp, &
        'model 2 energy invariant was not checked')

    call evaluate_gc_full_orbit_frequency(full_context, eta, -1, GC_ORBIT_PASSING, &
        period_estimate, reverse, status)
    call require(status == GC_FREQUENCY_SUCCESS, &
        'model 2 reverse frequency provider failed')
    call require(reverse%backend_id == GC_FREQUENCY_BACKEND_CYLINDRICAL, &
        'reverse model 2 frequency used legacy backend')
    call require(reverse%omega_b < 0.0_dp, 'reverse model 2 frequency lost sign')
    call require(abs(reverse%omega_prec - &
        (reverse%delta_phi + 2.0_dp*pi*full_context%cylindrical_backend%q_fieldline) / &
        reverse%period) < 1.0e-8_dp, &
        'model 2 return frequency is inconsistent')

    call evaluate_gc_full_orbit_phase_average(full_context, eta, 1, &
        GC_ORBIT_PASSING, period_estimate, full%omega_b, full%omega_phi, 0, 0, &
        constant_perturbation, average, status)
    call require(status == GC_FREQUENCY_SUCCESS, &
        'model 2 phase-average provider failed')
    call require(abs(real(average%perturbation_average, dp) - 1.0_dp) < 2.0e-7_dp, &
        'model 2 phase average did not use cylindrical trajectory')
    call require(abs(average%inverse_b_average) > 0.0_dp .and. &
        abs(average%b_average) > 0.0_dp, &
        'model 2 phase average omitted cylindrical trajectory averages')

    call get_gc_frequency_runtime_metadata(metadata)
    call require(trim(metadata%backend) == 'cylindrical_full_fow', &
        'runtime metadata did not identify cylindrical backend')
    call require(trim(metadata%coordinates) == 'R,Z,phi', &
        'runtime metadata did not identify cylindrical coordinates')
    call require(metadata%wall_certified .and. len_trim(metadata%wall_hash) == 64, &
        'runtime metadata did not certify the loaded polygon wall')
    call require(trim(metadata%wall_units) == 'cm' .and. &
        trim(metadata%wall_backend_units) == 'cm', &
        'runtime metadata wall units are not explicit')
    call require(metadata%cylindrical_entry_count > 0 .and. &
        metadata%legacy_entry_count == 0, &
        'model 2 runtime dispatch counters are contaminated by legacy entry')
    call require(.not. metadata%canonical_measure_certified .and. &
        .not. metadata%component_identity_certified .and. &
        metadata%nonlocal_transport_required, &
        'model 2 exposed an uncertified local torque measure')
    call evaluate_gc_cylindrical_nonlocal_measure(full_context%cylindrical_backend, eta, &
        speed, 1, GC_ORBIT_PASSING, measure, status)
    call require(status == GC_CYL_BACKEND_MEASURE_UNAVAILABLE .and. &
        .not. measure%canonical_measure_certified .and. &
        .not. measure%component_identity_certified, &
        'uncertified local measure was not fail-closed')

    call initialize_gc_frequency_context(surface, th0, 1.0_dp, 0.0_dp, &
        mass, qe, speed, full_context, status, &
        selected_frequency_model=FREQUENCY_MODEL_GC_FULL, &
        wall_file=narrow_wall_file, wall_units='cm')
    call require(status == GC_FREQUENCY_SUCCESS, 'narrow wall initialization')
    call evaluate_gc_full_orbit_frequency(full_context, eta, 1, GC_ORBIT_PASSING, &
        period_estimate, full, status)
    call require(status /= GC_FREQUENCY_SUCCESS .and. &
        full%orbit_status == GC_ORBIT_WALL_LOSS, &
        'model 2 did not report a physical polygon-wall event')

    call delete_file(large_wall_file)
    call delete_file(narrow_wall_file)

    call pass_test

contains

    subroutine constant_perturbation(position, bmod, amplitude, status)
        real(dp), intent(in) :: position(3), bmod
        complex(dp), intent(out) :: amplitude
        integer, intent(out) :: status

        associate (unused_position => position, unused_bmod => bmod)
        end associate
        amplitude = cmplx(1.0_dp, 0.0_dp, dp)
        status = 0
    end subroutine constant_perturbation

    subroutine require(condition, label)
        logical, intent(in) :: condition
        character(*), intent(in) :: label

        if (.not. condition) then
            write(*, '(a)') 'FAIL: '//trim(label)
            error stop 1
        end if
    end subroutine require

    subroutine write_rectangle(path, r_left, r_right, z_bottom, z_top)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: r_left, r_right, z_bottom, z_top
        integer :: unit

        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, *) r_left, z_bottom
        write (unit, *) r_right, z_bottom
        write (unit, *) r_right, z_top
        write (unit, *) r_left, z_top
        close (unit)
    end subroutine write_rectangle

    subroutine delete_file(path)
        character(len=*), intent(in) :: path
        integer :: unit, ios

        open (newunit=unit, file=path, status='old', iostat=ios)
        if (ios == 0) close (unit, status='delete')
    end subroutine delete_file

end program test_gc_cylindrical_backend_entry
