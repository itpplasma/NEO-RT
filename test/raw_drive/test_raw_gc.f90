program test_raw_gc
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neo_perturbation_field, only: UNITS_SI, PERTFIELD_OK
    use neort_raw_drive, only: raw_source_t, raw_harmonics, RAW_OK
    use neort_raw_orbit, only: raw_background_t, raw_trajectory_t, raw_species_t, &
        raw_gc_point, raw_gc_trajectory, raw_gc_actions
    use test_raw_fixture, only: speed_unit, time_unit, charge_si, energy_unit, &
        moment, mass_si, epsilon, chi_scale, gauge_omega, circular_background, &
        circular_perturbation, linear_gauge, read_reference, species_si
    implicit none

    real(dp), allocatable :: reference(:, :), time(:), x(:, :), velocity(:, :)
    real(dp), allocatable :: b(:, :)
    type(raw_trajectory_t) :: orbit
    type(raw_source_t) :: source
    type(raw_species_t) :: species
    real(dp) :: omega_b, omega_phi, action_b, pphi, closure(3)
    integer :: ierr, n
    character(len=1024) :: path

    call get_command_argument(1, path)
    if (len_trim(path) == 0) error stop "raw GC reference path required"
    call read_reference(trim(path), reference)
    species = species_si()
    n = size(reference, 2)
    allocate (time(n), x(3, n), velocity(3, n), b(3, n))
    time = reference(1, :)*time_unit
    x = reference(2:4, :)
    velocity = reference(6:8, :)*speed_unit
    b = reference(9:11, :)
    call check_pointwise_velocity()
    call check_integrated_orbit()
    call check_physical_source()
    call check_actual_orbit_gauge(0)
    call check_actual_orbit_gauge(1)
    print *, "PASS test_raw_gc"

contains

    subroutine check_pointwise_velocity()
        type(raw_background_t) :: background
        real(dp) :: actual(3), acceleration, bstar_parallel, error
        integer :: k, status

        error = 0.0_dp
        do k = 1, n
            call circular_background(x(:, k), background, status)
            if (status /= RAW_OK) error stop "analytic background rejected"
            call raw_gc_point(background, species, moment, &
                reference(5, k)*speed_unit, actual, acceleration, &
                bstar_parallel, status)
            if (status /= RAW_OK) error stop "full guiding-center velocity rejected"
            error = max(error, maxval(abs(actual/speed_unit - reference(6:8, k))))
        end do
        call check_error("independent physical full-GC velocity", error, 5.0e-12_dp)
    end subroutine check_pointwise_velocity

    subroutine check_integrated_orbit()
        real(dp) :: initial(4), scaled(4, n), energy_error, momentum_error
        real(dp) :: radii(n), minimum, maximum
        initial(1:3) = reference(2:4, 1)
        initial(4) = reference(5, 1)*speed_unit
        call raw_gc_trajectory(circular_background, species, moment, time, &
            initial, 1.0e-11_dp, [1.0e-12_dp, 1.0e-12_dp, 1.0e-12_dp, 1.0e-8_dp], &
            orbit, ierr)
        if (ierr /= RAW_OK) error stop "full guiding-center integration failed"
        scaled = orbit%state
        scaled(4, :) = scaled(4, :)/speed_unit
        call check_error("independent closed radial-GC trajectory", &
            maxval(abs(scaled - reference(2:5, :))), 2.0e-7_dp)
        energy_error = (maxval(orbit%energy) - minval(orbit%energy)) &
            /(0.001_dp*energy_unit)
        momentum_error = (maxval(orbit%pphi) - minval(orbit%pphi))/abs(orbit%pphi(1))
        call check_error("full-GC energy invariant", energy_error, 2.0e-7_dp)
        call check_error("canonical toroidal momentum invariant", &
            momentum_error, 2.0e-7_dp)
        radii = hypot(orbit%state(1, :) - 1.0_dp, orbit%state(3, :))
        minimum = minval(radii)
        maximum = maxval(radii)
        if (maximum - minimum < 0.005_dp) error stop "radial orbit was frozen"
        call raw_gc_actions(orbit, [2.0e-7_dp, 2.0e-7_dp, 0.02_dp], 1, action_b, &
            pphi, omega_b, omega_phi, closure, ierr)
        if (ierr /= RAW_OK) error stop "primitive closed-orbit action rejected"
        call check_action_quadrature()
        print '(a,6(es22.14,a))', 'RAW_GC_CSV,', omega_b, ',', omega_phi, ',', &
            action_b, ',', pphi, ',', minimum, ',', maximum, ','
    end subroutine check_integrated_orbit

    subroutine check_action_quadrature()
        real(dp) :: density(n), r, p, bmag, canonical_r, canonical_z, oracle
        integer :: k

        do k = 1, n
            r = hypot(x(1, k) - 1.0_dp, x(3, k))
            p = r/(1.5_dp + 8.0_dp*r*r)
            bmag = sqrt(1.0_dp + p*p)/x(1, k)
            canonical_r = mass_si*reference(5, k)*speed_unit*b(1, k)/bmag
            canonical_z = -charge_si*log(x(1, k)) &
                + mass_si*reference(5, k)*speed_unit*b(3, k)/bmag
            density(k) = canonical_r*velocity(1, k) + canonical_z*velocity(3, k)
        end do
        oracle = sum(0.5_dp*(time(2:) - time(:n - 1)) &
            *(density(2:) + density(:n - 1)))/(2.0_dp*acos(-1.0_dp))
        call check_error("independent reduced-action quadrature", &
            abs(action_b - oracle)/(charge_si*0.01_dp), 2.0e-7_dp)
        call check_error("independent canonical Pphi", &
            abs(pphi - charge_si*flux_at_initial() &
            - mass_si*reference(5, 1)*speed_unit*x(1, 1)*b(2, 1) &
            /norm2(b(:, 1)))/abs(pphi), 1.0e-12_dp)
    end subroutine check_action_quadrature

    real(dp) function flux_at_initial() result(flux)
        real(dp) :: r

        r = hypot(x(1, 1) - 1.0_dp, x(3, 1))
        flux = log(1.0_dp + 8.0_dp*r*r/1.5_dp)/16.0_dp
    end function flux_at_initial

    subroutine check_physical_source()
        complex(dp) :: h(1), a(3, 1), db(3, 1), oracle
        real(dp) :: error, scale
        integer :: k, status

        call source%field%init_analytic(0.65_dp, 1.40_dp, 128, -0.4_dp, 0.4_dp, &
            128, [1], circular_perturbation, UNITS_SI, status, with_potential=.true.)
        if (status /= PERTFIELD_OK) error stop "physical raw source init failed"
        call source%set_frequencies([0.0_dp], status)
        if (status /= RAW_OK) error stop "physical raw source frequency rejected"
        call source%field%eval_modes(x(1, 1), x(3, 1), a, db, status)
        if (status /= PERTFIELD_OK) error stop "rational mode sample rejected"
        call check_error("finite rational magnetic mode retained", &
            abs(a(2, 1) - epsilon/x(1, 1))/(epsilon/x(1, 1)), 5.0e-6_dp)
        scale = maxval(abs(cmplx(reference(12, :), reference(13, :), dp)))
        error = 0.0_dp
        do k = 1, n
            call source%evaluate(x(:, k), time(k), b(:, k), velocity(:, k), &
                moment, charge_si, h, status)
            if (status /= RAW_OK) error stop "actual physical raw source rejected"
            oracle = cmplx(reference(12, k), reference(13, k), dp)
            error = max(error, abs(h(1)/energy_unit - oracle)/scale)
        end do
        call check_error("independent full-orbit physical energy samples", &
            error, 5.0e-6_dp)
    end subroutine check_physical_source

    subroutine check_actual_orbit_gauge(mb)
        integer, intent(in) :: mb
        complex(dp) :: harmonic(1, 1), values(n), chi_hat, endpoint, oracle
        real(dp) :: detuning, scale, frequencies(1)
        integer :: status, modes(1)

        gauge_omega = 0.0_dp
        if (mb == 0) gauge_omega = omega_phi
        call source%field%init_analytic(0.6_dp, 1.5_dp, 8, -0.5_dp, 0.5_dp, &
            8, [1], linear_gauge, UNITS_SI, status, with_potential=.true.)
        if (status /= PERTFIELD_OK) error stop "full-GC pure gauge init failed"
        frequencies(1) = gauge_omega
        modes(1) = mb
        call source%set_frequencies(frequencies, status)
        if (status /= RAW_OK) error stop "full-GC paired time gauge rejected"
        detuning = mb*omega_b + omega_phi - gauge_omega
        call gauge_oracle(detuning, values, chi_hat, endpoint, oracle)
        call raw_harmonics(source, time, x, velocity, b, moment, charge_si, modes, &
            omega_b, omega_phi, harmonic, status)
        if (status /= RAW_OK) error stop "actual radial-GC gauge projection failed"
        scale = charge_si*omega_b*chi_scale
        call check_error("actual radial-GC endpoint plus detuning", &
            abs(harmonic(1, 1) - oracle)/scale, 2.0e-7_dp)
        if (mb == 0) call check_error("actual radial-GC resonant pure gauge zero", &
            abs(harmonic(1, 1))/scale, 2.0e-7_dp)
    end subroutine check_actual_orbit_gauge

    subroutine gauge_oracle(detuning, values, chi_hat, endpoint, oracle)
        real(dp), intent(in) :: detuning
        complex(dp), intent(out) :: values(n), chi_hat, endpoint, oracle
        real(dp) :: period
        integer :: k

        do k = 1, n
            values(k) = chi_scale*x(1, k)*exp(cmplx(0.0_dp, &
                x(2, k) - (gauge_omega + detuning)*time(k), dp))
        end do
        period = time(n) - time(1)
        chi_hat = cmplx(0.0_dp, 0.0_dp, dp)
        do k = 2, n
            chi_hat = chi_hat + 0.5_dp*(time(k) - time(k - 1)) &
                *(values(k) + values(k - 1))/period
        end do
        endpoint = (values(n) - values(1))/period
        oracle = -charge_si*endpoint - cmplx(0.0_dp, charge_si*detuning, dp)*chi_hat
    end subroutine gauge_oracle

    subroutine check_error(label, error, tolerance)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: error, tolerance

        print '(a,1x,es12.4)', label, error
        if (.not. ieee_is_finite(error)) error stop "nonfinite raw GC check"
        if (error > tolerance) error stop "raw GC disagrees with independent oracle"
    end subroutine check_error

end program test_raw_gc
