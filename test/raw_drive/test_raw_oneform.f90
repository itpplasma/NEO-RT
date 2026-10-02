program test_raw_oneform
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neo_perturbation_field, only: UNITS_GAUSSIAN, UNITS_SI, PERTFIELD_OK
    use neort_raw_drive, only: raw_source_t, raw_harmonics, RAW_OK
    implicit none

    integer, parameter :: npoint = 2049, mb = 1
    integer, parameter :: gauge_modes(1) = [1], bounce_modes(1) = [mb]
    real(dp), parameter :: pi = acos(-1.0_dp), light = 29979245800.0_dp
    real(dp), parameter :: charge = 4.803204712570263e-10_dp
    real(dp), parameter :: radius = 150.0_dp, excursion = 10.0_dp
    real(dp), parameter :: period = 1.0e-4_dp, beta = 0.3_dp, chi_scale = 2.0_dp
    real(dp), parameter :: om_b = 2.0_dp*pi/period, om_phi = 0.25_dp*om_b
    type(raw_source_t) :: source
    real(dp) :: time(npoint), position(3, npoint), velocity(3, npoint)
    real(dp) :: background(3, npoint), source_omega
    complex(dp) :: harmonic(1, 1), expected, chi_harmonic, endpoint
    integer :: ierr

    call check_energy_units()
    call check_si_vector_energy()
    call check_magnetic_energy()
    source_omega = mb*om_b + om_phi
    call initialize_gauge()
    call prescribed_orbit(period)
    call integrate_gauge()
    call check_close("closed on-resonance pure gauge", harmonic(1, 1), &
        cmplx(0.0_dp, 0.0_dp, dp), charge/light*om_b*chi_scale*radius, 1.0e-10_dp)

    source_omega = 0.0_dp
    call initialize_gauge()
    call integrate_gauge()
    chi_harmonic = chi_scale*(radius*bessel_jn(mb, beta) &
        + excursion/2.0_dp*(bessel_jn(mb - 1, beta) + bessel_jn(mb + 1, beta)))
    expected = -cmplx(0.0_dp, charge/light*(mb*om_b + om_phi), dp)*chi_harmonic
    call check_close("off-resonance gauge detuning", harmonic(1, 1), expected, &
        abs(expected), 1.0e-10_dp)
    if (abs(harmonic(1, 1)) < 0.5_dp*abs(expected)) then
        error stop "raw off-resonance gauge must retain parallel content"
    end if

    source_omega = mb*om_b + om_phi
    call initialize_gauge()
    call prescribed_orbit(0.75_dp*period)
    call integrate_gauge()
    endpoint = (gauge_value(npoint) - gauge_value(1))/time(npoint)
    expected = -charge/light*endpoint
    call check_close("open-segment gauge endpoint", harmonic(1, 1), expected, &
        abs(expected), 1.0e-5_dp)
    print *, "PASS test_raw_oneform"

contains

    subroutine initialize_gauge()
        real(dp) :: frequencies(1)

        call source%field%init_analytic(radius - 2.0_dp*excursion, &
            radius + 2.0_dp*excursion, 8, -2.0_dp*excursion, 2.0_dp*excursion, &
            8, gauge_modes, pure_gauge_amplitude, UNITS_GAUSSIAN, ierr, &
            with_potential=.true.)
        if (ierr /= PERTFIELD_OK) error stop "pure-gauge source initialization failed"
        frequencies(1) = source_omega
        call source%set_frequencies(frequencies, ierr)
        if (ierr /= RAW_OK) error stop "pure-gauge frequencies rejected"
    end subroutine initialize_gauge

    subroutine pure_gauge_amplitude(n, r_cyl, z_cyl, a, potential)
        integer, intent(in) :: n
        real(dp), intent(in) :: r_cyl, z_cyl
        complex(dp), intent(out) :: a(3), potential

        ! chi_n=chi_scale*R: chi=chi_scale*(x+i*y) in Cartesian space.
        ! Its physical cylindrical gradient is constant; the spline represents
        ! it exactly, and curl(grad chi)=0 independent of the sampled orbit.
        if (n /= 1) error stop "pure gauge has toroidal mode one"
        if (abs(z_cyl) > 2.0_dp*excursion) error stop "unexpected gauge grid"
        a(1) = cmplx(chi_scale, 0.0_dp, dp)
        a(2) = cmplx(0.0_dp, chi_scale, dp)
        a(3) = cmplx(0.0_dp, 0.0_dp, dp)
        potential = cmplx(0.0_dp, source_omega*chi_scale*r_cyl/light, dp)
    end subroutine pure_gauge_amplitude

    subroutine prescribed_orbit(span)
        real(dp), intent(in) :: span
        real(dp) :: angle
        integer :: k

        ! Declared kinematics, independent of any production orbit solver.
        do k = 1, npoint
            time(k) = span*real(k - 1, dp)/real(npoint - 1, dp)
            angle = om_b*time(k)
            position(1, k) = radius + excursion*cos(angle)
            position(2, k) = om_phi*time(k) + beta*sin(angle)
            position(3, k) = excursion*sin(angle)
            velocity(1, k) = -excursion*om_b*sin(angle)
            velocity(2, k) = position(1, k)*(om_phi + beta*om_b*cos(angle))
            velocity(3, k) = excursion*om_b*cos(angle)
            background(:, k) = [0.0_dp, 1.0_dp, 0.0_dp]
        end do
    end subroutine prescribed_orbit

    subroutine integrate_gauge()
        call raw_harmonics(source, time, position, velocity, background, &
            1.0e-8_dp, charge, bounce_modes, om_b, om_phi, harmonic, ierr)
        if (ierr /= RAW_OK) error stop "raw gauge quadrature failed"
    end subroutine integrate_gauge

    complex(dp) function gauge_value(k) result(chi)
        integer, intent(in) :: k

        chi = chi_scale*position(1, k)*exp(cmplx(0.0_dp, &
            position(2, k) - source_omega*time(k), dp))
    end function gauge_value

    subroutine check_energy_units()
        complex(dp) :: h(3)
        integer :: status

        call source%field%init_analytic(radius - excursion, radius + excursion, &
            8, -excursion, excursion, 8, [-2, 0, 1], electrostatic_amplitude, &
            UNITS_GAUSSIAN, status, with_potential=.true.)
        if (status /= PERTFIELD_OK) error stop "electrostatic source init failed"
        call source%set_frequencies([0.0_dp, 0.0_dp, 0.0_dp], status)
        if (status /= RAW_OK) error stop "signed mode frequencies rejected"
        call source%evaluate([radius, 0.0_dp, 0.0_dp], 0.0_dp, &
            [0.0_dp, 1.0_dp, 0.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], &
            0.0_dp, charge, h, status)
        if (status /= RAW_OK) error stop "electrostatic raw source failed"
        expected = cmplx(3.0_dp*charge, -2.0_dp*charge, dp)
        call check_close("Gaussian potential energy, negative n", h(1), &
            expected, abs(expected), 1.0e-12_dp)
        call check_close("Gaussian potential energy, n zero", h(2), &
            expected, abs(expected), 1.0e-12_dp)
        call check_close("Gaussian potential energy, positive n", h(3), &
            expected, abs(expected), 1.0e-12_dp)
        call check_signed_phase()
    end subroutine check_energy_units

    subroutine check_signed_phase()
        complex(dp) :: h(3), oracle
        real(dp), parameter :: angle = pi/3.0_dp
        integer :: status, k

        call source%evaluate([radius, angle, 0.0_dp], 0.0_dp, &
            [0.0_dp, 1.0_dp, 0.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], &
            0.0_dp, charge, h, status)
        if (status /= RAW_OK) error stop "signed-phase evaluation failed"
        do k = 1, 3
            oracle = charge*cmplx(3.0_dp, -2.0_dp, dp) &
                *exp(cmplx(0.0_dp, source%field%ntor(k)*angle, dp))
            call check_close("physical signed toroidal phase", h(k), oracle, &
                abs(oracle), 1.0e-12_dp)
        end do
    end subroutine check_signed_phase

    subroutine electrostatic_amplitude(n, r_cyl, z_cyl, a, potential)
        integer, intent(in) :: n
        real(dp), intent(in) :: r_cyl, z_cyl
        complex(dp), intent(out) :: a(3), potential

        if (abs(n) > 2) error stop "unexpected electrostatic mode"
        if (r_cyl <= 0.0_dp) error stop "invalid cylindrical radius"
        if (abs(z_cyl) > excursion) error stop "unexpected electrostatic grid"
        a = cmplx(0.0_dp, 0.0_dp, dp)
        potential = cmplx(3.0_dp, -2.0_dp, dp)
    end subroutine electrostatic_amplitude

    subroutine check_si_vector_energy()
        real(dp), parameter :: charge_si = 1.602176634e-19_dp
        complex(dp) :: h(1), oracle
        integer :: status

        call source%field%init_analytic(radius - excursion, radius + excursion, &
            8, -excursion, excursion, 8, [0], uniform_si_amplitude, &
            UNITS_SI, status, with_potential=.true.)
        if (status /= PERTFIELD_OK) error stop "SI source initialization failed"
        call source%set_frequencies([0.0_dp], status)
        if (status /= RAW_OK) error stop "SI mode frequencies rejected"
        call source%evaluate([radius, 0.0_dp, 0.0_dp], 0.0_dp, &
            [0.0_dp, 1.0_dp, 0.0_dp], [4.0_dp, 0.0_dp, 5.0_dp], &
            0.0_dp, charge_si, h, status)
        if (status /= RAW_OK) error stop "SI raw source evaluation failed"
        oracle = -charge_si*cmplx(8.0_dp, -15.0_dp, dp)
        call check_close("SI vector-potential energy in J", h(1), oracle, &
            abs(oracle), 1.0e-12_dp)
    end subroutine check_si_vector_energy

    subroutine uniform_si_amplitude(n, r_cyl, z_cyl, a, potential)
        integer, intent(in) :: n
        real(dp), intent(in) :: r_cyl, z_cyl
        complex(dp), intent(out) :: a(3), potential

        if (n /= 0) error stop "uniform SI control has n zero"
        if (r_cyl <= 0.0_dp) error stop "invalid SI control radius"
        if (abs(z_cyl) > excursion) error stop "unexpected SI grid"
        a = [cmplx(2.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 0.0_dp, dp), &
            cmplx(0.0_dp, -3.0_dp, dp)]
        potential = cmplx(0.0_dp, 0.0_dp, dp)
    end subroutine uniform_si_amplitude

    subroutine check_magnetic_energy()
        real(dp), parameter :: moment = 3.0e-8_dp
        complex(dp) :: h(1), oracle
        integer :: status

        call source%field%init_analytic(radius - excursion, radius + excursion, &
            8, -excursion, excursion, 8, [0], uniform_magnetic_amplitude, &
            UNITS_GAUSSIAN, status, with_potential=.true.)
        if (status /= PERTFIELD_OK) error stop "magnetic source initialization failed"
        call source%set_frequencies([0.0_dp], status)
        if (status /= RAW_OK) error stop "magnetic source frequency rejected"
        call source%evaluate([radius, 0.0_dp, 0.0_dp], 0.0_dp, &
            [0.0_dp, 0.0_dp, 9.0_dp], [0.0_dp, 0.0_dp, 0.0_dp], &
            moment, charge, h, status)
        if (status /= RAW_OK) error stop "magnetic source evaluation failed"
        oracle = moment*cmplx(5.0_dp, 2.0_dp, dp)
        call check_close("Gaussian magnetic-moment energy in erg", h(1), oracle, &
            abs(oracle), 1.0e-12_dp)
    end subroutine check_magnetic_energy

    subroutine uniform_magnetic_amplitude(n, r_cyl, z_cyl, a, potential)
        integer, intent(in) :: n
        real(dp), intent(in) :: r_cyl, z_cyl
        complex(dp), intent(out) :: a(3), potential

        if (n /= 0) error stop "uniform magnetic control has n zero"
        if (abs(z_cyl) > excursion) error stop "unexpected magnetic grid"
        a = cmplx(0.0_dp, 0.0_dp, dp)
        a(2) = 0.5_dp*r_cyl*cmplx(5.0_dp, 2.0_dp, dp)
        potential = cmplx(0.0_dp, 0.0_dp, dp)
    end subroutine uniform_magnetic_amplitude

    subroutine check_close(label, actual, oracle, scale, tolerance)
        character(len=*), intent(in) :: label
        complex(dp), intent(in) :: actual, oracle
        real(dp), intent(in) :: scale, tolerance
        real(dp) :: error

        if (.not. ieee_is_finite(real(actual, dp))) error stop "nonfinite raw drive"
        if (.not. ieee_is_finite(aimag(actual))) error stop "nonfinite raw drive"
        error = abs(actual - oracle)/scale
        print '(a,1x,es12.4)', label, error
        if (error > tolerance) then
            error stop "raw one-form disagrees with independent oracle"
        end if
    end subroutine check_close

end program test_raw_oneform
