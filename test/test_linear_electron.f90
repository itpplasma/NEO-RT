program test_linear_electron
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite, ieee_invalid, &
        ieee_set_halting_mode
    use do_magfie_mod, only: s, sign_theta, q, psi_pr, R0, iota, a
    use driftorbit, only: nonlin, dVds, B0, mph, mth, sign_vpar
    use neort_profiles, only: prepare_plasma_splines, init_plasma_at_s, &
        init_thermodynamic_forces, ni1, ni2, Ti1, Te, vth, dvthds, &
        Om_tE, A1, A2, prepare_profile_splines
    use neort_lib, only: config_t, neort_init, neort_setup_at_s
    use neort, only: compute_transport_harmonic, set_to_passing_region
    use neort_freq, only: Om_th, Om_tB
    use neort_resonance, only: resonance_value
    use neort_datatypes, only: transport_harmonic_t
    use neort_transport, only: D11int, D12int, Tphi_int
    use util, only: qi, qe, mi, mu, ev, c, pi
    implicit none

    integer, parameter :: n = 6
    real(dp), parameter :: electron_mass_u = 0.0005485861437678632_dp
    real(dp), parameter :: omega_mars = 8.0e3_dp
    real(dp) :: plasma(n, 6), rotation(n, 2), expected_vth, expected_a1
    real(dp) :: Dco(2), Dctr(2), Dt(2), Tco, Tctr, Tt, raw_D(2), scale
    real(dp) :: total_D(3)
    real(dp) :: eta_min, eta_max, eta, Omth, OmtB, dv, deta
    real(dp) :: ux_mid, residual, residual_derivative, frequency_scale
    integer :: i, j, direction
    character(32) :: mode
    type(config_t) :: config
    type(transport_harmonic_t) :: harmonic

    ! The old unused ion-collision initialization takes sqrt(negative density).
    call ieee_set_halting_mode(ieee_invalid, .true.)
    mph = 3
    call get_command_argument(1, mode)
    nonlin = trim(mode) == "nonlinear"
    do i = 1, n
        plasma(i, 1) = real(i - 1, dp) / real(n - 1, dp)
        plasma(i, 2) = 8.0e13_dp * (1.0_dp - 0.2_dp * plasma(i, 1))
        plasma(i, 3) = 0.0_dp
        plasma(i, 4) = 8.0e3_dp * (1.0_dp - 0.3_dp * plasma(i, 1))
        plasma(i, 5) = 1.0_dp
        plasma(i, 6) = plasma(i, 4)
    end do
    call prepare_plasma_splines(n, electron_mass_u, 2.0_dp, -1.0_dp, &
        1.0_dp, plasma)
    s = 0.4_dp
    call init_plasma_at_s()
    if (nonlin) error stop "nonlinear electron initialization was accepted"

    expected_vth = sqrt(2.0_dp * 8.0e3_dp * (1.0_dp - 0.3_dp * s) * ev &
        / (electron_mass_u * mu))
    call close("electron charge", qi, -qe)
    call close("electron mass", mi, electron_mass_u * mu)
    call close("electron thermal velocity", vth, expected_vth)
    call close("thermal velocity derivative", dvthds, &
        -0.15_dp * expected_vth / (1.0_dp - 0.3_dp * s))
    call close("physical density", ni1, 8.0e13_dp * (1.0_dp - 0.2_dp * s))
    if (ni2 /= 0.0_dp) error stop "electron fixture acquired a neutralizer"
    call close("electron temperature", Ti1, Te)
    Om_tE = -omega_mars
    call init_thermodynamic_forces(2.0e7_dp, 1.8_dp)
    expected_a1 = -0.2_dp / (1.0_dp - 0.2_dp * s) &
        - (-qe) / (Ti1 * ev) * sign_theta * 2.0e7_dp / (1.8_dp * c) &
        * (-omega_mars) + 0.45_dp / (1.0_dp - 0.3_dp * s)
    call close("signed density/electric force", A1, expected_a1)
    call close("temperature force", A2, -0.3_dp / (1.0_dp - 0.3_dp * s))

    if (trim(mode) /= "transport") then
        ! Nonzero local current-to-torque oracle, independent of orbit root finding.
        q = 1.8_dp
        psi_pr = 2.0e7_dp
        dVds = 1.2e7_dp
        do i = 1, 3
            raw_D = [D11int(0.5_dp * i, 2.0e-5_dp, 3.0e-15_dp), &
                D12int(0.5_dp * i, 2.0e-5_dp, 3.0e-15_dp)]
            call check_current_torque("local signed torque", raw_D, &
                Tphi_int(0.5_dp * i, 2.0e-5_dp, 3.0e-15_dp))
        end do
        print *, "Physical linear electron initialization and signed force/torque: PASS"
        stop
    end if

    config%s = s
    config%qs = -1.0_dp
    config%ms = electron_mass_u
    config%vth = expected_vth
    config%epsmn = 1.0e-3_dp
    config%mph = 3
    config%m0 = 0
    config%comptorque = .true.
    config%supban = .false.
    config%drift_phase = .false.
    config%fow_radial = .false.
    config%magdrift = .true.
    config%magdrift_passing = 1
    config%inp_swi = 8
    config%efac = -1.0_dp
    config%vsteps = 16
    config%log_level = -1
    call neort_init(config, "in_file")
    do i = 1, n
        rotation(i, 1) = plasma(i, 1)
        rotation(i, 2) = omega_mars * R0 / &
            sqrt(2.0_dp * plasma(i, 4) * ev / (electron_mass_u * mu))
    end do
    call prepare_profile_splines(rotation)
    call neort_setup_at_s(0.4_dp)
    call close("physical rotation mapping", Om_tE, -omega_mars)
    scale = pi * vth**3 / (16.0_dp * R0 * iota &
        * (qi * B0 / (mi * c))**2) * (2.0_dp / a * sqrt(s))**2
    total_D = 0.0_dp
    do j = -1, 1
        Dco = 0.0_dp
        Dctr = 0.0_dp
        Dt = 0.0_dp
        Tco = 0.0_dp
        Tctr = 0.0_dp
        Tt = 0.0_dp
        call compute_transport_harmonic(j, Dco, Dctr, Dt, Tco, Tctr, Tt, harmonic)
        call check_current_torque("co-passing", scale * Dco, Tco)
        call check_current_torque("counter-passing", scale * Dctr, Tctr)
        call check_current_torque("trapped", scale * Dt, Tt)
        total_D = total_D + abs([Dco(1), Dctr(1), Dt(1)])
    end do
    ! Low physical electron Mach number need not excite passing resonances
    ! in ell=-1..1. Manufacture an interior ell=0 passing root at a velocity
    ! midpoint: Omega_E = -q*Omega_theta-Omega_d. The signed current/torque
    ! oracle remains independent of this choice and must be nonzero for each
    ! passing direction; the physical-rotation trapped case above is retained.
    call set_to_passing_region(eta_min, eta_max)
    eta = 0.5_dp * (eta_min + eta_max)
    ux_mid = 1.0e-6_dp + 1.5_dp &
        * (config%vmax_over_vth - 1.0e-6_dp) / real(config%vsteps, dp)
    do direction = -1, 1, 2
        sign_vpar = real(direction, dp)
        call Om_th(ux_mid * vth, eta, Omth, dv, deta)
        call Om_tB(ux_mid * vth, eta, OmtB, dv, deta)
        Om_tE = -Omth / iota - OmtB
        print *, "Manufactured passing direction/eta/ux/Omega_E:", &
            direction, eta, ux_mid, Om_tE
        mth = 0
        call resonance_value(ux_mid * vth, eta, residual, residual_derivative)
        print *, "Native residual/eta derivative at predicted root:", &
            residual, residual_derivative
        frequency_scale = mph * max(abs(Om_tE), abs(Omth / iota), abs(OmtB))
        if (abs(residual) > 1.0e-12_dp * frequency_scale) &
            error stop "manufactured passing root prediction failed"
        if (.not. ieee_is_finite(residual)) &
            error stop "nonfinite passing resonance value"
        if (.not. ieee_is_finite(residual_derivative)) &
            error stop "nonfinite passing resonance derivative"
        if (residual_derivative == 0.0_dp) &
            error stop "manufactured passing root is not simple"
        call init_thermodynamic_forces(psi_pr, q)
        Dco = 0.0_dp
        Dctr = 0.0_dp
        Dt = 0.0_dp
        Tco = 0.0_dp
        Tctr = 0.0_dp
        Tt = 0.0_dp
        call compute_transport_harmonic(0, Dco, Dctr, Dt, Tco, Tctr, Tt, harmonic)
        call check_current_torque("manufactured co-passing", scale * Dco, Tco)
        call check_current_torque("manufactured counter-passing", scale * Dctr, Tctr)
        call check_current_torque("manufactured trapped", scale * Dt, Tt)
        if (direction == 1) then
            if (Dco(1) <= 0.0_dp) error stop "co-passing oracle is vacuous"
        else
            if (Dctr(1) <= 0.0_dp) error stop "counter-passing oracle is vacuous"
        end if
        total_D = total_D + abs([Dco(1), Dctr(1), Dt(1)])
    end do
    print *, "Population diffusion magnitudes:", total_D
    if (any(total_D <= 0.0_dp)) error stop "electron population oracle is vacuous"
    print *, "Native linear electron population current-to-torque identity: PASS"

contains

    subroutine check_current_torque(name, current_coefficients, actual_torque)
        character(*), intent(in) :: name
        real(dp), intent(in) :: current_coefficients(2), actual_torque
        real(dp) :: expected_torque

        expected_torque = sign(1.0_dp, psi_pr * q * sign_theta) * ni1 * qi &
            * dVds * abs(psi_pr) / (c * q) &
            * (A1 * current_coefficients(1) + A2 * current_coefficients(2))
        call close(name, actual_torque, expected_torque)
    end subroutine check_current_torque

    subroutine close(name, actual, expected)
        character(*), intent(in) :: name
        real(dp), intent(in) :: actual, expected

        if (.not. ieee_is_finite(actual)) error stop "nonfinite electron observable"
        if (.not. ieee_is_finite(expected)) error stop "nonfinite electron oracle"
        if (abs(actual - expected) > 1.0e-10_dp * max(abs(expected), 1.0e-100_dp)) then
            print *, trim(name), "actual/expected:", actual, expected
            error stop "electron oracle mismatch"
        end if
    end subroutine close

end program test_linear_electron
