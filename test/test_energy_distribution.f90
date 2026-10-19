program test_energy_distribution
    use iso_fortran_env, only: dp => real64
    use kinetic_energy_selection, only: maxwellian_speed_density, &
        speed_from_energy_ratio, energy_sample_count, energy_sample_speed, &
        energy_sample_weight
    use neort_transport, only: D11int, D11mono, D12int, Tphi_int
    use do_magfie_mod, only: q, psi_pr, sign_theta
    use neort_profiles, only: ni1, A1, A2
    use neort_magfie, only: dVds
    use driftorbit, only: mph, vth
    use util, only: c, qi
    use neort_config, only: config_t, set_config
    use neort, only: resolved_monoenergetic_x => monoenergetic_x

    implicit none

    integer, parameter :: steps = 200000
    real(dp), parameter :: upper_speed = 8.0_dp
    real(dp), parameter :: tolerance = 2.0e-9_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: expected, hamiltonian_squared, integral, speed, step, bounce_time
    integer :: index
    type(config_t) :: config

    step = upper_speed/real(steps, dp)
    integral = 0.0_dp
    do index = 1, steps
        speed = (real(index, dp) - 0.5_dp)*step
        integral = integral + step*maxwellian_speed_density(speed)
    end do

    if (abs(integral - 1.0_dp) > tolerance) &
        error stop "Maxwellian speed density is not normalized"
    if (abs(speed_from_energy_ratio(2.25_dp) - 1.5_dp) > epsilon(1.0_dp)) &
        error stop "x=E/T must map to u=sqrt(x)"

    mph = 3
    vth = 2.0_dp
    q = 1.25_dp
    qi = 4.0_dp
    psi_pr = -5.0_dp
    dVds = 7.0_dp
    speed = 1.5_dp
    bounce_time = 11.0_dp
    hamiltonian_squared = 13.0_dp
    expected = acos(-1.0_dp)**2/4.0_dp*mph**2*c**2*q*vth &
        /(qi**2*dVds*abs(psi_pr))*speed*bounce_time*hamiltonian_squared
    if (abs(D11mono(speed, bounce_time, hamiltonian_squared)/expected - 1.0_dp) &
        > 20.0_dp*epsilon(1.0_dp)) error stop "monoenergetic D11 kernel mismatch"
    ! Independent oracle: the literal thermal kernels of main before the
    ! monoenergetic factorization.
    expected = pi**1.5_dp*mph**2*c**2*q*vth/(qi**2*dVds*abs(psi_pr)) &
        *speed**3*exp(-speed**2)*bounce_time*hamiltonian_squared
    if (abs(D11int(speed, bounce_time, hamiltonian_squared)/expected - 1.0_dp) &
        > 20.0_dp*epsilon(1.0_dp)) error stop "thermal D11 kernel changed"
    if (abs(D12int(speed, bounce_time, hamiltonian_squared) &
        /(expected*speed**2) - 1.0_dp) > 20.0_dp*epsilon(1.0_dp)) &
        error stop "thermal D12 kernel changed"
    ni1 = 3.0e13_dp
    A1 = 0.7_dp
    A2 = -0.2_dp
    expected = sign(1.0_dp, psi_pr*q*sign_theta)*pi**1.5_dp*mph**2*ni1*c*vth/qi &
        *speed**3*exp(-speed**2)*bounce_time*hamiltonian_squared &
        *(A1 + A2*speed**2)
    if (abs(Tphi_int(speed, bounce_time, hamiltonian_squared)/expected - 1.0_dp) &
        > 20.0_dp*epsilon(1.0_dp)) error stop "thermal torque kernel changed"

    if (config%monoenergetic_x >= 0.0_dp) &
        error stop "thermal integration must remain the default"
    config%monoenergetic_x = 2.25_dp
    call set_config(config)
    if (resolved_monoenergetic_x /= config%monoenergetic_x) &
        error stop "monoenergetic energy was not applied"
    if (energy_sample_count(-1.0_dp, 8) /= 8) &
        error stop "thermal sample count changed"
    if (energy_sample_count(2.25_dp, 8) /= 1) &
        error stop "monoenergetic mode must use one sample"
    if (energy_sample_speed(2.25_dp, 1, 0.0_dp, 4.0_dp, 8) /= 1.5_dp) &
        error stop "monoenergetic sample is not at the requested energy"
    if (energy_sample_weight(2.25_dp, 1.5_dp, 0.5_dp) /= 1.0_dp) &
        error stop "delta-normalized sample weight must be one"
    expected = 0.5_dp*maxwellian_speed_density(1.5_dp)
    if (energy_sample_weight(-1.0_dp, 1.5_dp, 0.5_dp) /= expected) &
        error stop "thermal midpoint weight mismatch"
end program test_energy_distribution
