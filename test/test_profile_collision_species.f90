program test_profile_collision_species
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use collis_alp, only: efcolf, velrat, enrat, loacol_nbi
    use do_magfie_mod, only: s, sign_theta
    use neort_profiles, only: prepare_plasma_splines, init_plasma_at_s, &
        prepare_profile_splines, init_profile_at_s, init_thermodynamic_forces, &
        profile_species, ns, ts, dnsds, dtsds, vth, dvthds, A1, A2, Om_tE, dOm_tEds
    use util, only: ev, qe, mu, qi, mi, c
    implicit none

    integer, parameter :: n = 6
    real(dp), parameter :: electron_mass = 9.109382e-28_dp, pmass = 1.6726e-24_dp
    real(dp), parameter :: base(5) = [8.0e13_dp, 2.0e13_dp, 8.0e3_dp, &
        7.0e3_dp, 9.0e3_dp]
    real(dp), parameter :: slope(5) = [-1.6e13_dp, -2.0e12_dp, -2.4e3_dp, &
        -1.4e3_dp, -2.25e3_dp]
    real(dp), parameter :: rotation = -1.2e4_dp, radius = 620.0_dp
    real(dp), parameter :: psi_pr = 7.0e7_dp, safety_factor = 1.7_dp
    real(dp) :: plasma(n, 6), masses(2), charges(2)
    character(len=32) :: mode
    integer :: i, species

    if (profile_species /= 1) error stop "first ion must remain the default"
    do i = 1, n
        plasma(i, 1) = real(i - 1, dp) / real(n - 1, dp)
        plasma(i, 2:6) = base + slope * plasma(i, 1)
    end do
    s = 0.43_dp
    masses = [2.013553_dp, 3.015501_dp]
    charges = [1.0_dp, 1.0_dp]
    call get_command_argument(1, mode)
    if (mode == "invalid-density") plasma(:, 2) = 0.0_dp
    if (mode == "invalid-temperature") plasma(:, 4) = -1.0_dp
    if (mode == "nan-density") plasma(:, 2) = ieee_value(0.0_dp, ieee_quiet_nan)
    if (mode == "nan-temperature") plasma(:, 4) = ieee_value(0.0_dp, ieee_quiet_nan)
    if (mode == "infinite-mass") masses(1) = ieee_value(0.0_dp, ieee_positive_inf)
    if (mode == "infinite-charge") charges(1) = ieee_value(0.0_dp, ieee_positive_inf)
    call prepare_plasma_splines(n, masses(1), masses(2), charges(1), charges(2), plasma)
    if (len_trim(mode) > 0) then
        call init_plasma_at_s()
        stop
    end if
    call check_species(1)
    do species = 2, 0, -1
        profile_species = species
        call check_species(species)
    end do
    ! A multiply charged mixture checks charge-weighted neutrality and retains
    ! the configured-helium collision regression.
    masses = [4.002602_dp, 2.013553_dp]
    charges = [2.0_dp, 1.0_dp]
    call prepare_plasma_splines(n, masses(1), masses(2), charges(1), charges(2), plasma)
    do species = 0, 2
        profile_species = species
        call check_species(species)
    end do

contains

    subroutine check_species(selected)
        integer, intent(in) :: selected
        real(dp) :: expected(6), actual(6), old_collision(3, 3)

        if (selected == 0) then
            old_collision(:, 1) = efcolf
            old_collision(:, 2) = velrat
            old_collision(:, 3) = enrat
        end if
        call expected_species(selected, expected)
        call init_plasma_at_s()
        actual = [mi, qi, ns, ts, dnsds, dtsds]
        call assert_close("selected physical profiles", actual, expected)
        call assert_close("thermal speed", [vth, dvthds], &
            [sqrt(2.0_dp * expected(4) * ev / expected(1)), &
            0.5_dp * expected(6) / expected(4) * vth])
        call check_forces(expected)
        if (selected == 0) then
            call assert_close("electron bypass efcolf", efcolf, old_collision(:, 1))
            call assert_close("electron bypass velrat", velrat, old_collision(:, 2))
            call assert_close("electron bypass enrat", enrat, old_collision(:, 3))
        else
            call check_ion_collision(selected)
        end if
    end subroutine check_species

    subroutine expected_species(selected, expected)
        integer, intent(in) :: selected
        real(dp), intent(out) :: expected(6)
        real(dp) :: values(5)

        values = base + slope * s
        if (selected == 0) then
            expected = [electron_mass, -qe, dot_product(charges, values(1:2)), &
                values(5), dot_product(charges, slope(1:2)), slope(5)]
        else
            expected = [masses(selected) * mu, charges(selected) * qe, &
                values(selected), values(selected + 2), slope(selected), &
                slope(selected + 2)]
        end if
    end subroutine expected_species

    subroutine check_forces(expected)
        real(dp), intent(in) :: expected(6)
        real(dp) :: profile(n, 2)

        profile(:, 1) = plasma(:, 1)
        profile(:, 2) = radius * rotation / vth
        call prepare_profile_splines(profile)
        call init_profile_at_s(radius, 1.0_dp, 1.0_dp)
        call assert_close("fixed physical electric rotation", [Om_tE, dOm_tEds], &
            [rotation, 0.5_dp * rotation * expected(6) / expected(4)])
        call check_force_value(expected, rotation)
        Om_tE = 0.0_dp
        call check_force_value(expected, 0.0_dp)
    end subroutine check_forces

    subroutine check_force_value(expected, omega)
        real(dp), intent(in) :: expected(6), omega
        real(dp) :: potential_gradient, expected_a1, expected_a2

        ! Maxwellian force: d(log n)/ds + q*dPhi/ds/T - 3*d(log T)/2ds.
        potential_gradient = -sign_theta * psi_pr * omega / (safety_factor * c)
        expected_a2 = expected(6) / expected(4)
        expected_a1 = expected(5) / expected(3) + expected(2) * &
            potential_gradient / (expected(4) * ev) - 1.5_dp * expected_a2
        call init_thermodynamic_forces(psi_pr, safety_factor)
        call assert_close("signed thermodynamic force", [A1, A2], &
            [expected_a1, expected_a2])
    end subroutine check_force_value

    subroutine check_ion_collision(selected)
        integer, intent(in) :: selected
        real(dp) :: actual_efcolf(3), actual_velrat(3), actual_enrat(3), values(5)
        real(dp) :: amb, zb, ebeam, v0, dchichi, slowrate, dchichi_norm, slowrate_norm

        actual_efcolf = efcolf
        actual_velrat = velrat
        actual_enrat = enrat
        values = base + slope * s
        amb = masses(selected)
        zb = charges(selected)
        v0 = vth
        ebeam = amb * pmass * v0**2 / (2.0_dp * ev)
        call loacol_nbi(amb, masses(1), masses(2), zb, charges(1), charges(2), &
            values(1), values(2), values(3), values(4), values(5), ebeam, v0, &
            dchichi, slowrate, dchichi_norm, slowrate_norm)
        call assert_close("configured ion efcolf", actual_efcolf, efcolf)
        call assert_close("configured ion velrat", actual_velrat, velrat)
        call assert_close("configured ion enrat", actual_enrat, enrat)
    end subroutine check_ion_collision

    subroutine assert_close(name, actual, expected)
        character(*), intent(in) :: name
        real(dp), intent(in) :: actual(:), expected(:)
        real(dp) :: scale(size(expected))

        scale = max(abs(expected), tiny(1.0_dp))
        if (any(abs(actual - expected) > 1.0e-12_dp * scale)) then
            print *, trim(name)
            print *, "actual:  ", actual
            print *, "expected:", expected
            error stop 1
        end if
    end subroutine assert_close

end program test_profile_collision_species
