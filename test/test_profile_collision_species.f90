program test_profile_collision_species
    use iso_fortran_env, only: dp => real64
    use collis_alp, only: efcolf, velrat, enrat, loacol_nbi
    use do_magfie_mod, only: s
    use neort_profiles, only: prepare_plasma_splines, init_plasma_at_s, vth
    use util, only: ev
    implicit none

    integer, parameter :: n = 6
    real(dp), parameter :: am1 = 4.002602_dp, am2 = 2.013553_dp
    real(dp), parameter :: z1 = 2.0_dp, z2 = 1.0_dp
    real(dp), parameter :: pmass = 1.6726e-24_dp
    real(dp) :: plasma(n, 6), actual_efcolf(3), actual_velrat(3), actual_enrat(3)
    real(dp) :: amb, zb, ebeam, v0, dchichi, slowrate, dchichi_norm, slowrate_norm
    integer :: i

    do i = 1, n
        plasma(i, 1) = real(i - 1, dp) / real(n - 1, dp)
        plasma(i, 2) = 8.0e13_dp * (1.0_dp - 0.2_dp * plasma(i, 1))
        plasma(i, 3) = 2.0e13_dp * (1.0_dp - 0.1_dp * plasma(i, 1))
        plasma(i, 4) = 8.0e3_dp * (1.0_dp - 0.3_dp * plasma(i, 1))
        plasma(i, 5) = 7.0e3_dp * (1.0_dp - 0.2_dp * plasma(i, 1))
        plasma(i, 6) = 9.0e3_dp * (1.0_dp - 0.25_dp * plasma(i, 1))
    end do

    call prepare_plasma_splines(n, am1, am2, z1, z2, plasma)
    s = 0.43_dp
    call init_plasma_at_s()
    actual_efcolf = efcolf
    actual_velrat = velrat
    actual_enrat = enrat

    ! Independently initialize the legacy collision kernel with the configured
    ! first species. The public profile initializer must produce this contract.
    amb = am1
    zb = z1
    v0 = vth
    ebeam = amb * pmass * v0**2 / (2.0_dp * ev)
    call loacol_nbi(amb, am1, am2, zb, z1, z2, plasma_value(2), plasma_value(3), &
        plasma_value(4), plasma_value(5), plasma_value(6), ebeam, v0, &
        dchichi, slowrate, dchichi_norm, slowrate_norm)

    call assert_close("efcolf", actual_efcolf, efcolf)
    call assert_close("velrat", actual_velrat, velrat)
    call assert_close("enrat", actual_enrat, enrat)

contains

    function plasma_value(column) result(value)
        integer, intent(in) :: column
        real(dp) :: value
        real(dp) :: fraction

        fraction = s
        select case (column)
        case (2)
            value = 8.0e13_dp * (1.0_dp - 0.2_dp * fraction)
        case (3)
            value = 2.0e13_dp * (1.0_dp - 0.1_dp * fraction)
        case (4)
            value = 8.0e3_dp * (1.0_dp - 0.3_dp * fraction)
        case (5)
            value = 7.0e3_dp * (1.0_dp - 0.2_dp * fraction)
        case (6)
            value = 9.0e3_dp * (1.0_dp - 0.25_dp * fraction)
        case default
            error stop "invalid plasma column"
        end select
    end function plasma_value

    subroutine assert_close(name, actual, expected)
        character(*), intent(in) :: name
        real(dp), intent(in) :: actual(:), expected(:)
        real(dp) :: scale

        scale = max(maxval(abs(expected)), tiny(1.0_dp))
        if (maxval(abs(actual - expected)) > 1.0e-12_dp * scale) then
            print *, trim(name), " does not use the configured first species"
            print *, "actual:  ", actual
            print *, "expected:", expected
            error stop 1
        end if
    end subroutine assert_close

end program test_profile_collision_species
