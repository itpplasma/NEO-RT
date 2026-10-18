program test_linear_carbon_background
    use iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use collis_alp, only: efcolf
    use neort_lib, only: config_t, transport_data_t, neort_init, &
        neort_compute_at_s
    use neort_profiles, only: prepare_plasma_splines, prepare_profile_splines
    implicit none

    integer, parameter :: n = 6
    real(dp), parameter :: surfaces(2) = [0.06_dp, 0.12_dp]
    real(dp) :: plasma(n, 6), rotation(n, 2), collisions(3, 2)
    type(transport_data_t) :: reference, carbon
    integer :: k

    call initialize_library()
    call manufacture_profiles()
    do k = 1, size(surfaces)
        call compute_background(.false., surfaces(k), reference, collisions(:, 1))
        call compute_background(.true., surfaces(k), carbon, collisions(:, 2))
        if (maxval(abs(collisions(:, 1) - collisions(:, 2))) <= &
            1.0e-6_dp * maxval(abs(collisions(:, 1)))) &
            error stop "carbon did not change the native collision background"
        if (sum(abs(reference%summary%Dco)) + sum(abs(reference%summary%Dctr)) &
            + sum(abs(reference%summary%Dt)) <= tiny(1.0_dp)) &
            error stop "transport invariant is vacuous: zero diffusion"
        if (sum(abs([reference%torque%Tco, reference%torque%Tctr, &
            reference%torque%Tt])) <= tiny(1.0_dp)) &
            error stop "transport invariant is vacuous: zero torque"
        call compare_transport(reference, carbon)
        print *, "surface, D-only torque, carbon torque:", surfaces(k), &
            reference%torque%Tt, carbon%torque%Tt
        print *, "D-only collision coefficients:", collisions(:, 1)
        print *, "carbon collision coefficients:", collisions(:, 2)
    end do
    print *, "Linear ion transport is invariant under carbon-only background change"

contains

    subroutine initialize_library()
        type(config_t) :: config

        config%s = 0.06_dp
        config%qs = 1.0_dp
        config%ms = 2.014_dp
        config%vth = 4.0e7_dp
        config%epsmn = 0.001_dp
        config%mph = 3
        config%comptorque = .true.
        config%magdrift = .true.
        config%magdrift_passing = 1
        config%noshear = .true.
        config%inp_swi = 8
        config%vsteps = 32
        config%log_level = -1
        ! nonlin remains false: the production QL path skips collisions.
        call neort_init(config, "in_file")
    end subroutine initialize_library

    subroutine manufacture_profiles()
        integer :: i

        do i = 1, n
            plasma(i, 1) = real(i - 1, dp) / real(n - 1, dp)
            plasma(i, 2) = 8.0e13_dp * (1.0_dp - 0.2_dp * plasma(i, 1))
            plasma(i, 3) = 0.0_dp
            plasma(i, 4) = 8.0e3_dp * (1.0_dp - 0.3_dp * plasma(i, 1))
            plasma(i, 5) = 1.0_dp
            plasma(i, 6) = plasma(i, 4)
            rotation(i, :) = [plasma(i, 1), 0.036_dp]
        end do
        call prepare_profile_splines(rotation)
    end subroutine manufacture_profiles

    subroutine compute_background(add_carbon, surface, result, collision)
        logical, intent(in) :: add_carbon
        real(dp), intent(in) :: surface
        type(transport_data_t), intent(out) :: result
        real(dp), intent(out) :: collision(3)
        real(dp) :: carbon_mass, carbon_charge

        carbon_mass = 2.013562985_dp
        carbon_charge = 1.0_dp
        plasma(:, 3) = 0.0_dp
        plasma(:, 5) = 1.0_dp
        if (add_carbon) then
            carbon_mass = 12.0_dp
            carbon_charge = 6.0_dp
            plasma(:, 3) = plasma(:, 2) / 24.0_dp
            plasma(:, 5) = plasma(:, 4)
        end if
        call prepare_plasma_splines(n, 2.013562985_dp, carbon_mass, &
            1.0_dp, carbon_charge, plasma)
        call neort_compute_at_s(surface, result)
        collision = efcolf
    end subroutine compute_background

    subroutine compare_transport(expected, actual)
        type(transport_data_t), intent(in) :: expected, actual
        integer :: i

        call close(expected%summary%Dco, actual%summary%Dco)
        call close(expected%summary%Dctr, actual%summary%Dctr)
        call close(expected%summary%Dt, actual%summary%Dt)
        call close([expected%torque%Tco, expected%torque%Tctr, expected%torque%Tt], &
            [actual%torque%Tco, actual%torque%Tctr, actual%torque%Tt])
        if (size(expected%harmonics) /= size(actual%harmonics)) &
            error stop "carbon changed harmonic count"
        do i = 1, size(expected%harmonics)
            if (expected%harmonics(i)%mth /= actual%harmonics(i)%mth) &
                error stop "carbon changed harmonic index"
            call close(expected%harmonics(i)%Dresco, actual%harmonics(i)%Dresco)
            call close(expected%harmonics(i)%Dresctr, actual%harmonics(i)%Dresctr)
            call close(expected%harmonics(i)%Drest, actual%harmonics(i)%Drest)
            call close([expected%harmonics(i)%Tresco, expected%harmonics(i)%Tresctr, &
                expected%harmonics(i)%Trest], [actual%harmonics(i)%Tresco, &
                actual%harmonics(i)%Tresctr, actual%harmonics(i)%Trest])
        end do
    end subroutine compare_transport

    subroutine close(expected, actual)
        real(dp), intent(in) :: expected(:), actual(:)

        if (.not. all(ieee_is_finite(expected)) .or. &
            .not. all(ieee_is_finite(actual))) error stop "nonfinite transport"
        if (any(abs(actual - expected) > &
            1.0e-12_dp * max(abs(expected), tiny(1.0_dp)))) &
            error stop "carbon-only change altered linear ion transport"
    end subroutine close

end program test_linear_carbon_background
