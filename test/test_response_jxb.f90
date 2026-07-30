program test_response_jxb
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use response_jxb, only: cylindrical_toroidal_torque, integrate_mars_profile, &
        mars_half_mesh_torque, mars_surface_torque

    implicit none

    call test_cylindrical_cross_product
    call test_mars_harmonic_sum
    call test_mars_half_mesh_stagger
    call test_common_phase_invariance
    call test_linear_profile_integral

contains

    subroutine test_cylindrical_cross_product
        complex(dp) :: jr(2), jz(2), br(2), bz(2)
        complex(dp) :: phase
        real(dp), parameter :: pi = acos(-1.0_dp)
        integer, parameter :: number_of_phases = 128
        real(dp) :: radius(2), weight(2), actual, expected, angle
        integer :: index

        radius = [2.0_dp, 3.0_dp]
        weight = [0.25_dp, 0.75_dp]
        jr = [cmplx(1.0_dp, 2.0_dp, dp), cmplx(-2.0_dp, 1.0_dp, dp)]
        jz = [cmplx(3.0_dp, -1.0_dp, dp), cmplx(0.5_dp, 4.0_dp, dp)]
        br = [cmplx(-1.0_dp, 3.0_dp, dp), cmplx(2.0_dp, -2.0_dp, dp)]
        bz = [cmplx(2.0_dp, 1.0_dp, dp), cmplx(-3.0_dp, 0.5_dp, dp)]
        expected = 0.0_dp
        do index = 0, number_of_phases - 1
            angle = 2.0_dp*pi*real(index, dp)/real(number_of_phases, dp)
            phase = cmplx(cos(angle), sin(angle), dp)
            expected = expected + sum(weight*radius*( &
                real(jz*phase, dp)*real(br*phase, dp) - &
                real(jr*phase, dp)*real(bz*phase, dp)))
        end do
        expected = expected/real(number_of_phases, dp)

        actual = cylindrical_toroidal_torque(radius, weight, jr, jz, br, bz)
        call assert_close(actual, expected, "cylindrical cross-product oracle")
    end subroutine test_cylindrical_cross_product

    subroutine test_mars_harmonic_sum
        complex(dp) :: j1(2), j2(2), b1(2), b2(2)
        real(dp) :: actual, expected

        j1 = [cmplx(1.0_dp, 2.0_dp, dp), cmplx(-0.5_dp, 1.0_dp, dp)]
        j2 = [cmplx(3.0_dp, -1.0_dp, dp), cmplx(2.0_dp, 0.25_dp, dp)]
        b1 = [cmplx(-2.0_dp, 1.0_dp, dp), cmplx(1.5_dp, -3.0_dp, dp)]
        b2 = [cmplx(0.5_dp, 4.0_dp, dp), cmplx(-1.0_dp, 2.0_dp, dp)]
        expected = 7.875_dp

        actual = mars_surface_torque(j1, j2, b1, b2)
        call assert_close(actual, expected, "MARS Eq. (8) harmonic oracle")
    end subroutine test_mars_harmonic_sum

    subroutine test_mars_half_mesh_stagger
        complex(dp) :: j1_half(1), b2_half(1)
        complex(dp) :: j2_left(1), b1_left(1), j2_right(1), b1_right(1)
        real(dp) :: actual

        j1_half = cmplx(1.0_dp, 1.0_dp, dp)
        b2_half = cmplx(3.0_dp, 3.0_dp, dp)
        j2_left = cmplx(1.0_dp, 0.0_dp, dp)
        b1_left = cmplx(4.0_dp, 0.0_dp, dp)
        j2_right = cmplx(0.0_dp, 2.0_dp, dp)
        b1_right = cmplx(0.0_dp, 3.0_dp, dp)

        actual = mars_half_mesh_torque( &
            j1_half, b2_half, j2_left, b1_left, j2_right, b1_right)
        call assert_close(actual, 0.5_dp, "MARS KCTORQ=2 radial staggering")
    end subroutine test_mars_half_mesh_stagger

    subroutine test_common_phase_invariance
        complex(dp), parameter :: phase = cmplx(0.6_dp, 0.8_dp, dp)
        complex(dp) :: j1(1), j2(1), b1(1), b2(1)
        real(dp) :: reference, rotated

        j1 = cmplx(1.0_dp, -2.0_dp, dp)
        j2 = cmplx(0.5_dp, 3.0_dp, dp)
        b1 = cmplx(-4.0_dp, 0.25_dp, dp)
        b2 = cmplx(2.0_dp, 1.5_dp, dp)
        reference = mars_surface_torque(j1, j2, b1, b2)
        rotated = mars_surface_torque(phase*j1, phase*j2, phase*b1, phase*b2)

        call assert_close(rotated, reference, "common Fourier phase invariance")
    end subroutine test_common_phase_invariance

    subroutine test_linear_profile_integral
        real(dp), parameter :: pi = acos(-1.0_dp)
        real(dp), parameter :: scale = 2.5_dp
        real(dp) :: edges(4), midpoint(3), density(3), actual, expected

        edges = [0.1_dp, 0.2_dp, 0.6_dp, 0.9_dp]
        midpoint = 0.5_dp*(edges(:3) + edges(2:))
        density = 3.0_dp - 2.0_dp*midpoint
        expected = 4.0_dp*pi**2*scale*( &
            3.0_dp*(edges(4) - edges(1)) - edges(4)**2 + edges(1)**2)

        actual = integrate_mars_profile(edges, density, scale)
        call assert_close(actual, expected, "linear profile analytic integral")
    end subroutine test_linear_profile_integral

    subroutine assert_close(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label
        real(dp), parameter :: tolerance = 1.0e-13_dp

        if (abs(actual - expected) > tolerance) then
            print *, trim(label), ": expected", expected, "got", actual
            error stop
        end if
    end subroutine assert_close

end program test_response_jxb
