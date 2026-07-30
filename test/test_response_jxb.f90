program test_response_jxb
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use response_jxb, only: cylindrical_toroidal_torque, mars_surface_torque

    implicit none

    call test_cylindrical_cross_product
    call test_mars_harmonic_sum
    call test_common_phase_invariance

contains

    subroutine test_cylindrical_cross_product
        complex(dp) :: jr(2), jz(2), br(2), bz(2)
        real(dp) :: radius(2), weight(2), actual, expected

        radius = [2.0_dp, 3.0_dp]
        weight = [0.25_dp, 0.75_dp]
        jr = [cmplx(1.0_dp, 2.0_dp, dp), cmplx(-2.0_dp, 1.0_dp, dp)]
        jz = [cmplx(3.0_dp, -1.0_dp, dp), cmplx(0.5_dp, 4.0_dp, dp)]
        br = [cmplx(-1.0_dp, 3.0_dp, dp), cmplx(2.0_dp, -2.0_dp, dp)]
        bz = [cmplx(2.0_dp, 1.0_dp, dp), cmplx(-3.0_dp, 0.5_dp, dp)]
        expected = 0.5_dp*sum(weight*radius*real(jz*conjg(br) - jr*conjg(bz), dp))

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
        expected = 0.5_dp*sum(real(conjg(j1)*b2 - conjg(j2)*b1, dp))

        actual = mars_surface_torque(j1, j2, b1, b2)
        call assert_close(actual, expected, "MARS Eq. (8) harmonic oracle")
    end subroutine test_mars_harmonic_sum

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
