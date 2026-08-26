program test_electromagnetic_torque
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use electromagnetic_torque, only: cross_contraction, staggered_cross_contraction

    implicit none

    call test_hand_derived_mode_contraction
    call test_explicit_full_half_stagger
    call test_common_complex_phase

contains

    subroutine test_hand_derived_mode_contraction
        complex(dp) :: j1(2), j2(2), b1(2), b2(2)
        real(dp) :: actual

        j1 = [cmplx(1.0_dp, 2.0_dp, dp), cmplx(-0.5_dp, 1.0_dp, dp)]
        j2 = [cmplx(3.0_dp, -1.0_dp, dp), cmplx(2.0_dp, 0.25_dp, dp)]
        b1 = [cmplx(-2.0_dp, 1.0_dp, dp), cmplx(1.5_dp, -3.0_dp, dp)]
        b2 = [cmplx(0.5_dp, 4.0_dp, dp), cmplx(-1.0_dp, 2.0_dp, dp)]

        actual = sum(cross_contraction(j1, b2, j2, b1))
        call assert_close(actual, 7.875_dp, "hand-derived two-mode contraction")
    end subroutine test_hand_derived_mode_contraction

    subroutine test_explicit_full_half_stagger
        complex(dp) :: j1_half(1, 1), b2_half(1, 1)
        complex(dp) :: j2_full(2, 1), b1_full(2, 1)
        real(dp), allocatable :: actual(:)

        j1_half = cmplx(1.0_dp, 1.0_dp, dp)
        b2_half = cmplx(3.0_dp, 3.0_dp, dp)
        j2_full(:, 1) = [cmplx(1.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 2.0_dp, dp)]
        b1_full(:, 1) = [cmplx(4.0_dp, 0.0_dp, dp), cmplx(0.0_dp, 3.0_dp, dp)]

        actual = staggered_cross_contraction(j1_half, b2_half, j2_full, b1_full)
        call assert_close(actual(1), 0.5_dp, "explicit full/half-grid contraction")
    end subroutine test_explicit_full_half_stagger

    subroutine test_common_complex_phase
        complex(dp), parameter :: phase = cmplx(0.6_dp, 0.8_dp, dp)
        complex(dp) :: j1(1), j2(1), b1(1), b2(1)
        real(dp) :: reference, rotated

        j1 = cmplx(1.0_dp, -2.0_dp, dp)
        j2 = cmplx(0.5_dp, 3.0_dp, dp)
        b1 = cmplx(-4.0_dp, 0.25_dp, dp)
        b2 = cmplx(2.0_dp, 1.5_dp, dp)
        reference = sum(cross_contraction(j1, b2, j2, b1))
        rotated = sum(cross_contraction(phase*j1, phase*b2, phase*j2, phase*b1))

        call assert_close(rotated, reference, "common complex phase")
    end subroutine test_common_complex_phase

    subroutine assert_close(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label
        real(dp), parameter :: tolerance = 1.0e-13_dp

        if (abs(actual - expected) > tolerance) then
            print *, trim(label), ": expected", expected, "got", actual
            error stop
        end if
    end subroutine assert_close

end program test_electromagnetic_torque
