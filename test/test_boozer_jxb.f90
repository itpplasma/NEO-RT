program test_boozer_jxb
    ! Independent Cartesian oracle.  An oblique coordinate basis is used to
    ! convert chosen physical complex J and B amplitudes into the coordinate
    ! inputs.  The expected result is then obtained by explicitly sampling the
    ! real fields over their common phase and evaluating
    !
    !   Jcal e_phi . (J_real cross B_real)
    !
    ! in Cartesian components.  This does not reuse the coordinate contraction
    ! or metric-cofactor formula implemented by the generated kernel.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_jxb, only: boozer_local_toroidal_torque
    implicit none

    integer, parameter :: phase_count = 512
    real(dp), parameter :: tolerance = 3.0e-13_dp
    real(dp), parameter :: pi = acos(-1.0_dp)
    real(dp) :: basis(3, 3), metric(3, 3), jacobian
    real(dp) :: expected, phase
    complex(dp) :: physical_b(3), physical_j(3)
    complex(dp) :: b_covariant(3), j_contravariant(3)
    complex(dp) :: weighted_current(3), weighted_magnetic(3)
    real(dp) :: actual
    real(dp) :: real_b(3), real_j(3)
    integer :: sample

    basis(:, 1) = [1.2_dp, -0.3_dp, 0.4_dp]
    basis(:, 2) = [0.2_dp, 2.1_dp, -0.5_dp]
    basis(:, 3) = [-0.4_dp, 0.6_dp, 1.7_dp]
    metric = matmul(transpose(basis), basis)
    jacobian = determinant3(basis)
    physical_b = [ &
        cmplx(0.7_dp, -1.1_dp, dp), &
        cmplx(-0.2_dp, 0.9_dp, dp), &
        cmplx(1.4_dp, 0.3_dp, dp)]
    physical_j = [ &
        cmplx(-0.8_dp, 0.4_dp, dp), &
        cmplx(1.1_dp, -0.6_dp, dp), &
        cmplx(0.5_dp, 1.3_dp, dp)]
    b_covariant = matmul(transpose(basis), physical_b)
    j_contravariant = solve3(basis, physical_j)
    weighted_current = jacobian*j_contravariant

    call boozer_local_toroidal_torque( &
        metric, jacobian, b_covariant, weighted_current, &
        weighted_magnetic, actual)

    expected = 0.0_dp
    do sample = 0, phase_count - 1
        phase = 2.0_dp*pi*real(sample, dp)/real(phase_count, dp)
        real_b = real(physical_b*cmplx(cos(phase), sin(phase), dp), dp)
        real_j = real(physical_j*cmplx(cos(phase), sin(phase), dp), dp)
        expected = expected + jacobian*dot_product( &
            basis(:, 2), cross(real_j, real_b))
    end do
    expected = expected/real(phase_count, dp)

    call assert_close(actual, expected, "Cartesian phase average")
    j_contravariant = solve3(basis, physical_b)
    call assert_close( &
        abs(weighted_magnetic(1) - jacobian*j_contravariant(1)), 0.0_dp, &
        "radial metric raise")
    call assert_close( &
        abs(weighted_magnetic(3) - jacobian*j_contravariant(3)), 0.0_dp, &
        "poloidal metric raise")
    print *, "test_boozer_jxb: all checks passed"

contains

    pure function determinant3(matrix) result(value)
        real(dp), intent(in) :: matrix(3, 3)
        real(dp) :: value

        value = matrix(1, 1)*( &
            matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)) - &
            matrix(1, 2)*( &
            matrix(2, 1)*matrix(3, 3) - matrix(2, 3)*matrix(3, 1)) + &
            matrix(1, 3)*( &
            matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1))
    end function determinant3

    pure function solve3(matrix, right_hand_side) result(solution)
        real(dp), intent(in) :: matrix(3, 3)
        complex(dp), intent(in) :: right_hand_side(3)
        complex(dp) :: solution(3)
        real(dp) :: inverse(3, 3), determinant

        determinant = determinant3(matrix)
        inverse(1, 1) = matrix(2, 2)*matrix(3, 3) - matrix(2, 3)*matrix(3, 2)
        inverse(1, 2) = matrix(1, 3)*matrix(3, 2) - matrix(1, 2)*matrix(3, 3)
        inverse(1, 3) = matrix(1, 2)*matrix(2, 3) - matrix(1, 3)*matrix(2, 2)
        inverse(2, 1) = matrix(2, 3)*matrix(3, 1) - matrix(2, 1)*matrix(3, 3)
        inverse(2, 2) = matrix(1, 1)*matrix(3, 3) - matrix(1, 3)*matrix(3, 1)
        inverse(2, 3) = matrix(1, 3)*matrix(2, 1) - matrix(1, 1)*matrix(2, 3)
        inverse(3, 1) = matrix(2, 1)*matrix(3, 2) - matrix(2, 2)*matrix(3, 1)
        inverse(3, 2) = matrix(1, 2)*matrix(3, 1) - matrix(1, 1)*matrix(3, 2)
        inverse(3, 3) = matrix(1, 1)*matrix(2, 2) - matrix(1, 2)*matrix(2, 1)
        solution = matmul(inverse/determinant, right_hand_side)
    end function solve3

    pure function cross(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3)

        product = [ &
            left(2)*right(3) - left(3)*right(2), &
            left(3)*right(1) - left(1)*right(3), &
            left(1)*right(2) - left(2)*right(1)]
    end function cross

    subroutine assert_close(actual_value, expected_value, label)
        real(dp), intent(in) :: actual_value, expected_value
        character(*), intent(in) :: label

        if (abs(actual_value - expected_value) > tolerance) then
            print *, trim(label), ": expected", expected_value, "got", actual_value
            error stop
        end if
    end subroutine assert_close

end program test_boozer_jxb
