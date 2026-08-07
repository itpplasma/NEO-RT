program test_gc_eqdsk_cut_chart_symbolic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_cut_numerator_symbolic, only: &
        evaluate_neort_eqdsk_cut_numerator
    use neort_eqdsk_cut_r_chart_symbolic, only: &
        evaluate_neort_eqdsk_cut_r_chart
    use neort_eqdsk_cut_z_chart_symbolic, only: &
        evaluate_neort_eqdsk_cut_z_chart
    use util_for_test, only: pass_test
    implicit none

    real(dp) :: radius, z_position, b0, f0, psi_sep
    real(dp) :: numerator(3), numerator_repeat(3), expected(3)
    real(dp) :: r_chart(2), z_chart(2)
    real(dp) :: expected_r_chart(2), expected_z_chart(2)
    real(dp) :: oriented_plus, oriented_minus

    radius = 1.8_dp
    z_position = 0.42_dp
    b0 = 2.3_dp
    f0 = 1.1_dp
    psi_sep = 7.0_dp

    call evaluate_numerator(radius, z_position, b0, f0, psi_sep, numerator)
    call circular_n_oracle(radius, z_position, b0, f0, expected)
    call require_close('circular N', numerator(1), expected(1))
    call require_close('circular N_R', numerator(2), expected(2))
    call require_close('circular N_Z', numerator(3), expected(3))

    if (abs(numerator(2)) <= 1.0e-12_dp .or. &
            abs(numerator(3)) <= 1.0e-12_dp) then
        error stop 'chart oracle point is singular'
    end if
    call evaluate_neort_eqdsk_cut_r_chart(numerator(2), numerator(3), &
        r_chart(1), r_chart(2))
    call evaluate_neort_eqdsk_cut_z_chart(numerator(2), numerator(3), &
        z_chart(1), z_chart(2))
    expected_r_chart = [-expected(2)/expected(3), &
        sqrt(1.0_dp+(expected(2)/expected(3))**2)]
    expected_z_chart = [-expected(3)/expected(2), &
        sqrt(1.0_dp+(expected(3)/expected(2))**2)]
    call require_close('R-chart slope', r_chart(1), expected_r_chart(1))
    call require_close('R-chart arc derivative', r_chart(2), expected_r_chart(2))
    call require_close('Z-chart slope', z_chart(1), expected_z_chart(1))
    call require_close('Z-chart arc derivative', z_chart(2), expected_z_chart(2))
    call require_close('R-chart tangent orthogonality', &
        expected(2)+expected(3)*r_chart(1), 0.0_dp)
    call require_close('Z-chart tangent orthogonality', &
        expected(2)*z_chart(1)+expected(3), 0.0_dp)

    call evaluate_numerator(radius, z_position, b0, f0, psi_sep, numerator_repeat)
    call require_close('N is orientation independent', numerator_repeat(1), &
        numerator(1))
    oriented_plus = numerator(1)
    oriented_minus = -numerator(1)
    call require_close('positive orientation numerator', oriented_plus, &
        expected(1))
    call require_close('negative orientation numerator sign', oriented_minus, &
        -expected(1))

    write (*, '(A)') 'test_gc_eqdsk_cut_chart_symbolic OK'
    call pass_test

contains

    subroutine evaluate_numerator(radius_in, z_in, b0_in, f0_in, psi_sep_in, &
            result)
        real(dp), intent(in) :: radius_in, z_in, b0_in, f0_in, psi_sep_in
        real(dp), intent(out) :: result(3)
        real(dp) :: psi_r, psi_z, psi_rr, psi_rz, psi_zz

        psi_r = b0_in*radius_in
        psi_z = b0_in*z_in
        psi_rr = b0_in
        psi_rz = 0.0_dp
        psi_zz = b0_in
        call evaluate_neort_eqdsk_cut_numerator(radius_in, psi_r, psi_z, &
            psi_rr, psi_rz, psi_zz, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, f0_in, &
            0.0_dp, psi_sep_in, result(1), result(2), result(3))
    end subroutine evaluate_numerator

    subroutine circular_n_oracle(radius_in, z_in, b0_in, f0_in, result)
        real(dp), intent(in) :: radius_in, z_in, b0_in, f0_in
        real(dp), intent(out) :: result(3)
        real(dp) :: g

        g = b0_in**2*(radius_in**2+z_in**2)+f0_in**2
        result(1) = b0_in*z_in*g
        result(2) = 2.0_dp*b0_in**3*radius_in*z_in
        result(3) = b0_in*g+2.0_dp*b0_in**3*z_in**2
    end subroutine circular_n_oracle

    subroutine require_close(label, got, reference)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, reference
        real(dp), parameter :: tolerance = 2.0e-11_dp

        if (abs(got-reference) > tolerance*max(1.0_dp, abs(reference))) then
            write (*, '(A,1X,ES24.16,1X,ES24.16)') trim(label), got, reference
            error stop 'independent chart oracle mismatch'
        end if
    end subroutine require_close

end program test_gc_eqdsk_cut_chart_symbolic
