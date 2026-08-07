program test_gc_eqdsk_cut_jet_symbolic
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_cut_jet_symbolic, only: evaluate_neort_eqdsk_cut_jet
    use neort_eqdsk_quintic_cell_jet_symbolic, only: &
        evaluate_neort_eqdsk_quintic_cell_jet
    use neort_eqdsk_quintic_profile_jet_symbolic, only: &
        evaluate_neort_eqdsk_quintic_profile_jet
    use util_for_test, only: pass_test
    implicit none

    real(dp) :: cell_coefficient(0:5,0:5), cell_actual(10), cell_expected(10)
    real(dp) :: profile_coefficient(0:5), profile_actual(4)
    real(dp) :: profile_expected(4), cut_actual(7), cut_expected(7)
    real(dp) :: cut_expected_plus(7)
    real(dp) :: cut_unit(7), cut_double(7), cut_phi(7)
    real(dp) :: delta_r, delta_z, profile_delta
    real(dp) :: radius, z_position, b0, f0, psi_sep
    real(dp) :: dot_r, dot_phi, dot_z, field_scale
    integer :: i, j, orientation

    do i = 0, 5
        do j = 0, 5
            cell_coefficient(i,j) = real((i+1)*(j+2), dp)/17.0_dp
        end do
    end do
    delta_r = 0.37_dp
    delta_z = -0.21_dp
    call evaluate_cell_kernel(delta_r, delta_z, cell_coefficient, cell_actual)
    call cell_polynomial_oracle(delta_r, delta_z, cell_coefficient, cell_expected)
    do i = 1, 10
        call require_close('cell jet', cell_actual(i), cell_expected(i))
    end do

    profile_coefficient = [0.75_dp, -1.2_dp, 0.4_dp, 1.1_dp, -0.6_dp, 0.2_dp]
    profile_delta = -0.31_dp
    call evaluate_neort_eqdsk_quintic_profile_jet(profile_delta, &
        profile_coefficient(0), profile_coefficient(1), profile_coefficient(2), &
        profile_coefficient(3), profile_coefficient(4), profile_coefficient(5), &
        2.4_dp, 1.7_dp, profile_actual(1), profile_actual(2), &
        profile_actual(3), profile_actual(4))
    call profile_polynomial_oracle(profile_delta, profile_coefficient, &
        2.4_dp, 1.7_dp, profile_expected)
    do i = 1, 4
        call require_close('profile jet', profile_actual(i), profile_expected(i))
    end do

    radius = 1.8_dp
    z_position = 0.42_dp
    b0 = 2.3_dp
    f0 = 1.1_dp
    psi_sep = 7.0_dp
    dot_r = 0.13_dp
    dot_phi = -0.08_dp
    dot_z = -0.17_dp
    field_scale = 1.7_dp
    orientation = 1
    call evaluate_cut_kernel(radius, field_scale, orientation, radius, z_position, &
        b0, f0, psi_sep, dot_r, dot_phi, dot_z, cut_actual)
    call circular_cut_oracle(radius, z_position, b0, f0, field_scale, &
        orientation, dot_r, dot_z, cut_expected)
    do i = 1, 7
        call require_close('circular Eq13 cut', cut_actual(i), cut_expected(i))
    end do
    call require_close('zero arc-phi derivative', cut_actual(3), 0.0_dp)

    call evaluate_cut_kernel(radius, 1.0_dp, 1, radius, z_position, b0, f0, &
        psi_sep, dot_r, dot_phi, dot_z, cut_unit)
    call evaluate_cut_kernel(radius, 2.0_dp, 1, radius, z_position, b0, f0, &
        psi_sep, dot_r, dot_phi, dot_z, cut_double)
    do i = 1, 7
        call require_close('field scale squared', cut_double(i), &
            4.0_dp*cut_unit(i))
    end do

    call evaluate_cut_kernel(radius, field_scale, -1, radius, z_position, b0, f0, &
        psi_sep, dot_r, dot_phi, dot_z, cut_actual)
    call circular_cut_oracle(radius, z_position, b0, f0, field_scale, -1, &
        dot_r, dot_z, cut_expected)
    do i = 1, 7
        call require_close('orientation reversal', cut_actual(i), cut_expected(i))
    end do
    call require_close('orientation absolute rate', cut_actual(6), &
        abs(cut_unit(5))*field_scale**2)

    call circular_cut_oracle(radius, z_position, b0, f0, field_scale, 1, &
        dot_r, dot_z, cut_expected_plus)
    call evaluate_cut_kernel(radius, field_scale, 1, radius, z_position, b0, f0, &
        psi_sep, dot_r, 0.51_dp, dot_z, cut_phi)
    do i = 1, 7
        call require_close('arc-phi independence', cut_phi(i), &
            cut_expected_plus(i))
    end do
    write (*, '(A)') 'test_gc_eqdsk_cut_jet_symbolic OK'
    call pass_test

contains

    subroutine evaluate_cell_kernel(delta_r_in, delta_z_in, coefficient, result)
        real(dp), intent(in) :: delta_r_in, delta_z_in, coefficient(0:5,0:5)
        real(dp), intent(out) :: result(10)

        call evaluate_neort_eqdsk_quintic_cell_jet(delta_r_in, delta_z_in, &
            coefficient(0,0), coefficient(0,1), coefficient(0,2), &
            coefficient(0,3), coefficient(0,4), coefficient(0,5), &
            coefficient(1,0), coefficient(1,1), coefficient(1,2), &
            coefficient(1,3), coefficient(1,4), coefficient(1,5), &
            coefficient(2,0), coefficient(2,1), coefficient(2,2), &
            coefficient(2,3), coefficient(2,4), coefficient(2,5), &
            coefficient(3,0), coefficient(3,1), coefficient(3,2), &
            coefficient(3,3), coefficient(3,4), coefficient(3,5), &
            coefficient(4,0), coefficient(4,1), coefficient(4,2), &
            coefficient(4,3), coefficient(4,4), coefficient(4,5), &
            coefficient(5,0), coefficient(5,1), coefficient(5,2), &
            coefficient(5,3), coefficient(5,4), coefficient(5,5), &
            result(1), result(2), result(3), result(4), result(5), result(6), &
            result(7), result(8), result(9), result(10))
    end subroutine evaluate_cell_kernel

    subroutine cell_polynomial_oracle(delta_r_in, delta_z_in, coefficient, result)
        real(dp), intent(in) :: delta_r_in, delta_z_in, coefficient(0:5,0:5)
        real(dp), intent(out) :: result(10)
        integer :: i, j

        result = 0.0_dp
        do i = 0, 5
            do j = 0, 5
                result(1) = result(1) + coefficient(i,j)*delta_r_in**i* &
                    delta_z_in**j
                if (i >= 1) result(2) = result(2) + i*coefficient(i,j)* &
                    delta_r_in**(i-1)*delta_z_in**j
                if (j >= 1) result(3) = result(3) + j*coefficient(i,j)* &
                    delta_r_in**i*delta_z_in**(j-1)
                if (i >= 2) result(4) = result(4) + i*(i-1)*coefficient(i,j)* &
                    delta_r_in**(i-2)*delta_z_in**j
                if (i >= 1 .and. j >= 1) result(5) = result(5) + i*j* &
                    coefficient(i,j)*delta_r_in**(i-1)*delta_z_in**(j-1)
                if (j >= 2) result(6) = result(6) + j*(j-1)*coefficient(i,j)* &
                    delta_r_in**i*delta_z_in**(j-2)
                if (i >= 3) result(7) = result(7) + i*(i-1)*(i-2)* &
                    coefficient(i,j)*delta_r_in**(i-3)*delta_z_in**j
                if (i >= 2 .and. j >= 1) result(8) = result(8) + i*(i-1)*j* &
                    coefficient(i,j)*delta_r_in**(i-2)*delta_z_in**(j-1)
                if (i >= 1 .and. j >= 2) result(9) = result(9) + i*j*(j-1)* &
                    coefficient(i,j)*delta_r_in**(i-1)*delta_z_in**(j-2)
                if (j >= 3) result(10) = result(10) + j*(j-1)*(j-2)* &
                    coefficient(i,j)*delta_r_in**i*delta_z_in**(j-3)
            end do
        end do
    end subroutine cell_polynomial_oracle

    subroutine profile_polynomial_oracle(delta_in, coefficient, btf_in, rtf_in, &
            result)
        real(dp), intent(in) :: delta_in, coefficient(0:5), btf_in, rtf_in
        real(dp), intent(out) :: result(4)
        integer :: i

        result = 0.0_dp
        do i = 0, 5
            result(1) = result(1) + coefficient(i)*delta_in**i
            if (i >= 1) result(2) = result(2) + i*coefficient(i)*delta_in**(i-1)
            if (i >= 2) result(3) = result(3) + i*(i-1)*coefficient(i)* &
                delta_in**(i-2)
        end do
        result(4) = btf_in*rtf_in
    end subroutine profile_polynomial_oracle

    subroutine evaluate_cut_kernel(radius_in, field_scale_in, orientation_in, &
            r_position, z_position_in, b0_in, f0_in, psi_sep_in, dot_r_in, &
            dot_phi_in, dot_z_in, result)
        real(dp), intent(in) :: radius_in, field_scale_in, r_position
        integer, intent(in) :: orientation_in
        real(dp), intent(in) :: z_position_in, b0_in, f0_in, psi_sep_in
        real(dp), intent(in) :: dot_r_in, dot_phi_in, dot_z_in
        real(dp), intent(out) :: result(7)
        real(dp) :: psi, psi_r, psi_z, psi_rr, psi_rz, psi_zz
        real(dp) :: psi_rrr, psi_rrz, psi_rzz, psi_zzz

        psi = b0_in*(r_position**2+z_position_in**2)/2.0_dp
        psi_r = b0_in*r_position
        psi_z = b0_in*z_position_in
        psi_rr = b0_in
        psi_rz = 0.0_dp
        psi_zz = b0_in
        psi_rrr = 0.0_dp
        psi_rrz = 0.0_dp
        psi_rzz = 0.0_dp
        psi_zzz = 0.0_dp
        call evaluate_neort_eqdsk_cut_jet(radius_in, field_scale_in, &
            real(orientation_in,dp), psi, psi_r, psi_z, psi_rr, psi_rz, psi_zz, &
            psi_rrr, psi_rrz, psi_rzz, psi_zzz, f0_in, 0.0_dp, psi_sep_in, &
            dot_r_in, dot_phi_in, dot_z_in, result(1), result(2), result(3), &
            result(4), result(5), result(6), result(7))
    end subroutine evaluate_cut_kernel

    subroutine circular_cut_oracle(radius_in, z_in, b0_in, f0_in, scale_in, &
            orientation_in, dot_r_in, dot_z_in, result)
        real(dp), intent(in) :: radius_in, z_in, b0_in, f0_in, scale_in
        integer, intent(in) :: orientation_in
        real(dp), intent(in) :: dot_r_in, dot_z_in
        real(dp), intent(out) :: result(7)
        real(dp) :: s, c, c_r, c_z, c_dot, scale_squared

        s = sqrt(b0_in**2*(radius_in**2+z_in**2)+f0_in**2)
        c = b0_in*z_in*s/radius_in**3
        c_r = b0_in*z_in*(b0_in**2/(s*radius_in**2)-3.0_dp*s/radius_in**4)
        c_z = b0_in*(s**2+b0_in**2*z_in**2)/(s*radius_in**3)
        scale_squared = scale_in**2
        c_dot = scale_squared*real(orientation_in,dp)*(c_r*dot_r_in+ &
            c_z*dot_z_in)
        result(1) = scale_squared*real(orientation_in,dp)*c
        result(2) = scale_squared*real(orientation_in,dp)*c_r
        result(3) = 0.0_dp
        result(4) = scale_squared*real(orientation_in,dp)*c_z
        result(5) = c_dot
        result(6) = abs(c_dot)
        result(7) = c_dot
    end subroutine circular_cut_oracle

    subroutine require_close(label, actual, expected)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, expected
        real(dp), parameter :: tolerance = 2.0e-12_dp

        if (abs(actual-expected) > tolerance*max(1.0_dp, abs(expected))) then
            write (*, '(A,2(1X,ES24.16))') trim(label), actual, expected
            error stop 'EqDSK Fortsym kernel oracle failed'
        end if
    end subroutine require_close

end program test_gc_eqdsk_cut_jet_symbolic
