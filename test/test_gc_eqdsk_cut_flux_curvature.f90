program test_gc_eqdsk_cut_flux_curvature
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_quintic_cell_fourth_jet_symbolic, only: &
        evaluate_neort_eqdsk_quintic_cell_fourth_jet
    use neort_eqdsk_cut_numerator_hessian_symbolic, only: &
        evaluate_neort_eqdsk_cut_numerator_hessian
    use neort_eqdsk_cut_r_flux_curvature_symbolic, only: &
        evaluate_neort_eqdsk_cut_r_flux_curvature
    implicit none

    real(dp), parameter :: local_r = 0.23_dp
    real(dp), parameter :: local_z = -0.17_dp
    real(dp) :: coefficient(0:5,0:5), jet(15), fourth_jet(5)
    integer :: i, j

    do i = 0, 5
        do j = 0, 5
            coefficient(i,j) = real((-1)**(i+2*j)*(i+1)*(j+2), dp) &
                /real(80*(i+j+1), dp)
        end do
    end do
    call independent_cell_jet(local_r, local_z, coefficient, jet)
    call generated_fourth_jet(local_r, local_z, coefficient, fourth_jet)
    do i = 1, 5
        call require_close(fourth_jet(i), jet(i+10), 2.0e-13_dp, &
            'quintic fourth derivative')
    end do

    call check_profile_second_derivative_oracle
    call check_independent_graph_oracle
    write (*, '(a)') 'test_gc_eqdsk_cut_flux_curvature OK'

contains

    subroutine generated_fourth_jet(r, z, c, value)
        real(dp), intent(in) :: r, z, c(0:5,0:5)
        real(dp), intent(out) :: value(5)

        call evaluate_neort_eqdsk_quintic_cell_fourth_jet(r, z, &
            c(0,0), c(0,1), c(0,2), c(0,3), c(0,4), c(0,5), &
            c(1,0), c(1,1), c(1,2), c(1,3), c(1,4), c(1,5), &
            c(2,0), c(2,1), c(2,2), c(2,3), c(2,4), c(2,5), &
            c(3,0), c(3,1), c(3,2), c(3,3), c(3,4), c(3,5), &
            c(4,0), c(4,1), c(4,2), c(4,3), c(4,4), c(4,5), &
            c(5,0), c(5,1), c(5,2), c(5,3), c(5,4), c(5,5), &
            value(1), value(2), value(3), value(4), value(5))
    end subroutine generated_fourth_jet

    pure subroutine independent_cell_jet(r, z, c, value)
        real(dp), intent(in) :: r, z, c(0:5,0:5)
        real(dp), intent(out) :: value(15)
        integer :: derivative_r(15), derivative_z(15), k, p, q

        derivative_r = [0, 1, 0, 2, 1, 0, 3, 2, 1, 0, 4, 3, 2, 1, 0]
        derivative_z = [0, 0, 1, 0, 1, 2, 0, 1, 2, 3, 0, 1, 2, 3, 4]
        value = 0.0_dp
        do k = 1, 15
            do p = derivative_r(k), 5
                do q = derivative_z(k), 5
                    value(k) = value(k)+c(p,q) &
                        *falling_factorial(p, derivative_r(k)) &
                        *falling_factorial(q, derivative_z(k)) &
                        *r**(p-derivative_r(k))*z**(q-derivative_z(k))
                end do
            end do
        end do
    end subroutine independent_cell_jet

    subroutine check_profile_second_derivative_oracle
        real(dp) :: n_rr, n_rz, n_zz
        real(dp) :: slope, second, flux_first, flux_second

        ! Independent exact oracle: x=R-2,
        ! psi_hat=x+Z+3*x**2/4 at x=Z=0 and psi_sep=1.
        ! Direct differentiation of N=R*K+psi_Z*(|grad psi|**2+F**2)
        ! gives the constants below.  The paired profiles have identical F
        ! and F' at the point, so this also detects a missing F'' term.
        call evaluate_neort_eqdsk_cut_numerator_hessian(2.0_dp, &
            1.0_dp, 1.0_dp, 1.5_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp, 2.0_dp, 1.0_dp, n_rr, n_rz, n_zz)
        call require_close(n_rr, 9.0_dp, 2.0e-14_dp, &
            'quadratic-profile N_RR')
        call require_close(n_rz, 6.0_dp, 2.0e-14_dp, &
            'quadratic-profile N_RZ')
        call require_close(n_zz, 6.0_dp, 2.0e-14_dp, &
            'quadratic-profile N_ZZ')
        call evaluate_neort_eqdsk_cut_r_flux_curvature(-1.0_dp, 2.0_dp, &
            n_rr, n_rz, n_zz, 1.0_dp, 1.0_dp, 1.5_dp, 0.0_dp, 0.0_dp, &
            1.0_dp, slope, second, flux_first, flux_second)
        call require_close(slope, 0.5_dp, 2.0e-14_dp, &
            'quadratic-profile graph slope')
        call require_close(second, -33.0_dp/4.0_dp, 2.0e-14_dp, &
            'quadratic-profile graph curvature')
        call require_close(flux_first, 1.5_dp, 2.0e-14_dp, &
            'quadratic-profile flux derivative')
        call require_close(flux_second, -27.0_dp/4.0_dp, 2.0e-14_dp, &
            'quadratic-profile flux curvature')

        call evaluate_neort_eqdsk_cut_numerator_hessian(2.0_dp, &
            1.0_dp, 1.0_dp, 1.5_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            1.0_dp, 1.0_dp, 0.0_dp, 1.0_dp, n_rr, n_rz, n_zz)
        call require_close(n_rr, 5.0_dp, 2.0e-14_dp, &
            'affine-profile N_RR')
        call require_close(n_rz, 2.0_dp, 2.0e-14_dp, &
            'affine-profile N_RZ')
        call require_close(n_zz, 2.0_dp, 2.0e-14_dp, &
            'affine-profile N_ZZ')
        call evaluate_neort_eqdsk_cut_r_flux_curvature(-1.0_dp, 2.0_dp, &
            n_rr, n_rz, n_zz, 1.0_dp, 1.0_dp, 1.5_dp, 0.0_dp, 0.0_dp, &
            1.0_dp, slope, second, flux_first, flux_second)
        call require_close(second, -15.0_dp/4.0_dp, 2.0e-14_dp, &
            'affine-profile graph curvature')
        call require_close(flux_second, -9.0_dp/4.0_dp, 2.0e-14_dp, &
            'affine-profile flux curvature')
    end subroutine check_profile_second_derivative_oracle

    subroutine check_independent_graph_oracle
        real(dp), parameter :: alpha = 0.35_dp
        real(dp), parameter :: graph_r = 0.4_dp
        real(dp), parameter :: a = 1.7_dp
        real(dp), parameter :: b = -0.6_dp
        real(dp), parameter :: d = 0.9_dp
        real(dp), parameter :: separatrix = 2.3_dp
        real(dp) :: graph_z, slope, second, flux_first, flux_second
        real(dp) :: oracle_slope, oracle_second, oracle_flux_first
        real(dp) :: oracle_flux_second, psi_r, psi_z

        graph_z = alpha*graph_r**2
        psi_r = a*graph_r+b*graph_z
        psi_z = b*graph_r+d*graph_z
        call evaluate_neort_eqdsk_cut_r_flux_curvature( &
            -2.0_dp*alpha*graph_r, 1.0_dp, -2.0_dp*alpha, 0.0_dp, &
            0.0_dp, psi_r, psi_z, a, b, d, separatrix, slope, second, &
            flux_first, flux_second)
        oracle_slope = 2.0_dp*alpha*graph_r
        oracle_second = 2.0_dp*alpha
        oracle_flux_first = (a*graph_r+3.0_dp*b*alpha*graph_r**2 &
            +2.0_dp*d*alpha**2*graph_r**3)/separatrix
        oracle_flux_second = (a+6.0_dp*b*alpha*graph_r &
            +6.0_dp*d*alpha**2*graph_r**2)/separatrix
        call require_close(slope, oracle_slope, 2.0e-14_dp, &
            'implicit graph slope')
        call require_close(second, oracle_second, 2.0e-14_dp, &
            'implicit graph curvature')
        call require_close(flux_first, oracle_flux_first, 2.0e-14_dp, &
            'graph flux derivative')
        call require_close(flux_second, oracle_flux_second, 2.0e-14_dp, &
            'graph flux curvature')
    end subroutine check_independent_graph_oracle

    pure real(dp) function falling_factorial(n, order)
        integer, intent(in) :: n, order
        integer :: k

        falling_factorial = 1.0_dp
        do k = 0, order-1
            falling_factorial = falling_factorial*real(n-k, dp)
        end do
    end function falling_factorial

    subroutine require_close(actual, expected, tolerance, label)
        real(dp), intent(in) :: actual, expected, tolerance
        character(*), intent(in) :: label

        if (abs(actual-expected) > tolerance*max(1.0_dp, abs(expected))) then
            write (*, '(a,2es24.15)') trim(label)//' actual/oracle: ', &
                actual, expected
            error stop trim(label)//' mismatch'
        end if
    end subroutine require_close

end program test_gc_eqdsk_cut_flux_curvature
