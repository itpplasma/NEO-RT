program test_gc_eqdsk_cut_jet_fourth
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: icp, ipoint, nrad, nzet, psi_sep, rad, zet, splpsi
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        eqdsk_cylindrical_field_t, initialize_eqdsk_cylindrical_field
    use neort_gc_eqdsk_cut_jet, only: EQDSK_CUT_JET_SUCCESS, &
        eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
    use util_for_test, only: pass_test
    implicit none

    type :: taylor2_t
        real(dp) :: coefficient(0:2,0:2) = 0.0_dp
    end type taylor2_t

    type(eqdsk_cylindrical_field_t) :: field
    type(eqdsk_cut_jet_t) :: actual
    real(dp) :: coefficient(0:5,0:5), expected_jet(15)
    real(dp) :: expected_numerator(3), position(3)
    character(len=1024) :: path
    integer :: status, cell_R, cell_Z, cell_pointer, i, j

    call get_environment_variable('EQDSK_FILE', path)
    call require(len_trim(path) > 0, 'EQDSK_FILE is required')
    call initialize_eqdsk_cylindrical_field(trim(path), 1.0_dp, field, status)
    call require(status == 0, 'failed to initialize EQDSK fixture')
    call require(nrad >= 3, 'EQDSK R grid is too small')
    call require(nzet >= 3, 'EQDSK Z grid is too small')
    call require(icp > 0, 'EQDSK spline-cell table is empty')

    cell_R = max(1, min(nrad-1, nrad/2))
    cell_Z = max(1, min(nzet-1, nzet/2))
    cell_pointer = ipoint(cell_R, cell_Z)
    call require(cell_pointer >= 1, 'EQDSK cell pointer is invalid')
    call require(cell_pointer <= size(splpsi, 3), &
        'EQDSK cell pointer exceeds spline-cell table')
    do i = 0, 5
        do j = 0, 5
            coefficient(i,j) = splpsi(i+1,j+1,cell_pointer)
        end do
    end do

    position = [0.5_dp*(rad(cell_R)+rad(cell_R+1)), &
        0.5_dp*(zet(cell_Z)+zet(cell_Z+1)), 0.0_dp]
    call independent_cell_jet(position(1)-rad(cell_R), &
        position(2)-zet(cell_Z), coefficient, expected_jet)
    call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, &
        [0.0_dp, 0.0_dp, 0.0_dp], actual, status)
    call require(status == EQDSK_CUT_JET_SUCCESS, &
        'runtime EQDSK cut-jet evaluation failed')
    do i = 1, 15
        call require_close(actual%psi_jet(i), expected_jet(i), &
            'runtime tensor polynomial jet')
    end do

    call numerator_hessian_oracle(position(1), actual%psi_jet, actual%f_jet, &
        psi_sep, expected_numerator)
    call require_close(actual%d2_cut_numerator_d_R2, expected_numerator(1), &
        'runtime numerator N_RR')
    call require_close(actual%d2_cut_numerator_d_RdZ, expected_numerator(2), &
        'runtime numerator N_RZ')
    call require_close(actual%d2_cut_numerator_d_Z2, expected_numerator(3), &
        'runtime numerator N_ZZ')

    write (*, '(a)') 'test_gc_eqdsk_cut_jet_fourth OK'
    call pass_test

contains

    pure subroutine independent_cell_jet(delta_r, delta_z, coefficient, value)
        real(dp), intent(in) :: delta_r, delta_z, coefficient(0:5,0:5)
        real(dp), intent(out) :: value(15)
        integer, parameter :: derivative_r(15) = [0, 1, 0, 2, 1, 0, 3, 2, &
            1, 0, 4, 3, 2, 1, 0]
        integer, parameter :: derivative_z(15) = [0, 0, 1, 0, 1, 2, 0, 1, &
            2, 3, 0, 1, 2, 3, 4]
        integer :: k, p, q

        value = 0.0_dp
        do k = 1, 15
            do p = derivative_r(k), 5
                do q = derivative_z(k), 5
                    value(k) = value(k)+coefficient(p,q) &
                        *falling_factorial(p, derivative_r(k)) &
                        *falling_factorial(q, derivative_z(k)) &
                        *delta_r**(p-derivative_r(k)) &
                        *delta_z**(q-derivative_z(k))
                end do
            end do
        end do
    end subroutine independent_cell_jet

    pure function constant_taylor(value) result(result)
        real(dp), intent(in) :: value
        type(taylor2_t) :: result

        result = taylor2_t()
        result%coefficient(0,0) = value
    end function constant_taylor

    pure function derivative_taylor(value, first_r, first_z, second_rr, &
            second_rz, second_zz) result(result)
        real(dp), intent(in) :: value, first_r, first_z, second_rr, second_rz
        real(dp), intent(in) :: second_zz
        type(taylor2_t) :: result

        result = taylor2_t()
        result%coefficient(0,0) = value
        result%coefficient(1,0) = first_r
        result%coefficient(0,1) = first_z
        result%coefficient(2,0) = 0.5_dp*second_rr
        result%coefficient(1,1) = second_rz
        result%coefficient(0,2) = 0.5_dp*second_zz
    end function derivative_taylor

    pure function add_taylor(left, right) result(result)
        type(taylor2_t), intent(in) :: left, right
        type(taylor2_t) :: result

        result%coefficient = left%coefficient+right%coefficient
    end function add_taylor

    pure function multiply_taylor(left, right) result(result)
        type(taylor2_t), intent(in) :: left, right
        type(taylor2_t) :: result
        integer :: i, j, p, q

        result = taylor2_t()
        do i = 0, 2
            do j = 0, 2-i
                do p = 0, i
                    do q = 0, j
                        result%coefficient(i,j) = result%coefficient(i,j) &
                            +left%coefficient(p,q)* &
                            right%coefficient(i-p,j-q)
                    end do
                end do
            end do
        end do
    end function multiply_taylor

    pure function scale_taylor(factor, value) result(result)
        real(dp), intent(in) :: factor
        type(taylor2_t), intent(in) :: value
        type(taylor2_t) :: result

        result%coefficient = factor*value%coefficient
    end function scale_taylor

    subroutine numerator_hessian_oracle(radius, psi_jet, f_jet, separatrix, &
            value)
        real(dp), intent(in) :: radius, psi_jet(15), f_jet(3), separatrix
        real(dp), intent(out) :: value(3)
        type(taylor2_t) :: psi_R, psi_Z, psi_RR, psi_RZ, psi_ZZ
        type(taylor2_t) :: flux, flux_increment, f_squared, q, a, numerator
        type(taylor2_t) :: radius_taylor, difference, squared_gradient

        psi_R = derivative_taylor(psi_jet(2), psi_jet(4), psi_jet(5), &
            psi_jet(7), psi_jet(8), psi_jet(9))
        psi_Z = derivative_taylor(psi_jet(3), psi_jet(5), psi_jet(6), &
            psi_jet(8), psi_jet(9), psi_jet(10))
        psi_RR = derivative_taylor(psi_jet(4), psi_jet(7), psi_jet(8), &
            psi_jet(11), psi_jet(12), psi_jet(13))
        psi_RZ = derivative_taylor(psi_jet(5), psi_jet(8), psi_jet(9), &
            psi_jet(12), psi_jet(13), psi_jet(14))
        psi_ZZ = derivative_taylor(psi_jet(6), psi_jet(9), psi_jet(10), &
            psi_jet(13), psi_jet(14), psi_jet(15))

        difference = constant_taylor(psi_jet(1))
        difference%coefficient(0,0) = 0.0_dp
        difference%coefficient(1,0) = psi_jet(2)
        difference%coefficient(0,1) = psi_jet(3)
        difference%coefficient(2,0) = 0.5_dp*psi_jet(4)
        difference%coefficient(1,1) = psi_jet(5)
        difference%coefficient(0,2) = 0.5_dp*psi_jet(6)
        flux_increment = difference
        flux = add_taylor(constant_taylor(f_jet(1)), &
            add_taylor(scale_taylor(f_jet(2)/separatrix, flux_increment), &
            scale_taylor(0.5_dp*f_jet(3)/separatrix**2, &
            multiply_taylor(flux_increment, flux_increment))))

        f_squared = multiply_taylor(flux, flux)
        squared_gradient = add_taylor(multiply_taylor(psi_R, psi_R), &
            multiply_taylor(psi_Z, psi_Z))
        q = add_taylor(f_squared, squared_gradient)
        difference = add_taylor(psi_ZZ, scale_taylor(-1.0_dp, psi_RR))
        a = add_taylor(multiply_taylor(multiply_taylor(psi_R, psi_Z), difference), &
            multiply_taylor(psi_RZ, add_taylor( &
            multiply_taylor(psi_R, psi_R), &
            scale_taylor(-1.0_dp, multiply_taylor(psi_Z, psi_Z)))))
        radius_taylor = constant_taylor(radius)
        radius_taylor%coefficient(1,0) = 1.0_dp
        numerator = add_taylor(multiply_taylor(psi_Z, q), &
            multiply_taylor(radius_taylor, a))
        value = [2.0_dp*numerator%coefficient(2,0), &
            numerator%coefficient(1,1), 2.0_dp*numerator%coefficient(0,2)]
    end subroutine numerator_hessian_oracle

    pure real(dp) function falling_factorial(n, order)
        integer, intent(in) :: n, order
        integer :: k

        falling_factorial = 1.0_dp
        do k = 0, order-1
            falling_factorial = falling_factorial*real(n-k, dp)
        end do
    end function falling_factorial

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

    subroutine require_close(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label
        real(dp), parameter :: tolerance = 3.0e-11_dp

        if (abs(actual-expected) > tolerance*max(1.0_dp, abs(expected))) then
            write (*, '(a,2(1x,es24.16))') trim(label), actual, expected
            error stop trim(label)//' mismatch'
        end if
    end subroutine require_close

end program test_gc_eqdsk_cut_jet_fourth
