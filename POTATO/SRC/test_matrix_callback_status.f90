program test_matrix_callback_status
    use matrix_callback_status_mod, only: matrix_callback_bounce_failed

    implicit none

    integer :: ierr
    integer :: i
    double precision :: integral, exact
    external :: failing_class_matrix, valid_class_matrix
    external :: failing_outer_matrix, valid_outer_matrix
    external :: cusp_outer_matrix

    call configure_class_sampler
    call sample_matrix(failing_class_matrix, ierr)
    if (ierr /= matrix_callback_bounce_failed) then
        error stop 'sample_matrix did not propagate callback failure'
    endif

    call configure_class_sampler
    call sample_matrix(valid_class_matrix, ierr)
    if (ierr /= 0) error stop 'sample_matrix rejected valid callback'

    call configure_outer_sampler
    call sample_matrix_out(failing_outer_matrix, ierr)
    if (ierr /= matrix_callback_bounce_failed) then
        error stop 'sample_matrix_out did not propagate callback failure'
    endif

    call configure_outer_sampler
    call sample_matrix_out(valid_outer_matrix, ierr)
    if (ierr /= 0) error stop 'sample_matrix_out rejected valid callback'

    call configure_cusp_sampler
    call sample_matrix_out(cusp_outer_matrix, ierr)
    if (ierr /= 0) error stop 'sample_matrix_out did not converge a valid cusp'
    call integrate_outer_grid(integral)
    exact = (2.d0/3.d0)*(1.d0-0.37d0)**1.5d0
    if (abs(integral-exact) > 2.d-3*exact) then
        error stop 'sample_matrix_out cusp integral missed analytic oracle'
    endif

contains

    subroutine configure_class_sampler
        use sample_matrix_mod, only: nlagr, n1, n2, itermax, eps, xbeg, xend

        nlagr = 2
        n1 = 1
        n2 = 1
        itermax = 2
        eps = 1.0d6
        xbeg = 0.0d0
        xend = 1.0d0
    end subroutine configure_class_sampler

    subroutine configure_outer_sampler
        use sample_matrix_out_mod, only: nlagr, n1, n2, npoi, itermax, &
            eps, xbeg, xend

        nlagr = 2
        n1 = 1
        n2 = 1
        npoi = 5
        itermax = 2
        eps = 1.0d6
        xbeg = 0.0d0
        xend = 1.0d0
    end subroutine configure_outer_sampler

    subroutine configure_cusp_sampler
        use sample_matrix_out_mod, only: nlagr, n1, n2, npoi, itermax, &
            eps, xbeg, xend

        nlagr = 3
        n1 = 1
        n2 = 1
        npoi = 11
        itermax = 10
        eps = 1.d-3
        xbeg = 0.d0
        xend = 1.d0
    end subroutine configure_cusp_sampler

    subroutine integrate_outer_grid(result)
        use sample_matrix_out_mod, only: npoi, xarr, amat_arr
        double precision, intent(out) :: result

        result = 0.d0
        do i = 2, npoi
            result = result + 0.5d0*(amat_arr(1,1,i-1)+amat_arr(1,1,i)) &
                *(xarr(i)-xarr(i-1))
        enddo
    end subroutine integrate_outer_grid

end program test_matrix_callback_status

subroutine failing_class_matrix
    use sample_matrix_mod, only: x, amat
    use matrix_callback_status_mod, only: matrix_callback_bounce_failed, &
        set_matrix_callback_error

    amat(1,1) = cmplx(x, 0.0d0, kind=8)
    if (x > 0.4d0) call set_matrix_callback_error(matrix_callback_bounce_failed)
end subroutine failing_class_matrix

subroutine valid_class_matrix
    use sample_matrix_mod, only: x, amat

    amat(1,1) = cmplx(x, 0.0d0, kind=8)
end subroutine valid_class_matrix

subroutine failing_outer_matrix
    use sample_matrix_out_mod, only: x, amat
    use matrix_callback_status_mod, only: matrix_callback_bounce_failed, &
        set_matrix_callback_error

    amat(1,1) = cmplx(x, 0.0d0, kind=8)
    if (x > 0.4d0) call set_matrix_callback_error(matrix_callback_bounce_failed)
end subroutine failing_outer_matrix

subroutine valid_outer_matrix
    use sample_matrix_out_mod, only: x, amat

    amat(1,1) = cmplx(x, 0.0d0, kind=8)
end subroutine valid_outer_matrix

subroutine cusp_outer_matrix
    use sample_matrix_out_mod, only: x, amat

    amat(1,1) = sqrt(max(0.d0, x-0.37d0))
end subroutine cusp_outer_matrix
