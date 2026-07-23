program test_matrix_callback_status
    use matrix_callback_status_mod, only: matrix_callback_bounce_failed

    implicit none

    integer :: ierr
    external :: failing_class_matrix, valid_class_matrix
    external :: failing_outer_matrix, valid_outer_matrix

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
