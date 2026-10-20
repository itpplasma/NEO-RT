program test_curvature_tables
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: bfac, inp_swi, set_s
    use do_magfie_pert_mod, only: inp_swi_pert, read_boozer_pert_file, &
        read_boozer_curvature_file, init_magfie_pert_at_s, magfie_pert_thread_init, &
        do_magfie_pert_amp, do_magfie_pert_kappa_amp
    implicit none
    integer :: step, format, nm, j
    real(dp) :: x(3), scale
    complex(dp) :: b, k, expected

    inp_swi = 9
    bfac = 1.7_dp
    do step = 1, 6
        format = 9
        nm = 3
        if (step == 2) format = 8
        if (step == 4) nm = 2
        if (step >= 5) nm = 4
        if (step == 5) format = 8
        inp_swi_pert = format
        call write_table('pert.bc', format, nm, 3, 1.0_dp)
        call write_table('kappa.bc', format, nm, 4, real(step, dp))
        call read_boozer_pert_file('pert.bc')
        call read_boozer_curvature_file('kappa.bc')
        ! Reuse persistent OpenMP workers across column/mode-count changes.
        !$omp parallel private(j, x, scale, b, k, expected)
        if (step == 1) call magfie_pert_thread_init()
        call set_s(0.4_dp)
        call init_magfie_pert_at_s()
        do j = 1, 5
            x = [0.2_dp + 0.1_dp*j, 0.0_dp, -0.7_dp + 0.3_dp*j]
            scale = 1.0_dp + x(1)
            expected = 0.0_dp
            block
                integer :: m
                do m = 0, nm - 1
                    if (format == 8) then
                        expected = expected + scale*(m + 1)*cos(m*x(3))
                    else
                        expected = expected + scale*cmplx(m + 1.0_dp, -0.3_dp*m, dp) &
                            *exp(cmplx(0.0_dp, m*x(3), dp))
                    end if
                end do
            end block
            expected = bfac*expected
            call do_magfie_pert_amp(x, b)
            call do_magfie_pert_kappa_amp(x, k)
            if (abs(b - expected) > 1e-11_dp) error stop 'perturbation reload oracle'
            if (abs(k - step*expected) > 1e-11_dp) error stop 'curvature reload oracle'
        end do
        !$omp end parallel
    end do
    call delete_file('pert.bc')
    call delete_file('kappa.bc')
    print *, 'curvature table reload and parallel oracle PASSED'
contains
    subroutine write_table(path, format, nm, ns, factor)
        character(len=*), intent(in) :: path
        integer, intent(in) :: format, nm, ns
        real(dp), intent(in) :: factor
        integer :: u, i, m
        real(dp) :: s, c, sn
        open(newunit=u, file=path, status='replace')
        write(u, '(a)') 'CC', 'CC', 'CC', 'CC', 'header'
        write(u, *) nm - 1, 0, ns, 1, 1.0_dp, 1.0_dp, 1.0_dp
        do i = 1, ns
            s = 0.1_dp + 0.8_dp*real(i - 1, dp)/real(ns - 1, dp)
            write(u, '(a)') 'surface', 'units'
            write(u, *) s, 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write(u, '(a)') 'harmonics'
            do m = 0, nm - 1
                c = factor*1e-4_dp*(1 + s)*(m + 1)
                sn = factor*1e-4_dp*(1 + s)*0.3_dp*m
                if (format == 8) then
                    write(u, *) m, 3, 0.0_dp, 0.0_dp, 0.0_dp, c
                else
                    write(u, *) m, 3, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, c, sn
                end if
            end do
        end do
        close(u)
    end subroutine write_table

    subroutine delete_file(path)
        character(len=*), intent(in) :: path
        integer :: u
        open(newunit=u, file=path, status='old')
        close(u, status='delete')
    end subroutine delete_file
end program test_curvature_tables
