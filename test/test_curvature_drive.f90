program test_curvature_drive
    ! General ideal drive H/E = eta*B*dB_L/B + 2*(1 - eta*B)*kappa.xi_perp with
    ! kappa.xi_perp read from a second Boozer table (doc/hamiltonian_drive.md).
    ! Oracle: the manufactured Fourier coefficients evaluated here, independently
    ! of the NEO-RT readers, and the closed-form ratio to the default Boozer drive
    ! (2 - eta*B)*dB_L/B. Checks: (1) no table keeps the default drive; (2) the
    ! perturbation file itself as table reproduces it; (3) a zero table leaves the
    ! mirror term only; (4) an independent table gives the closed form at every
    ! sample; (5) a later initialisation without the flag disables the table.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_lib, only: config_t, neort_init
    use do_magfie_mod, only: do_magfie, init_magfie_at_s, magfie_thread_init, set_s
    use do_magfie_pert_mod, only: init_magfie_pert_at_s, curvature_drive
    use driftorbit, only: etatp
    use neort_transport, only: timestep_transport

    implicit none

    character(len=*), parameter :: eq_file = 'curv_equilibrium.bc', pert_file = 'curv_perturbation.bc'
    character(len=*), parameter :: zero_file = 'curv_zero.bc', other_file = 'curv_other.bc'
    real(dp), parameter :: pert_c(3) = [1.0e-4_dp, 0.5e-4_dp, 2.0e-4_dp]
    real(dp), parameter :: pert_s(3) = [0.3e-4_dp, 0.0_dp, 0.7e-4_dp]
    real(dp), parameter :: other_c(3) = [-0.4e-4_dp, 0.2e-4_dp, 0.0_dp]
    real(dp), parameter :: other_s(3) = [0.1e-4_dp, -0.3e-4_dp, 0.9e-4_dp]
    real(dp), parameter :: zeros(3) = 0.0_dp
    real(dp), parameter :: v = 1.0e8_dp, eta = 0.2_dp, tol = 1.0e-12_dp
    integer, parameter :: nsample = 7
    complex(dp) :: base(nsample), h(nsample), expect
    real(dp) :: theta(nsample), etab(nsample), worst
    integer :: k

    call write_equilibrium(eq_file)
    call write_table(pert_file, pert_c, pert_s)
    call write_table(zero_file, zeros, zeros)
    call write_table(other_file, other_c, other_s)
    do k = 1, nsample
        theta(k) = -3.0_dp + 6.0_dp*real(k - 1, dp)/real(nsample - 1, dp)
    end do

    call evaluate(.false., pert_file, base)
    if (curvature_drive) error stop 'curvature drive active without pertfile_curvature'

    call evaluate(.true., pert_file, h)
    if (.not. curvature_drive) error stop 'curvature table not activated'
    worst = maxval(abs(h - base)/abs(base))
    print *, 'table = perturbation file: max relative change', worst
    if (worst > tol) error stop 'kappa.xi = dB_L/B must reproduce the Boozer drive'

    call evaluate(.true., zero_file, h)
    worst = 0.0_dp
    do k = 1, nsample
        expect = base(k)*etab(k)/(2.0_dp - etab(k))
        worst = max(worst, abs(h(k) - expect)/abs(base(k)))
    end do
    print *, 'zero table: max deviation from eta*B/(2 - eta*B)', worst
    if (worst > tol) error stop 'zero curvature must leave the mirror term only'

    call evaluate(.true., other_file, h)
    worst = 0.0_dp
    do k = 1, nsample
        expect = base(k)/(2.0_dp - etab(k))*(etab(k) &
            + 2.0_dp*(1.0_dp - etab(k))*series(other_c, other_s, theta(k)) &
            /series(pert_c, pert_s, theta(k)))
        worst = max(worst, abs(h(k) - expect)/abs(base(k)))
    end do
    print *, 'independent table: max deviation from closed form', worst
    if (worst > tol) error stop 'general drive differs from eta*B*dB_L/B + 2(1-eta*B)kappa.xi'

    call evaluate(.false., other_file, h)
    if (curvature_drive) error stop 'reinitialisation without the flag must disable the table'
    if (maxval(abs(h - base)/abs(base)) > 0.0_dp) error stop 'default drive changed after reset'

    call delete_file(eq_file)
    call delete_file(pert_file)
    call delete_file(zero_file)
    call delete_file(other_file)
    print *, 'test_curvature_drive PASSED'

contains

    subroutine evaluate(use_table, table, values)
        logical, intent(in) :: use_table
        character(len=*), intent(in) :: table
        complex(dp), intent(out) :: values(nsample)
        type(config_t) :: config
        real(dp) :: y(8), ydot(8), x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        integer :: j

        config%s = 0.5_dp
        config%vth = 4.0e7_dp
        config%qs = 1.0_dp
        config%ms = 2.0_dp
        config%inp_swi = 9
        config%pertfile = .true.
        config%pertfile_curvature = use_table
        config%log_level = 0
        call neort_init(config, eq_file, pert_file, table)
        call magfie_thread_init()
        call set_s(0.5_dp)
        call init_magfie_at_s()
        call init_magfie_pert_at_s()
        etatp = 0.0_dp     ! sample as trapped: v_parallel is taken from y(2)
        do j = 1, nsample
            y = 0.0_dp
            y(1) = theta(j)
            y(2) = 0.3_dp*v
            call timestep_transport(v, eta, size(y), 0.0_dp, y, ydot)
            values(j) = cmplx(ydot(3), ydot(4), dp)
            x = [0.5_dp, 0.0_dp, theta(j)]
            call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            etab(j) = eta*bmod
        end do
    end subroutine evaluate

    complex(dp) function series(c, s, th)
        ! NEO convention: sum (bmnc - i bmns) exp(i m theta), m = -1, 0, 1.
        real(dp), intent(in) :: c(3), s(3), th
        integer :: m

        series = (0.0_dp, 0.0_dp)
        do m = -1, 1
            series = series + cmplx(c(m + 2), -s(m + 2), dp)*exp((0.0_dp, 1.0_dp)*real(m, dp)*th)
        end do
    end function series

    subroutine write_equilibrium(path)
        character(len=*), intent(in) :: path
        integer :: unit, surface
        real(dp), parameter :: surfaces(3) = [0.2_dp, 0.5_dp, 0.8_dp]

        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, '(a)') 'CC manufactured equilibrium field'
        write (unit, '(a)') 'CC'
        write (unit, '(a)') 'CC'
        write (unit, '(a)') 'CC'
        write (unit, '(a)') 'm0b n0b nsurf nper flux a R'
        write (unit, '(*(g0,1x))') 1, 0, 3, 1, 1.0_dp, 0.5_dp, 5.0_dp
        do surface = 1, size(surfaces)
            write (unit, '(a)') 's iota Jpol Itor pprime sqrtg00'
            write (unit, '(a)') 'units A A Pa m3'
            write (unit, '(*(g0,1x))') surfaces(surface), 0.4_dp, 100.0_dp, 10.0_dp, 0.0_dp, 1.0_dp
            write (unit, '(a)') 'm n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
            write (unit, '(*(g0,1x))') 0, 0, 5.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 2.0_dp, 0.0_dp
            write (unit, '(*(g0,1x))') 1, 0, 0.5_dp, 0.0_dp, 0.0_dp, 0.5_dp, 0.0_dp, 0.0_dp, -0.2_dp, 0.0_dp
        end do
        close (unit)
    end subroutine write_equilibrium

    subroutine write_table(path, c, s)
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: c(3), s(3)
        integer :: unit, surface, m
        real(dp), parameter :: surfaces(3) = [0.2_dp, 0.5_dp, 0.8_dp]

        open (newunit=unit, file=path, status='replace', action='write')
        write (unit, '(a)') 'CC manufactured perturbation table'
        write (unit, '(a)') 'CC'
        write (unit, '(a)') 'CC'
        write (unit, '(a)') 'CC'
        write (unit, '(a)') 'm0b n0b nsurf nper flux a R'
        write (unit, '(*(g0,1x))') 2, 0, 3, 1, 1.0_dp, 1.0_dp, 1.0_dp
        do surface = 1, size(surfaces)
            write (unit, '(a)') 's iota Jpol Itor pprime sqrtg'
            write (unit, '(a)') 'units A A Pa m3'
            write (unit, '(*(es24.16e3,1x))') surfaces(surface), 1.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp
            write (unit, '(a)') 'm n rmnc rmns zmnc zmns vmnc vmns bmnc bmns'
            do m = -1, 1
                write (unit, '(*(g0,1x))') m, 3, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
                    c(m + 2), s(m + 2)
            end do
        end do
        close (unit)
    end subroutine write_table

    subroutine delete_file(path)
        character(len=*), intent(in) :: path
        integer :: unit

        open (newunit=unit, file=path, status='old')
        close (unit, status='delete')
    end subroutine delete_file

end program test_curvature_drive
