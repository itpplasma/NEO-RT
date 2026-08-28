! Dumps reference axisymmetric-field values from the real do_magfie_mod at
! sample (s, theta) points, for cross-checking the C field port.
program field_ref
    use do_magfie_mod, only: set_s, read_boozer_file, init_magfie_at_s, do_magfie, &
                             inp_swi, bfac
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    integer, parameter :: ntheta = 33
    real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: s, theta, pi
    integer :: i
    character(len=256) :: path

    call get_command_argument(1, path)
    pi = 4.0d0 * atan(1.0d0)
    s = 0.5d0
    inp_swi = 9
    bfac = 1.0d0

    call read_boozer_file(trim(path))
    call set_s(s)
    call init_magfie_at_s()

    do i = 0, ntheta - 1
        theta = -pi + 2.0d0 * pi * real(i, dp) / real(ntheta - 1, dp)
        x = [s, 0.0d0, theta]
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        write (*, '(11(1X,ES24.16E3))') theta, bmod, sqrtg, bder(1), bder(3), &
            hcovar(2), hcovar(3), hctrvr(2), hctrvr(3), real(s, dp), 0.0d0
    end do
end program field_ref
