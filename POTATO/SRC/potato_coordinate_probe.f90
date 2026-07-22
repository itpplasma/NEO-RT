program potato_coordinate_probe
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: ierrfield, psi_axis, psi_sep
    use field_sub, only: psif

    implicit none

    character(len=1024) :: points_file
    integer :: input_unit, ios, sample, status
    real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(dp) :: b_cartesian(3), cphi, sphi

    interface
        subroutine magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            import dp
            real(dp), intent(in) :: x(3)
            real(dp), intent(out) :: bmod, sqrtg
            real(dp), intent(out) :: bder(3), hcovar(3), hctrvr(3), hcurl(3)
        end subroutine magfie
    end interface

    call get_command_argument(1, points_file, status=status)
    if (status /= 0 .or. len_trim(points_file) == 0) then
        error stop 'usage: potato_coordinate_probe.x POINTS_FILE'
    end if

    open(newunit=input_unit, file=trim(points_file), status='old', action='read')
    sample = 0
    do
        read(input_unit, *, iostat=ios) x
        if (ios < 0) exit
        if (ios > 0) error stop 'invalid coordinate-probe point'
        sample = sample + 1

        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        if (ierrfield /= 0) error stop 'coordinate-probe point outside field domain'

        cphi = cos(x(2))
        sphi = sin(x(2))
        b_cartesian = bmod*[ &
            hctrvr(1)*cphi - x(1)*hctrvr(2)*sphi, &
            hctrvr(1)*sphi + x(1)*hctrvr(2)*cphi, &
            hctrvr(3)]

        write(*, '(a,1x,i0,1x,*(es24.16e3,1x))') 'POTATO_COORD_SAMPLE', sample, &
            x, b_cartesian, bmod, sqrtg, psif, psi_axis, psi_sep, &
            bder, hcovar, hctrvr, hcurl
    end do
    close(input_unit)

end program potato_coordinate_probe
