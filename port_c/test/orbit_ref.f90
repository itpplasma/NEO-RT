! Reference dump of bounce time and bounce averages from the real Fortran orbit
! integration (DVODE), for cross-checking the C CVODE port on real orbits.
program orbit_ref
    use do_magfie_mod, only: set_s, read_boozer_file, init_magfie_at_s, do_magfie, &
                             inp_swi, bfac, R0
    use neort_magfie, only: init_flux_surface_average
    use neort_profiles, only: read_and_init_plasma_input, read_and_init_profile_input, vth
    use driftorbit, only: etatp, etadt, sign_vpar
    use neort_orbit, only: bounce, bounce_time
    use logger, only: set_log_level, LVL_SILENT
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3), s
    real(dp) :: vv(3), ee(3), taub, bavg(7), tt
    character(len=256) :: bpath, ppath, rpath
    integer :: i

    call get_command_argument(1, bpath)
    call get_command_argument(2, ppath)
    call get_command_argument(3, rpath)
    call set_log_level(LVL_SILENT)
    s = 0.5d0
    inp_swi = 9
    bfac = 1.0d0
    call read_boozer_file(trim(bpath))
    call set_s(s)
    call init_magfie_at_s()
    x = [s, 0.0d0, 0.0d0]
    call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    call init_flux_surface_average(s)
    call read_and_init_plasma_input(trim(ppath), s)
    call read_and_init_profile_input(trim(rpath), s, R0, 1.0d0, 1.0d0)

    sign_vpar = 1.0d0
    vv = [vth, vth, 1.5d0*vth]
    ee = [0.5d0*(etatp + etadt), 0.3d0*etatp, 0.8d0*etatp + 0.2d0*etadt]

    do i = 1, 3
        tt = bounce_time(vv(i), ee(i))
        call bounce(vv(i), ee(i), taub, bavg)
        write (*, '(6(1X,ES24.16E3))') tt, taub, bavg(1), bavg(2), bavg(3), bavg(6)
    end do
end program orbit_ref
