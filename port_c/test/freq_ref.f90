! Reference dump of canonical frequencies from the real Fortran neort_freq.
program freq_ref
    use do_magfie_mod, only: set_s, read_boozer_file, init_magfie_at_s, do_magfie, &
                             inp_swi, bfac, R0
    use neort_magfie, only: init_flux_surface_average
    use neort_profiles, only: read_and_init_plasma_input, read_and_init_profile_input, vth
    use driftorbit, only: etatp, etadt, sign_vpar, magdrift
    use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline, &
                          Om_th, Om_tB, Om_ph
    use logger, only: set_log_level, LVL_SILENT
    implicit none
    integer, parameter :: dp = kind(1.0d0)
    real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3), s, v
    real(dp) :: ee(4), eta, a1, a2, a3
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
    magdrift = .true.
    sign_vpar = 1.0d0

    call init_canon_freq_trapped_spline
    call init_canon_freq_passing_spline

    v = vth
    ee = [etatp + 0.3d0*(etadt - etatp), 0.7d0*etatp, &
          etatp*(1.0d0 + 1.0d-7), etatp*(1.0d0 - 1.0d-7)]
    do i = 1, 4
        eta = ee(i)
        call Om_th(v, eta, a1, a2, a3); write (*, '(3(1X,ES24.16E3))') a1, a2, a3
        call Om_tB(v, eta, a1, a2, a3); write (*, '(3(1X,ES24.16E3))') a1, a2, a3
        call Om_ph(v, eta, a1, a2, a3); write (*, '(3(1X,ES24.16E3))') a1, a2, a3
    end do
end program freq_ref
