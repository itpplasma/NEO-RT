program test_frequencies
    ! Canonical frequencies vs pitch at one flux surface and one kinetic energy.
    !
    ! Usage: test_frequencies.x <runname> [ux]
    !   runname : reads <runname>.in, in_file (, in_file_pert), plasma.in,
    !             profile.in from the current directory, like neo_rt.x
    !   ux      : v = ux*vth, i.e. E_kin = ux^2 * T_local (default 1.0)
    !
    ! Writes <runname>_ux<tag>_freq_vs_eta_{t,pco,pct}.dat (uniform eta grids)
    ! and    <runname>_ux<tag>_freq_vs_eta_tzoom.dat (trapped side, log-spaced
    ! distance to the trapped-passing boundary down to (eta-etatp)/etatp=1e-6,
    ! header carries etatp, etadt, ux, s). Columns: eta [1/G], omega_b [rad/s],
    ! Omega_tor [rad/s].
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: s
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_freq, only: Om_th, Om_ph, Om_tB
    use neort_profiles, only: Om_tE
    use driftorbit, only: vth, etadt, etatp, sign_vpar
    use util, only: files_exist

    implicit none

    integer, parameter :: neta = 1000    ! steps of the uniform grids
    integer, parameter :: nzoom = 400    ! steps of the log-spaced zoom grid
    real(dp), parameter :: distmin = 1.0e-6_dp  ! innermost (eta-etatp)/etatp

    character(1024) :: runname, arg2, tag
    integer :: i, fid
    real(dp) :: Omth, dOmthdv, dOmthdeta, Omph, dOmphdv, dOmphdeta
    real(dp) :: OmtB, dOmtBdv, dOmtBdeta
    real(dp) :: v, eta, ux, xlo, xhi, x

    call get_command_argument(1, runname)
    if (len_trim(runname) == 0) then
        print *, 'usage: test_frequencies.x <runname> [ux]'
        stop 1
    end if
    ux = 1.0_dp
    call get_command_argument(2, arg2)
    if (len_trim(arg2) > 0) read (arg2, *) ux

    call neort_init(trim(runname)//".in", "in_file", "in_file_pert")
    if (files_exist("plasma.in", "profile.in")) then
        call neort_prepare_splines("plasma.in", "profile.in")
    end if
    call neort_setup_at_s(s)

    v = ux*vth

    write (tag, '(F0.2)') ux
    do i = 1, len_trim(tag)
        if (tag(i:i) == '.') tag(i:i) = 'p'
    end do

    print *, 'test_frequencies: s =', s, ' ux =', ux, ' v =', v
    print *, '  etatp =', etatp, ' etadt =', etadt

    open (newunit=fid, file=trim(runname)//'_ux'//trim(tag)//'_freq_vs_eta_t.dat')
    write (fid, *) '%# eta [1/G], omega_b, Omega_tor, Om_tB, Om_tE [rad/s] trapped'
    sign_vpar = 1
    do i = 1, neta
        eta = etatp + i*(etadt - etatp)/neta
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
        write (fid, *) eta, Omth, Omph, OmtB, Om_tE
    end do
    close (fid)

    open (newunit=fid, file=trim(runname)//'_ux'//trim(tag)//'_freq_vs_eta_tzoom.dat')
    write (fid, '(A,ES23.15,A,ES23.15,A,F8.4,A,ES23.15)') &
        '%# etatp =', etatp, ' etadt =', etadt, ' ux =', ux, ' s =', s
    write (fid, *) '%# eta [1/G], omega_b, Omega_tor, Om_tB, Om_tE [rad/s] trapped zoom'
    sign_vpar = 1
    xlo = log10(distmin)
    xhi = log10(0.999_dp*(etadt - etatp)/etatp)
    do i = 0, nzoom
        x = xlo + (xhi - xlo)*i/nzoom
        eta = etatp*(1.0_dp + 10.0_dp**x)
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
        write (fid, *) eta, Omth, Omph, OmtB, Om_tE
    end do
    close (fid)

    open (newunit=fid, file=trim(runname)//'_ux'//trim(tag)//'_freq_vs_eta_pco.dat')
    write (fid, *) '%# eta [1/G], omega_b, Omega_tor, Om_tB, Om_tE [rad/s] co-passing'
    sign_vpar = 1
    do i = 0, neta - 1
        eta = i*etatp/neta
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
        write (fid, *) eta, Omth, Omph, OmtB, Om_tE
    end do
    close (fid)

    open (newunit=fid, file=trim(runname)//'_ux'//trim(tag)//'_freq_vs_eta_pct.dat')
    write (fid, *) '%# eta [1/G], omega_b, Omega_tor, Om_tB, Om_tE [rad/s] counter-passing'
    sign_vpar = -1
    do i = 0, neta - 1
        eta = i*etatp/neta
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
        write (fid, *) eta, -Omth, Omph, OmtB, Om_tE
    end do
    close (fid)

end program test_frequencies
