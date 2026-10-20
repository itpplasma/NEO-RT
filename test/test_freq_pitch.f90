program test_freq_pitch
    ! Thin-orbit canonical frequencies versus pitch at one flux surface and one
    ! kinetic energy, for the POTATO frequency comparison.
    !
    ! Usage: test_freq_pitch.x <runname> [lambda_file]
    !   runname     : reads <runname>.in, in_file, plasma.in, profile.in from the
    !                 current directory, like neo_rt.x; the energy is the thermal
    !                 energy vth of that input (v = vth)
    !   lambda_file : v_par/v at the outboard-midplane field minimum, one value
    !                 per line (default lambda.in)
    ! Writes <runname>_freq_pitch.dat with columns lambda, eta [1/G], omega_b,
    ! Omega_tor, Om_tB, Om_tE [rad/s] for sign_vpar = +1; the header carries
    ! etatp, etadt, s, iota, v.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: s, iota
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort_freq, only: Om_th, Om_ph, Om_tB
    use neort_profiles, only: Om_tE
    use driftorbit, only: vth, etadt, etatp, sign_vpar
    use util, only: files_exist

    implicit none

    character(1024) :: runname, lamfile
    integer :: fid, lid, ios
    real(dp) :: lam, eta, v, Omth, dOmthdv, dOmthdeta
    real(dp) :: Omph, dOmphdv, dOmphdeta, OmtB, dOmtBdv, dOmtBdeta

    call get_command_argument(1, runname)
    if (len_trim(runname) == 0) then
        print *, 'usage: test_freq_pitch.x <runname> [lambda_file]'
        stop 1
    end if
    call get_command_argument(2, lamfile)
    if (len_trim(lamfile) == 0) lamfile = 'lambda.in'

    call neort_init(trim(runname)//".in", "in_file", "in_file_pert")
    if (files_exist("plasma.in", "profile.in")) then
        call neort_prepare_splines("plasma.in", "profile.in")
    end if
    call neort_setup_at_s(s)

    v = vth
    sign_vpar = 1

    open (newunit=fid, file=trim(runname)//'_freq_pitch.dat')
    write (fid, '(A,5ES23.15)') '%# etatp etadt s iota v =', etatp, etadt, s, iota, v
    write (fid, '(A)') '%# lambda, eta [1/G], omega_b, Omega_tor, Om_tB, Om_tE [rad/s]'
    open (newunit=lid, file=trim(lamfile), status='old', action='read')
    do
        read (lid, *, iostat=ios) lam
        if (ios /= 0) exit
        eta = (1.0_dp - lam**2)*etadt
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
        write (fid, '(6ES23.15)') lam, eta, Omth, Omph, OmtB, Om_tE
    end do
    close (lid)
    close (fid)
end program test_freq_pitch
