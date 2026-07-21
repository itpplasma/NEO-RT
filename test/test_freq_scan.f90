program test_freq_scan
    ! Thin-orbit canonical frequencies vs flux surface for the rung-2 benchmark
    ! against POTATO. At each s the bounce frequency Om_th and the toroidal
    ! precession Om_ph are evaluated at a fixed pitch (eta = 0.84/Bmin, i.e.
    ! v_par/v is supplied as an optional second command-line argument) and a
    ! fixed energy. No E x B (Om_tE = 0) to match POTATO's
    ! zero-potential profile, so Om_ph is the magnetic precession alone.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: s, inp_swi, do_magfie_init, iota
    use neort, only: init
    use neort_freq, only: Om_th, Om_ph
    use driftorbit, only: mph, mth, vth, etadt, Om_tE, sign_vpar, magdrift

    implicit none

    integer, parameter :: ns = 60
    real(dp), parameter :: smin = 0.05_dp, smax = 0.95_dp
    real(dp), parameter :: v_10keV = 9.787e7_dp ! 10 keV deuteron, cm/s
    real(dp) :: lambda0

    integer :: i, fid
    character(len=32) :: input_switch_arg, lambda_arg
    real(dp) :: v, eta, Omth, dOmthdv, dOmthdeta, Omph, dOmphdv, dOmphdeta

    inp_swi = 9 ! default: ASDEX Upgrade Boozer format
    lambda0 = 0.4_dp
    call get_command_argument(1, input_switch_arg)
    if (len_trim(input_switch_arg) > 0) read(input_switch_arg, *) inp_swi
    call get_command_argument(2, lambda_arg)
    if (len_trim(lambda_arg) > 0) read(lambda_arg, *) lambda0
    if (abs(lambda0) >= 1.0_dp) error stop "midplane pitch must satisfy |lambda|<1"
    magdrift = .true.
    s = smin
    call do_magfie_init("in_file")

    mth = 1
    mph = 2
    vth = v_10keV
    v = v_10keV
    Om_tE = 0.0_dp ! no E x B, magnetic precession only
    sign_vpar = 1.0_dp

    open(newunit=fid, file='freq_scan_neort.dat')
    write(fid, '(A)') '# s rho_tor omega_b[1/s] omega_phi[1/s] iota'
    do i = 1, ns
        s = smin + (smax - smin)*real(i - 1, dp)/real(ns - 1, dp)
        call init
        eta = (1.0_dp - lambda0**2)*etadt ! = 0.84/Bmin, lambda=0.4 at midplane
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        write(fid, '(5ES16.7)') s, sqrt(s), Omth, Omph, iota
    end do
    close(fid)
    print *, 'freq_scan_neort.dat written, ns =', ns

end program test_freq_scan
