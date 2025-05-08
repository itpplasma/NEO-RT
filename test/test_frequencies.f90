program test_frequencies
  use util
  use do_magfie_mod, only: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds,&
  q, dqds, iota, R0, a, eps, inp_swi, do_magfie_init, do_magfie
  use neort, only: init
  use neort_orbit, only: Om_th, Om_ph
  use driftorbit, only: epsmn, mph, mth, vth, qi, mi, &
    etadt, etatp, Om_tE, sigv

  implicit none

  integer, parameter :: neta = 1000             ! Steps in eta
  real(8), parameter :: M_t = 0.129884652289d0  ! Unscaled Mach number

  integer :: i, fid
  real(8) :: Omth, dOmthdv, dOmthdeta, Omph, dOmphdv, dOmphdeta
  real(8), parameter :: scalfac_energy = 100.0d0
  real(8), parameter :: scalfac_efield = 100.0d0

  real(8) :: v, eta

  inp_swi = 9  ! ASDEX Upgrade format
  s = 0.17     ! Normalized toroidal flux

  !print *, s
  call do_magfie_init

  mth = 1                  ! Canonical poloidal harmonic
  mph = 2                  ! Canonical toroidal harmonic
  vth = 44299218.58052201  ! Unscaled thermal velocity
  Om_tE = vth*M_t/(R0*scalfac_efield)  ! ExB precession frequency
  vth = vth/sqrt(scalfac_energy)

  call init

  v = vth

  open(newunit=fid, file='canonical_freqs_vs_eta_t.dat')
  write(fid,*)  '%# eta [1/G], omega_b [rad/s], Omega_tor [rad/s] trapped'
  sigv = 1d0
  do i = 1,neta
    eta = etatp + i*(etadt-etatp)/neta
    call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
    call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
    write(fid,*) eta, Omth, Omph
  end do
  close(fid)

  open(newunit=fid, file='canonical_freqs_vs_eta_pco.dat')
  write(fid,*)  '%# eta [1/G], omega_b [rad/s], Omega_tor [rad/s] co-passing'
  sigv = 1d0
  do i = 0,neta-1
    eta = i*etatp/neta
    call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
    call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
    write(fid,*) eta, Omth, Omph
  end do
  close(fid)

  open(newunit=fid, file='canonical_freqs_vs_eta_pct.dat')
  write(fid,*)  '%# eta [1/G], omega_b [rad/s], Omega_tor [rad/s] co-passing'
  sigv = -1d0
  do i = 0,neta-1
    eta = i*etatp/neta
    call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
    call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
    write(fid,*) eta, -Omth, Omph
  end do
  close(fid)

end program test_frequencies
