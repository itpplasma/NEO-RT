program test_frequencies
  use common
  use do_magfie_mod, only: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds,&
  q, dqds, iota, R0, a, eps, inp_swi, do_magfie_init, do_magfie
  use driftorbit, only: init, epsmn, mph, mth, vth, qi, mi, Jperp, &
    eta, etadt, etatp, v, Om_th, Om_ph, Om_tE

  implicit none

  integer, parameter :: neta = 100     ! Steps in eta
  real(8), parameter :: M_t = 0.036d0  ! Mach number

  integer :: i
  real(8) :: Omth, dOmthdv, dOmthdeta, Omph, dOmphdv, dOmphdeta

  inp_swi = 9  ! ASDEX Upgrade format
  s = 0.5      ! Normalized toroidal flux

  !print *, s
  call do_magfie_init

  mth = 1             ! Canonical poloidal harmonic
  mph = 2             ! Canonical toroidal harmonic
  vth = 40000000.0d0  ! Thermal velocity
  Om_tE = vth*M_t/R0  ! ExB precession frequency

  call init

  v = vth
  do i = 1,neta
    eta = etatp + i*(etadt-etatp)/neta
    call Om_th(Omth, dOmthdv, dOmthdeta)
    call Om_ph(Omph, dOmphdv, dOmphdeta)
    write(99,*) Jperp(), Omth, Omph
  end do

end program test_frequencies
