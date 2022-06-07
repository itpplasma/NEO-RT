program test_frequencies
  use do_magfie_mod, only: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds,&
  q, dqds, iota, R0, a, eps, inp_swi, do_magfie_init
  use driftorbit, only: init

  implicit none

  inp_swi = 9
  s = 0.5

  call do_magfie_init

  call init

end program test_frequencies
