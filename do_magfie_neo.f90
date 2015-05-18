module do_magfie_mod
  use common
  use neo_exchange, only: nper,b_min,b_max, &
       theta_bmin,theta_bmax,phi_bmin,phi_bmax,rt0
  use magfie_mod, only: magfie, stevvo, magfie_deallocate
  use neo_magfie_mod, only: magfie_spline, magfie_sarray, boozer_curr_tor_hat,&
       boozer_curr_pol_hat, boozer_curr_tor_hat_s, boozer_curr_pol_hat_s,&
       boozer_psi_pr_hat, boozer_sqrtg11, boozer_isqrg, boozer_iota
  use nrtype, only: twopi
  USE partpa_mod,  ONLY : bmod0
  use neo_magfie_mod, only: boozer_curr_tor_hat, boozer_curr_pol_hat,&
       boozer_curr_tor_hat_s, boozer_curr_pol_hat_s, boozer_psi_pr_hat,&
       boozer_sqrtg11, boozer_isqrg
  
  implicit none
  save

  real(8) :: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds, q, iota, R0, eps

contains
  
  subroutine do_magfie_init
    real(8) :: bmod, bmod1, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)

    bmod0 = 1d-4
    magfie_spline = 1
    
    allocate(magfie_sarray(1))
    magfie_sarray = s

    x(1) = s
    x(2:3) = 0.0d0
    call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    R0 = 1d2*rt0
    x(3) = pi
    call do_magfie(x, bmod1, sqrtg, bder, hcovar, hctrvr, hcurl)
    eps = (bmod1/bmod-1d0)/(bmod1/bmod+1d0)
    
  end subroutine do_magfie_init
  
  subroutine do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    real(8), dimension(:),       intent(in)         :: x
    real(8),                     intent(out)        :: bmod
    real(8),                     intent(out)        :: sqrtg
    real(8), dimension(size(x)), intent(out)        :: bder
    real(8), dimension(size(x)), intent(out)        :: hcovar
    real(8), dimension(size(x)), intent(out)        :: hctrvr
    real(8), dimension(size(x)), intent(out)        :: hcurl
    call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    bmod      = 1d4*bmod
    sqrtg     = abs(sqrtg)
    psi_pr    = boozer_psi_pr_hat
    Bthcov    = boozer_curr_tor_hat
    Bphcov    = boozer_curr_pol_hat
    dBthcovds = boozer_curr_tor_hat_s
    dBphcovds = boozer_curr_pol_hat_s
    iota      = boozer_iota
    q         = 1d0/iota
    R0        = 1d2*rt0
  end subroutine do_magfie
  
end module do_magfie_mod
