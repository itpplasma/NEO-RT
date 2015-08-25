module do_magfie_mod
  use common
  use neo_exchange, only: nper,b_min,b_max, &
       theta_bmin,theta_bmax,phi_bmin,phi_bmax,rt0
  use magfie_mod, only: magfie, stevvo, magfie_deallocate
  use neo_magfie_mod, only: magfie_spline, magfie_sarray, boozer_curr_tor_hat,&
       boozer_curr_pol_hat, boozer_curr_tor_hat_s, boozer_curr_pol_hat_s,&
       boozer_psi_pr_hat, boozer_sqrtg11, boozer_isqrg, boozer_iota,&
       boozer_iota_s
  use nrtype, only: twopi
  use partpa_mod,  ONLY : bmod0
  use neo_magfie_mod, only: boozer_curr_tor_hat, boozer_curr_pol_hat,&
       boozer_curr_tor_hat_s, boozer_curr_pol_hat_s, boozer_psi_pr_hat,&
       boozer_sqrtg11, boozer_isqrg
  use neo_input, only: bmnc
  
  implicit none
  save

  real(8) :: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds, q, dqds,&
       iota, R0, eps
  ! B0h is the 0th theta harmonic of bmod on current flux surface
  ! and B00 the 0th theta harmonic of bmod on the innermost flux surface

  real(8), parameter :: a = 4.6d1 ! TODO 1: make minor radius a changeable

contains
  
  subroutine do_magfie_init
    real(8) :: bmod, bmod1, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)

    !! VERY DANGEROUS: CHANGING Q MANUALLY
    !boozer_iota      = 1d0/3d0
    
    bmod0 = 1d-4
    !print *, bmod0
    magfie_spline = 1

    if (.not. allocated(magfie_sarray)) then
       allocate(magfie_sarray(1))
       magfie_sarray = s
    end if

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
    !! VERY DANGEROUS: CHANGING Q MANUALLY
    boozer_iota      = 1d0/3d0
    call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
    bmod      = 1d4*bmod
    sqrtg     = abs(sqrtg)
    psi_pr    = boozer_psi_pr_hat
    Bthcov    = boozer_curr_tor_hat
    Bphcov    = boozer_curr_pol_hat
    dBthcovds = boozer_curr_tor_hat_s
    dBphcovds = boozer_curr_pol_hat_s
    iota      = boozer_iota
    !! VERY DANGEROUS: CHANGING Q MANUALLY
    !iota      = 1d0/3d0
    q         = 1d0/iota
    dqds      = -boozer_iota_s/iota**2
    R0        = 1d2*rt0
  end subroutine do_magfie
  
end module do_magfie_mod
