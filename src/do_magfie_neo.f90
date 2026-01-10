module do_magfie_mod
    use util
    use neo_exchange, only: nper, b_min, b_max, &
                            theta_bmin, theta_bmax, phi_bmin, phi_bmax, rt0
    use magfie_mod, only: magfie, stevvo, magfie_deallocate
    use neo_magfie, only: magfie_spline, magfie_sarray, boozer_curr_tor_hat, &
                              boozer_curr_pol_hat, boozer_curr_tor_hat_s, boozer_curr_pol_hat_s, &
                              boozer_psi_pr_hat, boozer_sqrtg11, boozer_isqrg, boozer_iota, &
                              boozer_iota_s, inp_swi
    use nrtype, only: twopi
    use partpa_mod, ONLY: bmod0
    use neo_magfie, only: boozer_curr_tor_hat, boozer_curr_pol_hat, &
                              boozer_curr_tor_hat_s, boozer_curr_pol_hat_s, boozer_psi_pr_hat, &
                              boozer_sqrtg11, boozer_isqrg
    use neo_input, only: bmnc

    implicit none

    real(8), parameter :: sign_theta = -1.0d0  ! negative for left-handed

    real(8) :: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds, q, dqds, &
               iota, R0, eps, bfac ! TODO: implement bfac
    ! B0h is the 0th theta harmonic of bmod on current flux surface
    ! and B00 the 0th theta harmonic of bmod on the innermost flux surface

    real(8), parameter :: a = 4.6d1 ! TODO 1: make minor radius changeable

    !$omp threadprivate (s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds)
    !$omp threadprivate (q, dqds, iota, R0, eps, bfac)

contains

    subroutine do_magfie_init
        real(8) :: bmod, bmod1, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)

        print *, "start do_magfie_init"

        bmod0 = 1d-4
        magfie_spline = 1

        if (.not. allocated(magfie_sarray)) then
            allocate (magfie_sarray(1))
            magfie_sarray = s
        end if

        x(1) = s
        x(2:3) = 0.0d0
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        R0 = 1d2*rt0
        print *, R0
        x(3) = pi
        call do_magfie(x, bmod1, sqrtg, bder, hcovar, hctrvr, hcurl)
        eps = (bmod1/bmod - 1d0)/(bmod1/bmod + 1d0) ! TODO: make this correct

    end subroutine do_magfie_init

    subroutine do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        real(8), dimension(:), intent(in) :: x
        real(8), intent(out) :: bmod
        real(8), intent(out) :: sqrtg
        real(8), dimension(size(x)), intent(out) :: bder
        real(8), dimension(size(x)), intent(out) :: hcovar
        real(8), dimension(size(x)), intent(out) :: hctrvr
        real(8), dimension(size(x)), intent(out) :: hcurl
        call magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        bmod = 1d4*bmod
        psi_pr = boozer_psi_pr_hat
        Bthcov = boozer_curr_tor_hat
        Bphcov = boozer_curr_pol_hat
        dBthcovds = boozer_curr_tor_hat_s
        dBphcovds = boozer_curr_pol_hat_s
        iota = boozer_iota
        q = 1d0/iota
        dqds = -boozer_iota_s/iota**2
        R0 = 1d2*rt0
        ! set B_r to zero for now
        hcovar(1) = 0
    end subroutine do_magfie

end module do_magfie_mod

module do_magfie_pert_mod
    use neo_magfie_perturbation, only: neo_read_pert, &
                                       neo_init_spline_pert, m_phi, &
                                       neo_magfie_pert_amp
    use neo2_ql, only: read_in_namelists, set_default_values, init

    integer :: mph

contains

    subroutine do_magfie_pert_init
        integer :: iunit

        call init
        open(newunit=iunit,file='neo2.in',status='old')
        call read_in_namelists(iunit)
        close(iunit)

        call neo_read_pert
        call neo_init_spline_pert
        mph = nint(m_phi)
    end subroutine do_magfie_pert_init

    !
  ! Compute the perturbation field amplitude
  ! for constant toroidal modenumber n = ixm_pert(1)
  SUBROUTINE neo_magfie_pert_amp( x, bn_hat_pert )
    !
    ! input / output
    !
    REAL(8), DIMENSION(:), INTENT(in) :: x
    COMPLEX(8), INTENT(out) :: bn_hat_pert
    !
    ! local definitions
    !
    COMPLEX(8), PARAMETER :: imun=(0.0_dp,1.0_dp)
    INTEGER :: swd
    INTEGER :: i, m
    REAL(8) :: yp, ypp, yppp
    REAL(8) :: bmnc_pert_val, bmns_pert_val
    COMPLEX(8) :: expv
    !
    ! read Boozer file and prepare spline routines (1st call)
    !
    IF (.NOT. ALLOCATED(es_pert)) THEN
       CALL neo_read_pert()
       CALL neo_init_spline_pert()
    END IF

    bn_hat_pert = (0.0_dp,0.0_dp)
    DO i = 1, mnmax_pert
       swd = 1
       CALL splint_horner3(es_pert, a_bmnc_pert(:,i), b_bmnc_pert(:,i),&
            c_bmnc_pert(:,i), d_bmnc_pert(:,i), swd, r_mhalf_pert(i),&
            x(1), tf, tfp, tfpp, tfppp, bmnc_pert_val, yp, ypp, yppp)

       ! Additional data from Boozer files without Stellarator symmetry
       IF (inp_swi .EQ. 9) THEN ! ASDEX-U (E. Strumberger)
          CALL splint_horner3(es_pert, a_bmns_pert(:,i), b_bmns_pert(:,i),&
               c_bmns_pert(:,i), d_bmns_pert(:,i), swd, r_mhalf_pert(i),&
               x(1), tf, tfp, tfpp, tfppp, bmns_pert_val, yp, ypp, yppp)
       END IF

       IF (inp_swi .EQ. 8) THEN ! NEW IPP TOKAMAK
          m = (-1)*ixm_pert(i)
             expv = EXP(imun*m*x(3))
             bn_hat_pert = bn_hat_pert + bmnc_pert_val * expv
       ELSEIF (inp_swi .EQ. 9) THEN ! ASDEX-U (E. Strumberger)
             expv = EXP(imun*m*x(3))
             bn_hat_pert = bn_hat_pert + (bmnc_pert_val-imun*bmns_pert_val)*expv
       END IF

    END DO
    bn_hat_pert = bn_hat_pert / bmod0
  END SUBROUTINE neo_magfie_pert_amp
end module do_magfie_pert_mod
