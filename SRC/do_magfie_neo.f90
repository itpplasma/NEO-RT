module do_magfie_mod
    use common
    use neo_exchange, only: nper, b_min, b_max, &
                            theta_bmin, theta_bmax, phi_bmin, phi_bmax, rt0
    use magfie_mod, only: magfie, stevvo, magfie_deallocate
    use neo_magfie_mod, only: magfie_spline, magfie_sarray, boozer_curr_tor_hat, &
                              boozer_curr_pol_hat, boozer_curr_tor_hat_s, boozer_curr_pol_hat_s, &
                              boozer_psi_pr_hat, boozer_sqrtg11, boozer_isqrg, boozer_iota, &
                              boozer_iota_s, inp_swi
    use nrtype, only: twopi
    use partpa_mod, ONLY: bmod0
    use neo_magfie_mod, only: boozer_curr_tor_hat, boozer_curr_pol_hat, &
                              boozer_curr_tor_hat_s, boozer_curr_pol_hat_s, boozer_psi_pr_hat, &
                              boozer_sqrtg11, boozer_isqrg
    use neo_input, only: bmnc

    implicit none
    save

    real(8) :: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds, q, dqds, &
               iota, R0, eps, bfac ! TODO: implement bfac
    ! B0h is the 0th theta harmonic of bmod on current flux surface
    ! and B00 the 0th theta harmonic of bmod on the innermost flux surface

    real(8), parameter :: a = 4.6d1 ! TODO 1: make minor radius changeable

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
        print *, "before do_magfie"
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        print *, "after do_magfie"
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
        sqrtg = abs(sqrtg)
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
        hctrvr(1) = 0
    end subroutine do_magfie

end module do_magfie_mod

module do_magfie_pert_mod
    use neo_magfie_perturbation, only: neo_read_pert_control, neo_read_pert, &
                                       neo_init_spline_pert, neo_magfie_pert_amp, m_phi

    real(8) :: mph

contains

    subroutine do_magfie_pert_init
        call neo_read_pert_control
        call neo_read_pert
        call neo_init_spline_pert
        mph = m_phi
    end subroutine do_magfie_pert_init

    subroutine do_magfie_pert_amp(x, bamp)
        real(8) :: x(3)
        complex(8) :: bamp
        call neo_magfie_pert_amp(x, bamp)
    end subroutine do_magfie_pert_amp
end module do_magfie_pert_mod
