program main
    use neort, only: read_control, check_magfie, init_profile, init_plasma, init_run, &
                     compute_transport_harmonic, runname, s, M_t
    use driftorbit
    implicit none

    character(1024) :: TEST_RUN = "driftorbit64.new"
    real(8) :: ds = 0.001d0

    runname = TEST_RUN

    call read_control

    call setup
    call test_omega_prime

    s = s + 0.5d0*ds
    call setup
    call test_omega_prime

    s = s - ds
    call setup
    call test_omega_prime

contains

    subroutine setup
        call do_magfie_init
        call init_profile
        call init_plasma
        call init_run(use_thermodynamic_profiles=.true.)
        call check_magfie

        mth = -3
        sigv = 1
        etamin = (1 + epst)*etatp
        etamax = (1 - epst)*etadt
    end subroutine setup

    subroutine test_omega_prime
        real(8) :: ux
        real(8) :: Om, dOmdv, dOmdeta, Ompr_old, Ompr_new, dOmdpph, Omth, dOmthdeta
        real(8) :: roots(nlev, 3), eta_res(2)
        integer :: nroots

        ux = 0.9d0
        v = ux*vth

        call driftorbit_coarse(etamin, etamax, roots, nroots)
        eta_res = driftorbit_root(1d-8*abs(Om_tE), roots(1, 1), roots(1, 2))

        call compute_frequencies(ux, eta_res, Omth, Om, dOmdv, dOmdeta, dOmdpph)
        Ompr_old = omega_prime(ux, Omth, dOmdv, dOmdeta, dOmdpph)
        Ompr_new = omega_prime_new(ux, Omth, dOmdv, dOmdeta, dOmdpph)

        print *, s, v, eta, Om, Ompr_old, Ompr_new

    end subroutine test_omega_prime

    subroutine compute_frequencies(ux, etax, Omth, Om, dOmdv, dOmdeta, dOmdpph)
        real(8), intent(in) :: ux, etax(2)
        real(8), intent(out) :: Omth, Om, dOmdv, dOmdeta, dOmdpph

        real(8) :: Omph, dOmphdv, dOmphdeta, dOmphds, dOmthds, dOmthdv, dOmthdeta

        v = ux*vth
        eta = etax(1)

        call Om_th(Omth, dOmthdv, dOmthdeta)
        taub = 2d0*pi/abs(Omth)
        call bounce_fast
        call Om_ph(Omph, dOmphdv, dOmphdeta)
        call d_Om_ds(dOmthds, dOmphds)
        Om = mth*Omth + mph*Omph
        dOmdv = mth*dOmthdv + mph*dOmphdv
        dOmdeta = mth*dOmthdeta + mph*dOmphdeta
        dOmdpph = -(qi/c*iota*psi_pr)**(-1)*(mth*dOmthds + mph*dOmphds)
    end subroutine compute_frequencies

end program main
