program test_timestep
    use util_for_test, only: print_test, pass_test, fail_test
    use util, only: qe, mu
    implicit none

    real(8), parameter :: tol = 1e-13
    real(8) :: v, eta

    call setup
    call test_timestep_poloidal

contains

    subroutine setup
        use neort, only: init
        use driftorbit, only: do_magfie_init, s, inp_swi, etamin, etamax, &
            Om_tE, dOm_tEds, etatp, etadt, epst, sign_vpar, vth, M_t, dM_tds
        use do_magfie_mod, only: R0

        inp_swi = 8
        s = 3.9063d-2
        vth = 1.0d8
        M_t = 0.1d0

        call do_magfie_init
        call init

        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0

        etamin = (1+epst)*etatp
        etamax = (1-epst)*etadt
        sign_vpar = 1

        v = vth
        eta = 0.5d0*(etamin + etamax)
    end subroutine setup

    subroutine test_timestep_poloidal()
        use neort_orbit, only: timestep_poloidal_motion
        use do_magfie_mod, only: do_magfie, s

        integer, parameter :: neq = 2
        real(8) :: y(neq), ydot(neq)
        real(8) :: t, dt
        real(8) :: bmod, sqrtg
        real(8), dimension(3) :: x, hder, hcovar, hctrvr, hcurl
        real(8) :: vpar_val

        call print_test("test_timestep_poloidal")

        x(1) = s
        x(2) = 0d0
        x(3) = 0d0
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        vpar_val = sqrt(2d0*v**2*eta*bmod)

        y(1) = 0d0
        y(2) = vpar_val
        t = 0d0
        dt = 1d-6

        call timestep_poloidal_motion(v, eta, neq, t, y, ydot)

        if (abs(ydot(1)) > 0d0 .and. abs(ydot(1)) < 1d10 .and. &
            abs(ydot(2)) < 1d20) then
            call pass_test()
        else
            call fail_test()
        end if

    end subroutine test_timestep_poloidal

end program test_timestep
