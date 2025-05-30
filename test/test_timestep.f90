program test_timestep
    use util_for_test, only: print_test, pass_test, fail_test
    use util, only: qe, mu
    implicit none

    real(8), parameter :: tol = 1e-13
    real(8) :: v, eta

    call setup
    call test_timestep_poloidal()

contains

    subroutine setup
        use neort, only: init
        use driftorbit, only: do_magfie_init, etamin, etamax, &
            Om_tE, dOm_tEds, etatp, etadt, epst, sign_vpar, vth, M_t, dM_tds
        use do_magfie_mod, only: R0

        call setup_control
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

    subroutine setup_control
        use driftorbit, only: s, M_t, qi, mi, vth, epsmn, m0, &
            mph, mth, magdrift, nopassing, pertfile, &
            nonlin, bfac, efac, inp_swi
        use neort_orbit, only: noshear
        real(8) :: qs, ms

        ! Use same configuration as test/ripple_plateau
        s = 3.9063d-2
        M_t = 0.1d0
        qs = 1.0d0
        ms = 2.014d0
        vth = 1.0d8
        epsmn = 1.0d-3
        m0 = 0
        mph = 18
        mth = -1
        magdrift = .false.
        nopassing = .false.
        noshear = .false.
        pertfile = .false.
        nonlin = .false.
        bfac = 1.0d0
        efac = 1.0d0
        inp_swi = 8

        M_t = M_t*efac/bfac
        qi = qs*qe
        mi = ms*mu
    end subroutine setup_control

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

        ! Evaluate magnetic field at current position
        x(1) = s
        x(2) = 0d0
        x(3) = 0d0  ! Start at theta = 0
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        ! Calculate parallel velocity
        vpar_val = sqrt(2d0*v**2*eta*bmod)

        ! Initialize starting conditions
        y(1) = 0d0      ! Starting poloidal angle theta
        y(2) = vpar_val ! Starting parallel velocity vpar
        t = 0d0         ! Starting time
        dt = 1d-6       ! Small timestep

        ! Test single timestep of timestep_poloidal_motion
        call timestep_poloidal_motion(v, eta, neq, t, y, ydot)

        ! Basic validation: check that derivatives are reasonable
        ! ydot(1) should be theta_dot (poloidal velocity)
        ! ydot(2) should be vpar_dot (parallel acceleration)
        if (abs(ydot(1)) > 0d0 .and. abs(ydot(1)) < 1d10 .and. &
            abs(ydot(2)) < 1d20) then
            call pass_test()
        else
            call fail_test()
        end if

    end subroutine test_timestep_poloidal

end program test_timestep
