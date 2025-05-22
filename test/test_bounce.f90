program test_bounce_program
    use util, only: qe, mu

    implicit none

    real(8), parameter :: tol = 1e-13
    real(8) :: v, eta

    call setup
    call test_bounce
    call test_bounce_fast

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

        s = 0.153d0
        M_t = 0.1d0
        qs = 1.0d0
        ms = 2.014d0
        vth = 37280978.0d0
        epsmn = 1e-3
        m0 = 0
        mph = 18
        mth = -1
        magdrift = .false.
        nopassing = .false.
        noshear = .true.
        pertfile = .false.
        nonlin = .false.
        bfac = 1.0d0
        efac = 1.0d0
        inp_swi = 8

        M_t = M_t*efac/bfac
        qi = qs*qe
        mi = ms*mu
    end subroutine setup_control

    subroutine test_bounce
        use neort_orbit, only: bounce, nvar

        real(8) :: taub, bounceavg(nvar), bounceavg_tmp(nvar)

        bounceavg = 0d0
        call bounce(v, eta, taub, bounceavg_tmp)

        bounceavg = 0d0
        call bounce(v, eta, taub, bounceavg)

        if (maxval(abs(bounceavg - bounceavg_tmp)) > tol) then
            print *, 'test_bounce failed'
            error stop
        end if

        print *, 'test_bounce OK'
    end subroutine test_bounce


    subroutine test_bounce_fast
        use neort_orbit, only: bounce, bounce_fast, nvar, timestep
        real(8) :: taub, bounceavg(nvar), bounceavg_tmp(nvar)

        bounceavg = 0d0

        bounceavg = 0d0
        call bounce(v, eta, taub, bounceavg_tmp, taub)

        bounceavg = 0d0
        call bounce_fast(v, eta, taub, bounceavg, timestep)

        if (maxval(abs((bounceavg - bounceavg_tmp)/(bounceavg+1d-6))) &
            > 1d-3) then
            print *, 'test_bounce_fast failed'
            error stop
        end if

        print *, 'test_bounce_fast OK'
    end subroutine test_bounce_fast
end program test_bounce_program
