program test_bounce_program
    use iso_fortran_env, only: dp => real64
    use util, only: qe, mu

    implicit none

    real(dp), parameter :: tol = 1.0e-13_dp
    real(dp) :: v, eta

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
        call do_magfie_init("in_file")
        call init

        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0

        etamin = (1+epst)*etatp
        etamax = (1-epst)*etadt
        sign_vpar = 1

        v = vth
        eta = 0.5_dp*(etamin + etamax)

    end subroutine setup

    subroutine setup_control
        use driftorbit, only: s, M_t, qi, mi, vth, epsmn, m0, &
            mph, mth, magdrift, nopassing, pertfile, &
            nonlin, bfac, efac, inp_swi
        use neort_orbit, only: noshear
        real(dp) :: qs, ms

        s = 0.153_dp
        M_t = 0.1_dp
        qs = 1.0_dp
        ms = 2.014_dp
        vth = 37280978.0_dp
        epsmn = 1.0e-3_dp
        m0 = 0
        mph = 18
        mth = -1
        magdrift = .false.
        nopassing = .false.
        noshear = .true.
        pertfile = .false.
        nonlin = .false.
        bfac = 1.0_dp
        efac = 1.0_dp
        inp_swi = 8

        M_t = M_t*efac/bfac
        qi = qs*qe
        mi = ms*mu
    end subroutine setup_control

    subroutine test_bounce
        use neort_orbit, only: bounce, nvar

        real(dp) :: taub, bounceavg(nvar), bounceavg_tmp(nvar)

        bounceavg = 0.0_dp
        call bounce(v, eta, taub, bounceavg_tmp)

        bounceavg = 0.0_dp
        call bounce(v, eta, taub, bounceavg)

        if (maxval(abs(bounceavg - bounceavg_tmp)) > tol) then
            print *, 'test_bounce failed'
            error stop
        end if

        print *, 'test_bounce OK'
    end subroutine test_bounce


    subroutine test_bounce_fast
        use neort_orbit, only: bounce, bounce_fast, nvar, timestep
        real(dp) :: taub, bounceavg(nvar), bounceavg_tmp(nvar)

        bounceavg = 0.0_dp

        bounceavg = 0.0_dp
        call bounce(v, eta, taub, bounceavg_tmp, taub)

        bounceavg = 0.0_dp
        call bounce_fast(v, eta, taub, bounceavg, timestep)

        if (maxval(abs((bounceavg - bounceavg_tmp)/(bounceavg+1.0e-6_dp))) &
            > 1.0e-3_dp) then
            print *, 'test_bounce_fast failed'
            error stop
        end if

        print *, 'test_bounce_fast OK'
    end subroutine test_bounce_fast
end program test_bounce_program
