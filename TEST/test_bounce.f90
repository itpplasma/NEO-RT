program test_bounce_program
    use common, only: qe, mu

    implicit none

    call setup
    call test_bounce

    contains

    subroutine setup
        use driftorbit, only: do_magfie_init, init, etamin, etamax, &
            Om_tE, dOm_tEds, etatp, etadt, epst, sigv, vth, M_t, dM_tds
        use do_magfie_mod, only: R0

        call setup_control
        call do_magfie_init
        call init

        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0

        etamin = (1+epst)*etatp
        etamax = (1-epst)*etadt
        sigv = 1

    end subroutine setup

    subroutine setup_control
        use driftorbit, only: s, M_t, qi, mi, vth, epsmn, m0, &
            mph, mth, supban, magdrift, nopassing, noshear, pertfile, &
            nonlin, bfac, efac, inp_swi, orbit_mode_avg, &
            orbit_mode_transp, vsteps, comptorque, intoutput

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
        supban = .false.
        magdrift = .false.
        nopassing = .false.
        noshear = .true.
        pertfile = .false.
        nonlin = .false.
        bfac = 1.0d0
        efac = 1.0d0
        inp_swi = 8
        orbit_mode_avg = 0
        orbit_mode_transp = 0
        vsteps = -1
        comptorque = .false.
        intoutput = .false.

        M_t = M_t*efac/bfac
        qi = qs*qe
        mi = ms*mu
    end subroutine setup_control

    subroutine test_bounce
        use driftorbit, only: bounce

        call bounce
    end subroutine test_bounce
end program test_bounce_program