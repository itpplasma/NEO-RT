program test_line_drive
    !! drive_form = 'line' on the examples/base Boozer equilibrium (noshear off)
    !! with the ideal_helical perturbation.  Oracles independent of the drive
    !! code: finite differences of do_magfie (collapse condition, div dB = 0),
    !! p_phi conservation (radial drift), NEO-RT's own precession Om_tB, the
    !! Boozer scalar route, and the exact by-parts value of a gauge term.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie_init, R0, s, q, do_magfie, Bphcov, iota, &
        psi_pr, sign_theta
    use driftorbit, only: efac, bfac, mth, m0, epsmn, etatp, etadt, sign_vpar
    use neort, only: init
    use neort_config, only: read_and_set_config
    use neort_orbit, only: noshear, vpar
    use neort_freq, only: Om_th, Om_tB
    use neort_profiles, only: init_profiles, read_and_init_plasma_input, &
        read_and_init_profile_input, vth
    use neort_line_drive, only: set_line_drive_options, line_bounce, &
        line_harmonics_t, exact_orbit_period, gauge_amp, collapse_eta, &
        local_field, local_field_t, drift_velocity, pert_eval, pert_point_t
    use do_magfie_pert_mod, only: mph
    use util, only: pi, imun, mi, qi, c

    implicit none

    integer :: nfail = 0
    complex(dp) :: offres_prediction = 0
    type(line_harmonics_t) :: rb, rg, rr, rgen

    call read_and_set_config("driftorbit.in")
    noshear = .false.
    m0 = 2
    epsmn = 1.0e-3_dp
    call set_line_drive_options('line', 'ideal_helical', .false., 8, noshear)
    call do_magfie_init("in_file")
    call init_profiles(R0)
    call read_and_init_plasma_input("plasma.in", s)
    call read_and_init_profile_input("profile.in", s, R0, efac, bfac)
    call init

    call check_pointwise_identities()

    ! trapped, bounce harmonic m_b = 0 on its resonance n (Om_tE + Om_tB) = 0
    call resonant_harmonics(0.6_dp, 0, 1.0_dp, 0.0_dp, rb, rg)
    call check('collapse |H_line-H_boozer|/|H_boozer|', rel(rb%h_line, rb%h_boozer), &
        1.0e-9_dp)
    call check('line vs flux-surface drive h_FS', rel(rb%h_line, rb%h_fs), 1.0e-9_dp)
    call check('gauge chi (|grad chi|~10|dA|) on resonance', &
        rel(rg%h_line, rb%h_line), 1.0e-6_dp)
    call check('negative: naive thin dA.Xdot is gauge dependent', &
        rel(rg%h_naive, rb%h_naive), 1.0_dp, .true.)
    ! double count: symplectic term kept and the full Boozer H_1 added on top
    call check('negative: double count detected', &
        rel(rb%h_line - rb%mu_dbe + rb%h_boozer, rb%h_boozer), 0.1_dp, .true.)
    call check('double count error = m vpar^2 dB_L/B + mu xi.grad B', &
        rel(rb%h_line - rb%mu_dbe, rb%p_par + rb%m_xi), 1.0e-9_dp)
    call check('negative: Eulerian dB_E in the Boozer weight detected', &
        rel(rb%h_boozer_euler, rb%h_boozer), 0.01_dp, .true.)
    call check('line reproduces the xi.grad|B| term of the Boozer weight', &
        abs(rb%h_line - rb%h_boozer)/abs(rb%h_boozer - rb%h_boozer_euler), &
        1.0e-9_dp)
    call check('negative: Lagrangian mu dB_L in the line form detected', &
        abs(rb%m_xi)/abs(rb%h_boozer), 0.01_dp, .true.)

    ! off resonance the raw path term shifts by exactly i (m.Omega) (e/c) chi_m
    call resonant_harmonics(0.6_dp, 0, 1.0_dp, 0.05_dp, rr, rg)
    call check('off-resonance gauge shift = i (m.Omega) chi_m', &
        rel(rg%path - rr%path, offres_prediction), 1.0e-6_dp)
    call check('gauge-fixed H_line unchanged off resonance', &
        rel(rg%h_line, rr%h_line), 1.0e-6_dp)

    ! a generic, non-Boozer displacement: line = h_FS, Boozer weight O(1) off
    collapse_eta = .false.
    call resonant_harmonics(0.6_dp, 0, 1.0_dp, 0.0_dp, rgen, rg)
    collapse_eta = .true.
    call check('generic xi: line vs h_FS', rel(rgen%h_line, rgen%h_fs), 1.0e-9_dp)
    call check('negative: generic xi with the Boozer weight is O(1) off', &
        rel(rgen%h_boozer, rgen%h_line), 0.05_dp, .true.)

    ! passing with a small path rotation and trapped m_b = 1: agreement up to
    ! terms of the order of (Om_E R / v)^2 that NEO-RT's orbits neglect
    call resonant_harmonics(0.4_dp, -4, 1.0_dp, 0.0_dp, rr, rg)
    call check('passing m_b=-4: |H_line-H_boozer|/|H_boozer|', &
        rel(rr%h_line, rr%h_boozer), 1.0e-4_dp)
    call check('passing m_b=-4 gauge on resonance', rel(rg%h_line, rr%h_line), &
               1.0e-5_dp)
    call resonant_harmonics(0.6_dp, 1, 1.0_dp, 0.0_dp, rr, rg)
    call check('trapped m_b=1: |H_line-H_boozer|/|H_boozer|', &
        rel(rr%h_line, rr%h_boozer), 1.0e-4_dp)

    if (nfail > 0) then
        print *, 'test_line_drive: FAILED checks:', nfail
        error stop 1
    end if
    print *, 'test_line_drive: all checks passed'

contains

    real(dp) function rel(a, b)
        complex(dp), intent(in) :: a, b
        rel = abs(a - b)/abs(b)
    end function rel

    subroutine check(label, value, tol, lower_bound)
        character(*), intent(in) :: label
        real(dp), intent(in) :: value, tol
        logical, intent(in), optional :: lower_bound
        logical :: ok

        ok = value <= tol
        if (present(lower_bound)) then
            if (lower_bound) ok = value >= tol
        end if
        if (.not. ok) nfail = nfail + 1
        print '(a,a,es10.3,a,es9.2,a)', merge('  ok   ', '  FAIL ', ok), label//': ', &
            value, ' (bound ', tol, ')'
    end subroutine check

    subroutine resonant_harmonics(frac, mb, sv, detune, res, resg)
        !! Harmonics without and with the gauge shift chi on the orbit at pitch
        !! fraction frac (>0.5 trapped), bounce harmonic mb, with the path
        !! rotation Om_tE + Om_tB on (detune = 0) or off the resonance.  Om_tE is
        !! set so that it and the computed Om_tB add up to the path rotation.
        real(dp), intent(in) :: frac, sv, detune
        integer, intent(in) :: mb
        type(line_harmonics_t), intent(out) :: res, resg
        real(dp) :: eta, v, taub, omth, om_path, om_te, d1, d2, mdotom
        integer :: istat

        sign_vpar = sv
        v = vth
        if (frac > 0.5_dp) then
            eta = etatp + (frac - 0.5_dp)*2.0_dp*(etadt - etatp)
        else
            eta = 2.0_dp*frac*etatp
        end if
        mth = mb
        call Om_th(v, eta, omth, d1, d2)
        call exact_orbit_period(v, eta, 2.0_dp*pi/abs(omth), taub, omth)
        mdotom = mb*omth
        if (eta <= etatp) mdotom = (mb + q*mph)*omth
        om_path = (-mdotom + detune*abs(omth))/mph
        gauge_amp = 0
        call line_bounce(v, eta, taub, omth, om_path, 0.0_dp, res, istat)
        om_te = om_path - res%om_drift
        call line_bounce(v, eta, taub, omth, om_path, om_te, res, istat)
        gauge_amp = (10.0_dp, 3.0_dp)*epsmn
        call line_bounce(v, eta, taub, omth, om_path, om_te, resg, istat)
        gauge_amp = 0
        offres_prediction = imun*(mdotom + mph*om_path)*resg%chi
    end subroutine resonant_harmonics

    subroutine check_pointwise_identities()
        real(dp), parameter :: h = 1.0e-4_dp
        type(local_field_t) :: lf, lsp, lsm, ltp, ltm
        type(line_harmonics_t) :: res
        real(dp) :: th, eta, v, vp, vpdot, vd(3), ddsdt, err_div, err_divb, err_ds
        real(dp) :: omtb, d1, d2, taub, omth, s0
        integer :: j, istat

        s0 = s
        err_div = 0
        err_divb = 0
        err_ds = 0
        v = vth
        eta = etatp + 0.3_dp*(etadt - etatp)
        do j = 0, 7
            th = -2.5_dp + 0.6_dp*j
            call local_field(s0, th, lf)
            err_div = max(err_div, abs(div_xi_fd(s0, th, h) &
                + 2.0_dp*(lf%dBE + lf%xigradB)/lf%bmod) &
                /abs(2.0_dp*(lf%dBE + lf%xigradB)/lf%bmod))
            call local_field(s0 + h, th, lsp)
            call local_field(s0 - h, th, lsm)
            call local_field(s0, th + h, ltp)
            call local_field(s0, th - h, ltm)
            err_divb = max(err_divb, abs((lsp%sgB(1) - lsm%sgB(1))/(2*h) &
                + (ltp%sgB(2) - ltm%sgB(2))/(2*h) &
                + imun*mph*lf%sgB(3))/abs(mph*lf%sgB(3)))
            call local_field(s0, th, lf)
            if (eta*lf%bmod >= 1.0_dp) cycle
            vp = vpar(v, eta, lf%bmod)
            vpdot = -0.5_dp*v**2*eta*lf%hth*lf%dB_dth
            vd = drift_velocity(vp, 0.5_dp*mi*v**2*eta, lf)
            ddsdt = mi*c*Bphcov/(qi*iota*sign_theta*psi_pr) &
                *(vpdot/lf%bmod - vp*lf%dB_dth*vp*lf%hth/lf%bmod**2)
            err_ds = max(err_ds, abs(ddsdt - vd(1))/abs(vd(1)))
        end do
        call check('collapse condition div xi = -2 dB_L/B (FD of sqrt g)', err_div, &
            1.0e-5_dp)
        call check('div dB = 0 (FD of sqrt(g) dB)', err_divb, 1.0e-6_dp)
        call check('radial drift = d/dt of p_phi excursion', err_ds, 1.0e-10_dp)

        sign_vpar = 1.0_dp
        mth = 0
        call Om_th(v, eta, omth, d1, d2)
        call exact_orbit_period(v, eta, 2.0_dp*pi/abs(omth), taub, omth)
        call line_bounce(v, eta, taub, omth, 0.0_dp, 0.0_dp, res, istat)
        call Om_tB(v, eta, omtb, d1, d2)
        call check('orbit-averaged drift vs NEO-RT Om_tB (spline)', &
            abs(res%om_drift/omtb - 1), 5.0e-3_dp)
    end subroutine check_pointwise_identities

    complex(dp) function div_xi_fd(s0, th, h) result(div)
        !! div xi of xi = xi^s e_s + eta e_theta from do_magfie's sqrt(g).
        real(dp), intent(in) :: s0, th, h
        type(pert_point_t) :: pp, pm, p0
        real(dp) :: sgp, sgm, sg0, b, bder(3), hcov(3), hctr(3), hcurl(3)
        complex(dp) :: ds_term, dth_term

        call do_magfie([s0 + h, 0.0_dp, th], b, sgp, bder, hcov, hctr, hcurl)
        call pert_eval(s0 + h, th, pp)
        call do_magfie([s0 - h, 0.0_dp, th], b, sgm, bder, hcov, hctr, hcurl)
        call pert_eval(s0 - h, th, pm)
        ds_term = (sgp*pp%xs - sgm*pm%xs)/(2*h)
        call do_magfie([s0, 0.0_dp, th + h], b, sgp, bder, hcov, hctr, hcurl)
        call pert_eval(s0, th + h, pp)
        call do_magfie([s0, 0.0_dp, th - h], b, sgm, bder, hcov, hctr, hcurl)
        call pert_eval(s0, th - h, pm)
        dth_term = (sgp*pp%eta - sgm*pm%eta)/(2*h)
        call do_magfie([s0, 0.0_dp, th], b, sg0, bder, hcov, hctr, hcurl)
        call pert_eval(s0, th, p0)
        div = ds_term + dth_term
        div = div/sg0
    end function div_xi_fd

end program test_line_drive
