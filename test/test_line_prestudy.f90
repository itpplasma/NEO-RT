program test_line_prestudy
    !! Thin-orbit reproduction of the rmp-proposal pre-study (neort-realspace,
    !! results.md section 4): circular model field, trapped orbit r0 = 0.2 R0
    !! with vpar/v = 0.4 at the outboard midplane, n = 1, bounce harmonic 0 on
    !! the resonance n Om_phi = 0, generic smooth displacement xi.  Pre-study
    !! (finite orbit width, rho* -> 0): line vs flux-surface drive agree to
    !! O(rho*^2); Boozer weight with this xi is off by 0.17-0.19; with
    !! dPhi_E = 0 the historical finite-width misalignment potential carries
    !! 9-11 % of |H_0|. This thin-path test certifies the complex potential
    !! channel independently; the historical finite-width range is descriptive.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, bfac, read_boozer_file, set_s, &
        init_magfie_at_s, q, psi_pr, booz_to_cyl, do_magfie, &
        Bphcov, iota, sign_theta
    use driftorbit, only: mth, epsmn, sign_vpar, etatp
    use do_magfie_pert_mod, only: set_mph, mph
    use neort_magfie, only: init_flux_surface_average
    use neort_orbit, only: noshear, bounce_time, th0
    use line_period_oracle, only: thin_period_oracle
    use prestudy_potential_oracle, only: prestudy_potential_harmonic
    use neort_profiles, only: vth
    use neort_line_drive, only: set_line_drive_options, line_bounce
    use neort_line_drive, only: line_harmonics_t, exact_orbit_period
    use neort_line_drive, only: collapse_eta, ideal_potential, circ_a
    use neort_line_drive, only: sigma_theta, sigma_phi, gauge_amp
    use neort_line_drive, only: local_field, local_field_t, drift_velocity
    use util, only: pi, qi, mi, qe, mu, c

    implicit none

    real(dp), parameter :: r0 = 0.2_dp, vpar_over_v = 0.4_dp
    character(len=1024) :: boozer
    real(dp) :: dpsi_e, s0, xcyl(3), bmod, sqrtg, bder(3), hcov(3), hctr(3), hcurl(3)
    real(dp) :: eta, taub, omth, d_fs, d_boozer, misalign, d_gauge, q_exact, rr
    real(dp) :: om_te, independent_period, quadrature_error, prefactor
    complex(dp) :: potential_expected
    type(line_harmonics_t) :: res, resphi, resg
    integer :: nfail = 0, istat, quadrature_order

    call get_environment_variable("BOOZER_PRESTUDY_FILE", boozer)
    if (len_trim(boozer) == 0) error stop "BOOZER_PRESTUDY_FILE must be set"
    qi = qe
    mi = 2.0_dp*mu
    vth = 1.0e7_dp
    bfac = 1.0_dp
    noshear = .false.
    inp_swi = 9
    call read_boozer_file(trim(boozer))
    ! toroidal flux per radian F (R0 - sqrt(R0^2 - r^2)), F = 1 T m = 1e6 G cm
    dpsi_e = abs(psi_pr)/1.0e8_dp
    circ_a = sqrt(1.0_dp - (1.0_dp - dpsi_e)**2)
    s0 = (1.0_dp - sqrt(1.0_dp - r0**2))/dpsi_e
    call set_s(s0)
    call init_magfie_at_s()
    call init_flux_surface_average(s0)
    rr = r0
    q_exact = (1.5_dp + 8.0_dp*rr**2)/sqrt(1.0_dp - rr**2)
    call check('q(s0) vs closed form', abs(abs(q)/q_exact - 1), 1.0e-3_dp)

    ! orientation of the Boozer angles: Z > 0 at theta_B = pi/2 means theta_B
    ! runs counterclockwise like the pre-study angle; phi is then oriented so
    ! that the field line has the pre-study's positive pitch.
    call booz_to_cyl([s0, 0.0_dp, 0.5_dp*pi], xcyl)
    sigma_theta = sign(1.0_dp, xcyl(3))
    sigma_phi = sigma_theta*sign(1.0_dp, q)
    call set_mph(nint(sigma_phi))
    epsmn = 1.0e-3_dp
    collapse_eta = .false.
    call set_line_drive_options('line', 'circ_prestudy', .false., 9, noshear)
    print '(a,2f5.1,a,f8.5,a,f8.4)', ' sigma_theta, sigma_phi', sigma_theta, &
        sigma_phi, '  s0 =', s0, '  q =', q

    call check_cartesian_field()

    call do_magfie([s0, 0.0_dp, 0.0_dp], bmod, sqrtg, bder, hcov, hctr, hcurl)
    eta = (1.0_dp - vpar_over_v**2)/bmod
    sign_vpar = 1.0_dp
    mth = 0
    call exact_orbit_period(vth, eta, bounce_time(vth, eta), taub, omth)
    call thin_period_oracle(s0, th0, vth, eta, eta > etatp, independent_period, &
        quadrature_error, quadrature_order)
    call check('period vs independent energy quadrature', &
        abs(taub/independent_period - 1.0_dp), 5.0e-10_dp)
    call line_bounce(vth, eta, taub, omth, 0.0_dp, 0.0_dp, res, istat)
    if (istat /= 2) error stop "line-drive oracle: integration failed"
    ! E x B rotation cancels the precession: n Om_phi = 0
    om_te = -res%om_drift
    ideal_potential = .true.
    call line_bounce(vth, eta, taub, omth, 0.0_dp, om_te, res, istat)
    if (istat /= 2) error stop "line-drive oracle: integration failed"
    ideal_potential = .false.
    call line_bounce(vth, eta, taub, omth, 0.0_dp, om_te, resphi, istat)
    if (istat /= 2) error stop "line-drive oracle: integration failed"
    ideal_potential = .true.
    gauge_amp = (5.0_dp, -2.0_dp)*epsmn
    call line_bounce(vth, eta, taub, omth, 0.0_dp, om_te, resg, istat)
    if (istat /= 2) error stop "line-drive oracle: integration failed"
    gauge_amp = 0

    d_fs = abs(res%h_line - res%h_fs)/abs(res%h_line)
    d_boozer = abs(res%h_boozer - res%h_line)/abs(res%h_line)
    misalign = abs(resphi%h_line - res%h_line)/abs(res%h_line)
    d_gauge = abs(resg%h_line - res%h_line)/abs(res%h_line)
    print '(a,es12.4)', ' |H_0|/(m v^2/2)                ', abs(res%h_line)
    call check('line vs flux surfaces (pre-study 1.5e-7..1.5e-5, FOW)', d_fs, 1.0e-8_dp)
    call check('other gauge (pre-study <= 4e-10)', d_gauge, 1.0e-6_dp)
    call check('Boozer weight, generic xi (pre-study 0.17-0.19), > 0.15', d_boozer, &
        0.15_dp, .true.)
    call check('Boozer weight, generic xi, < 0.21', d_boozer, 0.21_dp)
    prefactor = qi*(-iota*sign_theta*psi_pr*om_te/c)/(0.5_dp*mi*vth**2)
    call prestudy_potential_harmonic(r0, dpsi_e, epsmn, eta, vth, sigma_theta, &
        real(mph, dp)*q, prefactor, potential_expected, &
        quadrature_error)
    call check('Cartesian electrostatic quadrature convergence', quadrature_error, &
        1.0e-10_dp)
    call check('misalignment complex potential vs Cartesian time quadrature', &
        abs(resphi%h_line - res%h_line - potential_expected) &
        /abs(potential_expected), 5.0e-5_dp)
    print '(a,es12.4)', ' thin-path misalignment ratio (finite-width range 0.09-0.11): ', &
        misalign

    if (nfail > 0) then
        print *, 'test_line_prestudy: FAILED checks:', nfail
        error stop 1
    end if
    print *, 'test_line_prestudy: all checks passed'

contains

    subroutine check_cartesian_field()
        real(dp), parameter :: angles(4) = [0.3_dp, 0.8_dp, 1.4_dp, 2.0_dp]
        type(local_field_t) :: lf
        complex(dp) :: dbe_reference
        real(dp) :: vd_reference(3), vd_model(3), vdscale(3)
        real(dp) :: vp, mu_gc, err_b, err_v, bsign
        integer :: j

        bsign = sigma_phi*sign(1.0_dp, Bphcov)
        err_b = 0.0_dp
        err_v = 0.0_dp
        do j = 1, size(angles)
            call local_field(s0, angles(j), lf)
            vp = 0.4_dp*vth
            mu_gc = 0.5_dp*mi*vth**2*(1.0_dp - 0.4_dp**2)/lf%bmod
            call prestudy_cartesian_oracle(s0, angles(j), circ_a, &
                sigma_theta, sigma_phi, epsmn, bsign, &
                vp, mu_gc, dbe_reference, vd_reference)
            vd_model = drift_velocity(vp, mu_gc, lf)
            err_b = max(err_b, abs(lf%dBE - dbe_reference)/abs(dbe_reference))
            vdscale = max(abs(vd_reference), 1.0e-8_dp*maxval(abs(vd_reference)))
            err_v = max(err_v, maxval(abs(vd_model - vd_reference)/vdscale))
        end do
        call check('Cartesian curl dBE', err_b, 5.0e-5_dp)
        call check('Cartesian magnetic drift', err_v, 5.0e-5_dp)
    end subroutine check_cartesian_field

    ! Independent Cartesian oracle helpers for test_line_prestudy CONTAINS.
    ! Host imports dp, mi, qi, c. Geometry R0=1m, B0=1T; arguments use cgs.
    ! bsign = sigma_phi*sign(1.0_dp, Bphcov) aligns physical B with the fixture.
    ! Does not call local_field, pert_eval, drift_velocity, or do_magfie.

    subroutine prestudy_cartesian_oracle(s0, thb, edge_a, st, sp, amp, bsign, &
            vp, mu_gc, dbe, vd)
        real(dp), intent(in) :: s0, thb, edge_a, st, sp, amp, bsign, vp, mu_gc
        complex(dp), intent(out) :: dbe
        real(dp), intent(out) :: vd(3)
        real(dp), parameter :: h = 1.0e-4_dp
        real(dp) :: flux_edge, w, r, k, tg, rmaj, x(3), xp(3), xm(3)
        real(dp) :: b(3), bp(3), bm(3), bhat(3), dbhat(3, 3), gradb(3)
        real(dp) :: bmag, curlb(3), vc(3), er(3), eth(3), vr, vt, vf
        real(dp) :: dsdr, dtbdr, dtbdtg
        complex(dp) :: a(3), ap(3), am(3), da(3, 3), curl_a(3)
        integer :: j

        flux_edge = 1.0_dp - sqrt(1.0_dp - edge_a**2)
        w = 1.0_dp - s0*flux_edge
        r = sqrt(1.0_dp - w*w)
        k = sqrt((1.0_dp - r)/(1.0_dp + r))
        tg = 2.0_dp*atan2(sin(0.5_dp*st*thb), k*cos(0.5_dp*st*thb))
        rmaj = 1.0_dp + r*cos(tg)
        x = 100.0_dp*[rmaj, 0.0_dp, r*sin(tg)]
        call prestudy_cartesian_bundle(x, amp, bsign, b, a)
        bmag = sqrt(sum(b*b))
        bhat = b/bmag
        do j = 1, 3
            xp = x
            xm = x
            xp(j) = xp(j) + h
            xm(j) = xm(j) - h
            call prestudy_cartesian_bundle(xp, amp, bsign, bp, ap)
            call prestudy_cartesian_bundle(xm, amp, bsign, bm, am)
            da(:, j) = (ap - am)/(2.0_dp*h)
            dbhat(:, j) = (bp/sqrt(sum(bp*bp)) - bm/sqrt(sum(bm*bm)))/(2.0_dp*h)
            gradb(j) = (sqrt(sum(bp*bp)) - sqrt(sum(bm*bm)))/(2.0_dp*h)
        end do
        curl_a = [da(3, 2) - da(2, 3), da(1, 3) - da(3, 1), &
            da(2, 1) - da(1, 2)]
        curlb = [dbhat(3, 2) - dbhat(2, 3), dbhat(1, 3) - dbhat(3, 1), &
            dbhat(2, 1) - dbhat(1, 2)]
        dbe = sum(bhat*curl_a)
        vc = (mi*c*vp**2/qi)*curlb/bmag &
            + (c*mu_gc/qi)*prestudy_cross_rr(bhat, gradb)/bmag
        er = [cos(tg), 0.0_dp, sin(tg)]
        eth = [-sin(tg), 0.0_dp, cos(tg)]
        vr = sum(vc*er)/100.0_dp
        vt = sum(vc*eth)/100.0_dp
        vf = vc(2)/100.0_dp
        dsdr = r/(w*flux_edge)
        dtbdr = -sin(tg)/(w*rmaj)
        dtbdtg = w/rmaj
        vd = [dsdr*vr, st*(dtbdr*vr + dtbdtg*vt/r), sp*vf/rmaj]
    end subroutine prestudy_cartesian_oracle

    subroutine prestudy_cartesian_bundle(x, amp, bsign, b, a)
        real(dp), intent(in) :: x(3), amp, bsign
        real(dp), intent(out) :: b(3)
        complex(dp), intent(out) :: a(3)
        real(dp) :: xx(3), rmaj, r, tg, ph, psi_r
        real(dp) :: er(3), eth(3), eph(3)
        complex(dp) :: p1, p2, xr, xt, xf, xi(3)

        xx = x/100.0_dp
        rmaj = sqrt(xx(1)**2 + xx(2)**2)
        r = sqrt((rmaj - 1.0_dp)**2 + xx(3)**2)
        tg = atan2(xx(3), rmaj - 1.0_dp)
        ph = atan2(xx(2), xx(1))
        er = [cos(tg)*cos(ph), cos(tg)*sin(ph), sin(tg)]
        eth = [-sin(tg)*cos(ph), -sin(tg)*sin(ph), cos(tg)]
        eph = [-sin(ph), cos(ph), 0.0_dp]
        psi_r = r/(1.5_dp + 8.0_dp*r*r)
        b = 1.0e4_dp*bsign*(psi_r*eth + eph)/rmaj
        p1 = exp(-((r - 0.2_dp)/0.15_dp)**2) &
            *exp(cmplx(0.0_dp, -2.0_dp*tg + ph, dp))
        p2 = 0.5_dp*exp(-((r - 0.25_dp)/0.2_dp)**2) &
            *exp(cmplx(0.0_dp, -3.0_dp*tg + ph, dp))
        xr = p1 + p2
        xt = cmplx(0.0_dp, 0.4_dp, dp)*p1 - cmplx(0.0_dp, 0.2_dp, dp)*p2
        xf = 0.3_dp*p1
        xi = 100.0_dp*amp*(xr*er + xt*eth + xf*eph)
        a = [xi(2)*b(3) - xi(3)*b(2), xi(3)*b(1) - xi(1)*b(3), &
            xi(1)*b(2) - xi(2)*b(1)]
    end subroutine prestudy_cartesian_bundle

    pure function prestudy_cross_rr(a, b) result(cross_ab)
        real(dp), intent(in) :: a(3), b(3)
        real(dp) :: cross_ab(3)

        cross_ab = [a(2)*b(3) - a(3)*b(2), a(3)*b(1) - a(1)*b(3), &
            a(1)*b(2) - a(2)*b(1)]
    end function prestudy_cross_rr

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

end program test_line_prestudy
