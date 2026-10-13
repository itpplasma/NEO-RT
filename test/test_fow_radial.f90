program test_fow_radial
    ! First-order radial orbit width (fow_radial).
    ! 1. Size: the displacement at theta = 0 matches the large-aspect-ratio
    !    poloidal-gyroradius estimate Delta s = 2 q R0 v_par / (a**2 omega_c)
    !    (Delta r = q R0 v_par/(r omega_c), s = r**2/a**2) to O(eps).
    ! 2. Labelling: s is the p_phi label. The displacement averages to zero over
    !    trapped orbits; for passing orbits its mean is nonzero and odd in v_par.
    ! 3. Sign: d(Delta s)/dt along an integrated orbit (central differences)
    !    equals the radial magnetic drift v_d . grad s evaluated from the field,
    !    v_d = (m c / (q B^3)) (v_par^2 + v_perp^2/2) B x grad B, i.e.
    !    ds/dt = m c W B_phi dB/dtheta / (q B^3 J) with J = -sqrtg the Jacobian
    !    of x = (s, phi, theta). The orientation J = -sqrtg is fixed independently
    !    of fow_displacement by test_fow_gc_orbit, whose full guiding-centre
    !    precession reproduces Om_tB only for J = -sqrtg.
    ! 4. No side effects: a fow_radial timestep leaves the flux-surface module
    !    quantities of s (q, dqds, iota, B_theta, B_phi and their s-derivatives)
    !    and the orbit equations unchanged.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie_init, do_magfie, s, a, R0, q, dqds, iota, &
        Bthcov, Bphcov, dBthcovds, dBphcovds
    use do_magfie_pert_mod, only: do_magfie_pert_init
    use neort, only: init_profiles, init, set_to_trapped_region
    use neort_config, only: read_and_set_config
    use neort_profiles, only: read_and_init_plasma_input, read_and_init_profile_input
    use neort_orbit, only: nvar, bounce_fast, bounce_time
    use neort_freq, only: Om_th
    use neort_transport, only: timestep_transport, fow_displacement, Omth
    use util, only: pi, c, qi
    use driftorbit

    implicit none

    real(dp) :: v, d1, d2, mr_trap(2), mr_co(2), mr_ctr(2)
    integer :: failures

    failures = 0
    call read_and_set_config("driftorbit.in")
    call do_magfie_init("in_file")
    call init_profiles(R0)
    call read_and_init_plasma_input("plasma.in", s)
    call read_and_init_profile_input("profile.in", s, R0, efac, bfac)
    call init
    sign_vpar = 1
    mth = 0
    call set_to_trapped_region(etamin, etamax)
    v = 0.9_dp*vth

    call check_size(0.5_dp*(etamin + etamax))
    call check_mean("trapped", 0.5_dp*(etamin + etamax), mr_trap)
    if (abs(mr_trap(1)) > 1.0e-3_dp*mr_trap(2)) failures = failures + 1
    call check_mean("passing", 0.5_dp*etatp, mr_co)
    sign_vpar = -1
    call check_mean("passing counter", 0.5_dp*etatp, mr_ctr)
    sign_vpar = 1
    ! The passing mean is the label offset, O(rms), and odd in v_par
    if (.not. (abs(mr_co(1)) > 0.5_dp*mr_co(2))) failures = failures + 1
    if (abs(mr_co(1) + mr_ctr(1)) > 1.0e-3_dp*abs(mr_co(1))) failures = failures + 1

    call check_no_side_effects(0.5_dp*etatp)
    call check_no_side_effects(0.5_dp*(etamin + etamax))

    call check_drift("trapped co", 0.5_dp*(etamin + etamax), 1.0_dp)
    call check_drift("trapped counter", 0.5_dp*(etamin + etamax), -1.0_dp)
    call check_drift("passing co", 0.5_dp*etatp, 1.0_dp)
    call check_drift("passing counter", 0.5_dp*etatp, -1.0_dp)

    if (failures > 0) error stop "test_fow_radial failed"
    print *, "test_fow_radial passed"

contains

    subroutine check_size(eta)
        real(dp), intent(in) :: eta
        real(dp) :: x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: v_par, ds, ds_ref, omc

        x = [s, 0.0_dp, 0.0_dp]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        v_par = v*sqrt(1.0_dp - eta*bmod)
        omc = qi*bmod/(mi*c)
        ds = fow_displacement(v_par, hcovar(2))
        ds_ref = 2.0_dp*q*R0*v_par/(a**2*omc)
        print '(a,es12.4,a,es12.4,a,f8.4)', "trapped theta=0: Delta s =", ds, &
            "  large-aspect =", ds_ref, "  |ratio| =", abs(ds/ds_ref)
        if (abs(abs(ds/ds_ref) - 1.0_dp) > 0.25_dp) failures = failures + 1
    end subroutine check_size

    subroutine check_drift(label, eta, dir)
        ! RK4 orbit in (theta, v_par); compare dDelta s/dt by central differences
        ! with the independent drift formula at the same points
        character(*), intent(in) :: label
        real(dp), intent(in) :: eta, dir
        integer, parameter :: nstep = 400
        real(dp) :: y(2), yk(2, 4), dt, ds(0:nstep), dr(0:nstep), x(3)
        real(dp) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3), worst, fd
        integer :: k

        x = [s, 0.0_dp, 0.3_dp]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        y = [0.3_dp, dir*v*sqrt(1.0_dp - eta*bmod)]
        dt = 0.15_dp*q*R0/v/nstep
        do k = 0, nstep
            ds(k) = fow_displacement(y(2), hphi_at(y(1)))
            dr(k) = drift_s(y(1), y(2), eta)
            if (k == nstep) exit
            call rhs(y, eta, yk(:, 1))
            call rhs(y + 0.5_dp*dt*yk(:, 1), eta, yk(:, 2))
            call rhs(y + 0.5_dp*dt*yk(:, 2), eta, yk(:, 3))
            call rhs(y + dt*yk(:, 3), eta, yk(:, 4))
            y = y + dt/6.0_dp*(yk(:, 1) + 2.0_dp*yk(:, 2) + 2.0_dp*yk(:, 3) + yk(:, 4))
        end do
        worst = 0.0_dp
        do k = 1, nstep - 1
            fd = (ds(k + 1) - ds(k - 1))/(2.0_dp*dt)
            worst = max(worst, abs(fd - dr(k))/maxval(abs(dr)))
        end do
        print '(a,a,es12.4,a,es12.4,a,es10.2)', label, ": max|d(Delta s)/dt| =", &
            maxval(abs(dr)), "  max |fd - drift|/max =", worst
        if (.not. (worst < 1.0e-4_dp)) failures = failures + 1
    end subroutine check_drift

    subroutine rhs(y, eta, ydot)
        real(dp), intent(in) :: y(2), eta
        real(dp), intent(out) :: ydot(2)
        real(dp) :: yfull(7), yd(7)

        yfull = 0.0_dp
        yfull(1:2) = y
        call timestep_transport(v, eta, 7, 0.0_dp, yfull, yd)
        ydot = yd(1:2)
    end subroutine rhs

    function drift_s(theta, v_par, eta) result(dsdt)
        ! v_d . grad s for low beta from B, dB/dtheta and the Jacobian only;
        ! (B x grad B)^s = B_phi dB/dtheta / J with J = -sqrtg for (s, phi, theta)
        real(dp), intent(in) :: theta, v_par, eta
        real(dp) :: dsdt
        real(dp) :: x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: w

        x = [s, 0.0_dp, theta]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        w = v_par**2 + 0.5_dp*v**2*eta*bmod  ! v_par^2 + v_perp^2/2
        dsdt = mi*c*w*(bmod*hcovar(2))*(bmod*hder(3))/(qi*bmod**3*(-sqrtg))
    end function drift_s

    subroutine check_mean(label, eta, mean_rms)
        ! Orbit mean and rms of Delta s
        character(*), intent(in) :: label
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: mean_rms(2)
        real(dp) :: taub, ba(nvar), ba_rms(nvar)

        ! Exact orbit period (bounce_time locates co-directed orbits only), so
        ! that the orbit closes and the mean formula is checked without the
        ! frequency-spline error
        call Om_th(v, eta, Omth, d1, d2)
        taub = 2.0_dp*pi/abs(Omth)
        if (sign_vpar > 0) then
            taub = bounce_time(v, eta, taub)
            Omth = sign(2.0_dp*pi/taub, Omth)
        end if
        call bounce_fast(v, eta, taub, ba, ts_mean)
        call bounce_fast(v, eta, taub, ba_rms, ts_rms)
        print '(a,a,es12.4,a,es12.4)', label, ": <Delta s> =", ba(3), &
            "  rms Delta s =", sqrt(ba_rms(3))
        if (.not. (sqrt(ba_rms(3)) > 0.0_dp)) failures = failures + 1
        mean_rms = [ba(3), sqrt(ba_rms(3))]
    end subroutine check_mean

    subroutine check_no_side_effects(eta)
        real(dp), intent(in) :: eta
        real(dp) :: before(7), after(7), y(7), yd_thin(7), yd_fow(7), x(3)
        real(dp) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x = [s, 0.0_dp, 0.3_dp]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        y = 0.0_dp
        y(1:2) = [0.3_dp, v*sqrt(1.0_dp - eta*bmod)]
        call timestep_transport(v, eta, 7, 0.0_dp, y, yd_thin)
        before = [q, dqds, iota, Bthcov, Bphcov, dBthcovds, dBphcovds]
        fow_radial = .true.
        call timestep_transport(v, eta, 7, 0.0_dp, y, yd_fow)
        fow_radial = .false.
        after = [q, dqds, iota, Bthcov, Bphcov, dBthcovds, dBphcovds]
        print '(a,es10.2,a,es10.2)', "fow_radial timestep: max state change =", &
            maxval(abs(after - before)), "  max orbit rhs change =", &
            maxval(abs(yd_fow(1:2) - yd_thin(1:2)))
        if (any(after /= before)) failures = failures + 1
        if (any(yd_fow(1:2) /= yd_thin(1:2))) failures = failures + 1
    end subroutine check_no_side_effects

    subroutine ts_mean(v_, eta_, neq, t, y, ydot)
        real(dp), intent(in) :: v_, eta_
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq)
        real(dp), intent(out) :: ydot(neq)

        call timestep_transport(v_, eta_, neq, t, y, ydot)
        ydot(3) = fow_displacement(y(2), hphi_at(y(1)))
        ydot(4:) = 0.0_dp
    end subroutine ts_mean

    subroutine ts_rms(v_, eta_, neq, t, y, ydot)
        real(dp), intent(in) :: v_, eta_
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq)
        real(dp), intent(out) :: ydot(neq)

        call timestep_transport(v_, eta_, neq, t, y, ydot)
        ydot(3) = fow_displacement(y(2), hphi_at(y(1)))**2
        ydot(4:) = 0.0_dp
    end subroutine ts_rms

    function hphi_at(theta) result(hphi)
        real(dp), intent(in) :: theta
        real(dp) :: hphi
        real(dp) :: x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x = [s, 0.0_dp, theta]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        hphi = hcovar(2)
    end function hphi_at

end program test_fow_radial
