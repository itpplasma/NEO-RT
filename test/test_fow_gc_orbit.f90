program test_fow_gc_orbit
    ! Guiding-centre orbit oracle for the first-order radial orbit width.
    !
    ! The full guiding-centre (Littlejohn) orbit is integrated in the
    ! examples/base field with state (s, theta, phi, v_par), fixed v and
    ! eta = v_perp**2/(v**2 B):
    !   dX/dt      = (v_par B* + (c mu/e) h x grad B)/B*_par
    !   dv_par/dt  = -(mu/m) B* . grad B / B*_par,  B* = B + (m c v_par/e) curl h
    ! in coordinates x = (s, phi, theta) with Jacobian J = eps*sqrtg. Curl and
    ! cross products are built here from h_phi, h_theta and B(s, theta), not from
    ! NEO-RT's hcurl, and p_phi = m v_par h_phi + (e/c) A_phi(s) with
    ! A_phi' = J B^theta = eps*sign_theta*iota*psi_pr is checked to be conserved.
    !
    ! 1. Handedness: the bounce-averaged toroidal precession of a trapped orbit
    !    started at its bounce point on s must equal NEO-RT's Om_tB (shear term
    !    on) for exactly one eps (tolerance 0.1; the wrong eps flips the sign).
    ! 2. Displacement (tolerance 0.05 of max|Delta s|, declared before the run):
    !    s(t) - s_label against fow_displacement along the orbit, where s_label
    !    is the p_phi label (e/c) A_phi(s_label) = p_phi, and s(t) - <s>_t
    !    against fow_displacement minus its orbit mean, for a trapped and a co-
    !    and counter-passing orbit.
    ! 3. Labelling: NEO-RT's passing drift frequency Om_tB(s) contains the shear
    !    term q' <Delta s_raw thetadot>. The full-orbit Om_phi - q(s) Om_theta
    !    must match it within 0.1 for the orbit whose p_phi label is s; the orbit
    !    whose time-mean radius is s is reported for comparison and must deviate
    !    by more (it differs at first order by q' <Delta s> Om_theta).
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie_init, do_magfie, s, R0, q, iota, psi_pr, &
        sign_theta, Bthcov, Bphcov, dBthcovds, dBphcovds
    use neort, only: init_profiles, init
    use neort_config, only: read_and_set_config
    use neort_profiles, only: read_and_init_plasma_input, read_and_init_profile_input
    use neort_orbit, only: noshear
    use neort_freq, only: Om_tB, Om_th
    use neort_transport, only: fow_displacement
    use util, only: pi, c, qi
    use driftorbit

    implicit none

    integer, parameter :: nstep = 20000
    real(dp), parameter :: tol_prec = 0.1_dp, tol_ds = 0.05_dp, tol_freq = 0.1_dp
    real(dp), parameter :: vfrac = 0.25_dp
    real(dp) :: v, eta_t, eta_p, eps_found, ratio(2), eps_try(2)
    real(dp) :: ys(4, 0:nstep), period
    integer :: failures, ie

    failures = 0
    call read_and_set_config("driftorbit.in")
    noshear = .false.
    call do_magfie_init("in_file")
    call init_profiles(R0)
    call read_and_init_plasma_input("plasma.in", s)
    call read_and_init_profile_input("profile.in", s, R0, efac, bfac)
    call init
    v = vfrac*vth
    eta_t = 0.5_dp*(etatp + etadt)
    eta_p = 0.5_dp*etatp

    ! 1. Handedness from the trapped precession
    eps_try = [1.0_dp, -1.0_dp]
    do ie = 1, 2
        ratio(ie) = trapped_precession_ratio(eps_try(ie))
    end do
    if (count(abs(ratio - 1.0_dp) < tol_prec) /= 1) then
        print *, "FAIL: no unique handedness reproduces Om_tB"
        failures = failures + 1
        eps_found = 1.0_dp
    else
        eps_found = merge(eps_try(1), eps_try(2), abs(ratio(1) - 1.0_dp) < tol_prec)
    end if
    print '(a,f5.1,a)', "J(s,phi,theta) = ", eps_found, " * sqrtg"

    ! 2. Displacement against the orbit; 3. labelling of passing orbits
    call check_trapped_displacement(eps_found)
    call check_passing(eps_found, 1.0_dp, "co-passing")
    call check_passing(eps_found, -1.0_dp, "counter-passing")

    if (failures > 0) error stop "test_fow_gc_orbit failed"
    print *, "test_fow_gc_orbit passed"

contains

    function trapped_precession_ratio(eps) result(r)
        real(dp), intent(in) :: eps
        real(dp) :: r, y0(4), OmtB, d1, d2, Omph

        y0 = [s, bounce_theta(eta_t), 0.0_dp, 0.0_dp]
        call run_orbit(eps, eta_t, y0, .true., ys, period)
        Omph = (ys(3, nstep) - ys(3, 0))/period
        call Om_tB(v, eta_t, OmtB, d1, d2)
        r = Omph/OmtB
        print '(a,f5.1,a,es13.5,a,es13.5,a,f9.5)', "trapped eps=", eps, &
            ": full-orbit Om_phi =", Omph, "  NEO-RT Om_tB =", OmtB, "  ratio =", r
    end function trapped_precession_ratio

    subroutine check_trapped_displacement(eps)
        real(dp), intent(in) :: eps
        real(dp) :: y0(4)

        ! Started at the bounce point on s, so the p_phi label is exactly s
        y0 = [s, bounce_theta(eta_t), 0.0_dp, 0.0_dp]
        call run_orbit(eps, eta_t, y0, .true., ys, period)
        call compare_displacement("trapped", eta_t, s)
    end subroutine check_trapped_displacement

    subroutine check_passing(eps, dir, label)
        real(dp), intent(in) :: eps, dir
        character(*), intent(in) :: label
        real(dp) :: y0(4), s0, omb_label, omb_mean, OmtB, d1, d2, smean
        integer :: it

        sign_vpar = nint(dir)
        call Om_tB(v, eta_p, OmtB, d1, d2)

        ! Orbit whose p_phi label is s
        s0 = s
        do it = 1, 4
            y0 = passing_start(s0, dir)
            s0 = s0 + (s - p_phi_label(eps, y0))
        end do
        y0 = passing_start(s0, dir)
        call run_orbit(eps, eta_p, y0, .false., ys, period)
        call compare_displacement(label, eta_p, p_phi_label(eps, y0))
        omb_label = drift_freq()
        smean = sum(ys(1, 0:nstep - 1))/nstep

        ! Orbit whose time-mean radius is s
        do it = 1, 3
            s0 = s0 - (smean - s)
            y0 = passing_start(s0, dir)
            call run_orbit(eps, eta_p, y0, .false., ys, period)
            smean = sum(ys(1, 0:nstep - 1))/nstep
        end do
        omb_mean = drift_freq()

        print '(a,a,es13.5,a,es13.5,a,es13.5)', label, &
            ": Om_phi - q Om_th: label s =", omb_label, "  mean s =", omb_mean, &
            "  NEO-RT Om_tB(s) =", OmtB
        print '(a,a,f9.5,a,f9.5)', label, ": relative deviation label s =", &
            abs(omb_label/OmtB - 1.0_dp), "  mean s =", abs(omb_mean/OmtB - 1.0_dp)
        if (.not. (abs(omb_label/OmtB - 1.0_dp) < tol_freq)) failures = failures + 1
        if (.not. (abs(omb_mean/OmtB - 1.0_dp) > abs(omb_label/OmtB - 1.0_dp))) &
            failures = failures + 1
        sign_vpar = 1
    end subroutine check_passing

    function drift_freq() result(omb)
        ! Om_phi - q(s) Om_theta of the stored orbit
        real(dp) :: omb

        call restore_surface
        omb = ((ys(3, nstep) - ys(3, 0)) - q*(ys(2, nstep) - ys(2, 0)))/period
    end function drift_freq

    subroutine compare_displacement(label, eta, slabel)
        character(*), intent(in) :: label
        real(dp), intent(in) :: eta, slabel
        real(dp) :: ds_orb(0:nstep), ds_fow(0:nstep), err_label, err_mean
        real(dp) :: mo, mf
        integer :: k

        call restore_surface
        do k = 0, nstep
            ds_orb(k) = ys(1, k) - slabel
            ds_fow(k) = fow_displacement(ys(4, k), hphi_at(s, ys(2, k)))
        end do
        mo = sum(ds_orb(0:nstep - 1))/nstep
        mf = sum(ds_fow(0:nstep - 1))/nstep
        err_label = maxval(abs(ds_orb - ds_fow))/maxval(abs(ds_orb))
        err_mean = maxval(abs((ds_orb - mo) - (ds_fow - mf)))/maxval(abs(ds_orb - mo))
        print '(a,a,es12.4,a,es12.4,a,es12.4)', label, ": max|s - s_label| =", &
            maxval(abs(ds_orb)), "  <s> - s_label =", mo, "  <fow> =", mf
        print '(a,a,f9.5,a,f9.5)', label, ": error vs label =", err_label, &
            "  error vs mean =", err_mean
        if (.not. (err_label < tol_ds)) failures = failures + 1
        if (.not. (err_mean < tol_ds)) failures = failures + 1
    end subroutine compare_displacement

    subroutine run_orbit(eps, eta, y0, trapped, traj, tper)
        ! Find the orbit period with a first pass, then integrate exactly one
        ! period with nstep RK4 steps and store the trajectory
        real(dp), intent(in) :: eps, eta, y0(4)
        logical, intent(in) :: trapped
        real(dp), intent(out) :: traj(4, 0:nstep), tper
        real(dp) :: y(4), yold(4), dt, Omth_, d1, d2, frac, p0, dp_max
        integer :: k, nchange, sref

        call Om_th(v, eta, Omth_, d1, d2)
        dt = 2.0_dp*pi/abs(Omth_)/nstep
        y = y0
        nchange = 0
        sref = 0
        tper = -1.0_dp
        do k = 0, 3*nstep
            yold = y
            call rk4(eps, eta, y, dt)
            if (trapped) then
                if (sref == 0) then
                    sref = nint(sign(1.0_dp, y(4)))
                else if (nint(sign(1.0_dp, y(4))) /= sref) then
                    sref = -sref
                    nchange = nchange + 1
                    if (nchange == 2) then
                        frac = yold(4)/(yold(4) - y(4))
                        tper = (k + frac)*dt
                        exit
                    end if
                end if
            else if (abs(y(2) - y0(2)) >= 2.0_dp*pi) then
                frac = (2.0_dp*pi - abs(yold(2) - y0(2)))/abs(y(2) - yold(2))
                tper = (k + frac)*dt
                exit
            end if
        end do
        if (tper < 0.0_dp) error stop "run_orbit: orbit did not close"

        dt = tper/nstep
        traj(:, 0) = y0
        p0 = p_phi_rel(eps, y0, y0)
        dp_max = 0.0_dp
        do k = 1, nstep
            traj(:, k) = traj(:, k - 1)
            call rk4(eps, eta, traj(:, k), dt)
            dp_max = max(dp_max, abs(p_phi_rel(eps, traj(:, k), y0) - p0))
        end do
        ! p_phi conservation, relative to m v |h_phi|
        if (.not. (dp_max < 1.0e-6_dp)) then
            print '(a,es10.2)', "FAIL: p_phi not conserved, relative drift", dp_max
            failures = failures + 1
        end if
        call restore_surface
    end subroutine run_orbit

    subroutine rk4(eps, eta, y, dt)
        real(dp), intent(in) :: eps, eta, dt
        real(dp), intent(inout) :: y(4)
        real(dp) :: k1(4), k2(4), k3(4), k4(4)

        call gc_rhs(eps, eta, y, k1)
        call gc_rhs(eps, eta, y + 0.5_dp*dt*k1, k2)
        call gc_rhs(eps, eta, y + 0.5_dp*dt*k2, k3)
        call gc_rhs(eps, eta, y + dt*k3, k4)
        y = y + dt/6.0_dp*(k1 + 2.0_dp*k2 + 2.0_dp*k3 + k4)
    end subroutine rk4

    subroutine gc_rhs(eps, eta, y, ydot)
        ! y = (s, theta, phi, v_par); coordinates x = (s, phi, theta), J = eps*sqrtg
        real(dp), intent(in) :: eps, eta, y(4)
        real(dp), intent(out) :: ydot(4)
        real(dp) :: x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: dBds, dBdth, hph, hth, dhph_ds, dhth_ds, dhph_dth, jac
        real(dp) :: curl(3), hxg(3), bstar(3), bpar, k, gm

        x = [y(1), 0.0_dp, y(2)]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        dBds = bmod*hder(1)
        dBdth = bmod*hder(3)
        hph = Bphcov/bmod
        hth = Bthcov/bmod
        dhph_ds = dBphcovds/bmod - Bphcov*dBds/bmod**2
        dhth_ds = dBthcovds/bmod - Bthcov*dBds/bmod**2
        dhph_dth = -Bphcov*dBdth/bmod**2
        jac = eps*sqrtg
        ! contravariant (s, phi, theta) components, h_s = 0, d/dphi = 0
        curl = [-dhph_dth, -dhth_ds, dhph_ds]/jac
        hxg = [hph*dBdth, hth*dBds, -hph*dBds]/jac
        k = mi*c*y(4)/qi
        bstar = bmod*[0.0_dp, hctrvr(2), hctrvr(3)] + k*curl
        bpar = bmod + k*(hph*curl(2) + hth*curl(3))
        gm = c*mi*v**2*eta/(2.0_dp*qi)  ! c mu/e
        ydot(1) = (y(4)*bstar(1) + gm*hxg(1))/bpar
        ydot(3) = (y(4)*bstar(2) + gm*hxg(2))/bpar
        ydot(2) = (y(4)*bstar(3) + gm*hxg(3))/bpar
        ydot(4) = -0.5_dp*v**2*eta*(bstar(1)*dBds + bstar(3)*dBdth)/bpar
    end subroutine gc_rhs

    function p_phi_rel(eps, y, yref) result(p)
        ! (p_phi(y) - (e/c) A_phi(yref(1)))/(m v |h_phi(yref)|) with
        ! A_phi(s) - A_phi(sref) = eps sign_theta psi_pr int iota ds
        real(dp), intent(in) :: eps, y(4), yref(4)
        real(dp) :: p

        p = (mi*y(4)*hphi_at(y(1), y(2)) &
             + qi/c*eps*sign_theta*psi_pr*int_iota(yref(1), y(1))) &
            /(mi*v*abs(hphi_at(yref(1), yref(2))))
    end function p_phi_rel

    function p_phi_label(eps, y) result(sl)
        ! s_label with (e/c) A_phi(s_label) = p_phi(y)
        real(dp), intent(in) :: eps, y(4)
        real(dp) :: sl, rhs
        integer :: it

        rhs = mi*y(4)*hphi_at(y(1), y(2))*c/(qi*eps*sign_theta*psi_pr)
        sl = y(1) + rhs/iota_at(y(1))
        do it = 1, 8
            if (sl == y(1)) exit
            sl = y(1) + rhs*(sl - y(1))/int_iota(y(1), sl)
        end do
        call restore_surface
    end function p_phi_label

    function int_iota(a1, b1) result(r)
        ! Simpson rule, ample for |b1 - a1| ~ orbit width
        real(dp), intent(in) :: a1, b1
        real(dp) :: r

        r = (b1 - a1)/6.0_dp*(iota_at(a1) + 4.0_dp*iota_at(0.5_dp*(a1 + b1)) &
                              + iota_at(b1))
    end function int_iota

    function iota_at(sv) result(r)
        real(dp), intent(in) :: sv
        real(dp) :: r, x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x = [sv, 0.0_dp, 0.0_dp]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        r = iota
    end function iota_at

    function hphi_at(sv, theta) result(hphi)
        real(dp), intent(in) :: sv, theta
        real(dp) :: hphi, x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x = [sv, 0.0_dp, theta]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        hphi = hcovar(2)
    end function hphi_at

    function bmod_at(sv, theta) result(bmod)
        real(dp), intent(in) :: sv, theta
        real(dp) :: bmod, x(3), sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x = [sv, 0.0_dp, theta]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
    end function bmod_at

    function bounce_theta(eta) result(th)
        ! Bounce point on s: eta B(s, th) = 1 with th in (0, pi)
        real(dp), intent(in) :: eta
        real(dp) :: th, lo, hi
        integer :: it

        lo = 0.0_dp
        hi = pi
        do it = 1, 200
            th = 0.5_dp*(lo + hi)
            if (eta*bmod_at(s, th) < 1.0_dp) then
                lo = th
            else
                hi = th
            end if
        end do
        th = 0.5_dp*(lo + hi)
        call restore_surface
    end function bounce_theta

    function passing_start(s0, dir) result(y0)
        real(dp), intent(in) :: s0, dir
        real(dp) :: y0(4)

        y0 = [s0, 0.0_dp, 0.0_dp, dir*v*sqrt(1.0_dp - eta_p*bmod_at(s0, 0.0_dp))]
        call restore_surface
    end function passing_start

    subroutine restore_surface
        ! Reset the module flux-surface quantities (q, iota, ...) to s
        real(dp) :: x(3), bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x = [s, 0.0_dp, 0.0_dp]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
    end subroutine restore_surface

end program test_fow_gc_orbit
