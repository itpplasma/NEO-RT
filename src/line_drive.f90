module neort_line_drive
    !! Coordinate-free line-integral drive, drive_form = 'line'.
    !!
    !! The bounce harmonic of the perturbed guiding-centre one-form (monograph
    !! eq:coorddrive, sign of NEO-RT's Hamiltonian perturbation),
    !!   H_m = (1/T) oint [ e dPhi + mu d|B|_E - (e/c) dA . Xdot ] e^{-i m.theta} dt,
    !! along NEO-RT's thin orbit in Boozer coordinates (s, theta, phi), with Xdot
    !! the guiding-centre velocity vpar h + v_d + v_E at the orbit point.  One
    !! gauge only: the Hamiltonian carries mu d|B|_E, the Eulerian d|B|_E =
    !! h . curl dA, and no (2 - eta B) weight.  H_m/(m v^2/2) goes into the
    !! bounceavg(3:4) slots, so Hmn2 in neort_transport is |H_m|^2 unchanged.
    !!
    !! Thin-orbit gauge: the thin orbit stays on s while Xdot has a radial drift,
    !! so grad(chi) . Xdot is not a total time derivative along it.  The defect
    !! is (e/c) vpar (dX . grad) dA_par with dX the finite-orbit excursion, of
    !! the order of the drive whenever dA has a parallel part (h_naive shows it).
    !! dA is therefore first gauged to dA_par = 0 on the surface (magnetic
    !! differential equation, Fourier modes on s and s +- GAUGE_DS); then the
    !! thin evaluation is the leading order of the exact orbit integral and H_m
    !! is gauge invariant.  Resonant modes iota m + n = 0 cannot be gauged away
    !! and are left in dA_par (thin islands, not treated here).
    !!
    !! Perturbation sources (pert_model):
    !!   'scalar'        NEO-RT's d|B| (epsmn or pertfile); 'line' is rejected.
    !!   'ideal_helical' dA = xi x B0 with xi^s = epsmn 4 s (1-s) e^{i(m0 theta+n phi)}
    !!                   and the poloidal displacement eta = xi^theta - iota xi^phi
    !!                   fixed in closed form by the Boozer collapse condition
    !!                   div xi = -2 d|B|_L/B (monograph eq:coord-hFS), so that the
    !!                   Boozer scalar d|B|_L = d|B|_E + xi.grad|B| applies exactly;
    !!                   'boozer' then uses epsn = d|B|_L/B of the same xi.
    !!   'circ_prestudy' the generic displacement of the rmp-proposal pre-study on
    !!                   its circular model field (test/fixtures/prestudy).
    !! A libneo perturbation_field_t source (dA_R, dA_phi, dA_Z on an R-Z grid)
    !! would enter through local_field: covariant dA_i = dA . dX/dx^i from the
    !! Boozer-to-cylinder map, sqrt(g) dB^i from its curl.  Not wired yet.
    use iso_fortran_env, only: dp => real64
    use util, only: imun, pi, c, qi, mi
    use do_magfie_mod, only: do_magfie, s, iota, q, dqds, psi_pr, sign_theta, &
        Bthcov, Bphcov, dBthcovds, dBphcovds

    implicit none
    private

    integer, parameter, public :: PERT_SCALAR = 0, PERT_IDEAL_HELICAL = 1, &
        PERT_CIRC_PRESTUDY = 2
    integer, parameter, public :: NCOMP = 10

    character(len=8), public :: drive_form = 'boozer'
    integer, public :: pert_model_id = PERT_SCALAR

    ! Diagnostic switches, shared and read-only during a transport run.
    ! collapse_eta = .false. drops the Boozer-collapse poloidal displacement (a
    ! generic, non-Boozer xi); ideal_potential = .false. sets dPhi_E = 0 instead of
    ! the ideal dPhi_E = -xi.grad(Phi0), leaving the misalignment potential.
    logical, public :: collapse_eta = .true.
    logical, public :: ideal_potential = .true.
    ! Gauge shift dA -> dA + grad(chi),
    ! chi = gauge_amp psi_tor' s^2 e^{i(gauge_m theta + n phi)}.
    complex(dp), public :: gauge_amp = (0.0_dp, 0.0_dp)
    integer, public :: gauge_m = 1

    ! circ_prestudy geometry: boundary minor radius over R0 and the orientation
    ! of NEO-RT's Boozer angles relative to the geometric (theta, phi).
    real(dp), public :: circ_a = 0.4_dp
    real(dp), public :: sigma_theta = 1.0_dp, sigma_phi = 1.0_dp
    integer, parameter :: circ_m = -2, circ_n = 1

    type, public :: pert_point_t
        !! Ideal displacement at (s, theta), without the factor e^{i n phi}.
        complex(dp) :: xs = 0, dxs_ds = 0, dxs_dth = 0 ! xi^s and derivatives
        complex(dp) :: eta = 0, deta_dth = 0 ! eta = xi^theta - iota xi^phi
    end type pert_point_t

    type, public :: line_harmonics_t
        !! Bounce harmonics over m v^2/2, all on the same orbit and quadrature.
        complex(dp) :: h_line = 0 ! 'line' result: dA gauged to dA_par = 0 on s
        complex(dp) :: h_boozer = 0 ! (m v^2/2)(2 - eta B) dB_L/B + e dPhi_L
        complex(dp) :: h_fs = 0 ! flux-surface drive h_FS, eq:coord-hFS
        complex(dp) :: p_par = 0 ! m vpar^2 dB_L/B (parallel coupling)
        complex(dp) :: m_xi = 0 ! mu xi.grad|B| (Lagrangian minus Eulerian)
        complex(dp) :: chi = 0 ! (e/c) chi_m of the gauge function
        complex(dp) :: h_boozer_euler = 0 ! Boozer weight with Eulerian dB_E
        complex(dp) :: h_naive = 0 ! same as h_line without the dA_par gauge fix
        complex(dp) :: path = 0 ! (e/c) dA.Xdot_g along the path, raw gauge
        complex(dp) :: mu_dbe = 0 ! mu dB_E (the line form's Hamiltonian part)
        real(dp) :: om_drift = 0 ! orbit average of v_d.grad(alpha), = Om_tB
        real(dp) :: inv_b = 0, b = 0 ! orbit averages of 1/B and B
    end type line_harmonics_t

    integer, parameter :: NFOURIER = 64 ! poloidal modes of the gauge fix
    real(dp), parameter :: GAUGE_DS = 1.0e-4_dp ! s step for d(chi)/ds

    type, public :: gauge_fix_t
        !! chi_par with B.grad(chi_par) = -B dA_par on s - GAUGE_DS, s, s + GAUGE_DS,
        !! as poloidal Fourier coefficients (index -NFOURIER..NFOURIER-1).
        complex(dp) :: coef(-NFOURIER:NFOURIER - 1, 3) = 0
    end type gauge_fix_t

    type, public :: local_field_t
        !! Field and perturbation at one point of the field-line path.
        real(dp) :: bmod = 0, sqrtg = 0, hth = 0, dB_ds = 0, dB_dth = 0
        real(dp) :: beta = 0, dbeta_dth = 0 ! radial covariant B for circular oracle
        complex(dp) :: A(3) = 0 ! covariant dA_(s,theta,phi)
        complex(dp) :: sgB(3) = 0 ! sqrt(g) dB^(s,theta,phi)
        complex(dp) :: dBE = 0, xigradB = 0, divxi = 0, xs = 0, chi = 0
    end type local_field_t

    public :: set_line_drive_options, line_drive_enabled, pert_model_active, &
        pert_eval, line_scalar_eps, line_bounce, line_bounce_transport, &
        line_point, drift_frequency_part, drift_velocity, local_field, &
        exact_orbit_period

contains

    subroutine set_line_drive_options(form, model, pertfile, inp_swi, noshear)
        character(len=*), intent(in) :: form, model
        logical, intent(in) :: pertfile, noshear
        integer, intent(in) :: inp_swi

        select case (trim(model))
        case ('scalar')
            pert_model_id = PERT_SCALAR
        case ('ideal_helical')
            pert_model_id = PERT_IDEAL_HELICAL
        case ('circ_prestudy')
            pert_model_id = PERT_CIRC_PRESTUDY
        case default
            error stop "pert_model must be 'scalar', 'ideal_helical' or 'circ_prestudy'"
        end select
        if (pert_model_id /= PERT_SCALAR .and. pertfile) then
            error stop "pert_model other than 'scalar' requires pertfile = .false."
        end if
        if (pert_model_id /= PERT_SCALAR .and. inp_swi /= 8 .and. inp_swi /= 9) then
            error stop "pert_model other than 'scalar' needs a Boozer .bc (inp_swi 8/9)"
        end if

        select case (trim(form))
        case ('boozer')
        case ('line')
            if (pert_model_id == PERT_SCALAR) then
                error stop "drive_form = 'line' needs a vector-potential pert_model"
            end if
            if (noshear) then
                error stop "drive_form = 'line' requires noshear = .false."
            end if
        case default
            error stop "drive_form must be 'boozer' or 'line'"
        end select
        drive_form = form
    end subroutine set_line_drive_options

    logical function line_drive_enabled()
        line_drive_enabled = (drive_form == 'line')
    end function line_drive_enabled

    logical function pert_model_active()
        pert_model_active = (pert_model_id /= PERT_SCALAR)
    end function pert_model_active

    subroutine pert_eval(s_, th, pp)
        !! Ideal displacement of the active pert_model at (s_, th).  For
        !! 'ideal_helical' the profile data (B_theta, B_phi, iota and their
        !! s-derivatives) are those of the last do_magfie call, i.e. of s.
        use driftorbit, only: epsmn, m0
        use do_magfie_pert_mod, only: mph

        real(dp), intent(in) :: s_, th
        type(pert_point_t), intent(out) :: pp
        complex(dp) :: e, bcoef
        real(dp) :: a, da, GI, kfac, diota, denom

        select case (pert_model_id)
        case (PERT_IDEAL_HELICAL)
            e = epsmn*exp(imun*m0*th)
            a = 4.0_dp*s_*(1.0_dp - s_)
            da = 4.0_dp*(1.0_dp - 2.0_dp*s_)
            pp%xs = a*e
            pp%dxs_ds = da*e
            pp%dxs_dth = imun*m0*pp%xs
            if (collapse_eta) then
                ! div xi + 2 (dB_E + xi.grad B)/B = 0 with B_s = 0 and flux-function
                ! B_theta = I, B_phi = G reduces to a pointwise relation between
                ! the radial and poloidal amplitudes (psi_tor'' = 0).
                GI = Bphcov + iota*Bthcov
                diota = -dqds*iota**2
                kfac = (dBphcovds + iota*dBthcovds - diota*Bthcov)/GI
                denom = m0*(iota*Bthcov - Bphcov) + 2.0_dp*mph*Bthcov
                if (abs(denom) < 1.0e-12_dp*abs(GI)) then
                    error stop "ideal_helical: collapse displacement singular"
                end if
                bcoef = (da - a*kfac)*GI/(imun*denom)
                pp%eta = bcoef*e
                pp%deta_dth = imun*m0*pp%eta
            else
                pp%eta = 0
                pp%deta_dth = 0
            end if
        case (PERT_CIRC_PRESTUDY)
            call circ_prestudy_eval(s_, th, pp)
        case default
            error stop "pert_eval: no displacement for pert_model 'scalar'"
        end select
    end subroutine pert_eval

    subroutine circ_prestudy_eval(s_, th, pp)
        !! Displacement of the rmp-proposal pre-study (rsdrive_fields.py, _ideal_xi)
        !! on its circular model field (R0 = 1, F = 1, q = (1.5+8r^2)/sqrt(1-r^2)),
        !! mapped to Boozer contravariant components.  On concentric circles with
        !! constant F the Boozer angle is theta_B = 2 atan(k tan(theta/2)),
        !! k = sqrt((1-r)/(1+r)), and phi_B is the cylindrical angle.  Derivatives
        !! are central differences of the closed form.
        real(dp), intent(in) :: s_, th
        type(pert_point_t), intent(out) :: pp
        real(dp), parameter :: hs = 1.0e-5_dp, hth = 1.0e-5_dp
        complex(dp) :: xs0, eta0, xsp, xsm, etap, etam

        call circ_xi_boozer(s_, th, xs0, eta0)
        pp%xs = xs0
        pp%eta = eta0
        call circ_xi_boozer(s_ + hs, th, xsp, etap)
        call circ_xi_boozer(s_ - hs, th, xsm, etam)
        pp%dxs_ds = (xsp - xsm)/(2.0_dp*hs)
        call circ_xi_boozer(s_, th + hth, xsp, etap)
        call circ_xi_boozer(s_, th - hth, xsm, etam)
        pp%dxs_dth = (xsp - xsm)/(2.0_dp*hth)
        pp%deta_dth = (etap - etam)/(2.0_dp*hth)
        if (.not. collapse_eta) return
        error stop "circ_prestudy: collapse_eta must be .false. (generic xi)"
    end subroutine circ_prestudy_eval

    subroutine circ_xi_boozer(s_, thb, xs, eta)
        use driftorbit, only: epsmn
        use do_magfie_pert_mod, only: mph

        real(dp), intent(in) :: s_, thb
        complex(dp), intent(out) :: xs, eta
        real(dp) :: dpsi_e, w, r, rmaj, thg, k, dk, ds_dr, dthb_dr, dthb_dthg, hc, hs2
        complex(dp) :: ph1, ph2, xr, xth, xph, xthb, xphb

        if (mph /= circ_n*nint(sigma_phi)) then
            error stop "circ_prestudy: mph must equal n times the phi orientation"
        end if
        dpsi_e = 1.0_dp - sqrt(1.0_dp - circ_a**2)
        w = 1.0_dp - s_*dpsi_e
        r = sqrt(1.0_dp - w**2)
        k = sqrt((1.0_dp - r)/(1.0_dp + r))
        dk = -1.0_dp/(k*(1.0_dp + r)**2)
        thg = 2.0_dp*atan2(sin(0.5_dp*sigma_theta*thb), k*cos(0.5_dp*sigma_theta*thb))
        rmaj = 1.0_dp + r*cos(thg)
        hc = cos(0.5_dp*thg)
        hs2 = sin(0.5_dp*thg)
        ds_dr = r/(w*dpsi_e)
        dthb_dthg = sqrt(1.0_dp - r**2)/rmaj
        dthb_dr = dk*sin(thg)/(hc**2 + k**2*hs2**2)
        ph1 = exp(imun*circ_m*thg)
        ph2 = exp(imun*(circ_m - 1)*thg)
        ! physical components along e_r, e_theta, e_phi
        xr = exp(-((r - 0.2_dp)/0.15_dp)**2)*ph1 &
            + 0.5_dp*exp(-((r - 0.25_dp)/0.2_dp)**2)*ph2
        xth = exp(-((r - 0.2_dp)/0.15_dp)**2)*ph1*0.4_dp*imun &
            - 0.5_dp*exp(-((r - 0.25_dp)/0.2_dp)**2)*ph2*0.2_dp*imun
        xph = exp(-((r - 0.2_dp)/0.15_dp)**2)*ph1*0.3_dp
        xthb = sigma_theta*(dthb_dr*xr + dthb_dthg*xth/r)
        xphb = sigma_phi*xph/rmaj
        xs = epsmn*ds_dr*xr
        eta = epsmn*(xthb - iota*xphb)
    end subroutine circ_xi_boozer

    subroutine local_field(s_, th, lf)
        !! Unperturbed field and perturbation at (s, th) in the right-handed order
        !! (s, theta, phi): psi_tor' = sign_theta psi_pr, sqrt(g) = psi_tor'(G +
        !! iota I)/B^2, B_theta = I, B_phi = G, B_s = 0.  dA = xi x B0 + grad chi
        !! has covariant components (psi' eta, -psi' xi^s, iota psi' xi^s).
        use do_magfie_pert_mod, only: mph

        real(dp), intent(in) :: s_, th
        type(local_field_t), intent(out) :: lf
        real(dp) :: x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3), psit, diota, GI
        type(pert_point_t) :: pp
        complex(dp) :: gphase

        x = [s_, 0.0_dp, th]
        call do_magfie(x, lf%bmod, lf%sqrtg, bder, hcovar, hctrvr, hcurl)
        lf%hth = hctrvr(3)
        lf%dB_ds = bder(1)*lf%bmod
        lf%dB_dth = bder(3)*lf%bmod
        call pert_eval(s_, th, pp)
        psit = sign_theta*psi_pr
        diota = -dqds*iota**2
        GI = Bphcov + iota*Bthcov
        lf%xs = pp%xs
        lf%A = [psit*pp%eta, -psit*pp%xs, iota*psit*pp%xs]
        gphase = gauge_amp*psit*exp(imun*gauge_m*th)
        lf%chi = gphase*s_**2
        lf%A = lf%A + [gphase*2.0_dp*s_, imun*gauge_m*lf%chi, imun*mph*lf%chi]
        lf%sgB(1) = psit*(iota*pp%dxs_dth + imun*mph*pp%xs)
        lf%sgB(2) = imun*mph*psit*pp%eta - psit*(diota*pp%xs + iota*pp%dxs_ds)
        lf%sgB(3) = -psit*(pp%dxs_ds + pp%deta_dth)
        if (pert_model_id == PERT_CIRC_PRESTUDY) then
            call circ_covariant_radial(s_, th, lf%beta, lf%dbeta_dth)
        end if
        lf%dBE = (lf%beta*lf%sgB(1) + Bthcov*lf%sgB(2) + Bphcov*lf%sgB(3)) &
            /(lf%sqrtg*lf%bmod)
        lf%xigradB = pp%xs*lf%dB_ds + pp%eta*lf%dB_dth
        lf%divxi = pp%dxs_ds + pp%deta_dth - 2.0_dp*pp%eta*bder(3) &
            + pp%xs*((dBphcovds + diota*Bthcov + iota*dBthcovds)/GI &
            - 2.0_dp*bder(1))
    end subroutine local_field

    subroutine circ_covariant_radial(s_, thb, beta, dbeta_dth)
        !! The circular map has B_s = B.dot(dX/ds), even though B^s = 0.
        !! Retain it for the Cartesian pre-study oracle. The historical Boozer
        !! background approximates this component as zero for other sources.
        real(dp), intent(in) :: s_, thb
        real(dp), intent(out) :: beta, dbeta_dth
        real(dp) :: dpsi, w, r, k, thg, rmaj

        dpsi = 1.0_dp - sqrt(1.0_dp - circ_a**2)
        w = 1.0_dp - s_*dpsi
        r = sqrt(1.0_dp - w**2)
        k = sqrt((1.0_dp - r)/(1.0_dp + r))
        thg = 2.0_dp*atan2(sin(0.5_dp*sigma_theta*thb), &
            k*cos(0.5_dp*sigma_theta*thb))
        rmaj = 1.0_dp + r*cos(thg)
        beta = Bthcov*sigma_theta*dpsi*sin(thg)/(r*rmaj)
        dbeta_dth = Bthcov*dpsi*(cos(thg) + r)/(r*w*rmaj)
    end subroutine circ_covariant_radial

    subroutine line_point(v, eta_p, th, vpar, om_path, om_te, gf, lf, loc, drate)
        !! Integrand pieces at one point (without the harmonic phase), in the
        !! order of line_harmonics_t.  drate = v_d.grad(alpha) including the
        !! shear part q' ds thetadot of the p_phi excursion ds; its orbit average
        !! is Om_tB.  om_path is the toroidal angular velocity of the path beyond
        !! streaming (Om_tE + Om_tB), om_te its E x B part.
        real(dp), intent(in) :: v, eta_p, th, vpar, om_path, om_te
        type(gauge_fix_t), intent(in) :: gf
        type(local_field_t), intent(out) :: lf
        complex(dp), intent(out) :: loc(NCOMP)
        real(dp), intent(out) :: drate
        real(dp) :: psip, mu, thdot, phi0p, ds, wboozer, vd(3), xdot(3)
        complex(dp) :: dBL, dphiE, dphiL, hmag

        call local_field(s, th, lf)
        psip = iota*sign_theta*psi_pr
        mu = 0.5_dp*mi*v**2*eta_p
        thdot = vpar*lf%hth
        phi0p = -psip*om_te/c
        dphiE = 0
        if (ideal_potential) dphiE = -lf%xs*phi0p
        dphiL = dphiE + lf%xs*phi0p
        dBL = lf%dBE + lf%xigradB
        ds = mi*c*vpar*Bphcov/(qi*psip*lf%bmod)
        vd = drift_velocity(vpar, mu, lf)
        ! guiding-centre velocity vpar h + v_d + v_E (contravariant s, theta, phi)
        xdot = vd + [0.0_dp, thdot, q*thdot] &
            + c*phi0p*[0.0_dp, Bphcov, -Bthcov]/(lf%sqrtg*lf%bmod**2)
        hmag = mu*lf%dBE + qi*dphiE
        drate = vd(3) - q*vd(2) + dqds*ds*thdot
        wboozer = 0.5_dp*mi*v**2*(2.0_dp - eta_p*lf%bmod)
        loc(1) = hmag - (qi/c)*sum((lf%A + gauge_fix_eval(gf, th))*xdot)
        loc(2) = wboozer*dBL/lf%bmod + qi*dphiL
        loc(3) = (mu*lf%bmod - mi*vpar**2)*dBL/lf%bmod - mi*vpar**2*lf%divxi &
            + qi*dphiL
        loc(4) = mi*vpar**2*dBL/lf%bmod
        loc(5) = mu*lf%xigradB
        loc(6) = (qi/c)*lf%chi
        loc(7) = wboozer*lf%dBE/lf%bmod
        loc(8) = hmag - (qi/c)*sum(lf%A*xdot)
        loc(9) = (qi/c)*(lf%A(2)*thdot + lf%A(3)*(q*thdot + om_path))
        loc(10) = mu*lf%dBE
    end subroutine line_point

    subroutine build_gauge_fix(gf)
        !! Solve (iota d/dtheta + i n) chi = -(iota dA_theta + dA_phi), which is
        !! B.grad(chi) = -B dA_par for B_s = 0, mode by mode on three surfaces.
        !! Resonant modes |iota m + n| < 1e-8 are left out (dA_par stays there).
        use do_magfie_pert_mod, only: mph

        type(gauge_fix_t), intent(out) :: gf
        type(local_field_t) :: lf
        complex(dp) :: src(0:2*NFOURIER - 1)
        real(dp) :: th, s_k, iota_k, kpar
        integer :: k, j, m

        do k = 1, 3
            s_k = s + (k - 2)*GAUGE_DS
            do j = 0, 2*NFOURIER - 1
                th = pi*j/NFOURIER
                call local_field(s_k, th, lf)
                iota_k = iota
                src(j) = iota_k*lf%A(2) + lf%A(3)
            end do
            do m = -NFOURIER, NFOURIER - 1
                kpar = iota_k*m + mph
                if (abs(kpar) < 1.0e-8_dp) cycle
                gf%coef(m, k) = -sum(src*exp(-imun*m*pi*[(j, j=0, 2*NFOURIER - 1)] &
                    /NFOURIER))/(2*NFOURIER*imun*kpar)
            end do
        end do
        call local_field(s, 0.0_dp, lf)
    end subroutine build_gauge_fix

    function gauge_fix_eval(gf, th) result(dchi)
        !! Covariant grad(chi_par) at (s, th), without e^{i n phi}.
        use do_magfie_pert_mod, only: mph

        type(gauge_fix_t), intent(in) :: gf
        real(dp), intent(in) :: th
        complex(dp) :: dchi(3)
        complex(dp) :: ex(-NFOURIER:NFOURIER - 1), chi0
        integer :: m

        do m = -NFOURIER, NFOURIER - 1
            ex(m) = exp(imun*m*th)
        end do
        chi0 = sum(gf%coef(:, 2)*ex)
        dchi(1) = sum((gf%coef(:, 3) - gf%coef(:, 1))*ex)/(2.0_dp*GAUGE_DS)
        dchi(2) = sum(imun*[(m, m=-NFOURIER, NFOURIER - 1)]*gf%coef(:, 2)*ex)
        dchi(3) = imun*mph*chi0
    end function gauge_fix_eval

    pure function drift_velocity(vpar, mu, lf) result(vd)
        !! Contravariant (s, theta, phi) magnetic drift from the guiding-centre
        !! Lagrangian with A* = A + (m c/e) vpar h, B_s = 0, B*_par ~ B:
        !! v_d = [vpar (B* - B) + (c mu/e) h x grad B]/B.
        real(dp), intent(in) :: vpar, mu
        type(local_field_t), intent(in) :: lf
        real(dp) :: vd(3)
        real(dp) :: rhop, f, b2, gob_s, iob_s, gob_th, bob_th

        rhop = mi*c*vpar/qi
        f = 1.0_dp/(lf%sqrtg*lf%bmod)
        b2 = lf%bmod**2
        gob_s = dBphcovds/lf%bmod - Bphcov*lf%dB_ds/b2
        iob_s = dBthcovds/lf%bmod - Bthcov*lf%dB_ds/b2
        gob_th = -Bphcov*lf%dB_dth/b2
        bob_th = lf%dbeta_dth/lf%bmod - lf%beta*lf%dB_dth/b2
        vd(1) = f*(vpar*rhop*gob_th - (c*mu/qi)*Bphcov*lf%dB_dth/lf%bmod)
        vd(2) = f*(-vpar*rhop*gob_s + (c*mu/qi)*Bphcov*lf%dB_ds/lf%bmod)
        vd(3) = f*(vpar*rhop*(iob_s - bob_th) &
            + (c*mu/qi)*(lf%beta*lf%dB_dth - Bthcov*lf%dB_ds)/lf%bmod)
    end function drift_velocity

    subroutine line_bounce(v, eta_p, taub, omth, om_path, om_te, res, istate)
        !! Integrate the drive over one bounce or transit period taub, with the
        !! harmonic phase of timestep_transport (mth from driftorbit).
        use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK
        use neort_orbit, only: th0, evaluate_bfield_local, vpar
        use driftorbit, only: sign_vpar, sign_vpar_htheta

        real(dp), intent(in) :: v, eta_p, taub, omth, om_path, om_te
        type(line_harmonics_t), intent(out) :: res
        integer, intent(out) :: istate
        integer, parameter :: neq = 3 + 2*NCOMP + 2
        real(dp), parameter :: rtol = 1.0e-10_dp
        real(dp) :: y0(neq), atol(neq), bmod, htheta
        real(dp), allocatable :: yend(:)
        type(vode_state_t) :: vstate
        type(fortnum_status_t) :: status
        type(gauge_fix_t) :: gf

        call build_gauge_fix(gf)
        call evaluate_bfield_local(bmod, htheta)
        sign_vpar_htheta = sign(1.0_dp, htheta)*sign_vpar
        y0 = 0.0_dp
        y0(1) = th0
        y0(2) = sign_vpar_htheta*vpar(v, eta_p, bmod)
        atol = 1.0e-14_dp
        atol(2) = 1.0e-12_dp*v
        call vode_init(vstate, neq, 0.0_dp, y0)
        call vode_integrate_to(rhs, vstate, taub, rtol, atol, yend, status)
        istate = merge(2, -1, status%code == FORTNUM_OK)
        if (status%code /= FORTNUM_OK) then
            res = line_harmonics_t()
            return
        end if
        call unpack_harmonics(yend, res)

    contains

        subroutine rhs(t_, y_, dydt_, ctx_)
            real(dp), intent(in) :: t_
            real(dp), intent(in) :: y_(:)
            real(dp), intent(out) :: dydt_(:)
            class(*), intent(in), optional :: ctx_

            associate (dummy => ctx_)
            end associate
            call line_rhs(v, eta_p, taub, omth, om_path, om_te, gf, t_, y_, dydt_)
        end subroutine rhs
    end subroutine line_bounce

    subroutine line_rhs(v, eta_p, taub, omth, om_path, om_te, gf, t, y, ydot)
        use driftorbit, only: mth, etatp
        use do_magfie_pert_mod, only: mph

        real(dp), intent(in) :: v, eta_p, taub, omth, om_path, om_te, t
        type(gauge_fix_t), intent(in) :: gf
        real(dp), intent(in) :: y(:)
        real(dp), intent(out) :: ydot(:)
        type(local_field_t) :: lf
        complex(dp) :: loc(NCOMP), ph
        real(dp) :: drate, phase, wnorm
        integer :: j

        call line_point(v, eta_p, y(1), y(2), om_path, om_te, gf, lf, loc, drate)
        if (eta_p > etatp) then
            phase = mph*q*y(1) - mth*omth*t
        else
            phase = mph*q*y(1) - (mth + q*mph)*omth*t
        end if
        wnorm = 1.0_dp/(0.5_dp*mi*v**2*taub)
        ph = exp(imun*phase)*wnorm
        ydot(1) = y(2)*lf%hth
        ydot(2) = -0.5_dp*v**2*eta_p*lf%hth*lf%dB_dth
        ydot(3) = drate/taub
        do j = 1, NCOMP
            ydot(2 + 2*j) = real(loc(j)*ph)
            ydot(3 + 2*j) = aimag(loc(j)*ph)
        end do
        ydot(4 + 2*NCOMP) = 1.0_dp/(lf%bmod*taub)
        ydot(5 + 2*NCOMP) = lf%bmod/taub
    end subroutine line_rhs

    subroutine unpack_harmonics(y, res)
        real(dp), intent(in) :: y(:)
        type(line_harmonics_t), intent(out) :: res
        complex(dp) :: h(NCOMP)
        integer :: j

        do j = 1, NCOMP
            h(j) = cmplx(y(2 + 2*j), y(3 + 2*j), dp)
        end do
        res%h_line = h(1)
        res%h_boozer = h(2)
        res%h_fs = h(3)
        res%p_par = h(4)
        res%m_xi = h(5)
        res%chi = h(6)
        res%h_boozer_euler = h(7)
        res%h_naive = h(8)
        res%path = h(9)
        res%mu_dbe = h(10)
        res%om_drift = y(3)
        res%inv_b = y(4 + 2*NCOMP)
        res%b = y(5 + 2*NCOMP)
    end subroutine unpack_harmonics

    subroutine line_bounce_transport(v, eta_p, taub, omth, bounceavg, istate)
        !! Drop-in for bounce_fast(..., timestep_transport, ...) in drive_form
        !! 'line': fills bounceavg(3:4) with H_m/(m v^2/2) and (5:6) for nonlin.
        use neort_profiles, only: Om_tE
        use driftorbit, only: nonlin

        real(dp), intent(in) :: v, eta_p, taub, omth
        real(dp), intent(out) :: bounceavg(:)
        integer, intent(out) :: istate
        type(line_harmonics_t) :: res

        call line_bounce(v, eta_p, taub, omth, drift_frequency_part(v, eta_p, omth), &
            Om_tE, res, istate)
        bounceavg = 0.0_dp
        bounceavg(3) = real(res%h_line)
        bounceavg(4) = aimag(res%h_line)
        if (nonlin) then
            bounceavg(5) = res%inv_b
            bounceavg(6) = res%b
        end if
    end subroutine line_bounce_transport

    real(dp) function drift_frequency_part(v, eta_p, omth) result(om)
        !! Toroidal angular velocity of the thin path beyond field-line streaming:
        !! Om_tE + Om_tB, as in NEO-RT's canonical toroidal frequency.
        use neort_freq, only: Om_ph
        use driftorbit, only: etatp

        real(dp), intent(in) :: v, eta_p, omth
        real(dp) :: domdv, domdeta

        call Om_ph(v, eta_p, om, domdv, domdeta)
        if (eta_p <= etatp) om = om - omth/iota
    end function drift_frequency_part

    complex(dp) function line_scalar_eps(th) result(epsn)
        !! Relative Boozer scalar dB_L/B of the active pert_model at (s, th).
        real(dp), intent(in) :: th
        type(local_field_t) :: lf

        call local_field(s, th, lf)
        epsn = (lf%dBE + lf%xigradB)/lf%bmod
    end function line_scalar_eps

    subroutine exact_orbit_period(v, eta_p, taub_est, taub, omth)
        !! Refine the bounce/transit period inside the spline estimate's bracket.
        !! Fresh VODE shots have a finite closure floor; certify the bracket,
        !! endpoint state and velocity sign instead of a sub-noise Newton step.
        use neort_orbit, only: th0, evaluate_bfield_local, vpar
        use driftorbit, only: sign_vpar, sign_vpar_htheta, etatp

        real(dp), intent(in) :: v, eta_p, taub_est
        real(dp), intent(out) :: taub, omth
        real(dp) :: y(2), yinitial(2), target, bmod, htheta, lower, upper
        real(dp) :: f_lower, f_upper, residual
        integer :: it

        call evaluate_bfield_local(bmod, htheta)
        sign_vpar_htheta = sign(1.0_dp, htheta)*sign_vpar
        if (taub_est <= 0.0_dp) error stop "exact_orbit_period: nonpositive estimate"
        target = th0
        if (eta_p <= etatp) target = th0 + sign(2.0_dp*pi, sign_vpar_htheta)
        yinitial = [th0, sign_vpar_htheta*vpar(v, eta_p, bmod)]
        lower = 0.7_dp*taub_est
        upper = 1.3_dp*taub_est
        y = yinitial
        call poloidal_state_at(v, eta_p, lower, y)
        f_lower = y(1) - target
        y = yinitial
        call poloidal_state_at(v, eta_p, upper, y)
        f_upper = y(1) - target
        if (f_lower*f_upper >= 0.0_dp) then
            error stop "exact_orbit_period: period not bracketed within 30% of estimate"
        end if
        do it = 1, 80
            taub = 0.5_dp*(lower + upper)
            y = yinitial
            call poloidal_state_at(v, eta_p, taub, y)
            residual = y(1) - target
            if (f_lower*residual <= 0.0_dp) then
                upper = taub
            else
                lower = taub
                f_lower = residual
            end if
            if (upper - lower <= 1.0e-12_dp*taub_est) exit
        end do
        if (it > 80) error stop "exact_orbit_period: bracket did not converge"
        if (abs(residual) > 1.0e-10_dp*(1.0_dp + abs(target))) then
            error stop "exact_orbit_period: orbit angle did not close"
        end if
        if (y(2)*yinitial(2) <= 0.0_dp) then
            error stop "exact_orbit_period: wrong bounce branch"
        end if
        if (abs(y(2) - yinitial(2)) > 1.0e-10_dp*v) then
            error stop "exact_orbit_period: parallel velocity did not close"
        end if
        omth = sign(2.0_dp*pi/taub, sign_vpar_htheta)
    end subroutine exact_orbit_period

    subroutine poloidal_state_at(v, eta_p, t_end, y)
        use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK

        real(dp), intent(in) :: v, eta_p, t_end
        real(dp), intent(inout) :: y(2)
        real(dp), allocatable :: yend(:)
        type(vode_state_t) :: vstate
        type(fortnum_status_t) :: status

        call vode_init(vstate, 2, 0.0_dp, y)
        call vode_integrate_to(rhs, vstate, t_end, 1.0e-13_dp, &
            [1.0e-15_dp, 1.0e-13_dp*v], yend, status)
        if (status%code /= FORTNUM_OK) then
            error stop "exact_orbit_period: orbit integration failed"
        end if
        y = yend
    contains
        subroutine rhs(t_, y_, dydt_, ctx_)
            real(dp), intent(in) :: t_
            real(dp), intent(in) :: y_(:)
            real(dp), intent(out) :: dydt_(:)
            class(*), intent(in), optional :: ctx_

            associate (dummy => ctx_, dummy_t => t_)
            end associate
            call poloidal_rate(v, eta_p, y_, dydt_)
        end subroutine rhs
    end subroutine poloidal_state_at

    subroutine poloidal_rate(v, eta_p, y, ydot)
        real(dp), intent(in) :: v, eta_p, y(:)
        real(dp), intent(out) :: ydot(:)
        real(dp) :: bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

        call do_magfie([s, 0.0_dp, y(1)], bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        ydot(1) = y(2)*hctrvr(3)
        ydot(2) = -0.5_dp*v**2*eta_p*hctrvr(3)*bder(3)*bmod
    end subroutine poloidal_rate

end module neort_line_drive
