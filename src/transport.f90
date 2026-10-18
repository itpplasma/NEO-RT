module neort_transport
    use neort_mixture_collisions, only: finite_krook_active, effective_krook_rate
    use iso_fortran_env, only: dp => real64, error_unit
    use util, only: imun, pi, c, qi
    use logger, only: trace, debug, warning, error
    use do_magfie_mod, only: do_magfie, s, a, R0, iota, q, psi_pr, eps, &
        bphcov, dbthcovds, dbphcovds, q, dqds, sign_theta, Bthcov
    use do_magfie_pert_mod, only: do_magfie_pert_amp, do_magfie_pert_kappa_amp, curvature_drive
    use neort_magfie, only: dVds, B0
    use neort_profiles, only: ns, Om_tE
    use neort_species, only: selected_species
    use neort_nonlin, only: nonlinear_attenuation
    use neort_freq, only: Om_th, Om_ph, Om_tB
    use neort_orbit, only: bounce_fast, nvar, poloidal_velocity, Om_tB_local_v2, &
        bounce_rtol, bounce_atol, th0, vpar
    use neort_resonance, only: driftorbit_coarse, driftorbit_root, resonance_value
    use neort_krook, only: krook_collisions, krook_nu, krook_model, krook_harmonic, &
        krook_geometry_trapped, krook_geometry_passing, krook_eta_steps, &
        krook_max_depth, krook_rtol, krook_atol, krook_rate, integrate_krook
    use driftorbit, only: vth, mth, mph, mi, B0, Bmin, Bmax, comptorque, epsmn, &
        pertfile_scale, &
        etamin, etamax, A1, A2, nlev, pertfile, nonlin, m0, etatp, etadt, &
        sign_vpar_htheta, sign_vpar, drift_phase, magdrift, magdrift_passing, &
        fow_radial

  implicit none

  ! These values remain live while bounce_fast invokes timestep_transport.
  ! Each transport worker needs an independent callback state.
  real(dp) :: Omth, dOmthdv, dOmthdeta
  !$omp threadprivate (Omth, dOmthdv, dOmthdeta)

  abstract interface
    function drift_osc_i(t) result(Om_osc)
      import :: dp
      real(dp), intent(in) :: t
      real(dp) :: Om_osc
    end function drift_osc_i
  end interface

  ! Test seam: when associated (and drift_phase is set), replaces the local
  ! oscillating toroidal drift Om_tB(theta) - <Om_tB> [rad/s] by a prescribed
  ! function of orbit time, so tests can use analytic Bessel sidebands.
  procedure(drift_osc_i), pointer :: drift_osc_override => null()

contains

  pure function fmt_dbg(msg1, v1, msg2, v2, msg3, v3, msg4, v4) result(s)
    ! Helper to compose a short debug line
    character(*), intent(in) :: msg1, msg2
    character(*), intent(in), optional :: msg3, msg4
    real(dp), intent(in) :: v1, v2
    real(dp), intent(in), optional :: v3, v4
    character(len=256) :: s
    character(len=64) :: a1, a2, a3, a4
    a3 = ''; a4 = ''
    write(a1,'(ES12.5)') v1
    write(a2,'(ES12.5)') v2
    if (present(v3)) write(a3,'(ES12.5)') v3
    if (present(v4)) write(a4,'(ES12.5)') v4
    if (present(msg4)) then
      s = trim(msg1)//trim(a1)//' '//trim(msg2)//trim(a2)//' '//trim(msg3)//trim(a3)//' '//trim(msg4)//trim(a4)
    else if (present(msg3)) then
      s = trim(msg1)//trim(a1)//' '//trim(msg2)//trim(a2)//' '//trim(msg3)//trim(a3)
    else
      s = trim(msg1)//trim(a1)//' '//trim(msg2)//trim(a2)
    end if
  end function fmt_dbg

! original contains follows

    pure function D11int(ux, taub, Hmn2)
        real(dp) :: D11int
        real(dp), intent(in) :: ux, taub, Hmn2

        D11int = pi**(3.0_dp / 2.0_dp) * mph**2 * c**2 * q * vth &
                 / (qi**2 * dVds * abs(psi_pr)) * ux**3 * exp(-ux**2) * taub * Hmn2
    end function D11int

    pure function D12int(ux, taub, Hmn2)
        real(dp) :: D12int
        real(dp), intent(in) :: ux, taub, Hmn2

        D12int = D11int(ux, taub, Hmn2) * ux**2
    end function D12int

    pure function Tphi_int(ux, taub, Hmn2)
        real(dp) :: Tphi_int
        real(dp), intent(in) :: ux, taub, Hmn2

        Tphi_int = sign(1.0_dp, psi_pr * q * sign_theta) * pi**(3.0_dp / 2.0_dp) * mph**2 * ns * &
                   c * vth / qi &
                   * ux**3 * exp(-ux**2) * taub * Hmn2 * (A1 + A2 * ux**2)
    end function Tphi_int

    subroutine compute_transport_integral(vmin, vmax, vsteps, D, T)
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(dp), intent(out) :: D(2), T
        if (finite_krook_active()) then
            call compute_transport_krook(vmin, vmax, vsteps, D, T)
        else
            call compute_transport_collisionless(vmin, vmax, vsteps, D, T)
        end if
    end subroutine compute_transport_integral

    subroutine compute_transport_collisionless(vmin, vmax, vsteps, D, T)
        ! compute transport integral via midpoint rule
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(dp), intent(out) :: D(2), T  ! Transport coefficients D and torque density T
        real(dp) :: D_plateau, dsdreff  ! Plateau diffusion coefficient and ds/dreff=<|grad s|>
        real(dp) :: ux, du, dD11, dD12, dT, v, eta
        real(dp) :: eta_res(2)
        real(dp) :: taub, bounceavg(nvar)
        integer :: istate_dv
        real(dp) :: Hmn2, attenuation_factor
        real(dp) :: roots(nlev, 3)
        integer :: nroots, kr, ku

        call debug(fmt_dbg('compute_transport_integral: vmin=', vmin, ' vmax=', vmax, ' vsteps=', dble(vsteps)))

        D = 0.0_dp
        T = 0.0_dp
        du = (vmax - vmin) / (vsteps * vth)
        ux = vmin / vth + du / 2.0_dp

        do ku = 1, vsteps
            v = ux * vth
            call driftorbit_coarse(v, etamin, etamax, roots, nroots)
            ! No explicit nroots==0 guard: the do-loop below is empty when
            ! nroots==0, and a bare `cycle` here would skip the trailing
            ! `ux = ux + du` velocity-grid increment and stall the sweep.
            do kr = 1, nroots
                eta_res = driftorbit_root(v, 1.0e-8_dp * abs(Om_tE), roots(kr, 1), roots(kr, 2))
                if (eta_res(1) < 0.0_dp) cycle  ! bracket-failure sentinel
                eta = eta_res(1)

                call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)

                taub = 2.0_dp * pi / abs(Omth)
                call transport_bounce(v, eta, taub, bounceavg, istate_dv)
                if (istate_dv == -1) then
                    call error(fmt_dbg('VODE MXSTEP: mth=', dble(mth), ' ux=', ux, ' eta=', eta, ' taub=', taub))
                else if (istate_dv /= 2) then
                    call warning(fmt_dbg('dvode istate=', dble(istate_dv), ' at mth=', dble(mth), ' ux=', ux, ' eta=', eta))
                else
                    if (abs(eta - etatp) < 1.0e-8_dp*etatp) then
                        call trace(fmt_dbg('near etatp: mth=', dble(mth), ' ux=', ux, ' eta=', eta, ' taub=', taub))
                    end if
                end if
                Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2) * (mi * (ux * vth)**2 / 2.0_dp)**2
                attenuation_factor = nonlinear_attenuation(ux, eta, bounceavg, Omth, &
                                                           dOmthdv, dOmthdeta, Hmn2)

                dD11 = du * D11int(ux, taub, Hmn2) / abs(eta_res(2))
                dD12 = du * D12int(ux, taub, Hmn2) / abs(eta_res(2))
                D(1) = D(1) + dD11 * attenuation_factor
                D(2) = D(2) + dD12 * attenuation_factor

                if (comptorque) then
                    dT = du * Tphi_int(ux, taub, Hmn2) / abs(eta_res(2))
                    T = T + dT * attenuation_factor
                end if
            end do
            ux = ux + du
        end do

        D_plateau = pi * vth**3 / (16.0_dp * R0 * iota * (qi * B0 / (mi * c))**2)
        dsdreff = 2.0_dp / a * sqrt(s)  ! TODO: Use exact value instead of this approximation
        D = dsdreff**(-2) * D / D_plateau

        call debug("compute_transport_integral complete")

    end subroutine compute_transport_collisionless

    subroutine compute_transport_krook(vmin, vmax, vsteps, D, T)
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(dp), intent(out) :: D(2), T
        real(dp) :: du, ux, values(3), integral(3), D_plateau, dsdreff
        integer :: ku
        integral = 0.0_dp
        du = (vmax-vmin)/(vsteps*vth)
        ux = vmin/vth + du/2.0_dp
        do ku = 1, vsteps
            call compute_krook_velocity(ux, values)
            integral = integral + du*values
            ux = ux + du
        end do
        D_plateau = pi*vth**3/(16.0_dp*R0*iota*(qi*B0/(mi*c))**2)
        dsdreff = 2.0_dp/a*sqrt(s)
        D = dsdreff**(-2)*integral(1:2)/D_plateau
        T = integral(3)
    end subroutine compute_transport_krook

    subroutine compute_krook_velocity(ux, values)
        real(dp), intent(in) :: ux
        real(dp), intent(out) :: values(3)
        real(dp) :: effective_harmonic, geometry, nu
        integer :: status
        effective_harmonic = real(mth, dp)
        geometry = krook_geometry_trapped
        if (etamax < etatp) then
            geometry = krook_geometry_passing
            if (krook_harmonic == 1) effective_harmonic = effective_harmonic + q*mph
        end if
        nu = effective_krook_rate(ux, effective_harmonic, geometry, ux*vth*0.01_dp)
        call integrate_krook(residual, weight, etamin, etamax, nu, krook_eta_steps, &
            krook_rtol, krook_atol, values, status, krook_max_depth)
        if (status /= 0) call report_krook_failure(status, ux, nu, &
            effective_harmonic, geometry)
    contains
        subroutine residual(eta, res, derivative)
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: res, derivative
            call resonance_value(ux*vth, eta, res, derivative)
        end subroutine residual
        subroutine weight(eta, result_)
            real(dp), intent(in) :: eta
            real(dp), intent(out) :: result_(:)
            call compute_transport_weight(ux, eta, result_)
        end subroutine weight
    end subroutine compute_krook_velocity

    subroutine report_krook_failure(status, ux, nu, effective_harmonic, geometry)
        integer, intent(in) :: status
        real(dp), intent(in) :: ux, nu, effective_harmonic, geometry
        ! Preserve the returned status and the consuming thread's physical state.
        ! The eta values are integration bounds, not a guessed failing node.
        !$omp critical(krook_failure_context)
        write(error_unit, '(a,6(a,i0),*(a,es25.17e3))') &
            'KROOK_QUADRATURE_FAILURE', ' status=', status, &
            ' species_index=', selected_species, ' mth=', mth, ' mph=', mph, &
            ' eta_steps=', krook_eta_steps, ' max_depth=', krook_max_depth, &
            ' s=', s, ' mass_g=', mi, ' charge_statc=', qi, &
            ' density_cm-3=', ns, ' omega_E_s-1=', Om_tE, ' ux=', ux, &
            ' vth_cm_s=', vth, ' nu_s-1=', nu, ' eta_min=', etamin, &
            ' eta_max=', etamax, ' eta_tp=', etatp, &
            ' sign_vpar=', sign_vpar, &
            ' effective_harmonic=', effective_harmonic, ' geometry=', geometry, &
            ' rtol=', krook_rtol, ' atol=', krook_atol
        flush(error_unit)
        !$omp end critical(krook_failure_context)
        error stop 'finite-width Krook quadrature did not converge'
    end subroutine report_krook_failure

    subroutine compute_transport_weight(ux, eta, values)
        real(dp), intent(in) :: ux, eta
        real(dp), intent(out) :: values(:)
        real(dp) :: taub, bounceavg(nvar), Hmn2
        integer :: istate
        call Om_th(ux*vth, eta, Omth, dOmthdv, dOmthdeta)
        taub = 2.0_dp*pi/abs(Omth)
        call transport_bounce(ux*vth, eta, taub, bounceavg, istate)
        if (istate /= 2) error stop 'orbit integration failed in Krook quadrature'
        Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2)*(mi*(ux*vth)**2/2.0_dp)**2
        values(1) = D11int(ux, taub, Hmn2)
        values(2) = D12int(ux, taub, Hmn2)
        values(3) = 0.0_dp
        if (comptorque) values(3) = Tphi_int(ux, taub, Hmn2)
    end subroutine compute_transport_weight

    subroutine transport_bounce(v, eta, taub, bounceavg, istate)
        ! Bounce averages of the transport weight, H_mn**2 = |<H_n>|**2.
        ! int H_n dt, with the dimensionless H_n = drive*exp(i*phase) (transport_rates),
        ! is an oscillatory, often strongly cancelling integral of size up to
        ! taub*max|H_n|. The fixed bounce_atol = 1e-10 can exceed that scale
        ! (small eps_n, short taub) and then leaves these components without
        ! error control: H_mn**2 jumped by O(1) between neighbouring eta as the
        ! DVODE step sequence changed. Their absolute tolerance is therefore
        ! bounce_rtol relative to taub*max|H_n|, which makes the weight (to rounding)
        ! quadratic in the perturbation amplitude.
        integer, parameter :: samples = 8
        real(dp), intent(in) :: v, eta, taub
        real(dp), intent(out) :: bounceavg(nvar)
        integer, intent(out) :: istate
        real(dp) :: atol(nvar), y(nvar), ydot(nvar), scale
        integer :: k
        atol = bounce_atol
        scale = 0.0_dp
        y = 0.0_dp
        do k = 0, samples - 1
            y(1) = th0 + 2.0_dp*pi*real(k, dp)/real(samples, dp)
            call timestep_transport(v, eta, nvar, 0.0_dp, y, ydot)
            scale = max(scale, abs(cmplx(ydot(3), ydot(4), dp)))
        end do
        scale = bounce_rtol*taub*scale
        if (scale > 0.0_dp) atol(3:4) = min(bounce_atol, scale)
        if (eta > etatp) then
            call trapped_bounce(v, eta, taub, bounceavg, istate, atol)
        else
            call bounce_fast(v, eta, taub, bounceavg, timestep_transport, istate, atol)
        end if
    end subroutine transport_bounce

    subroutine trapped_bounce(v, eta, taub, bounceavg, istate, atol_in)
        ! Bounce averages of a trapped orbit over its closed orbit.
        ! Time integration of (theta, v_parallel) accumulates an energy error
        ! which near the trapped-passing boundary is amplified by
        ! 1/(eta*B_max - 1) through the time spent near the X-point; the weight
        ! then became noisy at the 1e-4 level and the Krook pitch quadrature
        ! did not converge. Here the orbit is parametrised by the angle chi,
        !   theta = theta_m - direction*Delta*cos(chi),
        ! between the turning points theta_m -+ Delta, and t(chi) is integrated
        ! together with the bounce integrands. v_parallel follows from the
        ! energy invariant, and dt/dchi is regular at the turning points.
        use fortnum_quadrature, only: gauss_legendre
        use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
            FORTNUM_CONVERGENCE_ERROR
        ! Below this 1 - eta*B is integrated from the turning point: the direct
        ! difference has relative rounding error eps/(1 - eta*B), which would
        ! break the ODE step control next to the turning points.
        real(dp), parameter :: gap_switch = 1.0e-6_dp
        real(dp), intent(in) :: v, eta, taub, atol_in(nvar)
        real(dp), intent(out) :: bounceavg(nvar)
        integer, intent(out) :: istate
        real(dp) :: lower, upper, middle, half, direction, sign_htheta, chi0, saved_omth
        integer :: pass
        logical :: averages
        real(dp) :: nodes(8), weights(8), y0(nvar), atol(nvar), position(3)
        real(dp) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp), allocatable :: yend(:)
        type(vode_state_t) :: vstate
        type(fortnum_status_t) :: status

        call gauss_legendre(size(nodes), nodes, weights)
        direction = sign(1.0_dp, sign_vpar)
        position = [s, 0.0_dp, th0]
        call do_magfie(position, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        sign_htheta = sign(1.0_dp, hctrvr(3))
        lower = turning_point(-1.0_dp)
        upper = turning_point(1.0_dp)
        middle = (lower + upper)/2.0_dp
        half = (upper - lower)/2.0_dp
        chi0 = acos(max(-1.0_dp, min(1.0_dp, direction*(middle - th0)/half)))
        atol = atol_in
        atol(1) = bounce_rtol*taub
        bounceavg = 0.0_dp
        ! The bounce harmonics exp(-i*mth*Omth*t) must close on the orbit:
        ! a first pass gives its period, which replaces the interpolated
        ! bounce frequency in the phase during the second pass.
        saved_omth = Omth
        do pass = 1, 2
            averages = pass == 2
            y0 = 0.0_dp
            call vode_init(vstate, nvar, chi0, y0)
            call vode_integrate_to(chi_rhs, vstate, chi0 + 2.0_dp*pi, bounce_rtol, &
                                   atol, yend, status)
            if (status%code == FORTNUM_OK) then
                istate = 2
            else if (status%code == FORTNUM_CONVERGENCE_ERROR) then
                istate = -1
            else
                istate = 0
            end if
            if (istate /= 2) exit
            Omth = sign(2.0_dp*pi/yend(1), saved_omth)
        end do
        Omth = saved_omth
        if (istate /= 2) return
        bounceavg(3:) = yend(3:)/yend(1)

    contains

        subroutine chi_rhs(chi, y, dydchi, rhs_ctx)
            real(dp), intent(in) :: chi
            real(dp), intent(in) :: y(:)
            real(dp), intent(out) :: dydchi(:)
            class(*), intent(in), optional :: rhs_ctx
            real(dp) :: theta, gap, dtdchi, distance, slope, orbit(nvar), rate(nvar)
            real(dp) :: bmod_, sqrtg_, hder_(3), hcovar_(3), hctrvr_(3), hcurl_(3), x(3)
            associate (dummy => rhs_ctx)
            end associate
            theta = middle - direction*half*cos(chi)
            x = [s, 0.0_dp, theta]
            call do_magfie(x, bmod_, sqrtg_, hder_, hcovar_, hctrvr_, hcurl_)
            gap = 1.0_dp - eta*bmod_
            if (gap > gap_switch) then
                dtdchi = half*abs(sin(chi))/(v*sqrt(gap)*abs(hctrvr_(3)))
            else if (cos(chi) > 0.0_dp) then
                ! Turning point at chi = 0, distance 2*Delta*sin(chi/2)**2.
                distance = 2.0_dp*half*sin(chi/2.0_dp)**2
                slope = mean_slope(theta, middle - direction*half)
                gap = eta*slope*distance
                dtdchi = sqrt(2.0_dp*half/(eta*slope))*abs(cos(chi/2.0_dp))/ &
                    (v*abs(hctrvr_(3)))
            else
                distance = 2.0_dp*half*cos(chi/2.0_dp)**2
                slope = mean_slope(theta, middle + direction*half)
                gap = eta*slope*distance
                dtdchi = sqrt(2.0_dp*half/(eta*slope))*abs(sin(chi/2.0_dp))/ &
                    (v*abs(hctrvr_(3)))
            end if
            dydchi = 0.0_dp
            dydchi(1) = dtdchi
            if (.not. averages) return
            orbit = 0.0_dp
            orbit(1) = theta
            orbit(2) = direction*sign(1.0_dp, sin(chi))*sign_htheta*v*sqrt(gap)
            orbit(3:) = y(3:)
            call transport_rates(v, eta, nvar, y(1), orbit, bmod_, hder_, hcovar_, &
                                 hctrvr_, rate)
            dydchi(3:) = rate(3:)*dtdchi
        end subroutine chi_rhs

        real(dp) function mean_slope(theta, turning)
            ! |dB/dtheta| averaged between theta and the turning point, so
            ! that 1 - eta*B(theta) = eta*mean_slope*|turning - theta|.
            real(dp), intent(in) :: theta, turning
            real(dp) :: x(3), bmod_, sqrtg_, hder_(3), hcovar_(3), hctrvr_(3), hcurl_(3)
            integer :: k
            mean_slope = 0.0_dp
            x(1:2) = [s, 0.0_dp]
            do k = 1, size(nodes)
                x(3) = theta + (turning - theta)*(nodes(k) + 1.0_dp)/2.0_dp
                call do_magfie(x, bmod_, sqrtg_, hder_, hcovar_, hctrvr_, hcurl_)
                mean_slope = mean_slope + weights(k)*bmod_*hder_(3)
            end do
            mean_slope = abs(mean_slope)/2.0_dp
        end function mean_slope

        real(dp) function turning_point(side)
            ! First root of 1 - eta*B from th0 towards side. A grid scan finds
            ! a sign change or a local maximum of B that reflects the orbit;
            ! bisection then resolves the root to rounding.
            integer, parameter :: steps = 256
            real(dp), intent(in) :: side
            real(dp) :: inside, outside, theta, gap, slope, previous_slope
            integer :: k
            inside = th0
            call gap_and_slope(inside, side, gap, previous_slope)
            do k = 1, steps
                theta = th0 + side*2.0_dp*pi*real(k, dp)/real(steps, dp)
                call gap_and_slope(theta, side, gap, slope)
                if (gap < 0.0_dp) then
                    outside = theta
                    exit
                end if
                if (previous_slope > 0.0_dp .and. slope <= 0.0_dp) then
                    outside = local_maximum(inside, theta, side)
                    call gap_and_slope(outside, side, gap, slope)
                    if (gap < 0.0_dp) exit
                end if
                inside = theta
                previous_slope = slope
            end do
            if (k > steps) error stop 'trapped orbit without turning point'
            do
                theta = (inside + outside)/2.0_dp
                if (theta == inside .or. theta == outside) exit
                call gap_and_slope(theta, side, gap, slope)
                if (gap < 0.0_dp) then
                    outside = theta
                else
                    inside = theta
                end if
            end do
            turning_point = inside
        end function turning_point

        real(dp) function local_maximum(a, b, side)
            real(dp), intent(in) :: a, b, side
            real(dp) :: left, right, theta, gap, slope
            left = a
            right = b
            do
                theta = (left + right)/2.0_dp
                if (theta == left .or. theta == right) exit
                call gap_and_slope(theta, side, gap, slope)
                if (slope > 0.0_dp) then
                    left = theta
                else
                    right = theta
                end if
            end do
            local_maximum = left
        end function local_maximum

        subroutine gap_and_slope(theta, side, gap, slope)
            ! 1 - eta*B and dB/dtheta along side at theta.
            real(dp), intent(in) :: theta, side
            real(dp), intent(out) :: gap, slope
            real(dp) :: x(3), bmod_, sqrtg_, hder_(3), hcovar_(3), hctrvr_(3), hcurl_(3)
            x = [s, 0.0_dp, theta]
            call do_magfie(x, bmod_, sqrtg_, hder_, hcovar_, hctrvr_, hcurl_)
            gap = 1.0_dp - eta*bmod_
            slope = side*bmod_*hder_(3)
        end subroutine gap_and_slope
    end subroutine trapped_bounce

    subroutine timestep_transport(v, eta, neq, t, y, ydot)
        !
        !  Timestep function for orbit integration.
        !  Includes poloidal angle theta and parallel velocity.
        !  More integrands may be added starting from y(3)
        !

        real(dp), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq)
        real(dp), intent(out) :: ydot(neq)
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        call transport_rates(v, eta, neq, t, y, bmod, hder, hcovar, hctrvr, ydot)
    end subroutine timestep_transport

    subroutine transport_rates(v, eta, neq, t, y, bmod, hder, hcovar, hctrvr, ydot)
        ! timestep_transport for the field at theta = y(1).
        real(dp), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq), bmod, hder(3), hcovar(3), hctrvr(3)
        real(dp), intent(out) :: ydot(neq)

        ! BEGIN TODO: remove all of this after refactoring and re-use routine in orbit
        ! for y(1:3)
        real(dp) :: x(3), Om_tB_v
        real(dp) :: t0, Om_osc, OmtB_avg, dOmtBdv, dOmtBdeta, Omph_avg
        real(dp) :: dOmphdv, dOmphdeta
        real(dp) :: phase_drift, bmod1, vparallel, x1(3), sqrtg1, hder1(3), hcovar1(3), &
                    hctrvr1(3), hcurl1(3)
        complex(dp) :: epsn, Hn  ! relative amplitude of perturbation field epsn=Bn/B0
        ! and Hamiltonian Hn = (H - H0)_n
        complex(dp) :: kappan, phase_factor  ! kappa.xi_perp and the phase of Hn
        real(dp) :: w_lagrangian, w_curvature

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)
        vparallel = y(2)
        if (eta < etatp) then
            ! A passing orbit never reverses. Enforce its energy invariant:
            ! evolving v_parallel accumulates an ODE energy error which is
            ! amplified by 1/(1-eta*B) near the trapped-passing boundary.
            vparallel = sign(vpar(v, eta, bmod), y(2))
        end if
        call poloidal_velocity(v, eta, bmod, hctrvr(3), hder(3), vparallel, ydot)

        ! evaluate orbit averages of Hamiltonian perturbation, optionally at the
        ! first-order guiding-centre radius (finite radial orbit width). Only the
        ! perturbation amplitude and bmod are taken there; everything else below
        ! (Hn phase, drift phase) stays on the label surface s.
        x1 = x
        bmod1 = bmod
        if (fow_radial) then
            x1(1) = s + fow_displacement(vparallel, hcovar(2))
            call do_magfie_displaced(x1, bmod1, sqrtg1, hder1, hcovar1, hctrvr1, hcurl1)
        end if
        if (pertfile) then
            call do_magfie_pert_amp(x1, epsn)
            epsn = pertfile_scale * epsn / bmod1
        else
            epsn = epsmn * exp(imun * m0 * y(1))
        end if

        ! Oscillating part of the toroidal drift phase,
        ! int (Om_tB(theta) - <Om_tB>) dt, accumulated in y(7).  Its toroidal
        ! phase is mph times this integral; <Om_tB> is the value that enters
        ! the resonance condition (Om_ph - Om_tE for trapped orbits).
        Om_osc = 0.0_dp
        phase_drift = 0.0_dp
        if (drift_phase) then
            if (associated(drift_osc_override)) then
                Om_osc = drift_osc_override(t)
            else if (magdrift .and. (eta > etatp .or. magdrift_passing > 0)) then
                if (eta > etatp) then
                    call Om_ph(v, eta, Omph_avg, dOmphdv, dOmphdeta)
                    OmtB_avg = Omph_avg - Om_tE
                else
                    call Om_tB(v, eta, OmtB_avg, dOmtBdv, dOmtBdeta)
                end if
                Om_tB_v = Om_tB_local_v2(eta, bmod, hder(1), hctrvr(3))
                Om_osc = Om_tB_v * v**2 - OmtB_avg
            end if
            phase_drift = mph * y(7)
        end if

        if (abs(t) <= 0.0_dp) then
            ! Local Hamiltonian samples have no bounce-time phase and need no
            ! frequency installed by a previous orbit integration.
            phase_factor = exp(imun * (q * mph * y(1) + phase_drift))
        else if (eta > etatp) then
            !t0 = 0.25*2*pi/Omth ! Different starting position in orbit
            t0 = 0.0_dp
            phase_factor = exp(imun * (q * mph * (y(1)) - mth * (t - t0) * Omth + phase_drift))
        else
            phase_factor = exp(imun * (q * mph * (y(1)) - (mth + q * mph) * t * Omth + phase_drift))
        end if

        ! First-order guiding-centre Hamiltonian per unit energy (doc/hamiltonian_drive.md).
        ! Default (Boozer form, perturbed field Boozer-representable): (2 - eta B) dB_L/B.
        ! With a curvature table the general ideal drive for deltaA = xi_perp x B0,
        !   H/E = eta B dB_L/B + 2 (1 - eta B) kappa.xi_perp,
        ! which reduces to the default for kappa.xi_perp = dB_L/B.
        if (curvature_drive) then
            call do_magfie_pert_kappa_amp(x1, kappan)
            kappan = pertfile_scale * kappan / bmod1
            call drive_weight_leaf(eta * bmod1, w_lagrangian, w_curvature)
            Hn = (w_lagrangian * epsn + w_curvature * kappan) * phase_factor
        else
            Hn = (2.0_dp - eta * bmod1) * epsn * phase_factor
        end if
        ydot(3) = real(Hn)
        ydot(4) = aimag(Hn)

        ! evaluate orbit averages for nonlinear attenuation
        if (nonlin) then
            ydot(5) = 1.0_dp / bmod
            ydot(6) = bmod
        else
            ydot(5:6) = 0.0_dp
        end if

        ! y(7) integrates the oscillating drift phase (zero unless drift_phase).
        ! Zero any trailing integrands this routine does not compute. The DOP853 solver carries every component through
        ! its stage combinations, so an uninitialised derivative leaves a
        ! denormal in the state that trips the underflow FPE trap.
        ydot(7) = Om_osc
        ydot(8:) = 0.0_dp
    end subroutine transport_rates

    function fow_displacement(v_par, hphi) result(ds)
        ! First-order radial guiding-centre displacement s(t) - s from the thin
        ! orbit on s. Conservation of p_phi = m v_par h_phi + (e/c) A_phi(s) with
        ! A_phi' = -sign_theta*iota*psi_pr (x = (s, phi, theta) has Jacobian
        ! -sqrtg, test_fow_gc_orbit) gives
        !   Delta s = m c q v_par h_phi / (e sign_theta psi_pr).
        ! s is the p_phi label, (e/c) A_phi(s) = p_phi, for trapped and passing
        ! orbits alike. This is the label of NEO-RT's frequencies: the shear term
        ! of Om_tB is q' <Delta s thetadot> with this Delta s. The orbit mean is
        ! therefore not removed; for passing orbits it is <Delta s> ~ +-rho_pol.
        real(dp), intent(in) :: v_par, hphi
        real(dp) :: ds

        ds = mi * c * q * v_par * hphi / (qi * sign_theta * psi_pr)
    end function fow_displacement

    subroutine do_magfie_displaced(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        ! do_magfie at a displaced radius without side effects: do_magfie
        ! overwrites the flux-surface module quantities of the label surface s,
        ! which the caller and drift_phase (Om_tB_local_v2) still need. (The
        ! standalone B0h is also overwritten; nothing reads it after init.)
        real(dp), intent(in) :: x(3)
        real(dp), intent(out) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: saved(7)

        saved = [Bthcov, Bphcov, dBthcovds, dBphcovds, q, dqds, iota]
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        Bthcov = saved(1)
        Bphcov = saved(2)
        dBthcovds = saved(3)
        dBphcovds = saved(4)
        q = saved(5)
        dqds = saved(6)
        iota = saved(7)
    end subroutine do_magfie_displaced

    include 'generated/drive_leaf.inc'

end module neort_transport
