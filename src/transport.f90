module neort_transport
    use iso_fortran_env, only: dp => real64
    use util, only: imun, pi, c, qi
    use logger, only: trace, debug, warning, error
    use do_magfie_mod, only: do_magfie, s, a, R0, iota, q, psi_pr, eps, &
        bphcov, dbthcovds, dbphcovds, q, dqds, sign_theta, Bthcov
    use do_magfie_pert_mod, only: do_magfie_pert_amp
    use neort_magfie, only: dVds, B0
    use neort_profiles, only: ni1, Om_tE
    use neort_nonlin, only: nonlinear_attenuation
    use neort_freq, only: Om_th, Om_ph
    use neort_orbit, only: bounce_fast, nvar, noshear, poloidal_velocity
    use neort_resonance, only: driftorbit_coarse, driftorbit_root
    use driftorbit, only: vth, mth, mph, mi, B0, Bmin, Bmax, comptorque, epsmn, &
        etamin, etamax, A1, A2, nlev, pertfile, nonlin, m0, etatp, etadt, &
        sign_vpar_htheta, sign_vpar

  implicit none

  real(dp) :: Omth, dOmthdv, dOmthdeta

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

        Tphi_int = sign(1.0_dp, psi_pr * q * sign_theta) * pi**(3.0_dp / 2.0_dp) * mph**2 * ni1 * &
                   c * vth / qi &
                   * ux**3 * exp(-ux**2) * taub * Hmn2 * (A1 + A2 * ux**2)
    end function Tphi_int

    subroutine compute_transport_integral(vmin, vmax, vsteps, D, T)
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
            if (nroots == 0) continue
            do kr = 1, nroots
                eta_res = driftorbit_root(v, 1.0e-8_dp * abs(Om_tE), roots(kr, 1), roots(kr, 2))
                eta = eta_res(1)

                call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)

                taub = 2.0_dp * pi / abs(Omth)
                call bounce_fast(v, eta, taub, bounceavg, timestep_transport, istate_dv)
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

    end subroutine compute_transport_integral

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

        ! BEGIN TODO: remove all of this after refactoring and re-use routine in orbit
        ! for y(1:3)
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3), Om_tB_v
        real(dp) :: t0
        complex(dp) :: epsn, Hn  ! relative amplitude of perturbation field epsn=Bn/B0
        ! and Hamiltonian Hn = (H - H0)_n

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        call poloidal_velocity(v, eta, bmod, hctrvr(3), hder(3), y(2), ydot)

        ! evaluate orbit averages of Hamiltonian perturbation
        if (pertfile) then
            call do_magfie_pert_amp(x, epsn)
            epsn = epsmn * epsn / bmod
        else
            epsn = epsmn * exp(imun * m0 * y(1))
        end if

        if (eta > etatp) then
            !t0 = 0.25*2*pi/Omth ! Different starting position in orbit
            t0 = 0.0_dp
            Hn = (2.0_dp - eta * bmod) * epsn * exp(imun * (q * mph * (y(1)) - mth * (t - &
                                                                                      t0) * Omth))
        else
            Hn = (2.0_dp - eta * bmod) * epsn * exp(imun * (q * mph * (y(1)) - (mth + q * mph) &
                                                            * t * Omth))
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
    end subroutine timestep_transport

end module neort_transport
