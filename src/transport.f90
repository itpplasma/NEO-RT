module neort_transport
    use util, only: imun, pi, c, qi
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
        etamin, etamax, A1, A2, nlev, pertfile, init_done, nonlin, m0, etatp, etadt, &
        sign_vpar_htheta, sign_vpar

    implicit none

    real(8) :: Omth, dOmthdv, dOmthdeta

contains

    pure function D11int(ux, taub, Hmn2)
        real(8) :: D11int
        real(8), intent(in) :: ux, taub, Hmn2

        D11int = pi**(3d0/2d0)*mph**2*c**2*q*vth &
                 /(qi**2*dVds*sign_theta*psi_pr)*ux**3*exp(-ux**2) &
                 *taub*Hmn2
    end function D11int

    pure function D12int(ux, taub, Hmn2)
        real(8) :: D12int
        real(8), intent(in) :: ux, taub, Hmn2

        D12int = D11int(ux, taub, Hmn2)*ux**2
    end function D12int

    pure function Tphi_int(ux, taub, Hmn2)
        real(8) :: Tphi_int
        real(8), intent(in) :: ux, taub, Hmn2

        Tphi_int = -pi**(3d0/2d0)*mph**2*ni1*c*vth/qi*ux**3*exp(-ux**2)*taub &
                   *Hmn2*(A1 + A2*ux**2)
    end function Tphi_int

    subroutine compute_transport_integral(vmin, vmax, vsteps, D, T)
        ! compute transport integral via midpoint rule
        real(8), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps
        real(8), intent(out) :: D(2), T ! Transport coefficients D and torque density T
        real(8) :: Dp, dsdreff ! Plateau diffusion coefficient and ds/dreff=<|grad s|>
        real(8) :: ux, du, dD11, dD12, dT, v, eta
        real(8) :: eta_res(2)
        real(8) :: taub, bounceavg(nvar)
        real(8) :: Hmn2, attenuation_factor
        real(8) :: roots(nlev, 3)
        integer :: nroots, kr, ku

        D = 0d0
        T = 0d0
        du = (vmax - vmin)/(vsteps*vth)
        ux = vmin/vth + du/2d0

        do ku = 1, vsteps
            v = ux*vth
            call driftorbit_coarse(v, etamin, etamax, roots, nroots)
            if (nroots == 0) continue
            do kr = 1, nroots
                eta_res = driftorbit_root(v, 1d-8*abs(Om_tE), roots(kr, 1), roots(kr, 2))
                eta = eta_res(1)

                call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)

                taub = 2d0*pi/abs(Omth)
                call bounce_fast(v, eta, taub, bounceavg, timestep_transport)
                Hmn2 = (bounceavg(3)**2 + bounceavg(4)**2)*(mi*(ux*vth)**2/2d0)**2
                attenuation_factor = nonlinear_attenuation(ux, eta, bounceavg, Omth, &
                                                           dOmthdv, dOmthdeta, Hmn2)

                dD11 = du*D11int(ux, taub, Hmn2)/abs(eta_res(2))
                dD12 = du*D12int(ux, taub, Hmn2)/abs(eta_res(2))
                D(1) = D(1) + dD11*attenuation_factor
                D(2) = D(2) + dD12*attenuation_factor

                if (comptorque) then
                    dT = du*Tphi_int(ux, taub, Hmn2)/abs(eta_res(2))
                    T = T + dT*attenuation_factor
                end if
            end do
            ux = ux + du
        end do

        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
        dsdreff = 2d0/a*sqrt(s)  ! TODO: Use exact value instead of this approximation
        D = dsdreff**(-2)*D/Dp

    end subroutine compute_transport_integral

    subroutine timestep_transport(v, eta, neq, t, y, ydot)
        !
        !  Timestep function for orbit integration.
        !  Includes poloidal angle theta and parallel velocity.
        !  More integrands may be added starting from y(3)
        !

        real(8), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(8), intent(in) :: t
        real(8), intent(in) :: y(neq)
        real(8), intent(out) :: ydot(neq)

        ! BEGIN TODO: remove all of this after refactoring and re-use routine in orbit
        ! for y(1:3)
        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3), Om_tB_v
        real(8) :: t0
        complex(8) :: epsn, Hn ! relative amplitude of perturbation field epsn=Bn/B0
        ! and Hamiltonian Hn = (H - H0)_n

        x(1) = s
        x(2) = 0d0
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        call poloidal_velocity(v, eta, bmod, hctrvr(3), hder(3), y(2), ydot)

        ! evaluate orbit averages of Hamiltonian perturbation
        if (pertfile) then
            call do_magfie_pert_amp(x, epsn)
            epsn = epsmn*epsn/bmod
        else
            epsn = epsmn*exp(imun*m0*y(1))
        end if

        if (eta > etatp) then
            !t0 = 0.25*2*pi/Omth ! Different starting position in orbit
            t0 = 0d0
            Hn = (2d0 - eta*bmod)*epsn*exp(imun*(q*mph*(y(1)) - mth*(t - t0)*Omth))
        else
            Hn = (2d0 - eta*bmod)*epsn*exp(imun*(q*mph*(y(1)) - (mth + q*mph)*t*Omth))
        end if
        ydot(3) = real(Hn)
        ydot(4) = aimag(Hn)

        ! evaluate orbit averages for nonlinear attenuation
        if (nonlin) then
            ydot(5) = 1d0/bmod
            ydot(6) = bmod
        else
            ydot(5:6) = 0d0
        end if
    end subroutine timestep_transport

end module neort_transport
