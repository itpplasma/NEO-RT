module neort_transport
    use util, only: pi, c, qi
    use do_magfie_mod, only: s, a, R0, iota, q, psi_pr, eps
    use neort_magfie, only: dVds, B0
    use neort_profiles, only: ni1, Om_tE
    use neort_nonlin, only: nonlinear_attenuation
    use neort_orbit, only: Om_th, Om_ph, bounce_fast, nvar
    use neort_resonance, only: driftorbit_coarse, driftorbit_root
    use driftorbit, only: vth, mth, mph, mi, B0, Bmin, Bmax, comptorque, &
        etamin, etamax, A1, A2, nlev

    implicit none

contains

    pure function D11int(ux, taub, Hmn2)
        real(8) :: D11int
        real(8), intent(in) :: ux, taub, Hmn2

        D11int = pi**(3d0/2d0)*mph**2*c**2*q*vth &
                 /(qi**2*dVds*psi_pr)*ux**3*exp(-ux**2) &
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
        real(8) :: Omth, dOmthdv, dOmthdeta
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
                call bounce_fast(v, eta, taub, bounceavg)
                Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2
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

end module neort_transport
