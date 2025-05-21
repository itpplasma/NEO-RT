module neort_freq
    use util, only: pi
    use spline, only: spline_coeff, spline_val_0
    use neort_orbit, only: nvar, bounce
    use neort_profiles, only: vth, Om_tE, dOm_tEds
    use driftorbit, only: etamin, etamax, etatp, etadt, epsst_spl, epst_spl, magdrift, epssp_spl, epsp_spl, sigv, sigv_theta
    use do_magfie_mod, only: iota, s, Bthcov
    implicit none

    ! For splining in the trapped eta region
    integer, parameter :: netaspl = 100
    real(8) :: OmtB_spl_coeff(netaspl - 1, 5)
    real(8) :: Omth_spl_coeff(netaspl - 1, 5)
    real(8) :: vres_spl_coeff(netaspl - 1, 5)

    ! For splining in the passing eta region
    integer, parameter :: netaspl_pass = 100
    real(8) :: OmtB_pass_spl_coeff(netaspl_pass - 1, 5)
    real(8) :: Omth_pass_spl_coeff(netaspl_pass - 1, 5)
    real(8) :: vres_pass_spl_coeff(netaspl - 1, 5)

    real(8) :: k_taub_p=0d0, d_taub_p=0d0, k_taub_t=0d0, d_taub_t=0d0 ! extrapolation at tp bound
    real(8) :: k_OmtB_p=0d0, d_Omtb_p=0d0, k_Omtb_t=0d0, d_Omtb_t=0d0 ! extrapolation at tp bound

contains
    subroutine init_Om_spl
        ! Initialise splines for canonical frequencies of trapped orbits

        real(8) :: etarange(netaspl), Om_tB_v(netaspl), Omth_v(netaspl)
        integer :: k
        real(8) :: aa, b
        real(8) :: taub0, taub1, leta0, leta1, OmtB0, OmtB1
        real(8) :: v, eta, taub, bounceavg(nvar)

        print *, 'init_Om_spl'

        taub0 = 0d0
        taub1 = 0d0
        leta0 = 0d0
        leta1 = 0d0
        OmtB0 = 0d0
        OmtB1 = 0d0

        v = vth
        etamin = etatp
        etamax = etatp + (etadt - etatp)*(1d0 - epsst_spl)

        ! logspace for eta
        b = log(epst_spl)
        aa = 1d0/(netaspl - 1d0)*(log(etamax/etamin - 1d0) - b)

        do k = netaspl - 1, 0, -1
            eta = etamin*(1d0 + exp(aa*k + b))
            etarange(k + 1) = eta
            if (k == netaspl - 1) then
                call bounce(v, eta, taub, bounceavg)
            else
                call bounce(v, eta, taub, bounceavg, taub)
            end if
            if (magdrift) Om_tB_v(k + 1) = bounceavg(3)
            Omth_v(k + 1) = 2*pi/(v*taub)
            if (k == 0) then
                leta0 = log(eta - etatp)
                taub0 = v*taub
                if (magdrift) OmtB0 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
            if (k == 1) then
                leta1 = log(eta - etatp)
                taub1 = v*taub
                if (magdrift) OmtB1 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
        end do

        k_taub_t = (taub1 - taub0)/(leta1 - leta0)
        d_taub_t = taub0 - k_taub_t*leta0
        Omth_spl_coeff = spline_coeff(etarange, Omth_v)

        if (magdrift) then
            k_OmtB_t = (OmtB1 - OmtB0)/(leta1 - leta0)
            d_OmtB_t = OmtB0 - k_OmtB_t*leta0
            OmtB_spl_coeff = spline_coeff(etarange, Om_tB_v)
        end if

    end subroutine init_Om_spl

    subroutine init_Om_pass_spl
        ! Initialise splines for canonical frequencies of passing orbits

        real(8) :: etarange(netaspl_pass), Om_tB_v(netaspl_pass), Omth_v(netaspl_pass)
        real(8) :: aa, b
        integer :: k
        real(8) :: leta0, leta1, taub0, taub1, OmtB0, OmtB1
        real(8) :: v, eta, taub, bounceavg(nvar)

        print *, 'init_Om_pass_spl'

        taub0 = 0d0
        taub1 = 0d0
        leta0 = 0d0
        leta1 = 0d0
        OmtB0 = 0d0
        OmtB1 = 0d0

        v = vth
        etamin = etatp*epssp_spl
        etamax = etatp

        b = log((etamax - etamin)/etamax)
        aa = 1d0/(netaspl_pass - 1d0)*(log(epsp_spl) - b)

        do k = netaspl_pass - 1, 0, -1
            eta = etamax*(1d0 - exp(aa*k + b))
            etarange(k + 1) = eta
            if (k == netaspl_pass - 1) then
                call bounce(v, eta, taub, bounceavg)
            else
                call bounce(v, eta, taub, bounceavg, taub)
            end if
            if (magdrift) Om_tB_v(k + 1) = bounceavg(3)
            Omth_v(k + 1) = 2*pi/(v*taub)*sign(1d0, Bthcov)
            if (k == netaspl_pass - 2) then
                leta0 = log(etatp - eta)
                taub0 = v*taub*sign(1d0, Bthcov)
                if (magdrift) OmtB0 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
            if (k == netaspl_pass - 1) then
                leta1 = log(etatp - eta)
                taub1 = v*taub*sign(1d0, Bthcov)
                if (magdrift) OmtB1 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
        end do

        k_taub_p = (taub1 - taub0)/(leta1 - leta0)
        d_taub_p = taub0 - k_taub_p*leta0
        Omth_pass_spl_coeff = spline_coeff(etarange, Omth_v)

        if (magdrift) then
            k_OmtB_p = (OmtB1 - OmtB0)/(leta1 - leta0)
            d_OmtB_p = OmtB0 - k_OmtB_p*leta0
            OmtB_pass_spl_coeff = spline_coeff(etarange, Om_tB_v)
        end if
    end subroutine init_Om_pass_spl

    subroutine Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
        ! returns bounce averaged toroidal magnetic drift frequency
        ! and derivatives w.r.t. v and eta
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: OmtB, dOmtBdv, dOmtBdeta
        real(8) :: splineval(3)
        real(8) :: Omth, dOmthdv, dOmthdeta
        if (eta > etatp) then
            if (eta > etatp*(1 + epst_spl)) then
                splineval = spline_val_0(OmtB_spl_coeff, eta)
            else ! extrapolation
                call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
                splineval(1) = sigv*(k_OmtB_t*log(eta - etatp) + d_OmtB_t)*Omth/v
                splineval(2) = sigv*(Omth/v*k_OmtB_t/(eta - etatp) + &
                                dOmthdeta/v*(k_OmtB_t*log(eta - etatp) + d_OmtB_t))
            end if
        else
            if (eta < etatp*(1 - epsp_spl)) then
                splineval = spline_val_0(OmtB_pass_spl_coeff, eta)
            else ! extrapolation
                call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
                splineval(1) = sigv*(k_OmtB_p*log(etatp - eta) + d_OmtB_p)*Omth/v
                splineval(2) = sigv*(Omth/v*k_OmtB_p/(eta - etatp) + &
                                dOmthdeta/v*(k_OmtB_p*log(etatp - eta) + d_OmtB_p))
            end if
        end if
        OmtB = splineval(1)*v**2
        dOmtBdv = 2d0*splineval(1)*v
        dOmtBdeta = splineval(2)*v**2
    end subroutine Om_tB

    subroutine Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        ! returns canonical toroidal frequency
        ! and derivatives w.r.t. v and eta
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: Omph, dOmphdv, dOmphdeta
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: OmtB, dOmtBdv, dOmtBdeta

        if (eta > etatp) then
            Omph = Om_tE
            dOmphdv = 0d0
            dOmphdeta = 0d0
            if (magdrift) then
                call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
                Omph = Omph + OmtB
                dOmphdv = dOmphdv + dOmtBdv
                dOmphdeta = dOmphdeta + dOmtBdeta
            end if
        else
            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
            Omph = Om_tE + Omth/iota
            dOmphdv = dOmthdv/iota
            dOmphdeta = dOmthdeta/iota
            if (magdrift) then
                call Om_tB(v, eta, OmtB, dOmtBdv, dOmtBdeta)
                Omph = Omph + OmtB
                dOmphdv = dOmphdv + dOmtBdv
                dOmphdeta = dOmphdeta + dOmtBdeta
            end if
        end if
    end subroutine Om_ph

    subroutine Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        ! returns canonical poloidal frequency
        ! and derivatives w.r.t. v and eta
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: Omth, dOmthdv, dOmthdeta
        real(8) :: splineval(3)

        if (eta > etatp) then
            if (eta > etatp*(1 + epst_spl)) then
                splineval = spline_val_0(Omth_spl_coeff, eta)
            else ! extrapolation
                splineval(1) = 2d0*pi/(k_taub_t*log(eta - etatp) + d_taub_t)
                splineval(2) = -splineval(1)**2/(2d0*pi)*k_taub_t/(eta - etatp)
            end if
        else
            if (eta < etatp*(1 - epsp_spl)) then
                splineval = spline_val_0(Omth_pass_spl_coeff, eta)
            else ! extrapolation
                splineval(1) = 2d0*pi/(k_taub_p*log(etatp - eta) + d_taub_p)
                splineval(2) = -splineval(1)**2/(2d0*pi)*k_taub_p/(eta - etatp)
            end if
        end if
        Omth = sigv*splineval(1)*v
        dOmthdv = sigv*splineval(1)
        dOmthdeta = sigv*splineval(2)*v
    end subroutine Om_th

    subroutine d_Om_ds(v, eta, dOmthds, dOmphds)
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: dOmthds, dOmphds
        real(8) :: s0, ds, bounceavg(nvar)
        real(8) :: taub, Omth, Omph_noE
        ! store current flux surface values
        s0 = s

        ds = 2d-8
        s = s0 - ds/2d0
        call bounce(v, eta, taub, bounceavg)
        Omth = sigv_theta*2d0*pi/taub
        if (magdrift) then
            if (eta > etatp) then
                Omph_noE = bounceavg(3)*v**2
            else
                Omph_noE = bounceavg(3)*v**2 + Omth/iota
            end if
        else
            if (eta > etatp) then
                Omph_noE = 0d0
            else
                Omph_noE = Omth/iota
            end if
        end if
        s = s0 + ds/2d0
        call bounce(v, eta, taub, bounceavg, taub)
        dOmthds = sigv_theta*(2d0*pi/taub - sigv_theta*Omth)/ds
        if (magdrift) then
            if (eta > etatp) then
                dOmphds = dOm_tEds + (bounceavg(3)*v**2 - Omph_noE)/ds
            else
                dOmphds = dOm_tEds + (bounceavg(3)*v**2 + (2d0*pi/taub)/iota - Omph_noE)/ds
            end if
        else
            if (eta > etatp) then
                dOmphds = dOm_tEds
            else
                dOmphds = dOm_tEds + ((2d0*pi/taub)/iota - Omph_noE)/ds
            end if
        end if

        ! re-set current flux surface values
        s = s0
    end subroutine d_Om_ds
end module neort_freq
