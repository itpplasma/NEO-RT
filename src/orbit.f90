module neort_orbit
    use util, only: imun, pi, mi, qi, c
    use spline, only: spline_coeff, spline_val_0
    use do_magfie_mod, only: do_magfie, s, iota, R0, eps, psi_pr, &
        bphcov, dbthcovds, dbphcovds, q, dqds
    use do_magfie_pert_mod, only: do_magfie_pert_amp
    use neort_profiles, only: vth, Om_tE, dOm_tEds
    use driftorbit, only: etatp, etadt, etamin, etamax, epsmn, mth, mph, m0, mth, &
        init_done, pertfile, magdrift, nonlin, sigv, epsst_spl, epssp_spl, epst_spl, &
        epsp_spl

    implicit none

    integer, parameter :: nvar = 7
    real(8) :: th0

    logical :: noshear = .false.      ! neglect magnetic shear

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


    real(8) :: k_taub_p, d_taub_p, k_taub_t, d_taub_t ! extrapolation at tp bound
    real(8) :: k_OmtB_p, d_Omtb_p, k_Omtb_t, d_Omtb_t ! extrapolation at tp bound

contains

    subroutine bounce(v, eta, taub, bounceavg, taub_estimate)
        ! calculate all bounce averages
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: taub, bounceavg(nvar)
        real(8), optional :: taub_estimate  ! estimated bounce time (user input)
        real(8) :: findroot_res(nvar + 1)
        real(8) :: bmod
        real(8) :: y0(nvar)

        ! Initialize bounce-averated quantities y0. Their meaning
        ! is defined inside subroutine timestep (thin orbit integration)
        call evaluate_bfield_local(bmod)
        y0 = 1d-15
        y0(1) = th0         ! poloidal angle theta
        y0(2) = vpar(v, eta, bmod)  ! parallel velocity vpar
        y0(3) = 0d0         ! toroidal velocity v_ph for drift frequency Om_ph
        y0(4) = 0d0         ! perturbed Hamiltonian real part
        y0(5) = 0d0         ! perturbed Hamiltonian imaginary part
        y0(6) = 0d0         ! 1/abs(B)
        ! y0(7) = 0d0       ! abs(B)

        ! If bounce time estimate exists (elliptic integrals),
        ! initialize taub with it, owtherwise estimate here.
        if (present(taub_estimate)) then
            taub = taub_estimate
        else
            taub = 2.0*pi/abs(vperp(v, eta, bmod)*iota/R0*sqrt(eps/2d0))
        end if

        ! Look for exactly one orbit turn via root-finding.
        ! Start by looking for 5 points per turn.
        findroot_res = bounce_integral(v, eta, nvar, y0, taub/5d0, timestep)

        taub = findroot_res(1)
        bounceavg = findroot_res(2:)/taub
    end subroutine bounce

    subroutine evaluate_bfield_local(bmod)
        real(8), intent(out) :: bmod
        real(8) :: sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0d0
        x(3) = th0
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
    end subroutine evaluate_bfield_local

    pure function vpar(v, eta, bmod)
        !   parallel velocity
        real(8) :: vpar
        real(8), intent(in) :: v, eta, bmod
        vpar = v*sqrt(1d0 - eta*bmod)
        if (isnan(vpar)) then
            vpar = 0d0
        end if
    end function vpar

    pure function vperp(v, eta, bmod)
        !   perpendicular velocity
        real(8) :: vperp
        real(8), intent(in) :: v, eta, bmod
        vperp = v*sqrt(eta*bmod)
        if (isnan(vperp)) then
            vperp = 0d0
        end if
    end function vperp

    subroutine bounce_fast(v, eta, taub, bounceavg)
        use dvode_f90_m

        real(8), intent(in) :: v, eta, taub
        real(8), intent(out) :: bounceavg(nvar)

        real(8) :: t1, t2, bmod
        real(8) :: y(nvar)
        real(8) :: atol(nvar), rtol
        integer :: neq, itask, istate
        type(vode_opts) :: options

        t1 = 0d0
        t2 = taub

        call evaluate_bfield_local(bmod)
        y = 1d-15
        y(1) = th0
        y(2) = vpar(v, eta, bmod)
        y(3:6) = 0d0

        neq = nvar
        rtol = 1d-9
        atol = 1d-10
        itask = 1
        istate = 1
        options = set_normal_opts(abserr_vector=atol, relerr=rtol)

        call dvode_f90(timestep2, neq, y, t1, t2, itask, istate, options)

        bounceavg = y/taub

    contains

        subroutine timestep2(neq_, t_, y_, ydot_)
            ! Wrapper routine for timestep to work with VODE
            integer, intent(in) :: neq_
            real(8), intent(in) :: t_
            real(8), intent(in) :: y_(neq_)
            real(8), intent(out) :: ydot_(neq_)

            call timestep(v, eta, neq_, t_, y_, ydot_)
        end subroutine timestep2
    end subroutine bounce_fast

    function bounce_time(v, eta, taub_estimate) result(taub)
        use dvode_f90_m, only: dvode_f90, set_normal_opts, vode_opts

        real(8), intent(in) :: v, eta
        real(8), intent(in), optional :: taub_estimate
        real(8) :: taub

        integer, parameter :: neq = 2
        real(8) :: y0(neq), roots(neq+1)
        real(8) :: bmod

        call evaluate_bfield_local(bmod)

        y0(1) = th0         ! poloidal angle theta
        y0(2) = vpar(v, eta, bmod)  ! parallel velocity vpar

        if (present(taub_estimate)) then
            taub = taub_estimate
        else
            taub = 2.0*pi/abs(vperp(v, eta, bmod)*iota/R0*sqrt(eps/2d0))
        end if

        roots = bounce_integral(v, eta, neq, y0, taub, timestep_poloidal_motion)
        taub = roots(1)

    end function bounce_time

    subroutine timestep_poloidal_motion(v, eta, neq, t, y, ydot)
        real(8), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(8), intent(in) :: t, y(neq)
        real(8), intent(out) :: ydot(neq)

        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0d0
        x(3) = y(1)

        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        call timestep_poloidal_internal(v, eta, bmod, hctrvr(3), hder(3), neq, t, y, ydot)
    end subroutine timestep_poloidal_motion

    pure subroutine timestep_poloidal_internal(v, eta, bmod, hthctr, hderth, neq, t, y, ydot)
        real(8), intent(in) :: v, eta, bmod, hthctr, hderth
        integer, intent(in) :: neq
        real(8), intent(in) :: t, y(neq)
        real(8), intent(out) :: ydot(neq)

        ydot(1) = y(2)*hthctr                                   ! theta
        ydot(2) = -v**2*eta/2d0*hthctr*hderth*bmod              ! v_par
    end subroutine timestep_poloidal_internal

    function bounce_integral(v, eta, neq, y0, dt, ts)
        !
        !  Finds the root of an orbit after the first turn
        !
        use dvode_f90_m

        real(8) :: bounce_integral(neq + 1)
        real(8), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(8), intent(in) :: y0(neq), dt
        external ts

        integer :: n

        integer :: k, state, rootstate
        real(8) :: ti, told
        real(8) :: y(neq), yold(neq)

        logical :: passing

        real(8) :: atol(neq), rtol, tout
        integer :: itask, istate
        type(vode_opts) :: options

        rtol = 1d-9
        atol = 1d-10
        itask = 1
        istate = 1
        options = set_normal_opts(abserr_vector=atol, relerr=rtol, nevents=2)

        ! check for passing orbit
        passing = .false.
        if (eta < etatp) then
            passing = .true.
        end if

        n = 500
        rootstate = -1

        y = y0
        yold = y0
        ti = 0d0
        state = 1
        do k = 2, n
            yold = y
            told = ti

            tout = ti + dt
            call dvode_f90(timestep2, neq, y, ti, tout, itask, istate, options, &
                        g_fcn=bounceroots)
            if (istate == 3) then
                if (passing .or. (yold(1) - th0) < 0) then
                    exit
                end if

            end if

            istate = 2
        end do
        if (istate /= 3) then
            write (0, *) "ERROR: bounce_integral did not converge after 500 iterations"
            write (0, *) eta, etamin, etamax, y(1)
        end if

        bounce_integral(1) = ti
        bounce_integral(2:) = y

    contains

        subroutine timestep2(neq_, t_, y_, ydot_)
            ! Wrapper routine for timestep to work with VODE
            integer, intent(in) :: neq_
            real(8), intent(in) :: t_
            real(8), intent(in) :: y_(neq_)
            real(8), intent(out) :: ydot_(neq_)

            call ts(v, eta, neq_, t_, y_, ydot_)
        end subroutine timestep2
    end function bounce_integral

    subroutine bounceroots(NEQ, T, Y, NG, GOUT)
        integer, intent(in) :: NEQ, NG
        real(8), intent(in) :: T, Y(neq)
        real(8), intent(out) :: GOUT(ng)
        associate (dummy => T)
        end associate
        GOUT(1) = Y(1) - th0
        GOUT(2) = 2d0*pi - (Y(1) - th0)
        return
    end subroutine bounceroots

    subroutine timestep(v, eta, neq, t, y, ydot)
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

        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(8) :: Om_tB_v
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: t0
        real(8) :: shearterm
        complex(8) :: epsn, Hn ! relative amplitude of perturbation field epsn=Bn/B0
        ! and Hamiltonian Hn = (H - H0)_n

        x(1) = s
        x(2) = 0d0
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        shearterm = Bphcov*dqds
        if (noshear) then
            shearterm = 0
        end if

        Om_tB_v = mi*c*q/(2d0*qi*psi_pr*bmod)*( &      ! Om_tB/v**2
                  -(2d0 - eta*bmod)*bmod*hder(1) &
                  + 2d0*(1d0 - eta*bmod)*hctrvr(3)* &
                  (dBthcovds + q*dBphcovds + shearterm))

        ydot(1) = y(2)*hctrvr(3)                                    ! theta
        ydot(2) = -v**2*eta/2d0*hctrvr(3)*hder(3)*bmod              ! v_par
        ydot(3) = Om_tB_v                                           ! v_ph

        if (init_done) then
            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
            Omth = abs(Omth)

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
            ydot(4) = real(Hn)
            ydot(5) = aimag(Hn)

            ! evaluate orbit averages for nonlinear attenuation
            if (nonlin) then
                ydot(6) = 1d0/bmod
                ydot(7) = bmod
            else
                ydot(6:7) = 0d0
            end if
        else
            ydot(4:7) = 0d0
        end if
    end subroutine timestep

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
            Omth_v(k + 1) = 2*pi/(v*taub)
            if (k == netaspl_pass - 2) then
                leta0 = log(etatp - eta)
                taub0 = v*taub
                if (magdrift) OmtB0 = Om_tB_v(k + 1)/Omth_v(k + 1)
            end if
            if (k == netaspl_pass - 1) then
                leta1 = log(etatp - eta)
                taub1 = v*taub
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
        Omth = 2d0*pi/taub
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
        dOmthds = (2d0*pi/taub - Omth)/ds
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

end module neort_orbit
