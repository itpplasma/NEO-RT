! Resonant transport regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, since 2015

module driftorbit
    use util
    use do_magfie_mod
    use do_magfie_pert_mod, only: do_magfie_pert_amp, mph
    use spline
    use neort_profiles, only: vth, M_t, ni1, Om_tE, A1, A2, dM_tds, dni1ds, dOm_tEds
    use collis_alp, only: coleff

    implicit none

    integer, parameter :: nvar = 7

    ! Normalization factor for radial electric field
    real(8) :: efac = 1d0

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

    ! Harmonics TODO: make dynamic, multiple harmonics
    ! Default values are overridden by config file in driftorbit_test:read_control
    real(8) :: epsmn = 1d0            ! perturbation amplitude B1/B0
    integer :: m0 = 1                 ! Boozer poloidal perturbation mode
    integer :: mth = 1                ! canonical poloidal mode
    logical :: magdrift = .true.      ! consider magnetic drift
    logical :: nopassing = .false.    ! neglect passing particles
    logical :: noshear = .false.      ! neglect magnetic shear
    logical :: pertfile = .false.     ! read perturbation from file with neo_magfie_pert
    logical :: comptorque = .true.    ! compute torque

    ! Flux surface TODO: make a dynamic, multiple flux surfaces support
    real(8) :: dVds, etadt, etatp
    real(8) :: etamin, etamax

    ! Check if init is done
    logical :: init_done

    ! TODO: better B0 calculation (average magnetic field on flux surface)
    real(8) :: B0
    real(8) :: Bmin, Bmax, th0

    real(8), parameter :: epst_spl = 1d-6, epsp_spl = 1d-6   ! dist to tpb for spline
    real(8), parameter :: epsst_spl = 1d-3, epssp_spl = 1d-3 ! dist to deep for spline
    real(8), parameter :: epst = 1d-8, epsp = 1d-8 ! smallest eta distance to tp bound
    real(8) :: k_taub_p, d_taub_p, k_taub_t, d_taub_t ! extrapolation at tp bound
    real(8) :: k_OmtB_p, d_Omtb_p, k_Omtb_t, d_Omtb_t ! extrapolation at tp bound

    ! Number of levels for coarse root finding
    integer, parameter :: nlev = 100
    real(8) :: sigv = 1d0

    ! Nonlinear calculation switch
    logical :: nonlin = .false.
contains

    subroutine init
        init_done = .false.
        call init_fsa
        call init_misc
        call init_Om_spl       ! frequencies of trapped orbits
        if (.not. nopassing) call init_Om_pass_spl  ! frequencies of passing orbits
        sigv = 1
        call set_to_trapped_region(etamin, etamax)
        init_done = .true.
    end subroutine init

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

    subroutine init_fsa
        ! Calculate the flux surface areas for normalization
        integer, parameter :: nth = 1000
        integer :: k
        real(8) :: thrange(nth), dth
        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        print *, "       s: ", s

        thrange = -pi + (/(k*2*pi/nth, k=1, nth)/)

        dth = thrange(2) - thrange(1)
        x(1) = s
        x(2) = 0d0
        x(3) = 0d0

        dVds = 0d0
        B0 = 0d0
        print *, " eps orig: ", eps
        eps = 0d0

        Bmin = -1d0
        Bmax = 0d0

        do k = 1, nth
            x(3) = thrange(k)
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            dVds = dVds + sqrtg*dth
            B0 = B0 + bmod*dth
            eps = eps - cos(x(3))*bmod*dth

            ! TODO: do fine search
            if ((Bmin < 0) .or. (bmod < Bmin)) then
                Bmin = bmod
                th0 = x(3)
            end if
            if (bmod > Bmax) Bmax = bmod
        end do

        dVds = 2d0*pi*dVds
        B0 = B0/(2d0*pi)
        eps = eps/(B0*pi)

        print *, " eps calc: ", eps
        print *, "      th0: ", th0
        print *, "     dVds: ", dVds
        print *, "       B0: ", B0
        print *, "Bmin,Bmax: ", Bmin, Bmax
        print *, "        x: ", x
        print *, "     bmod: ", bmod
        print *, "    sqrtg: ", sqrtg
        print *, "     hder: ", hder
        print *, "   hcovar: ", hcovar
        print *, "   hctrvr: ", hctrvr
        print *, "    hcurl: ", hcurl

        call disp('init_fsa: iota       = ', iota)

    end subroutine init_fsa

    subroutine init_misc
        ! TODO: fine search for minima and maxima
        etatp = 1d0/Bmax
        etadt = 1d0/Bmin
    end subroutine init_misc

    pure function Jperp(v, eta)
        real(8) :: Jperp
        real(8), intent(in) :: v, eta
        Jperp = 0.5d0*mi*v**2*mi*c/qi*eta
    end function Jperp

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
            taub = 2.0*pi/(vperp(v, eta, bmod)*iota/R0*sqrt(eps/2d0))
        end if

        ! Look for exactly one orbit turn via root-finding.
        ! Start by looking for 5 points per turn.
        findroot_res = bounce_integral(v, eta, nvar, y0, taub/5d0, timestep)

        taub = findroot_res(1)
        bounceavg = findroot_res(2:)/taub
    end subroutine bounce

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
            taub = 2.0*pi/(vperp(v, eta, bmod)*iota/R0*sqrt(eps/2d0))
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

    subroutine evaluate_bfield_local(bmod)
        real(8), intent(out) :: bmod
        real(8) :: sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0d0
        x(3) = th0
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
    end subroutine evaluate_bfield_local

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

    subroutine driftorbit_coarse(v, eta_min, eta_max, roots, nroots)
        real(8), intent(in) :: v, eta_min, eta_max
        real(8), intent(out) :: roots(:, :)
        integer, intent(out) :: nroots
        real(8) :: deta
        real(8) :: Omph, dOmphdv, dOmphdeta
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: res, dresdv, dresdeta
        real(8) :: resold, dresdvold, dresdetaold
        real(8) :: eta
        integer :: k, ninterv

        ninterv = size(roots, 1)

        deta = (eta_max - eta_min)*1d0/ninterv
        nroots = 0

        do k = 0, ninterv
            eta = eta_min + k*deta
            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
            call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
            res = mth*Omth + mph*Omph
            dresdv = mth*dOmthdv + mph*dOmphdv
            dresdeta = mth*dOmthdeta + mph*dOmphdeta
            if (k > 0) then
                if (sign(1d0, res) /= sign(1d0, resold)) then
                    nroots = nroots + 1
                    roots(nroots, 1) = eta - deta
                    roots(nroots, 2) = eta
                end if
            end if
            resold = res
            dresdvold = dresdv
            dresdetaold = dresdeta
        end do
    end subroutine driftorbit_coarse

    function driftorbit_nroot(v, eta_min, eta_max)
        integer :: driftorbit_nroot
        real(8), intent(in) :: v
        real(8), intent(in) :: eta_min, eta_max
        real(8) :: etamin2, etamax2
        real(8) :: roots(nlev, 3)

        call driftorbit_coarse(v, eta_min, eta_max, roots, driftorbit_nroot)
    end function driftorbit_nroot

    function driftorbit_root(v, tol, eta_min, eta_max)
        real(8) :: driftorbit_root(2)
        real(8), intent(in) :: v, tol, eta_min, eta_max
        real(8) :: res, res_old, eta0, eta_old
        real(8) :: Omph, dOmphdv, dOmphdeta
        real(8) :: Omth, dOmthdv, dOmthdeta
        integer :: maxit, k, state
        real(8) :: etamin2, etamax2
        logical :: slope_pos
        real(8) :: resmin, resmax
        real(8) :: eta

        maxit = 100
        state = -2
        eta0 = eta
        eta_old = 0d0
        res = 0d0

        etamin2 = eta_min
        etamax2 = eta_max

        eta = etamin2
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        res = mph*Omph + mth*Omth
        resmin = res

        eta = etamax2
        call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
        call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
        resmax = mph*Omph + mth*Omth
        if (resmax - resmin > 0) then
            slope_pos = .true.
        else
            slope_pos = .false.
        end if

        if (driftorbit_nroot(v, etamin2, etamax2) == 0) then
            print *, "ERROR: driftorbit_root couldn't bracket 0 for v/vth = ", v/vth
            print *, "ERROR: etamin = ", etamin2, " etamax = ", etamax2

            return
        end if

        do k = 1, maxit
            res_old = res
            call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
            call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
            res = mph*Omph + mth*Omth

            driftorbit_root(1) = eta

            if (abs(res) < tol) then
                state = 1
                driftorbit_root(2) = mph*dOmphdeta + mth*dOmthdeta
                exit
            elseif ((slope_pos .and. res > 0) .or. &
                    ((.not. slope_pos) .and. res < 0)) then
                etamax2 = eta
                eta_old = eta
                eta = (eta + etamin2)/2d0
            else
                etamin2 = eta
                eta_old = eta
                eta = (eta + etamax2)/2d0
            end if
        end do
        if (state < 0) then
            driftorbit_root(2) = mph*dOmphdeta + mth*dOmthdeta
            print *, "ERROR: driftorbit_root did not converge in 100 iterations"
            print *, "v/vth  = ", v/vth, "mth    = ", mth, "sigv= ", sigv
            print *, "etamin = ", eta_min, "etamax = ", eta_max, "eta = ", eta
            print *, "resmin = ", resmin, "resmax = ", resmax, "res = ", res
            print *, "resold = ", res_old, "res    = ", res
        end if
        eta = eta0
    end function driftorbit_root

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

    function nonlinear_attenuation(ux, eta, bounceavg, Omth, dOmthdv, dOmthdeta, Hmn2)
        real(8), intent(in) :: ux, eta, bounceavg(nvar), Omth, dOmthdv, dOmthdeta, Hmn2
        real(8) :: nonlinear_attenuation

        real(8) :: dpp, dhh, fpeff, dres, dnorm, Omph, dOmphdv, dOmphdeta, dOmdv, &
                   dOmdeta, Ompr, dOmphds, dOmthds, dOmdpph, v

        if (.not. nonlin) then
            nonlinear_attenuation = 1d0
            return
        end if

        v = ux*vth

        if (nonlin) then
            call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
            call d_Om_ds(v, eta, dOmthds, dOmphds)
            dOmdv = mth*dOmthdv + mph*dOmphdv
            dOmdeta = mth*dOmthdeta + mph*dOmphdeta
            dOmdpph = -(qi/c*iota*psi_pr)**(-1)*(mth*dOmthds + mph*dOmphds)
            Ompr = omega_prime(ux, eta, bounceavg, Omth, dOmdv, dOmdeta, dOmdpph)
            call coleff(ux, dpp, dhh, fpeff)
            dhh = vth*dhh
            dpp = vth**3*dpp
            dres = dpp*(dOmdv/Ompr)**2 + dhh*eta*(bounceavg(6) - eta)*(dOmdeta/Ompr)**2
            dnorm = dres*sqrt(abs(Ompr))/sqrt(abs(Hmn2))**(3d0/2d0)
            call attenuation_factor(dnorm, nonlinear_attenuation)
        end if
    end function nonlinear_attenuation

    pure function omega_prime(ux, eta, bounceavg, Omth, dOmdv, dOmdeta, dOmdpph)
        real(8), intent(in) :: ux, eta, bounceavg(nvar), Omth, dOmdv, dOmdeta, dOmdpph
        real(8) :: omega_prime

        real(8) :: ma, mb, mc, md, me, mf, dvdJ, detadJ

        ma = mi*(ux*vth)*mi*c/qi*eta
        mb = mi*(ux*vth)**2/2*mi*c/qi
        mc = mi/(2d0*Omth)*(ux*vth)*(1d0 - eta*bounceavg(7))
        md = mi*(ux*vth)**2/2d0*Omth
        me = -mth/mph
        mf = 1d0/mph

        dvdJ = mb*me/(ma*md*mf - mb*mc*mf)
        detadJ = ma*me/(mb*mc*mf - ma*md*mf)

        omega_prime = dOmdv*dvdJ + dOmdeta*detadJ + mph*dOmdpph
    end function omega_prime

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

    pure subroutine set_to_trapped_region(eta_min, eta_max)
        real(8), intent(out) :: eta_min, eta_max
        eta_min = (1 + epst)*etatp
        eta_max = (1 - epst)*etadt
    end subroutine set_to_trapped_region

    pure subroutine set_to_passing_region(eta_min, eta_max)
        real(8), intent(out) :: eta_min, eta_max
        eta_min = epsp*etatp
        eta_max = (1 - epsp)*etatp
    end subroutine set_to_passing_region

end module driftorbit
