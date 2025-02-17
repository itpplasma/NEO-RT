! Resonant transport regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, since 2015

module driftorbit
    use do_magfie_mod
    use do_magfie_pert_mod, only: do_magfie_pert_amp, mph
    use collis_alp
    use spline

    implicit none
    save

    integer, parameter :: nvar = 7

    ! Orbit parameters
    real(8) :: v, eta
    real(8) :: taub, bounceavg(nvar)

    ! Plasma parameters
    real(8) :: Om_tE, dOm_tEds, vth, M_t, dM_tds, efac

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
    real(8) :: epsmn = 1d-3           ! perturbation amplitude B1/B0
    integer :: m0 = 1                 ! Boozer poloidal perturbation mode
    integer :: mth = 1                ! canonical poloidal mode
    logical :: supban = .False.     ! calculate superbanana plateau
    logical :: magdrift = .True.      ! consider magnetic drift
    logical :: nopassing = .False.    ! neglect passing particles
    logical :: noshear = .False.      ! neglect magnetic shear
    logical :: pertfile = .False.     ! read perturbation from file with neo_magfie_pert

    logical :: odeint = .False.       ! use ode integrator for orbit integration

    ! Flux surface TODO: make a dynamic, multiple flux surfaces support
    real(8) :: dVds, etadt, etatp
    real(8) :: etamin, etamax
    real(8) :: vmin2, vmax2 ! permanent vmin vmax for integration bracketing

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
    logical :: nonlin

    ! Orbit modes: 0 ... zeroth order, 1 ... first order, 2 ... full orbit
    integer :: orbit_mode_avg, orbit_mode_transp

    ! Number of integration steps in v, set 0 for adaptive integration by quadpack
    integer :: vsteps

    ! Box boundaries of s levels where to check intersections of orbits
    real(8), allocatable :: sbox(:)
    real(8) :: s1, snext, sprev ! first order correction to s

    ! Integrals for radial boxes
    real(8), allocatable :: fluxint_box(:, :)
    real(8), allocatable :: torque_int_box(:)
    real(8), allocatable :: taubins(:)

    ! Plasma parameters
    real(8) ni1, ni2, Ti1, Ti2, Te, dni1ds, dni2ds, dTi1ds, dTi2ds, dTeds

    ! Thermodynamic forces in radial variable s
    real(8) A1, A2

    ! Output integral quantities and resonance line
    logical :: intoutput

contains

    subroutine init
        init_done = .false.
        call init_fsa
        call init_misc
        call init_Om_spl       ! frequencies of trapped orbits
        if (.not. nopassing) call init_Om_pass_spl  ! frequencies of passing orbits
        init_done = .true.
    end subroutine init

    subroutine init_Om_spl
        ! Initialise splines for canonical frequencies of trapped orbits

        real(8) :: etarange(netaspl), Om_tB_v(netaspl), Omth_v(netaspl)
        integer :: k
        real(8) :: aa, b
        real(8) :: taub0, taub1, leta0, leta1, OmtB0, OmtB1

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
            !eta = etamin * (1d0 + (k/(netaspl-1d0))**4*(etamax/etamin-1d0))
            eta = etamin*(1d0 + exp(aa*k + b))
            etarange(k + 1) = eta
            if (k == netaspl - 1) then
                call bounce
            else
                call bounce(taub)
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
                call bounce
            else
                call bounce(taub)
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

    subroutine init_misc
        ! TODO: fine search for minima and maxima
        etatp = 1d0/Bmax
        etadt = 1d0/Bmin
    end subroutine init_misc

    function Jperp()
        real(8) :: Jperp
        Jperp = mi**2*c/(2d0*qi)*eta*v**2
    end function Jperp

    function vpar(bmod)
        !   parallel velocity
        real(8) :: bmod, vpar
        vpar = v*sqrt(1d0 - eta*bmod)
        if (isnan(vpar)) then
            vpar = 0d0
        end if
    end function vpar

    function vperp(bmod)
        !   perpendicular velocity
        real(8) :: bmod, vperp
        vperp = v*sqrt(eta*bmod)
        if (isnan(vperp)) then
            vperp = 0d0
        end if
    end function vperp

    subroutine jac
    end subroutine jac

    subroutine timestep2(neq, t, y, ydot)
        ! Wrapper routine for timestep to work with VODE
        integer, intent(in) :: neq
        real(8), intent(in) :: t
        real(8), intent(in) :: y(neq)
        real(8), intent(out) :: ydot(neq)

        call timestep(t, y, ydot)
    end subroutine timestep2

    subroutine timestep(t, y, ydot)
        !
        !  Timestep function for orbit integration.
        !  Includes poloidal angle theta and parallel velocity.
        !  More integrands may be added starting from y(3)
        !

        real(8), intent(in) :: t
        real(8), intent(in) :: y(*)
        real(8), intent(out) :: ydot(*)

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
            call Om_th(Omth, dOmthdv, dOmthdeta)
            Omth = abs(Omth)
            if (orbit_mode_avg == 1) then
                ! use first order radial orbit width
                x(1) = s + sigv*c*mi*y(2)*hcovar(2)*q/(qi*psi_pr)
                call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            end if

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
            if (orbit_mode_avg == 1) then
                ! reset original values
                x(1) = s
                call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            end if
        else
            ydot(4:7) = 0d0
        end if
    end subroutine timestep

    function findroot2(y0, dt)
        !
        !  Finds the root of an orbit after the first turn
        !
        use dvode_f90_m

        real(8) :: y0(nvar), dt, findroot2(nvar + 1)
        integer :: n

        integer :: k, state, rootstate
        real(8) :: ti, told
        real(8) :: y(nvar), yold(nvar)

        logical :: passing

        real(8) :: atol(nvar), rtol, tout
        integer :: neq, itask, istate
        type(vode_opts) :: options

        neq = nvar
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
            !print *, ti
            if (istate == 3) then
                if (passing .or. (yold(1) - th0) < 0) then
                    exit
                end if

            end if

            istate = 2
        end do
        if (istate /= 3) then
            write (0, *) "ERROR: findroot2 did not converge after 500 iterations"
            write (0, *) eta, etamin, etamax, y(1)
        end if

        findroot2(1) = ti
        findroot2(2:) = y
    end function findroot2

    subroutine bounceroots(NEQ, T, Y, NG, GOUT)
        integer, intent(in) :: NEQ, NG
        real(8), intent(in) :: T, Y(neq)
        real(8), intent(out) :: GOUT(ng)
        associate(dummy => T)
        end associate
        GOUT(1) = Y(1) - th0
        GOUT(2) = 2d0*pi - (Y(1) - th0)
        return
    end subroutine bounceroots

    subroutine bounce(taub_estimate)
        ! calculate all bounce averages
        real(8) :: findroot_res(nvar + 1)
        real(8), optional :: taub_estimate  ! estimated bounce time (user input)
        real(8) :: taub_est                 ! estimated bounce time (computed)
        real(8) :: bmod
        real(8) :: y0(nvar)

        ! Initialize bounce-averated quantities y0. Their meaning
        ! is defined inside subroutine timestep (thin orbit integration)
        call evaluate_bfield_local(bmod)
        y0 = 1d-15
        y0(1) = th0         ! poloidal angle theta
        y0(2) = vpar(bmod)  ! parallel velocity vpar
        y0(3) = 0d0         ! toroidal velocity v_ph for drift frequency Om_ph
        y0(4) = 0d0         ! perturbed Hamiltonian real part
        y0(5) = 0d0         ! perturbed Hamiltonian imaginary part
        y0(6) = 0d0         ! 1/abs(B)
        ! y0(7) = 0d0       ! abs(B)

        ! If bounce time estimate exists (elliptic integrals),
        ! initialize taub with it, owtherwise estimate here.
        if (present(taub_estimate)) then
            taub_est = taub_estimate
        else
            taub_est = 2.0*pi/(vperp(bmod)*iota/R0*sqrt(eps/2d0))
        end if

        ! Look for exactly one orbit turn via root-finding.
        ! Start by looking for 5 points per turn.
        findroot_res = findroot2(y0, taub_est/5d0)
        !print *, taub
        taub = findroot_res(1)
        bounceavg = findroot_res(2:)/taub

        call evaluate_bfield_local(bmod)

    end subroutine bounce

    subroutine bounce_fast
        use dvode_f90_m

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
        y(2) = vpar(bmod)
        y(3:6) = 0d0

        neq = nvar
        rtol = 1d-9
        atol = 1d-10
        itask = 1
        istate = 1
        options = set_normal_opts(abserr_vector=atol, relerr=rtol)

        call dvode_f90(timestep2, neq, y, t1, t2, itask, istate, options)

        bounceavg = y/taub
    end subroutine bounce_fast

    subroutine evaluate_bfield_local(bmod)
        real(8), intent(out) :: bmod
        real(8) :: sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0d0
        x(3) = th0
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
    end subroutine evaluate_bfield_local

    subroutine Om_tB(OmtB, dOmtBdv, dOmtBdeta)
        ! returns bounce averaged toroidal magnetic drift frequency
        ! and derivatives w.r.t. v and eta
        real(8), intent(out) :: OmtB, dOmtBdv, dOmtBdeta
        real(8) :: splineval(3)
        real(8) :: Omth, dOmthdv, dOmthdeta
        if (eta > etatp) then
            if (eta > etatp*(1 + epst_spl)) then
                splineval = spline_val_0(OmtB_spl_coeff, eta)
            else ! extrapolation
                call Om_th(Omth, dOmthdv, dOmthdeta)
                splineval(1) = sigv*(k_OmtB_t*log(eta - etatp) + d_OmtB_t)*Omth/v
                splineval(2) = sigv*(Omth/v*k_OmtB_t/(eta - etatp) + &
                                     dOmthdeta/v*(k_OmtB_t*log(eta - etatp) + d_OmtB_t))
            end if
        else
            if (eta < etatp*(1 - epsp_spl)) then
                splineval = spline_val_0(OmtB_pass_spl_coeff, eta)
            else ! extrapolation
                call Om_th(Omth, dOmthdv, dOmthdeta)
                splineval(1) = sigv*(k_OmtB_p*log(etatp - eta) + d_OmtB_p)*Omth/v
                splineval(2) = sigv*(Omth/v*k_OmtB_p/(eta - etatp) + &
                                     dOmthdeta/v*(k_OmtB_p*log(etatp - eta) + d_OmtB_p))
            end if
        end if
        OmtB = splineval(1)*v**2
        dOmtBdv = 2d0*splineval(1)*v
        dOmtBdeta = splineval(2)*v**2
    end subroutine Om_tB

    subroutine Om_ph(Omph, dOmphdv, dOmphdeta)
        ! returns canonical toroidal frequency
        ! and derivatives w.r.t. v and eta
        real(8), intent(out) :: Omph, dOmphdv, dOmphdeta
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: OmtB, dOmtBdv, dOmtBdeta

        if (eta > etatp) then
            Omph = Om_tE
            dOmphdv = 0d0
            dOmphdeta = 0d0
            if (magdrift) then
                call Om_tB(OmtB, dOmtBdv, dOmtBdeta)
                Omph = Omph + OmtB
                dOmphdv = dOmphdv + dOmtBdv
                dOmphdeta = dOmphdeta + dOmtBdeta
            end if
        else
            call Om_th(Omth, dOmthdv, dOmthdeta)
            Omph = Om_tE + Omth/iota
            dOmphdv = dOmthdv/iota
            dOmphdeta = dOmthdeta/iota
            if (magdrift) then
                call Om_tB(OmtB, dOmtBdv, dOmtBdeta)
                Omph = Omph + OmtB
                dOmphdv = dOmphdv + dOmtBdv
                dOmphdeta = dOmphdeta + dOmtBdeta
            end if
        end if
    end subroutine Om_ph

    subroutine Om_th(Omth, dOmthdv, dOmthdeta)
        ! returns canonical poloidal frequency
        ! and derivatives w.r.t. v and eta
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

    subroutine d_Om_ds(dOmthds, dOmphds)
        real(8) :: s0, taub0, bounceavg0(nvar), ds
        real(8) :: dOmthds, dOmphds, Omth, Omph_noE
        ! store current flux surface values
        s0 = s
        taub0 = taub
        bounceavg0 = bounceavg

!    ds = 1d-7*s0
        ds = 1d-7
        s = s0 - ds/2d0
        call bounce_fast
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
        call bounce_fast
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
        taub = taub0
        bounceavg = bounceavg0
    end subroutine d_Om_ds

    subroutine driftorbit_coarse(eta_min, eta_max, roots, nroots)
        real(8), intent(in) :: eta_min, eta_max
        real(8), intent(out) :: roots(:, :)
        integer, intent(out) :: nroots
        real(8) :: deta
        real(8) :: Omph, dOmphdv, dOmphdeta
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: res, dresdv, dresdeta
        real(8) :: resold, dresdvold, dresdetaold
        integer :: k, ninterv

        ninterv = size(roots, 1)

        deta = (eta_max - eta_min)*1d0/ninterv
        nroots = 0

        do k = 0, ninterv
            eta = eta_min + k*deta
            call Om_th(Omth, dOmthdv, dOmthdeta)
            call Om_ph(Omph, dOmphdv, dOmphdeta)
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

    function driftorbit_nroot(eta_min, eta_max)
        integer :: driftorbit_nroot
        real(8), optional :: eta_min, eta_max
        real(8) :: etamin2, etamax2
        real(8) :: Omph_etamin, Omph_etamax, dummy, dummy2, &
                   Omth_etamin, Omth_etamax, res_etamin, res_etamax

        if (present(eta_min) .and. present(eta_max)) then
            etamin2 = eta_min
            etamax2 = eta_max
        else
            ! default behavior for trapped particles
            etamin2 = etatp*(1d0 + epst)
            etamax2 = etadt*(1d0 - epst)
        end if

        driftorbit_nroot = 0
        ! TODO: return number of possible roots instead of 0 and 1
        eta = etamax2
        call Om_ph(Omph_etamin, dummy, dummy2)
        call Om_th(Omth_etamin, dummy, dummy2)
        eta = etamin2
        call Om_ph(Omph_etamax, dummy, dummy2)
        call Om_th(Omth_etamax, dummy, dummy2)

        res_etamin = mph*Omph_etamin + mth*Omth_etamin
        res_etamax = mph*Omph_etamax + mth*Omth_etamax
        if (sign(1d0, res_etamin) /= sign(1d0, res_etamax)) then
            driftorbit_nroot = 1
        end if
        if (isnan(res_etamin) .or. isnan(res_etamax)) then
            print *, "ERROR: driftorbit_nroot found NaN value in Om_ph_ba"
            return
        end if
    end function driftorbit_nroot

    function driftorbit_root(tol, eta_min, eta_max)
        real(8) :: driftorbit_root(2)
        real(8) :: tol, res, res_old, eta0, eta_old
        real(8) :: Omph, dOmphdv, dOmphdeta
        real(8) :: Omth, dOmthdv, dOmthdeta
        integer :: maxit, k, state
        real(8), intent(in) :: eta_min, eta_max
        real(8) :: etamin2, etamax2
        logical :: slope_pos
        real(8) :: resmin, resmax

        maxit = 100
        state = -2
        eta0 = eta
        eta_old = 0d0
        res = 0d0

        etamin2 = eta_min
        etamax2 = eta_max

        eta = etamin2
        call Om_ph(Omph, dOmphdv, dOmphdeta)
        call Om_th(Omth, dOmthdv, dOmthdeta)
        res = mph*Omph + mth*Omth
        resmin = res

        eta = etamax2
        call Om_ph(Omph, dOmphdv, dOmphdeta)
        call Om_th(Omth, dOmthdv, dOmthdeta)
        resmax = mph*Omph + mth*Omth
        if (resmax - resmin > 0) then
            slope_pos = .true.
        else
            slope_pos = .false.
        end if

        if (driftorbit_nroot(etamin2, etamax2) == 0) then
            print *, "ERROR: driftorbit_root couldn't bracket 0 for v/vth = ", v/vth
            print *, "ERROR: etamin = ", etamin2, " etamax = ", etamax2

            return
        end if

        do k = 1, maxit
            res_old = res
            call Om_ph(Omph, dOmphdv, dOmphdeta)
            call Om_th(Omth, dOmthdv, dOmthdeta)
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

    function find_vmin(vmin0, vmax0)
        real(8) :: find_vmin, tol, vmin0, vmax0, vmin, vmax, v0
        integer, parameter :: nit = 100
        integer :: k
        ! Bisection search for smallest possible v
        v0 = v
        tol = 1d-12*vth
        vmin = vmin0
        vmax = vmax0

        do k = 1, nit
            if (abs(vmax - vmin) < tol) then
                exit
            end if
            v = (vmax + vmin)/2d0
            if (driftorbit_nroot() /= 0) then
                vmax = v
            else
                vmin = v
            end if
        end do
        v = v0
        find_vmin = vmax
    end function find_vmin

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

        !th0 = 0d0 ! TODO remove this

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
        !call disp('init_fsa: fsa/psi_pr = ', fsa/psi_pr)

    end subroutine init_fsa

    function D11int(ux, etax)
        real(8) :: D11int
        real(8) :: ux, etax
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: dummy, dummy2
        real(8) :: OmtB
        real(8) :: Hmn2

        v = ux*vth
        eta = etax
        call Om_th(Omth, dOmthdv, dOmthdeta)
        call Om_tB(OmtB, dummy, dummy2)
        taub = 2d0*pi/abs(Omth)
        call bounce_fast
        Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2

        D11int = pi**(3d0/2d0)*mph**2*c**2*q*vth &
                 /(qi**2*dVds*psi_pr)*ux**3*exp(-ux**2) &
                 *taub*Hmn2 &
                 *nonlinear_attenuation(ux, dOmthdv, dOmthdeta, Hmn2)
    end function D11int

    function D12int(ux, etax)
        real(8) :: D12int
        real(8) :: ux, etax
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: dummy, dummy2
        real(8) :: OmtB
        real(8) :: Hmn2

        v = ux*vth
        eta = etax
        call Om_th(Omth, dOmthdv, dOmthdeta)
        call Om_tB(OmtB, dummy, dummy2)
        taub = 2d0*pi/abs(Omth)
        call bounce_fast
        Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2

        D12int = pi**(3d0/2d0)*mph**2*c**2*q*vth &
                 /(qi**2*dVds*psi_pr)*ux**3*exp(-ux**2) &
                 *taub*Hmn2*ux**2 &
                 *nonlinear_attenuation(ux, dOmthdv, dOmthdeta, Hmn2)
    end function D12int

    function nonlinear_attenuation(ux, dOmthdv, dOmthdeta, Hmn2)
        real(8), intent(in) :: ux, dOmthdv, dOmthdeta, Hmn2
        real(8) :: nonlinear_attenuation

        real(8) :: dpp, dhh, fpeff, dres, dnorm, Omph, dOmphdv, dOmphdeta, dOmdv, &
                   dOmdeta, Ompr, dOmphds, dOmthds, dOmdpph

        if (.not. nonlin) then
            nonlinear_attenuation = 1d0
            return
        end if

        call Om_ph(Omph, dOmphdv, dOmphdeta)
        call d_Om_ds(dOmthds, dOmphds)
        dOmdv = mth*dOmthdv + mph*dOmphdv
        dOmdeta = mth*dOmthdeta + mph*dOmphdeta
        dOmdpph = -(qi/c*iota*psi_pr)**(-1)*(mth*dOmthds + mph*dOmphds)
        Ompr = mth*(eta*dOmdeta - ux*vth/2*dOmdv)/(mi*(ux*vth)**2/(2d0*Omth)) + dOmdpph
        call coleff(ux, dpp, dhh, fpeff)
        dhh = vth*dhh
        dpp = vth**3*dpp
        dres = dpp*(dOmdv/Ompr)**2 + dhh*eta*(bounceavg(6) - eta)*(dOmdeta/Ompr)**2
        dnorm = dres*sqrt(abs(Ompr))/sqrt(abs(Hmn2))**(3d0/2d0)
        call attenuation_factor(dnorm, nonlinear_attenuation)
    end function

    function D11int_u(ux)
        real(8) :: ux
        real(8) :: D11int_u
        real(8) :: eta_res(2)
        real(8) :: roots(nlev, 3)
        integer :: nroots, kr

        v = ux*vth
        D11int_u = 0d0

        call driftorbit_coarse(etamin, etamax, roots, nroots)
        if (nroots == 0) return

        do kr = 1, nroots
            eta_res = driftorbit_root(1d-8*abs(Om_tE), roots(kr, 1), roots(kr, 2))
            eta = eta_res(1)
            D11int_u = D11int_u + D11int(ux, eta_res(1))*1d0/abs(eta_res(2))
        end do

    end function D11int_u

    function D12int_u(ux)
        real(8) :: ux
        real(8) :: D12int_u
        real(8) :: eta_res(2)
        real(8) :: roots(nlev, 3)
        integer :: nroots, kr
        v = ux*vth
        D12int_u = 0d0

        call driftorbit_coarse(etamin, etamax, roots, nroots)
        if (nroots == 0) return

        do kr = 1, nroots
            eta_res = driftorbit_root(1d-8*abs(Om_tE), roots(kr, 1), roots(kr, 2))
            eta = eta_res(1)
            D12int_u = D12int_u + D12int(ux, eta_res(1))*1d0/abs(eta_res(2))
        end do
    end function D12int_u

    function Tphi_int(ux, etax)
        real(8) :: Tphi_int
        real(8) :: ux, etax(2)
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: dummy, dummy2
        real(8) :: OmtB
        real(8) :: Hmn2
        real(8) :: dpp, dhh, fpeff, dres, dnorm, thatt, & ! for nonlin
                   Omph, dOmphdv, dOmphdeta, dOmdv, dOmdeta, Ompr, dOmphds, dOmthds, &
                   dOmdpph
        real(8) :: ma, mb, mc, md, me, mf
        real(8) :: dvdJ, detadJ

        v = ux*vth
        eta = etax(1)
        call Om_th(Omth, dOmthdv, dOmthdeta)
        call Om_tB(OmtB, dummy, dummy2)
        taub = 2d0*pi/abs(Omth)
        call bounce_fast
        Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2

        thatt = 1d0
        if (intoutput .or. nonlin) then
            call Om_ph(Omph, dOmphdv, dOmphdeta)
            call d_Om_ds(dOmthds, dOmphds)
            dOmdv = mth*dOmthdv + mph*dOmphdv
            dOmdeta = mth*dOmthdeta + mph*dOmphdeta
            dOmdpph = -(qi/c*iota*psi_pr)**(-1)*(mth*dOmthds + mph*dOmphds)

            ma = mi*(ux*vth)*mi*c/qi*eta
            mb = mi*(ux*vth)**2/2*mi*c/qi
            mc = mi/(2d0*Omth)*(ux*vth)*(1d0 - eta*bounceavg(7))
            md = mi*(ux*vth)**2/2d0*Omth
            me = -mth/mph
            mf = 1d0/mph

            dvdJ = mb*me/(ma*md*mf - mb*mc*mf)
            detadJ = ma*me/(mb*mc*mf - ma*md*mf)

            !Ompr=mth*(eta*dOmdeta-ux*vth/2*dOmdv)/(mi*(ux*vth)**2/(2d0*Omth))+dOmdpph

            Ompr = dOmdv*dvdJ + dOmdeta*detadJ + mph*dOmdpph

            if (intoutput) then
                ! 0:n, 1:l, 2:Eth, 3:Jperp_tp, 4:drphi/dpphi, 5:E/Eth, 6:Jperp/Jperp_tp, 7:rphi,
                ! 8:|Hmn|, 9:Omth, 10:Omph, 11:Ombarprime, 12:dOmdv, 13:dOmdeta, 14:dOmdpphi, 15:sigma
                ! 16:iota=1/q, 17: Om_tE, 18: dOmthds, 19: dOmphds, 20:Jperp_dt
                write (11, *) mth, mph, mi*vth**2/2d0, mi*(ux*vth)**2/2d0*mi*c/qi*etatp, &
                    -(qi/c*iota*psi_pr)**(-1), ux**2, eta/etatp, s, sqrt(Hmn2), Omth, Omph, &
                    Ompr, dOmdv, dOmdeta, dOmdpph, sigv, iota, Om_tE, dOmthds, dOmphds, &
                    mi*(ux*vth)**2/2d0*mi*c/qi*etadt
            end if
            if (nonlin) then
                call coleff(ux, dpp, dhh, fpeff)
                dhh = vth*dhh
                dpp = vth**3*dpp
                dres = dpp*(dOmdv/Ompr)**2 + dhh*eta*(bounceavg(6) - eta)*(dOmdeta/Ompr)**2
                dnorm = dres*sqrt(abs(Ompr))/sqrt(abs(Hmn2))**(3d0/2d0)
                call attenuation_factor(dnorm, thatt)
            end if
        end if

        Tphi_int = -pi**(3d0/2d0)*mph**2*ni1*c*vth/qi*ux**3*exp(-ux**2)*taub*Hmn2*thatt*(A1 + A2*ux**2)*1d0/abs(etax(2))
    end function Tphi_int

    function Tphi_int_u(ux)
        real(8) :: ux
        real(8) :: Tphi_int_u, dTphi_int_u
        real(8) :: eta_res(2)
        real(8) :: roots(nlev, 3)
        integer :: nroots, kr

        v = ux*vth
        Tphi_int_u = 0d0

        call driftorbit_coarse(etamin, etamax, roots, nroots)
        if (nroots == 0) return

        do kr = 1, nroots
            eta_res = driftorbit_root(1d-8*abs(Om_tE), roots(kr, 1), roots(kr, 2))
            eta = eta_res(1)
            dTphi_int_u = Tphi_int(ux, eta_res)
            Tphi_int_u = Tphi_int_u + dTphi_int_u

            if (orbit_mode_transp > 0) then
                call taurel
                torque_int_box = torque_int_box + dTphi_int_u*taubins
            end if
        end do

    end function Tphi_int_u

    function torque_integral_mid(vmin, vmax)
        ! compute flux integral via midpoint rule
        real(8) :: vmin, vmax
        real(8) :: torque_integral_mid
        real(8) :: ux, du
        integer :: ku

        if (orbit_mode_transp > 0) then
            torque_int_box = 0d0
        end if

        torque_integral_mid = 0d0
        du = (vmax - vmin)/(vsteps*vth)
        ux = vmin/vth + du/2d0

        do ku = 1, vsteps
            torque_integral_mid = torque_integral_mid + Tphi_int_u(ux)*du
            ux = ux + du
        end do

        if (orbit_mode_transp > 0) then
            torque_int_box = torque_int_box*du
        end if

    end function torque_integral_mid

    function flux_integral(vmin, vmax, tol)
        real(8) :: vmin, vmax
        real(8) :: tol(2)
        real(8) :: flux_integral(2), err
        real(8) :: Dp, dsdreff
        integer :: neval, ier

        flux_integral = 0d0

        call qag(D11int_u, &
                 (vmin + (vmax - vmin)*1d-10)/vth, (vmax - (vmax - vmin)*1d-10)/vth, &
                 tol(1), 1d-3, 6, flux_integral(1), err, neval, ier)

        call qag(D12int_u, &
                 (vmin + (vmax - vmin)*1d-10)/vth, (vmax - (vmax - vmin)*1d-10)/vth, &
                 tol(2), 1d-3, 6, flux_integral(2), err, neval, ier)

        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
        dsdreff = 2d0/a*sqrt(s)
        flux_integral = dsdreff**(-2)*flux_integral/Dp

    end function flux_integral

    subroutine taurel
        ! test box counting
        use dvode_f90_m2

        integer :: n

        integer :: k
        real(8) :: ti
        real(8) :: y(2), yold(2)

        real(8) :: atol(nvar), rtol, tout, rstats(22)
        integer :: neq, itask, istate, istats(31), numevents
        type(vode_opts) :: options

        real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(8) :: s1old, told
        integer :: sind, sind0 ! s index

        integer :: jroots(2)

        x(1) = s
        x(2) = 0d0
        x(3) = 0d0
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

        call bounce

        neq = 2
        rtol = 1d-12
        atol = 1d-13
        itask = 1
        istate = 1
        numevents = 2
        options = set_normal_opts(abserr_vector=atol, relerr=rtol, nevents=numevents)

        n = 3*size(sbox)
        taubins = 0d0

        y(1) = 0d0
        y(2) = sigv*vpar(bmod)
        ti = 0d0

        s1 = s + c*mi*y(2)*hcovar(2)*q/(qi*psi_pr)
        sind = size(sbox) + 1
        sprev = -1e5
        snext = 1e5
        do k = 1, size(sbox)
            if (sbox(k) > s1) then
                sind = k
                snext = sbox(k)
                if (k > 1) sprev = sbox(k - 1)
                exit
            end if
        end do
        sind0 = sind

        told = 0d0
        do k = 2, n
            yold = y
            s1old = s1
            tout = taub
            call dvode_f90(tsorb, neq, y, ti, tout, itask, istate, options, &
                           g_fcn=sroots)
            if (istate == 2) exit
            if (istate == 3) then
                taubins(sind) = taubins(sind) + ti - told
                told = ti
                call get_stats(rstats, istats, numevents, jroots)
                if (jroots(2) .ne. 0) then
                    sind = sind + 1 ! moving outwards
                    sprev = snext
                    if (sind == size(sbox) + 1) then
                        snext = 1e5
                    else
                        snext = sbox(sind)
                    end if
                end if
                if (jroots(1) .ne. 0) then
                    sind = sind - 1 ! moving inwards
                    snext = sprev
                    if (sind == 1) then
                        sprev = -1e5
                    else
                        sprev = sbox(sind - 1)
                    end if
                end if
            end if
        end do

        taubins(sind) = taubins(sind) + taub - told

        taubins = taubins/taub
    end subroutine taurel

    function flux_integral_mid(vmin, vmax)
        ! compute flux integral via midpoint rule
        real(8) :: vmin, vmax
        real(8) :: flux_integral_mid(2)
        real(8) :: Dp, dsdreff
        real(8) :: ux, du
        integer :: ku
        real(8) :: dD11, dD12

        if (orbit_mode_transp > 0) then
            fluxint_box = 0d0
        end if

        flux_integral_mid = 0d0
        du = (vmax - vmin)/(vsteps*vth)
        ux = vmin/vth + du/2d0

        do ku = 1, vsteps
            dD11 = D11int_u(ux)*du
            dD12 = D12int_u(ux)*du
            flux_integral_mid(1) = flux_integral_mid(1) + dD11
            flux_integral_mid(2) = flux_integral_mid(2) + dD12
            if (orbit_mode_transp > 0) then
                call taurel
                fluxint_box(1, :) = fluxint_box(1, :) + taubins*dD11
                fluxint_box(2, :) = fluxint_box(2, :) + taubins*dD12
            end if
            ux = ux + du
        end do

        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
        dsdreff = 2d0/a*sqrt(s)
        flux_integral_mid = dsdreff**(-2)*flux_integral_mid/Dp

        if (orbit_mode_transp > 0) then
            fluxint_box = dsdreff**(-2)*fluxint_box/Dp
        end if

    end function flux_integral_mid

    subroutine tsorb(neq, t, y, ydot)
    !
    !  Timestep function for orbit only.
    !  used for radial boxes
    !

        integer, intent(in) :: neq
        real(8), intent(in) :: t
        real(8), intent(in) :: y(neq)
        real(8), intent(out) :: ydot(neq)

        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        associate(dummy => t)
        end associate

        x(1) = s
        x(2) = 0d0
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        ydot(1) = y(2)*hctrvr(3)                                    ! theta
        ydot(2) = -v**2*eta/2d0*hctrvr(3)*hder(3)*bmod              ! v_par

        !s1 = s+c*mi*y(2)*hcovar(2)*q/(qi*psi_pr)
        !print *, ti/taub, y, s1
    end subroutine tsorb

    subroutine sroots(neq, t, y, ng, gout)
        integer, intent(in) :: neq, ng
        real(8), intent(in) :: t, y(neq)
        real(8), intent(out) :: gout(ng)

        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        associate(dummy => t)
        end associate

        x(1) = s
        x(2) = 0d0
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        s1 = s + c*mi*y(2)*hcovar(2)*q/(qi*psi_pr)

        gout(1) = s1 - sprev
        gout(2) = s1 - snext
    end subroutine sroots

end module driftorbit

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
! WIP: integration with ODEs                                                   !
!------------------------------------------------------------------------------!
module lineint
    use driftorbit

contains

    subroutine intstep(neq, t, y, ydot)
        !
        !  Timestep function for orbit integration.
        !  Includes poloidal angle theta and parallel velocity.
        !  More integrands may be added starting from y(3)
        !
        integer, intent(in) :: neq
        real(8), intent(in) :: t
        real(8), intent(in) :: y(neq)
        real(8), intent(out) :: ydot(neq)

        real(8) :: Omph, dOmphdv, dOmphdeta
        real(8) :: Omth, dOmthdv, dOmthdeta
        real(8) :: G, dGdv, dGdeta, absgradG
        real(8) :: deltatp
        real(8) :: Sveta

        associate(dummy => t)
        end associate

        ! check if trapped or passing region
        if (etamin > etatp) then
            deltatp = 0d0
        else
            deltatp = 1d0
        end if

        Sveta = (vmax2 - vmin2)*(etamax - etamin)

        v = vmin2 + y(1)*(vmax2 - vmin2)
        eta = etamin + y(2)*(etamax - etamin)

        call Om_ph(Omph, dOmphdv, dOmphdeta)
        call Om_th(Omth, dOmthdv, dOmthdeta)
        G = (mph*q*deltatp + mth)*Omth + mph*Omph
        dGdv = ((mph*q*deltatp + mth)*dOmthdv + mph*dOmphdv)*(vmax2 - vmin2)
        dGdeta = ((mph*q*deltatp + mth)*dOmthdeta + mph*dOmphdeta)*(etamax - etamin)
        absgradG = sqrt(dGdv**2 + dGdeta**2)

        ydot(1) = dGdeta           ! vbar
        ydot(2) = -dGdv            ! etabar
        taub = 2d0*pi/abs(Omth)
        call bounce_fast
        ydot(3) = (vmax2 - vmin2)*(etamax - etamin)*D11_ode()/vth            ! D11
        ydot(4) = (vmax2 - vmin2)*(etamax - etamin)*D11_ode()*(v/vth)**2/vth ! D12

        ydot = ydot/absgradG

        ! always go in positive velocity direction
        if (etamin > etatp) then
            ydot(1) = -ydot(1)
            ydot(2) = -ydot(2)
        end if
    end subroutine intstep

    function flux_integral_ode(vmin, vmax)
        use dvode_f90_m2
        real(8) :: flux_integral_ode(2)
        real(8) :: vmin, vmax

        real(8) :: sk, ds
        real(8) :: Dp, dsdreff, eta_res(2)
        integer :: ks, nums

        real(8) :: y(4), t, tout, rtol(4), atol(4), rwork(1024)
        integer :: neq, itask, istate, iopt, itol, iwork(1024), ng, jt
        type(VODE_OPTS) :: OPTIONS

        vmin2 = vmin + 1d-8*(vmax - vmin)
        vmax2 = vmax - 1d-8*(vmax - vmin)

        neq = 4
        ng = 4
        y(1) = 1d0 - 1d-15
        y(2) = 1d-15
        y(3) = 0d0
        y(4) = 0d0
        itol = 4
        t = 0.0d0
        rtol(1) = 1d-9
        rtol(2) = 1d-9
        rtol(3) = 1d-9
        rtol(4) = 1d-9
        atol(1) = 1d-9
        atol(2) = 1d-9
        atol(3) = 1d-9
        atol(4) = 1d-9
        itask = 1
        istate = 1
        iopt = 1
        jt = 2
        rwork = 0d0
        iwork = 0
        iwork(6) = 10000

        nums = 1000
        ds = 4d0/nums
        sk = 0d0

        dsdreff = 2d0/a*sqrt(s)
        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)

        OPTIONS = SET_OPTS(DENSE_J=.true., RELERR_VECTOR=RTOL, &
                           ABSERR_VECTOR=ATOL, NEVENTS=NG)

        v = vmax2
        if (mth /= 0) then
            if (etamin < etatp) then
                y(2) = 1d0 - 1d-15
            else
                y(2) = 1d-15
            end if
        else
            eta_res = driftorbit_root(max(1d-9*abs(Om_tE), 1d-12), etamin, etamax)
            y(2) = (eta_res(1) - etamin)/(etamax - etamin)
        end if

        do ks = 1, nums
            print *, ks, y(2)
            tout = t + ds
            call dvode_f90(intstep, neq, y, t, tout, itask, istate, options, g_fcn=bbox)
            if (istate == 3) then
                exit
            end if
            print *, ks
        end do

        flux_integral_ode(1) = dsdreff**(-2)*y(3)/Dp
        flux_integral_ode(2) = dsdreff**(-2)*y(4)/Dp
        print *, flux_integral_ode
    end function flux_integral_ode

    subroutine dummyjac
    end subroutine dummyjac

    subroutine bbox(NEQ, T, Y, NG, GOUT)
        integer, intent(in) :: NEQ, NG
        real(8), intent(in) :: T, Y(neq)
        real(8), intent(out) :: GOUT(ng)

        associate(dummy => T)
        end associate

        GOUT(1) = 1d0
        GOUT(2) = 1d0 - Y(1)
        GOUT(3) = Y(2)
        GOUT(4) = 1d0 - Y(2)
        return
    end subroutine bbox

    function D11_ode()
        real(8) :: D11_ode
        real(8) :: ux
        real(8) :: Hmn2

        ux = v/vth

        Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2

        D11_ode = pi**(3d0/2d0)*mph**2*c**2*q*vth/ &
                  (qi**2*dVds*psi_pr)*ux**3*exp(-ux**2)* &
                  taub*Hmn2
    end function D11_ode

end module lineint
