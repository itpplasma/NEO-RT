! Resonant regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, 2015

module driftorbit
  use do_magfie_mod
  use neo_magfie_perturbation, only: neo_read_pert_control, neo_read_pert,&
       neo_init_spline_pert, neo_magfie_pert_amp, m_phi
  use spline
  
  implicit none
  save

  integer, parameter :: nvar = 5

  ! Orbit parameters
  real(8) :: v, eta
  real(8) :: taub, bounceavg(nvar)

  ! Plasma parameters
  real(8) :: Om_tE, vth, M_t, n0

  ! For splining in the trapped eta region
  integer, parameter :: netaspl = 100
  real(8) :: OmtB_spl_coeff(netaspl-1, 5)
  real(8) :: Omth_spl_coeff(netaspl-1, 5)
  real(8) :: vres_spl_coeff(netaspl-1, 5)

  ! For splining in the passing eta region
  integer, parameter :: netaspl_pass = 100
  real(8) :: OmtB_pass_spl_coeff(netaspl_pass-1, 5)
  real(8) :: Omth_pass_spl_coeff(netaspl_pass-1, 5)
  real(8) :: vres_pass_spl_coeff(netaspl-1, 5)

  ! Harmonics TODO: make dynamic, multiple harmonics
  real(8) :: epsmn        ! perturbation amplitude B1/B0
  integer :: m0, mph      ! Boozer toroidal and poloidal perturbation mode
  integer :: mth          ! canonical poloidal mode
  logical :: supban       ! calculate superbanana plateau
  logical :: magdrift     ! consider magnetic drift 
  logical :: nopassing    ! neglect passing particles
  logical :: noshear      ! neglect magnetic shear
  logical :: pertfile     ! read perturbation from file with neo_magfie_pert

  ! Flux surface TODO: make a dynamic, multiple flux surfaces support
  real(8) :: dVds, etadt, etatp
  real(8) :: etamin, etamax
  real(8) :: vmin2, vmax2 ! permanent vmin vmax for integration bracketing

  ! Check if init is done
  logical :: init_done

  ! TODO: better B0 calculation (average magnetic field on flux surface)
  real(8) :: B0
  real(8) :: Bmin, Bmax, th0

  real(8), parameter :: epst_spl=1d-6, epsp_spl=1d-6   ! dist to tpb for spline
  real(8), parameter :: epsst_spl=1d-3, epssp_spl=1d-3 ! dist to deep for spline
  real(8), parameter :: epst=1d-8, epsp=1d-8 ! smallest eta distance to tp bound
  real(8) :: k_taub_p, d_taub_p, k_taub_t, d_taub_t ! extrapolation at tp bound
  real(8) :: k_OmtB_p, d_Omtb_p, k_Omtb_t, d_Omtb_t ! extrapolation at tp bound

  ! Number of levels for coarse root finding
  integer, parameter :: nlev = 100
  real(8) :: sigv

  ! Nonlinear calculation switch
  logical :: nonlin
contains

  subroutine init
    init_done = .false.
    call do_magfie_init
    if (pertfile) then
       call neo_read_pert_control
       call neo_read_pert
       call neo_init_spline_pert
       mph = m_phi
    end if
    call init_fsa
    call init_misc
    call init_Om_spl
    call init_Om_pass_spl
    init_done = .true.
  end subroutine init

  subroutine init_Om_spl
    ! Initialise splines for canonical frequencies of trapped orbits
    
    real(8) :: etarange(netaspl), Om_tB_v(netaspl), Omth_v(netaspl)
    integer :: k
    real(8) :: aa, b
    real(8) :: taub0, taub1, leta0, leta1, OmtB0, OmtB1
    
    taub0 = 0d0
    taub1 = 0d0
    leta0 = 0d0
    leta1 = 0d0
    OmtB0 = 0d0
    OmtB1 = 0d0
    
    v = vth    
    etamin = etatp
    etamax = etatp + (etadt-etatp)*(1d0-epsst_spl)

    b = log(epst_spl)
    aa = 1d0/(netaspl-1d0)*(log(etamax/etamin - 1d0) - b)
    
    do k = netaspl-1, 0, -1
       !eta = etamin * (1d0 + (k/(netaspl-1d0))**4*(etamax/etamin-1d0))
       eta = etamin*(1d0 + exp(aa*k+b))
       etarange(k+1) = eta
       if (k == netaspl-1) then
          call bounce
       else
          call bounce(taub)
       end if
       if (magdrift) Om_tB_v(k+1) = bounceavg(3)
       Omth_v(k+1) = 2*pi/(v*taub)
       if (k==0) then
          leta0 = log(eta-etatp)
          taub0 = v*taub
          if (magdrift) OmtB0 = Om_tB_v(k+1)/Omth_v(k+1)
       end if
       if (k==1) then
          leta1 = log(eta-etatp)
          taub1 = v*taub
          if (magdrift) OmtB1 = Om_tB_v(k+1)/Omth_v(k+1)
       end if
    end do
    
    k_taub_t = (taub1-taub0)/(leta1-leta0)
    d_taub_t = taub0 - k_taub_t*leta0
    Omth_spl_coeff = spline_coeff(etarange, Omth_v)
        
    if (magdrift) then
       k_OmtB_t = (OmtB1-OmtB0)/(leta1-leta0)
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
    taub0 = 0d0
    taub1 = 0d0
    leta0 = 0d0
    leta1 = 0d0
    OmtB0 = 0d0
    OmtB1 = 0d0
    
    v = vth
    etamin = etatp*epssp_spl
    etamax = etatp
    
    b = log((etamax-etamin)/etamax)
    aa = 1d0/(netaspl_pass-1d0)*(log(epsp_spl) - b)
    
    do k = netaspl_pass-1, 0, -1
       eta = etamax * (1d0 - exp(aa*k+b))
       etarange(k+1) = eta
       if (k == netaspl_pass-1) then
          call bounce
       else
          call bounce(taub)
       end if
       if (magdrift) Om_tB_v(k+1) = bounceavg(3)
       Omth_v(k+1) = 2*pi/(v*taub)
       if (k==netaspl_pass-2) then
          leta0 = log(etatp-eta)
          taub0 = v*taub
          if (magdrift) OmtB0 = Om_tB_v(k+1)/Omth_v(k+1)
       end if
       if (k==netaspl_pass-1) then
          leta1 = log(etatp-eta)
          taub1 = v*taub
          if (magdrift) OmtB1 = Om_tB_v(k+1)/Omth_v(k+1)
       end if
    end do

    k_taub_p = (taub1-taub0)/(leta1-leta0)
    d_taub_p = taub0 - k_taub_p*leta0
    Omth_pass_spl_coeff = spline_coeff(etarange, Omth_v)
    
    if (magdrift) then
       k_OmtB_p = (OmtB1-OmtB0)/(leta1-leta0)
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
    vpar = v*sqrt(1 - eta*bmod)
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
    integer, intent (in) :: neq
    real(8), intent (in) :: t
    real(8), intent (in) :: y(neq)
    real(8), intent (out) :: ydot(neq)

    call timestep(t, y, ydot)
  end subroutine timestep2
  
  subroutine timestep(t, y, ydot)
    !
    !  Timestep function for orbit integration.
    !  Includes poloidal angle theta and parallel velocity.
    !  More integrands may be added starting from y(3)
    !

    real(8), intent (in) :: t
    real(8), intent (in) :: y(*)
    real(8), intent (out) :: ydot(*)
    
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
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    call Om_th(Omth, dOmthdv, dOmthdeta)
    Omth = abs(Omth)
    
    if (pertfile) then
       call neo_magfie_pert_amp( x, epsn )
       epsn = epsn/bmod
    else
       epsn = epsmn*exp(imun*m0*y(1))
    end if

    shearterm = Bphcov*dqds
    if (noshear) then
       shearterm = 0
    end if

    Om_tB_v = mi*c*q/(2d0*qi*psi_pr*bmod)*(&      ! Om_tB/v**2
         -(2d0-eta*bmod)*bmod*hder(1)&
         +2d0*(1d0-eta*bmod)*hctrvr(3)*&
         (dBthcovds+q*dBphcovds+shearterm))
   
    ydot(1) = y(2)*hctrvr(3)                                    ! theta
    ydot(2) = -v**2*eta/2d0*hctrvr(3)*hder(3)*bmod              ! v_par
    ydot(3) = Om_tB_v                                           ! v_ph

    if (init_done) then
       if (eta > etatp) then
          !t0 = 0.25*2*pi/Omth ! Different starting position in orbit
          t0 = 0d0
          Hn = (2d0-eta*bmod)*epsn*exp(imun*(q*mph*(y(1))-mth*(t-t0)*Omth))
       else
          Hn = (2d0-eta*bmod)*epsn*exp(imun*(q*mph*(y(1))-(mth+q*mph)*t*Omth))
       end if
       ydot(4) = real(Hn)
       ydot(5) = aimag(Hn)
    else
       ydot(4:5) = 0d0
    end if
  end subroutine timestep
  
  function findroot2(y0, dt, tol)
    !
    !  Finds the root of an orbit after the first turn
    !
    use dvode_f90_m
    
    real(8) :: y0(nvar), dt, tol, findroot2(nvar+1)
    integer :: n

    integer :: k, state, rootstate
    real(8) :: ti, told
    real(8) :: y(nvar), yold(nvar)

    logical :: passing
    
    real(8) :: atol(nvar), rtol, tout
    integer :: neq, itask, istate
    type (vode_opts) :: options

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
    do k = 2,n
       yold = y
       told = ti

       tout = ti+dt
       call dvode_f90(timestep2, neq, y, ti, tout, itask, istate, options,&
            g_fcn = bounceroots)
       !print *, ti
       if (istate == 3) then
          if (passing .or. (yold(1)-th0)<0) then
             exit
          end if
          
       end if

       istate = 2
    end do
    if (istate /= 3) then
       print *, "ERROR: findroot2 did not converge after 500 iterations"
       print *, eta, etamin, etamax, y(1)
    end if
    
    findroot2(1)  = ti
    findroot2(2:) = y
  end function findroot2
  
  subroutine bounceroots (NEQ, T, Y, NG, GOUT)
     integer, intent(in) :: NEQ, NG
     real(8), intent(in) :: T, Y(neq)
     real(8), intent(out) :: GOUT(ng)
     GOUT(1) = Y(1) - th0
     GOUT(2) = 2d0*pi - (Y(1) - th0)
     return
   end subroutine bounceroots

  subroutine bounce(taub_estimate)
    ! calculate all bounce averages
    real(8) :: findroot_res(nvar+1)
    real(8), optional :: taub_estimate
    real(8) :: taub_est
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(8) :: y0(nvar)

    x(1) = s
    x(2) = 0d0
    x(3) = th0
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    y0 = 1d-15
    y0(1) = th0
    y0(2) = vpar(bmod)
    y0(3) = 0d0
    y0(4) = 0d0
    y0(5) = 0d0

    if (present(taub_estimate)) then
       taub_est = taub_estimate
    else
       taub_est = 2.0*pi/(vperp(bmod)*iota/R0*sqrt(eps/2d0))
    end if
       
    findroot_res = findroot2(y0, taub_est/5d0, 1d-10)
    !print *, taub
    taub = findroot_res(1)
    bounceavg = findroot_res(2:)/taub
  end subroutine bounce
  
  subroutine Om_tB(OmtB, dOmtBdv, dOmtBdeta)
    ! returns bounce averaged toroidal magnetic drift frequency
    ! and derivatives w.r.t. v and eta
    real(8), intent(out) :: OmtB, dOmtBdv, dOmtBdeta
    real(8) :: splineval(3)
    real(8) :: Omth, dOmthdv, dOmthdeta
    if (eta > etatp) then
       if (eta > etatp*(1+epst_spl)) then
          splineval = spline_val_0(OmtB_spl_coeff, eta)
       else ! extrapolation
          call Om_th(Omth, dOmthdv, dOmthdeta)
          splineval(1) = sigv*(k_OmtB_t*log(eta-etatp) + d_OmtB_t)*Omth/v
          splineval(2) = sigv*(Omth/v*k_OmtB_t/(eta-etatp) +&
               dOmthdeta/v*(k_OmtB_t*log(eta-etatp) + d_OmtB_t))
       end if
    else
       if (eta < etatp*(1-epsp_spl)) then
          splineval = spline_val_0(OmtB_pass_spl_coeff, eta)
       else ! extrapolation
          call Om_th(Omth, dOmthdv, dOmthdeta)
          splineval(1) = sigv*(k_OmtB_p*log(etatp-eta) + d_OmtB_p)*Omth/v
          splineval(2) = sigv*(Omth/v*k_OmtB_p/(eta-etatp) +&
               dOmthdeta/v*(k_OmtB_p*log(etatp-eta) + d_OmtB_p))
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
       if (eta > etatp*(1+epst_spl)) then
          splineval = spline_val_0(Omth_spl_coeff, eta)
       else ! extrapolation
          splineval(1) = 2d0*pi/(k_taub_t*log(eta-etatp) + d_taub_t)
          splineval(2) = -splineval(1)**2/(2d0*pi) * k_taub_t/(eta-etatp)
       end if
    else
       if (eta < etatp*(1-epsp_spl)) then
          splineval = spline_val_0(Omth_pass_spl_coeff, eta)
       else ! extrapolation
          splineval(1) = 2d0*pi/(k_taub_p*log(etatp-eta) + d_taub_p)
          splineval(2) = -splineval(1)**2/(2d0*pi) * k_taub_p/(eta-etatp)
       end if
    end if
    Omth = sigv*splineval(1)*v
    dOmthdv = sigv*splineval(1)
    dOmthdeta = sigv*splineval(2)*v
  end subroutine Om_th

  subroutine driftorbit_coarse(eta_min, eta_max, roots, nroots)
    real(8), intent(in) :: eta_min, eta_max
    real(8), intent(out) :: roots(:,:)
    real(8) :: deta
    real(8) :: Omph, dOmphdv, dOmphdeta
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: res, dresdv, dresdeta
    real(8) :: resold, dresdvold, dresdetaold
    integer :: k, ninterv, nroots

    ninterv = size(roots,1)
    
    deta = (eta_max-eta_min)*1d0/ninterv
    nroots = 0
    
    do k = 0,ninterv
       eta = eta_min + k*deta
       call Om_th(Omth, dOmthdv, dOmthdeta)
       call Om_ph(Omph, dOmphdv, dOmphdeta)
       res = mth*Omth + mph*Omph
       dresdv = mth*dOmthdv + mph*dOmphdv
       dresdeta = mth*dOmthdeta + mph*dOmphdeta
       if (k>0) then
          if (sign(1d0,res) /= sign(1d0,resold)) then
             nroots = nroots+1
             roots(nroots, 1) = eta - deta
             roots(nroots, 2) = eta
          end if
       end if
       resold = res
       dresdvold = dresdv
       dresdetaold = dresdeta
    end do
  end subroutine driftorbit_coarse
  
  function driftorbit_nroot(nlev, eta_min, eta_max)
    integer :: driftorbit_nroot, nlev
    real(8), optional :: eta_min, eta_max
    real(8) :: etamin2, etamax2
    real(8) :: Omph_etamin, Omph_etamax, dummy, dummy2, &
         Omth_etamin, Omth_etamax, res_etamin, res_etamax
    
    if (present(eta_min) .and. present(eta_max)) then
       etamin2 = eta_min
       etamax2 = eta_max
    else
       ! default behavior for trapped particles
       etamin2 = etatp*(1d0+epst)
       etamax2 = etadt*(1d0-epst)
    end if

    driftorbit_nroot = 0
    ! TODO: return number of possible roots instead of 0 and 1
    eta = etamax2
    call Om_ph(Omph_etamin, dummy, dummy2)
    call Om_th(Omth_etamin, dummy, dummy2)
    eta = etamin2
    call Om_ph(Omph_etamax, dummy, dummy2)
    call Om_th(Omth_etamax, dummy, dummy2)

    res_etamin = mph*Omph_etamin+mth*Omth_etamin
    res_etamax = mph*Omph_etamax+mth*Omth_etamax
    if(sign(1d0,res_etamin) /= sign(1d0,res_etamax)) then
       driftorbit_nroot = 1
    end if
    if(isnan(res_etamin) .or. isnan(res_etamax)) then
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
    if(resmax - resmin > 0) then
       slope_pos = .true.
    else
       slope_pos = .false.
    end if
   
    if(driftorbit_nroot(1, etamin2, etamax2) == 0) then
       print *, "ERROR: driftorbit_root couldn't bracket 0 for v/vth = ", v/vth
       print *, "ERROR: etamin = ", etamin2, " etamax = ", etamax2
       
       return
    end if
    
    do k = 1,maxit
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
          eta = (eta+etamin2)/2d0
       else
          etamin2 = eta
          eta_old = eta
          eta = (eta+etamax2)/2d0
       end if
    end do
    if (state < 0) then
       driftorbit_root(2) = mph*dOmphdeta + mth*dOmthdeta
       print *, "ERROR: driftorbit_root did not converge in 100 iterations"
       print *, "v/vth  = ", v/vth,   "mth    = ", mth,     "sigv= ", sigv
       print *, "etamin = ", eta_min, "etamax = ", eta_max, "eta = ", eta
       print *, "resmin = ", resmin,  "resmax = ", resmax,  "res = ", res
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
    
    do k=1,nit
       if(abs(vmax-vmin) < tol) then
          exit
       end if
       v = (vmax + vmin)/2d0
       if (driftorbit_nroot(1) /= 0) then
          vmax = v
       else
          vmin = v
       end if
    end do
    v = v0
    find_vmin = vmax
  end function find_vmin
  
  subroutine init_fsa
    ! Calculate the flux surface areas
    integer, parameter :: nth = 1000
    integer :: k
    real(8) :: thrange(nth), dth
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

    thrange = -pi + (/(k*2*pi/nth,k=1,nth)/)
    
    dth = thrange(2) - thrange(1) 
    x(1) = s
    x(2) = 0d0
    x(3) = 0d0
    
    dVds = 0d0
    B0  = 0d0
    print *, "eps orig: ", eps
    eps = 0d0

    Bmin = -1d0
    Bmax = 0d0
    
    do k = 1, nth
       x(3) = thrange(k)
       call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
       dVds = dVds + sqrtg*dth
       B0   = B0  + bmod*dth
       eps  = eps - cos(x(3))*bmod*dth

       ! TODO: do fine search
       if ((Bmin < 0) .or. (bmod < Bmin)) then
          Bmin = bmod
          th0 = x(3)
       end if
       if (bmod > Bmax) Bmax = bmod
    end do

    print *, th0
    !th0 = 0d0 ! TODO remove this

    dVds = 2d0*pi*dVds
    B0   = B0/(2d0*pi)
    eps  = eps/(B0*pi)
    print *, "eps calc:  ", eps
    print *, "Bmin,Bmax: ", Bmin, Bmax

   ! call disp('init_fsa: iota       = ', iota)
   ! call disp('init_fsa: fsa/psi_pr = ', fsa/psi_pr)

  end subroutine init_fsa
  
  function D11int(ux, etax)
    real(8) :: D11int
    real(8) :: ux, etax
    real(8) :: Omth, dummy, dummy2
    real(8) :: OmtB
    real(8) :: Hmn2
    
    v = ux*vth
    eta = etax
    call Om_th(Omth, dummy, dummy2)
    call Om_tB(OmtB, dummy, dummy2)
    call bounce(2d0*pi/abs(Omth))
    Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2
    D11int = pi**(3d0/2d0)*mph**2*c**2*q*vth/&
         (qi**2*dVds*psi_pr)*ux**3*exp(-ux**2)*&
         taub*Hmn2
  end function D11int
  
  function D12int(ux, etax)
    real(8) :: D12int
    real(8) :: ux, etax
    real(8) :: Omth, dummy, dummy2
    real(8) :: OmtB
    real(8) :: Hmn2
    
    v = ux*vth
    eta = etax
    call Om_th(Omth, dummy, dummy2)
    call Om_tB(OmtB, dummy, dummy2)
    call bounce(2d0*pi/abs(Omth))
    Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2     
    D12int = pi**(3d0/2d0)*mph**2*c**2*q*vth/&
         (qi**2*dVds*psi_pr)*ux**3*exp(-ux**2)*&
         taub*Hmn2*ux**2
  end function D12int
  
  function D11int_u(ux)
    real(8) :: ux
    real(8) :: D11int_u
    real(8) :: eta_res(2)
    real(8) :: roots(nlev, 3)
    integer :: nroots, kr
    
    v = ux*vth
    D11int_u = 0d0

    call driftorbit_coarse(etamin, etamax, roots, nroots)
    if(nroots == 0) return
    
    do kr = 1,nroots
       eta_res = driftorbit_root(1d-8*abs(Om_tE), roots(kr,1), roots(kr,2))
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
    if(nroots == 0) return
    
    do kr = 1,nroots       
       eta_res = driftorbit_root(1d-8*abs(Om_tE), roots(kr,1), roots(kr,2))
       eta = eta_res(1)
       D12int_u = D12int_u + D12int(ux, eta_res(1))*1d0/abs(eta_res(2))
    end do
  end function D12int_u
  
  function flux_integral(vmin, vmax, tol)
    real(8) :: vmin, vmax
    real(8) :: tol(2)
    real(8) :: flux_integral(2), err
    real(8) :: Dp, dsdreff
    integer :: neval, ier
      
    flux_integral = 0d0
      
    call qag(D11int_u ,&
         (vmin+(vmax-vmin)*1d-10)/vth, (vmax-(vmax-vmin)*1d-10)/vth,&
         tol(1), 1d-3, 6, flux_integral(1), err, neval, ier)       

    call qag(D12int_u ,&
         (vmin+(vmax-vmin)*1d-10)/vth, (vmax-(vmax-vmin)*1d-10)/vth,&
         tol(2), 1d-3, 6, flux_integral(2), err, neval, ier)   
      
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    dsdreff = 2d0/a*sqrt(s)
    flux_integral = dsdreff**(-2)*flux_integral/Dp

  end function flux_integral

!------------------------------------------------------------------------------!
!------------------------------------------------------------------------------!
! WIP: integration with ODEs                                                   !
!------------------------------------------------------------------------------!
end module driftorbit
module lineint
  use driftorbit

  contains
  
  subroutine intstep(neq, t, y, ydot)
    !
    !  Timestep function for orbit integration.
    !  Includes poloidal angle theta and parallel velocity.
    !  More integrands may be added starting from y(3)
    !
    integer, intent (in) :: neq
    real(8), intent (in) :: t
    real(8), intent (in) :: y(neq)
    real(8), intent (out) :: ydot(neq)

    real(8) :: Omph, dOmphdv, dOmphdeta
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: G, dGdv, dGdeta, absgradG
    real(8) :: deltatp
    real(8) :: Sveta

    ! check if trapped or passing region
    if (etamin > etatp) then
       deltatp = 0d0
    else
       deltatp = 1d0
    end if

    Sveta = (vmax2-vmin2)*(etamax-etamin)

    v = vmin2 + y(1)*(vmax2-vmin2)
    eta = etamin + y(2)*(etamax-etamin)
    
    call Om_ph(Omph, dOmphdv, dOmphdeta)
    call Om_th(Omth, dOmthdv, dOmthdeta)
    G = (mph*q*deltatp + mth)*Omth + mph*Omph
    dGdv = ((mph*q*deltatp + mth)*dOmthdv + mph*dOmphdv)*(vmax2-vmin2)
    dGdeta = ((mph*q*deltatp + mth)*dOmthdeta + mph*dOmphdeta)*(etamax-etamin)
    absgradG = sqrt(dGdv**2 + dGdeta**2)
    
    ydot(1) = dGdeta           ! vbar
    ydot(2) = -dGdv            ! etabar
    call bounce(2d0*pi/Omth)
    ydot(3) = (vmax2-vmin2)*(etamax-etamin)*D11_ode()/vth            ! D11
    ydot(4) = (vmax2-vmin2)*(etamax-etamin)*D11_ode()*(v/vth)**2/vth ! D12
    
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
    integer :: ks, ns

    real(8) :: y(4), t, tout, rtol(4), atol(4), rwork(1024)
    integer :: neq, itask, istate, iopt, itol, iwork(1024), ng, jt
    type(VODE_OPTS) :: OPTIONS

    vmin2 = vmin + 1d-8*(vmax-vmin)
    vmax2 = vmax - 1d-8*(vmax-vmin)

    neq = 4
    ng = 4
    y(1) = 1d0-1d-15
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

    ns = 1000
    ds = 4d0/ns
    sk = 0d0

    dsdreff = 2d0/a*sqrt(s)
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    
    OPTIONS = SET_OPTS(DENSE_J=.true.,RELERR_VECTOR=RTOL, &
         ABSERR_VECTOR=ATOL,NEVENTS=NG)

    v = vmax2
    if (mth /= 0) then
       if (etamin < etatp) then
          y(2) = 1d0 - 1d-15
       else
          y(2) = 1d-15
       end if
    else
       eta_res = driftorbit_root(max(1d-9*abs(Om_tE),1d-12), etamin, etamax)
       y(2) = (eta_res(1)-etamin)/(etamax-etamin)
    end if    
    
    do ks = 1, ns
       print *, ks, y(2)
       tout = t+ds
       call dvode_f90(intstep,neq,y,t,tout,itask,istate,options,g_fcn=bbox)
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

  subroutine bbox (NEQ, T, Y, NG, GOUT)
    integer, intent(in) :: NEQ, NG
    real(8), intent(in) :: T, Y(neq)
    real(8), intent(out) :: GOUT(ng)
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
    
    D11_ode = pi**(3d0/2d0)*mph**2*c**2*q*vth/&
         (qi**2*dVds*psi_pr)*ux**3*exp(-ux**2)*&
         taub*Hmn2
  end function D11_ode

end module lineint
