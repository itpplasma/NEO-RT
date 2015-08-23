! Resonant regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, 2015

module driftorbit
  use do_magfie_mod
  use spline
  use dvode_f90_m
  !use quadpack
  
  implicit none
  save

  integer, parameter :: nvar = 6

  ! Orbit parameters
  real(8) :: v, eta
  real(8) :: taub, bounceavg(nvar)

  ! Plasma parameters
  real(8) :: Om_tE, vth, M_t, n0

  ! For splining in the trapped eta region
  integer, parameter :: netaspl = 300
  real(8) :: Omph_spl_coeff(netaspl-1, 5)
  real(8) :: Omth_spl_coeff(netaspl-1, 5)

  ! For splining in the passing eta region
  integer, parameter :: netaspl_pass = 300
  real(8) :: Omth_pass_spl_coeff(netaspl_pass-1, 5)

  ! Harmonics TODO: make dynamic, multiple harmonics
  real(8) :: epsmn        ! perturbation amplitude B1/B0
  integer :: m0, mph      ! Boozer toroidal and poloidal perturbation mode
  integer :: mth          ! canonical poloidal mode
  logical :: nobdrift     ! neglect magnetic drift (for drift-orbit res.)
  logical :: nopassing    ! neglect passing particles
  logical :: noshear      ! neglect magnetic shear
  logical :: calcflux     ! calculation based on flux instead of torque

  ! Flux surface TODO: make a dynamic, multiple flux surfaces support
  real(8) :: fsa, etadt, etatp
  real(8) :: etamin, etamax
  real(8) :: vmin2, vmax2 ! permanent vmin vmax for integration bracketing

  ! Check if init is done
  logical :: init_done
contains

  subroutine init
    init_done = .false.
    call do_magfie_init
    call init_misc
    call init_Om_spl
    call init_Om_pass_spl
    call init_fsa
    init_done = .true.
  end subroutine init

  subroutine init_Om_spl
    ! Initialise splines for canonical frequencies of trapped orbits
    
    real(8) :: etarange(netaspl), Om_tB_v(netaspl), Omth_v(netaspl)
    integer :: k
    real(8) :: aa, b, delta
    
    v = vth    
    etamin = etatp
    etamax = etadt*(1d0-1d-9)

    delta = 1d-8   ! smallest relative distance to etamin
    b = log(delta)
    aa = 1d0/(netaspl-1d0)*(log(etamax/etamin - 1d0) - b)
    
    do k = netaspl-1, 0, -1
       !eta = etamin * (1d0 + (k/(netaspl-1d0))**4*(etamax/etamin-1d0))
       eta = etamin*(1d0 + exp(aa*k+b))
       etarange(k+1) = eta
       if (k == netaspl-1) then
          call bounce
          !print *, k, eta, taub, bounceavg(3)*v**2, 2*pi/taub, '*'
       else
          call bounce(taub)
          !print *, k, eta, taub, bounceavg(3)*v**2, 2*pi/taub
       end if
       Om_tB_v(k+1) = bounceavg(3)
       Omth_v(k+1) = 2*pi/(v*taub)
    end do
    
    Omph_spl_coeff = spline_coeff(etarange, Om_tB_v)
    Omth_spl_coeff = spline_coeff(etarange, Omth_v)
  end subroutine init_Om_spl

  
  subroutine init_Om_pass_spl
    ! Initialise splines for canonical frequencies of passing orbits
    
    real(8) :: etarange(netaspl_pass), Omth_v(netaspl_pass)
    real(8) :: delta, aa, b
    integer :: k
    
    v = vth    
    !etamin = etatp*1d-9
    !etamax = etatp*(1d0-1d-9)
    etamin = etatp*1d-9
    etamax = etatp*(1d0-1d-9)
    
    delta = 1d-9   ! smallest relative distance to etamin
    b = log((etamax-etamin)/etamax)
    aa = 1d0/(netaspl_pass-1d0)*(log(delta) - b)
    
    do k = netaspl_pass-1, 0, -1
       !eta = etamin * (1d0 + k/(netaspl_pass-1d0)*(etamax/etamin-1d0))
       eta = etamax * (1d0 - exp(aa*k+b))
       etarange(k+1) = eta
       if (k == netaspl_pass-1) then
          call bounce
          !print *, k, eta, taub, '*'
       else
          call bounce(taub)
          !print *, k, eta, taub
       end if
       Omth_v(k+1) = 2*pi/(v*taub)
    end do
    
    Omth_pass_spl_coeff = spline_coeff(etarange, Omth_v)
  end subroutine init_Om_pass_spl

  subroutine init_misc
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    x(1) = s
    x(2) = 0d0
    x(3) = 0d0
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    ! B0 = bmod
    ! a = 4.6d1
    ! eta deeply trapped and trapped passing limits
    etadt = 1d0/bmod
    x(3) = pi
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    etatp = 1d0/bmod

    !B0 = 1d0/2d0*(1d0/etadt + 1d0/etatp)
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

    !integer, intent (in) :: neq
    real(8), intent (in) :: t
    real(8), intent (in) :: y(*)
    real(8), intent (out) :: ydot(*)
    
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(8) :: Om_tB_v
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: t0
    real(8) :: shearterm

    x(1) = s
    x(2) = 0d0
    x(3) = y(1)
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    !call Om_th(Omth, dOmthdv, dOmthdeta)

    shearterm = Bphcov*dqds
    if (noshear) then
       shearterm = 0
    end if

    Om_tB_v = mi*c/(2d0*qi*sqrtg*hctrvr(3)*bmod**2)*(&      ! Om_tB/v**2
         -(2d0-eta*bmod)*bmod*hder(1)&
         +2d0*(1d0-eta*bmod)*hctrvr(3)*&
           (dBthcovds+q*dBphcovds+shearterm))

    ydot(1) = y(2)*hctrvr(3)                                    ! theta
    ydot(2) = -v**2*eta/2d0*hctrvr(3)*hder(3)*bmod              ! v_par
    ydot(3) = Om_tB_v                                           ! v_ph

    if (init_done) then
       if (eta > etatp) then
          !t0 = 0.25*2*pi/Omth ! Quarter of bounce time to set phase 0 at tips
          t0 = 0d0
!!$          ydot(4) = (2d0 - eta*bmod)*epsmn*&
!!$               cos((m0+q*mph)*y(1)-mth*(t-t0)*Omth) ! Re Hmn
!!$          ydot(5) = (2d0 - eta*bmod)*epsmn*&
!!$               sin((m0+q*mph)*y(1)-mth*(t-t0)*Omth) ! Im Hmn
          ydot(4) = (2d0 - eta*bmod)*epsmn*&
               cos((m0+q*mph)*y(1)+mph*Om_tE*(t-t0)) ! Re Hmn
          ydot(5) = (2d0 - eta*bmod)*epsmn*&
               sin((m0+q*mph)*y(1)+mph*Om_tE*(t-t0)) ! Im Hmn
          ydot(6) = 0d0
       else
!!$          ydot(4) = (2d0 - eta*bmod)*epsmn*&
!!$               cos((m0+q*mph)*y(1)-(mth+q*mph)*t*Omth) ! Re Hmn
!!$          ydot(5) = (2d0 - eta*bmod)*epsmn*&
!!$               sin((m0+q*mph)*y(1)-(mth+q*mph)*t*Omth) ! Im Hmn
          ydot(4) = (2d0 - eta*bmod)*epsmn*&
               cos((m0+q*mph)*y(1)+mph*Om_tE*t) ! Re Hmn
          ydot(5) = (2d0 - eta*bmod)*epsmn*&
               sin((m0+q*mph)*y(1)+mph*Om_tE*t) ! Im Hmn
          ydot(6) = 1d0/(bmod*y(2))                    ! Passing term
       end if
    else
       ydot(4:6) = 0d0
    end if
  end subroutine timestep
  
  function findroot2(y0, dt, tol)
    !
    !  Finds the root of an orbit after the first turn
    !
    ! TODO: might use Brent's method
    
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

    n = 100
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
       !print *, istate

       if (istate == 3) then
          if (passing .or. yold(1)<0) then
             !print *, k
             exit
          end if
          
       end if

       istate = 2
    end do
    if (istate /= 3) then
       print *, "ERROR: findroot2 did not converge after 100 iterations"
       print *, eta, etamin, etamax
       !print *, ti
       !print *, y(1)
    end if
    findroot2(1)  = ti
    findroot2(2:) = y
  end function findroot2
  
  subroutine bounceroots (NEQ, T, Y, NG, GOUT)
     integer, intent(in) :: NEQ, NG
     real(8), intent(in) :: T, Y(neq)
     real(8), intent(out) :: GOUT(ng)
     GOUT(1) = Y(1)
     GOUT(2) = 2d0*pi - Y(1)
     !print *, gout
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
    x(3) = 0d0
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    y0 = 1d-15
    y0(2) = vpar(bmod)

    if (present(taub_estimate)) then
       taub_est = taub_estimate
    else
       taub_est = 2.0*pi/(vperp(bmod)*iota/R0*sqrt(eps/2d0))
    end if
    !print *, 2*pi/taub_est
       
    findroot_res = findroot2(y0, taub_est/5d0, 1d-10)
    taub = findroot_res(1)
    bounceavg = findroot_res(2:)/taub
  end subroutine bounce
  
  subroutine Om_tB(OmtB, dOmtBdv, dOmtBdeta)
    ! returns canonical toroidal frequency
    ! and derivatives w.r.t. v and eta
    real(8), intent(out) :: OmtB, dOmtBdv, dOmtBdeta
    real(8) :: splineval(3)
    splineval = spline_val_0(Omph_spl_coeff, eta)*v**2
    OmtB = splineval(1)
    dOmtBdv = 2*splineval(1)/v
    dOmtBdeta = splineval(2)
  end subroutine Om_tB
    
  subroutine Om_ph(Omph, dOmphdv, dOmphdeta)
    ! returns canonical toroidal frequency
    ! and derivatives w.r.t. v and eta
    real(8), intent(out) :: Omph, dOmphdv, dOmphdeta
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: splineval(3)
    if (nobdrift) then      ! magnetic drift switch PoP 2014 paper
       if (eta > etatp) then
          Omph = Om_tE
          dOmphdv = 0d0
          dOmphdeta = 0d0
       else
          call Om_th(Omth, dOmthdv, dOmthdeta)
          Omph = Om_tE + Omth/iota
          dOmphdv = dOmthdv/iota
          dOmphdeta = dOmthdeta/iota
       end if
    else
       splineval = spline_val_0(Omph_spl_coeff, eta)
       Omph = splineval(1)*v**2 + Om_tE
       dOmphdv = 2d0*splineval(1)*v
       dOmphdeta = splineval(2)*v**2
    end if
  end subroutine Om_ph

  subroutine Om_th(Omth, dOmthdv, dOmthdeta)
    ! returns canonical poloidal frequency
    ! and derivatives w.r.t. v and eta
    real(8), intent(out) :: Omth, dOmthdv, dOmthdeta
    real(8) :: splineval(3)

    if (eta > etatp) then
       splineval = spline_val_0(Omth_spl_coeff, eta)
    else
       splineval = spline_val_0(Omth_pass_spl_coeff, eta)
    end if
    Omth = splineval(1)*v
    dOmthdv = splineval(1)
    dOmthdeta = splineval(2)*v
  end subroutine Om_th
  
  function driftorbit_nroot(maxlevel, eta_min, eta_max)
    integer :: driftorbit_nroot, maxlevel
    real(8), optional :: eta_min, eta_max
    real(8) :: etamin2, etamax2
    real(8) :: Omph_etamin, Omph_etamax, dummy, dummy2, &
         Omth_etamin, Omth_etamax, res_etamin, res_etamax
    
    if (present(eta_min) .and. present(eta_max)) then
       etamin2 = eta_min
       etamax2 = eta_max
    else
       ! default behavior for trapped particles
       etamin2 = etatp*(1d0+1d-8)
       etamax2 = etadt*(1d0-1d-8)
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
    !print *, v/vth, etamin2, etamax2, res_etamin, res_etamax
    if(sign(1d0,res_etamin) /= sign(1d0,res_etamax)) then
       driftorbit_nroot = 1
    !else
       !print *, 'driftorbit_nroot: res_etamin = '
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
    real(8), optional:: eta_min, eta_max
    real(8) :: etamin2, etamax2
    maxit = 100
    state = -2
    eta0 = eta
    eta_old = 0d0
    res = 0d0
    if (present(eta_min) .and. present(eta_max)) then
       etamin2 = eta_min
       etamax2 = eta_max
    else
       ! default behavior for trapped particles
       etamin2 = etatp*(1d0+1d-8)
       etamax2 = etadt*(1d0-1d-8)
    end if

    eta = etamin2
    !if (.not. nobdrift) then
    !   eta = etamax2
    !end if

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

       !if (k>0) then
       !   print *, "driftorbit_root: ", k, res, tol, etamin, eta
       !end if
       
       driftorbit_root(1) = eta

       if (abs(res) < tol) then
          state = 1
          driftorbit_root(2) = mph*dOmphdeta + mth*dOmthdeta
          !eta = eta*(1+1d-10)
          !call Om_ph(Omph, dOmphdv, dOmphdeta)
          !call Om_th(Omth, dOmthdv, dOmthdeta)
          !driftorbit_root(2) = (mph*Omph + mth*Omth - res)/(eta*1d-10)
          exit
       ! TODO: better condition based on slope of function
       elseif ((mth >= 0 .and. res > 0 .and. eta > etatp)      .or.&
               (mth >= -mph*q .and. res < 0 .and. eta < etatp) .or.&
               (mth < 0 .and. res < 0 .and. eta > etatp)       .or.&
               (mth < -mph*q .and. res > 0 .and. eta < etatp)) then
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
       driftorbit_root(1) = 0
       driftorbit_root(2) = 0
       print *, "ERROR: driftorbit_root did not converge in 100 iterations"  
    end if
    eta = eta0
  end function driftorbit_root
  
  function find_vmin(vmin0, vmax0)
    real(8) :: find_vmin, tol, vmin0, vmax0, vmin, vmax, v0
    integer, parameter :: nit = 100
    integer :: k
    ! Bisection search for smallest possible v
    v0 = v
    tol = 1d-10*vth
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
  
  subroutine find_vlim(vmin, vmax)
    real(8) :: vmin, vmax, eta0
    real(8) :: Omth, dOmthdv, dOmthdeta

    if (.not. nobdrift) then
       vmin = find_vmin(vmin, vmax)
       return
    end if

    ! if nobdrift is set for drift-orbit resonances,
    ! there is both a minimum and a maximum frequency
    eta0 = eta
    
    etamin = etatp*(1d0+1d-7)
    etamax = etadt*(1d0-1d-7)

    ! trapped orbits
    eta = etamax
    call Om_th(Omth, dOmthdv, dOmthdeta)
    vmin = max(vmin,-(mph*Om_tE)/(mth*Omth/v))
    !print *, mph, mth, Om_tE, Omth, v, vmin

    eta = etamin
    call Om_th(Omth, dOmthdv, dOmthdeta)
    vmax = min(vmax,-(mph*Om_tE)/(mth*Omth/v))
    !print *, mph, mth, Om_tE, Omth, v, vmax

    ! passing orbits
    etamin = etatp*1d-7
    etamax = etatp*(1d0-1d-7)
    
    eta = etamin
    call Om_th(Omth, dOmthdv, dOmthdeta)
    if (-(mph*Om_tE)/((q*mph+mth)*Omth/v) > 0) then
       vmin = min(vmin,-(mph*Om_tE)/((q*mph+mth)*Omth/v))
    end if
    !print *, mph, mth, Om_tE, Omth, v, vmin
    
    eta = etamax
    call Om_th(Omth, dOmthdv, dOmthdeta)
    vmax = max(vmax,-(mph*Om_tE)/((q*mph+mth)*Omth/v))
    !print *, mph, mth, Om_tE, Omth, v, vmax
        
    eta = eta0
  end subroutine find_vlim
  
 subroutine find_vlim_p(vmin, vmax)
    real(8) :: vmin, vmax, eta0
    real(8) :: Omth, dOmthdv, dOmthdeta

    if (.not. nobdrift) then
       vmin = find_vmin(vmin, vmax)
       return
    end if

    ! if nobdrift is set for drift-orbit resonances,
    ! there is both a minimum and a maximum frequency
    eta0 = eta

    ! passing orbits
    etamin = etatp*1d-7
    etamax = etatp*(1d0-1d-7)
    
    eta = etamin
    call Om_th(Omth, dOmthdv, dOmthdeta)
    vmin = max(vmin,-(mph*Om_tE)/((q*mph+mth)*Omth/v))
    
    eta = etamax
    call Om_th(Omth, dOmthdv, dOmthdeta)
    vmax = min(vmax,-(mph*Om_tE)/((q*mph+mth)*Omth/v))
        
    eta = eta0
  end subroutine find_vlim_p
  
  
  subroutine find_vlim_t(vmin, vmax)
    real(8) :: vmin, vmax, eta0
    real(8) :: Omth, dOmthdv, dOmthdeta

    if (mth == 0) then
       return
    end if
    
    if (.not. nobdrift) then
       vmin = find_vmin(vmin, vmax)
       return
    end if

    ! if nobdrift is set for drift-orbit resonances,
    ! there is both a minimum and a maximum frequency
    eta0 = eta
    
    etamin = etatp*(1d0+1d-7)
    etamax = etadt*(1d0-1d-7)

    ! trapped orbits
    eta = etamax
    call Om_th(Omth, dOmthdv, dOmthdeta)
    vmin = max(vmin,-(mph*Om_tE)/(mth*Omth/v))

    eta = etamin
    call Om_th(Omth, dOmthdv, dOmthdeta)
    vmax = min(vmax,-(mph*Om_tE)/(mth*Omth/v))
        
    eta = eta0
  end subroutine find_vlim_t
  
  function flux_integral(vr)
    real(8) :: vr
    real(8) :: flux_integral
    real(8) :: eta_res(2), Dp, dpsidreff, dresdeta
    
    flux_integral = 0d0

    !v0 = v
    !eta0 = eta
    
    ! Integral contribution
    v = vr
    if (abs(M_t) > 1d-12) then
       eta_res = driftorbit_root(1d-8*abs(Om_tE), etamin, etamax)
    else
       eta_res=driftorbit_root(1d-8*abs(c*mi*vth**2/(2*qi*psi_pr)),etamin,etamax)
    end if
    eta = eta_res(1)
    dresdeta = eta_res(2)

    flux_integral = integrand(vr, eta_res(1))*&
         1d0/abs(dresdeta)       

    ! Normalisation
    Dp = pi*q*vth**3/(16d0*R0*(qi*B0/(mi*c))**2)
    dpsidreff = 2d0/a*sqrt(s)
    flux_integral = dpsidreff**(-2)*flux_integral/Dp

    !call disp("flux_integral: ds/dr = ", dpsidreff)
    !call disp("flux_integral: Dp    = ", Dp)
    
    !v = v0
    !eta = eta0
  end function flux_integral

  
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
    
    fsa = 0d0
    do k = 1,nth
       x(3) = thrange(k)
       call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
       fsa = fsa + sqrtg*dth
       !fsa = fsa + psi_pr*(iota*hcovar(3)+hcovar(2))/bmod*dth
    end do

    fsa = 2*pi*fsa

   ! call disp('init_fsa: iota       = ', iota)
   ! call disp('init_fsa: fsa/psi_pr = ', fsa/psi_pr)

  end subroutine init_fsa
  
  function maxwellian(vrange)
    real(8) :: vrange(:)
    real(8) :: maxwellian(size(vrange))
    real(8) :: T

    T = mi*vth**2/2d0
    maxwellian = n0/(2d0*pi*mi*T)**(3d0/2d0)*exp(-(vrange/vth)**2)
  end function maxwellian

  function integrand(vx, etax)
    real(8) :: integrand
    real(8) :: vx, etax
    real(8) :: vr(1), maxres(1)
    real(8) :: Omth, dummy, dummy2
    real(8) :: Hmn2
    
    v = vx
    vr(1) = vx
    eta = etax
    maxres = maxwellian(vr)
    call Om_th(Omth, dummy, dummy2)
    call bounce(2d0*pi/Omth)
    Hmn2 = bounceavg(4)**2 + bounceavg(5)**2
    !print *, bounceavg(4), bounceavg(5)
    !Hmn2 = bounceavg(4)**2
!    integrand =  pi/(2d0*n0)*(2*pi)**3/fsa*mph**2*c/(qi*iota*psi_pr)*&
!         taub*mi**3*v**3*c/(qi*2d0*pi)*(mi*v**2/2d0)**2*abs(bounceavg(4))**2*&
!         maxres(1)

    if (calcflux) then
       integrand = (pi*mi*c/(2d0*qi))**2*pi/(iota*psi_pr*fsa)&
            *mph*v**7*taub*Hmn2*maxres(1)*mi**3/n0&
            *(mph - (mth+q*mph)*Bphcov*Omth*bounceavg(6))
    else          
       integrand = (pi*mi*c/(2d0*qi))**2*pi/(iota*psi_pr*fsa)&
            *mph**2*v**7*taub*Hmn2*maxres(1)*mi**3/n0
    end if
    ! TODO: factor of 2
  end function integrand

!!$  function integrand_v(vx)
!!$    real(8) :: vx
!!$    real(8) :: integrand_v
!!$    real(8) :: eta_res(2)
!!$    v = vx
!!$    if (abs(M_t) > 1d-12) then
!!$       eta_res = driftorbit_root(1d-8*abs(Om_tE), etamin, etamax)
!!$    else
!!$       eta_res=driftorbit_root(1d-8*abs(c*mi*vth**2/(2*qi*psi_pr)),&
!!$            etamin, etamax)
!!$    end if
!!$    eta = eta_res(1)
!!$    integrand_v = integrand(v, eta_res(1))*1d0/abs(eta_res(2))
!!$  end function integrand_v
  
!!$  function flux_integral2(vmin, vmax)
!!$    real(8) :: vmin, vmax
!!$    real(8) :: flux_integral2, err
!!$    real(8) :: Dp, dpsidreff
!!$    integer :: neval, ier
!!$    
!!$    flux_integral2 = 0d0
!!$
!!$    !v0 = v
!!$    !eta0 = eta
!!$    
!!$    ! Integral contribution
!!$        call qag(integrand_v,&
!!$         vmin+(vmax-vmin)*1d-10, vmax-(vmax-vmin)*1d-10,&
!!$         1d0, 1d-5, 6, flux_integral2, err, neval, ier)       
!!$
!!$    ! Normalisation
!!$    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
!!$    dpsidreff = 2d0/a*sqrt(s)
!!$    flux_integral2 = dpsidreff**(-2)*flux_integral2/Dp
!!$
!!$    !call disp("flux_integral: ds/dr = ", dpsidreff)
!!$    !call disp("flux_integral: Dp    = ", Dp)
!!$    
!!$    !v = v0
!!$    !eta = eta0
!!$  end function flux_integral2

  function D11int(ux, etax)
    real(8) :: D11int
    real(8) :: ux, etax
    real(8) :: Omth, dummy, dummy2
    real(8) :: Hmn2
    
    v = ux*vth
    eta = etax
    call Om_th(Omth, dummy, dummy2)
    call bounce(2d0*pi/Omth)
    Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2

    !if (calcflux) then
    !   integrand = (pi*mi*c/(2d0*qi))**2*pi/(iota*psi_pr*fsa)&
    !        *mph*v**7*taub*Hmn2*maxres(1)*mi**3/n0&
    !        *(mph - (mth+q*mph)*Bphcov*Omth*bounceavg(6))
    !else          
    D11int = pi**(3d0/2d0)*mph**2*c**2*q*vth/&
         (qi**2*fsa*psi_pr)*ux**3*exp(-ux**2)*&
         taub*Hmn2
    !end if
  end function D11int
  
  function integrand_u(ux)
    real(8) :: ux
    real(8) :: integrand_u
    real(8) :: eta_res(2)
    v = ux*vth
    if (abs(M_t) > 1d-12) then
       eta_res = driftorbit_root(1d-8*abs(Om_tE), etamin, etamax)
    else
       eta_res=driftorbit_root(1d-8*abs(c*mi*vth**2/(2*qi*psi_pr)),&
            etamin, etamax)
    end if
    eta = eta_res(1)
    integrand_u = D11int(ux, eta_res(1))*1d0/abs(eta_res(2))
  end function integrand_u
  
  function flux_integral3(vmin, vmax)
    real(8) :: vmin, vmax
    real(8) :: flux_integral3, err
    real(8) :: Dp, dpsidreff
    integer :: neval, ier
    
    flux_integral3 = 0d0

    !v0 = v
    !eta0 = eta
    
    ! Integral contribution
        call qag(integrand_u ,&
         (vmin+(vmax-vmin)*1d-10)/vth, (vmax-(vmax-vmin)*1d-10)/vth,&
         1d-9, 1d-3, 6, flux_integral3, err, neval, ier)       

    ! Normalisation
        
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    dpsidreff = 2d0/a*sqrt(s)
    flux_integral3 = dpsidreff**(-2)*flux_integral3/Dp

    !call disp("flux_integral: ds/dr = ", dpsidreff)
    !call disp("flux_integral: Dp    = ", Dp)
    
    !v = v0
    !eta = eta0
  end function flux_integral3
  
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

    !print *, 'intstep start'

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
    ydot(3) = D11_ode()        ! integral

    ydot = ydot/absgradG
    !print *, ydot
    
    !print *, ydot(1), ydot(2), dGdv, dGdeta, G

    ! always go in positive velocity direction
    if (etamin > etatp) then
       ydot(1) = -ydot(1)
       ydot(2) = -ydot(2)
    end if
    !if (mph*M_t < 0 .OR. mth == 0) then
    !   ydot(1) = -ydot(1)
    !   ydot(2) = -ydot(2)
    !end if
  end subroutine intstep

  !subroutine flux_integral_ode(vmin2, vmax2)
  !  real(8) :: atol(3), rtol, t, tout, y(3), rstats(22)
  !  integer :: neq, itask, istate
  !  type (vode_opts) :: options
            
  !  neq = 3
  !  y(1) = 1d0
  !  y(2) = 0d0
  !  y(3) = 0d0
  !  t = 0.0d0
  !  rtol = 1d-15
  !  atol(1) = 1d-16
  !  atol(2) = 1d-16
  !  atol(3) = 1d-16
  !  itask = 1
  !  istate = 1
  !  options = set_normal_opts(abserr_vector=atol, relerr=rtol)    
    
   ! flux_integral_ode = pi**(3d0/2d0)*mph**2*c**2*q*vth/(qi**2*fsa*psi_pr)*y(3)
  !end subroutine flux_integral_ode
  
  function D11_ode()
    real(8) :: D11_ode
    real(8) :: ux
    real(8) :: Hmn2
    
    ux = v/vth
    
    Hmn2 = (bounceavg(4)**2 + bounceavg(5)**2)*(mi*(ux*vth)**2/2d0)**2
    
    D11_ode = pi**(3d0/2d0)*mph**2*c**2*q*vth/&
         (qi**2*fsa*psi_pr)*ux**3*exp(-ux**2)*&
         taub*Hmn2
  end function D11_ode
  
!!$  function flux_integral3(vmin, vmax)
!!$    real(8) :: vmin, vmax
!!$    real(8) :: flux_integral3, err
!!$    real(8) :: Dp, dpsidreff
!!$    real(8) :: y(3)
!!$    integer :: neval, ier
!!$    
!!$    flux_integral3 = 0d0
!!$    
!!$
!!$    !v0 = v
!!$    !eta0 = eta
!!$    
!!$    ! Integral contribution
!!$    call qag(integrand_v,&
!!$         vmin+(vmax-vmin)*1d-10, vmax-(vmax-vmin)*1d-10,&
!!$         1d0, 1d-5, 6, flux_integral2, err, neval, ier)       
!!$
!!$    ! Normalisation
!!$    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
!!$    dpsidreff = 2d0/a*sqrt(s)
!!$    flux_integral2 = dpsidreff**(-2)*flux_integral2/Dp
!!$
!!$    !call disp("flux_integral: ds/dr = ", dpsidreff)
!!$    !call disp("flux_integral: Dp    = ", Dp)
!!$    
!!$    !v = v0
!!$    !eta = eta0
!!$  end function flux_integral3
end module driftorbit
