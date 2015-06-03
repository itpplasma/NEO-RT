! Resonant regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, 2015

module driftorbit
  use do_magfie_mod
  use spline
  use rkf45
  !use quadpack
  
  implicit none
  save

  integer, parameter :: nvar = 5

  ! Orbit parameters
  real(8) :: v, eta
  real(8) :: taub, bounceavg(nvar)

  ! Plasma parameters
  real(8) :: Om_tE, vth, M_t, n0

  ! For splining in the trapped eta region
  integer, parameter :: netaspl = 200
  real(8) :: Omph_spl_coeff(netaspl-1, 5)
  real(8) :: Omth_spl_coeff(netaspl-1, 5)

  ! For splining in the passing eta region
  integer, parameter :: netaspl_pass = 200
  real(8) :: Omth_pass_spl_coeff(netaspl_pass-1, 5)

  ! Harmonics TODO: make dynamic, multiple harmonics
  integer, parameter :: m0 = 0, mph = 18
  integer :: mth
  logical, parameter :: nobdrift = .true.
  logical, parameter :: nopassing = .false.

  ! Flux surface TODO: make a dynamic, multiple flux surfaces support
  real(8) :: fsa, B0, a, etadt, etatp
  real(8) :: etamin, etamax

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
    real(8) :: a, b, delta
    
    v = vth    
    etamin = etatp
    etamax = etadt*(1d0-1d-9)

    delta = 1d-9   ! smallest relative distance to etamin
    b = log(delta)
    a = 1d0/(netaspl-1d0)*(log(etamax/etamin - 1d0) - b)
    
    do k = netaspl-1, 0, -1
       !eta = etamin * (1d0 + (k/(netaspl-1d0))**4*(etamax/etamin-1d0))
       eta = etamin*(1d0 + exp(a*k+b))
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
    
    real(8) :: etarange(netaspl), Omth_v(netaspl)
    real(8) :: delta, a, b
    integer :: k
    
    v = vth    
    !etamin = etatp*1d-9
    !etamax = etatp*(1d0-1d-9)
    etamin = 1d-9
    etamax = etatp
    
    delta = 1d-9   ! smallest relative distance to etamin
    b = log((etamax-etamin)/etamax)
    a = 1d0/(netaspl-1d0)*(log(delta) - b)
    
    do k = netaspl_pass-1, 0, -1
       !eta = etamin * (1d0 + k/(netaspl_pass-1d0)*(etamax/etamin-1d0))
       eta = etamax * (1d0 - exp(a*k+b))
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
    B0 = bmod
    a = 4.6d1
    ! eta deeply trapped and trapped passing limits
    etadt = 1d0/bmod
    x(3) = pi
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    etatp = 1d0/bmod
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

  subroutine timestep(t, y, ydot)
    !
    !  Timestep function for orbit integration.
    !  Includes poloidal angle theta and parallel velocity.
    !  More integrands may be added starting from y(3)
    !

    real(8) :: t, y(*)
    real(8), intent(out) :: ydot(*)
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(8) :: Om_tB_v
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: t0

    x(1) = s
    x(2) = 0d0
    x(3) = y(1)
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    call Om_th(Omth, dOmthdv, dOmthdeta)

    Om_tB_v = mi*c/(2d0*qi*sqrtg*hctrvr(3)*bmod**2)*(&      ! Om_tB/v**2
         -(2d0-eta*bmod)*bmod*hder(1)&
         +2d0*(1d0-eta*bmod)*hctrvr(3)*&
           (dBthcovds+q*dBphcovds)&
     !         +Bphcov*dqds)& TODO: add extra term!
    )

    ydot(1) = y(2)*hctrvr(3)                                    ! theta
    ydot(2) = -v**2*eta/2d0*hctrvr(3)*hder(3)*bmod              ! v_par
    ydot(3) = Om_tB_v                                           ! v_ph

    if (init_done) then
       if (eta > etatp) then
          !t0 = 0.25*2*pi/Omth ! Quarter of bounce time to set phase 0 at tips
          t0 = 0d0
          ydot(4) = (2d0 - eta*bmod)*&
               cos((m0+mph/iota)*y(1)-mth*(t-t0)*Omth) ! Re Hmn
          ydot(5) = (2d0 - eta*bmod)*&
               sin((m0+mph/iota)*y(1)-mth*(t-t0)*Omth) ! Im Hmn
       else
          ydot(4) = (2d0 - eta*bmod)*&
               cos((m0+mph/iota)*y(1)-(mth+mph/iota)*t*Omth) ! Re Hmn
          ydot(5) = (2d0 - eta*bmod)*&
               sin((m0+mph/iota)*y(1)-(mth+mph/iota)*t*Omth) ! Im Hmn
       end if
    else
       ydot(4:5) = 0d0
    end if
  end subroutine timestep

  function findroot(y0, dt, tol)
    !
    !  Finds the root of an orbit after the first turn
    !
    ! TODO: might use Brent's method
    
    real(8) :: y0(nvar), dt, tol, findroot(nvar+1)
    integer :: n

    integer :: k, state, rootstate
    real(8) :: ti, told
    real(8) :: y(nvar), yold(nvar), yp(nvar)
    real(8) :: relerr, abserr

    logical :: pass

    ! check for passing orbit
    pass = .false.
    if (eta < etatp) then
       pass = .true.
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
       relerr = 1d-9
       abserr = 1d-10
       call r8_rkf45 ( timestep, nvar, y, yp, ti, ti+dt, relerr, abserr, state )
       !ti = ti+dt
       !if (k<5) then
       !   print *, k, ti, told+dt, y(1), state
       !end if

       ! check for full turn completed
       if (((.not. pass) .and. (abs(y(1)) < tol)) .or. &
          (pass .and. (abs(y(1)-2d0*pi) < tol)) .or. &
          (pass .and. (abs(y(1)+2d0*pi) < tol))) then
          !print *, "ROOT", k, ti, yold(1), y(1)
          rootstate = 0
          exit
       end if
          
       ! if orbit end is bracketed, do regula falsi search
       if ((.not. pass) .and. (yold(1)<0 .and. y(1)>0)) then
          dt = -yold(1)/(y(1)-yold(1))*dt
          y = yold
          ti = told
          !state = 1
       elseif (pass .and. (yold(1)<2d0*pi .and. y(1)>2d0*pi)) then
          !print *, "HIT", k, yold(1), y(1)
          dt = -(yold(1)-2d0*pi)/(y(1)-yold(1))*dt
          !dt = dt/2d0
          y = yold
          ti = told
          !state = 1
       elseif (pass .and. (yold(1)>-2d0*pi .and. y(1)<-2d0*pi)) then
          dt = -(yold(1)+2d0*pi)/(y(1)-yold(1))*dt
          !dt = dt/2d0
          y = yold
          ti = told
          !state = 1
       end if
    end do
    if (rootstate < 0) then
       print *, "ERROR: findroot did not converge in 500 iterations"
       print *, eta
       !print *, ti
       !print *, y(1)
    end if
    findroot(1)  = ti
    findroot(2:) = y
  end function findroot

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
    y0 = 0d0
    y0(2) = vpar(bmod)

    if (present(taub_estimate)) then
       taub_est = taub_estimate
    else
       taub_est = 2.0*pi/(vperp(bmod)*iota/R0*sqrt(eps/2d0))
    end if
    !print *, 2*pi/taub_est
       
    findroot_res = findroot(y0, taub_est/5d0, 1d-10)
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
    real(8) :: Omph_etamin, Omph_etamax, dummy,&
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
    call Om_ph(Omph_etamin, dummy, dummy)
    call Om_th(Omth_etamin, dummy, dummy)
    eta = etamin2
    call Om_ph(Omph_etamax, dummy, dummy)
    call Om_th(Omth_etamax, dummy, dummy)

    res_etamin = mph*Omph_etamin+mth*Omth_etamin
    res_etamax = mph*Omph_etamax+mth*Omth_etamax
    !print *, v/vth, etamin2, etamax2, res_etamin, res_etamax
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
       !   print *, "driftorbit_root: ", k, res, tol, eta
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
       elseif ((mth >= 0 .and.(((res > 0) .and. (eta > etatp)) .or.&
            ((res < 0) .and. (eta < etatp)))) .or.&
            (mth < 0 .and. (((res < 0) .and. (eta > etatp)) .or.&
            (res > 0) .and. (eta < etatp)))) then
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

    eta = etamin
    call Om_th(Omth, dOmthdv, dOmthdeta)
    vmax = min(vmax,-(mph*Om_tE)/(mth*Omth/v))

    ! passing orbits
    etamin = etatp*1d-7
    etamax = etatp*(1d0-1d-7)
    
    eta = etamin
    call Om_th(Omth, dOmthdv, dOmthdeta)
    if (-(mph*Om_tE)/((q*mph+mth)*Omth/v) > 0) then
       vmin = min(vmin,-(mph*Om_tE)/((q*mph+mth)*Omth/v))
    end if
    
    eta = etamax
    call Om_th(Omth, dOmthdv, dOmthdeta)
    vmax = max(vmax,-(mph*Om_tE)/((q*mph+mth)*Omth/v))
        
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
    real(8) :: vr, v0, eta0
    real(8) :: flux_integral
    real(8) :: eta_res(2), Dp, dpsidreff, dresdeta
    
    flux_integral = 0d0

    v0 = v
    eta0 = eta
    
    ! Integral contribution
    v = vr
    eta_res = driftorbit_root(1d-8*abs(Om_tE), etamin, etamax)
    eta = eta_res(1)
    dresdeta = eta_res(2)

    flux_integral = integrand(vr, eta_res(1))*&
         1d0/abs(dresdeta)       

    ! Normalisation
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    dpsidreff = 2d0/a*sqrt(s)
    flux_integral = dpsidreff**(-2)*flux_integral/Dp

    !call disp("flux_integral: ds/dr = ", dpsidreff)
    !call disp("flux_integral: Dp    = ", Dp)
    
    v = v0
    eta = eta0
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

    call disp('init_fsa: iota       = ', iota)
    call disp('init_fsa: fsa/psi_pr = ', fsa/psi_pr)

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
    real(8) :: Omth, dummy
    real(8) :: Hmn2
    
    v = vx
    vr(1) = vx
    eta = etax
    maxres = maxwellian(vr)
    call Om_th(Omth, dummy, dummy)
    call bounce(2d0*pi/Omth)
    Hmn2 = bounceavg(4)**2 + bounceavg(5)**2
    !print *, bounceavg(4), bounceavg(5)
    !Hmn2 = bounceavg(4)**2
!    integrand =  pi/(2d0*n0)*(2*pi)**3/fsa*mph**2*c/(qi*iota*psi_pr)*&
!         taub*mi**3*v**3*c/(qi*2d0*pi)*(mi*v**2/2d0)**2*abs(bounceavg(4))**2*&
!         maxres(1)
    integrand = (pi*mi*c/(2d0*qi))**2*pi/(iota*psi_pr*fsa)&
         *mph**2*v**7*taub*Hmn2*maxres(1)*mi**3/n0
    ! TODO: factor of 2
  end function integrand

  function integrand_v(vx)
    real(8) :: vx
    real(8) :: integrand_v
    real(8) :: eta_res(2)
    v = vx
    eta_res = driftorbit_root(1d-8*abs(Om_tE), etamin, etamax)
    eta = eta_res(1)
    integrand_v = integrand(v, eta_res(1))*1d0/abs(eta_res(2))  
  end function integrand_v
  
  function flux_integral2(vmin, vmax)
    real(8) :: vmin, vmax, v0, eta0
    real(8) :: flux_integral2, err
    real(8) :: Dp, dpsidreff
    integer :: neval, ier
    
    flux_integral2 = 0d0

    v0 = v
    eta0 = eta
    
    ! Integral contribution
    v = vmin + (vmax-vmin)*1d-10
    v = vmax - (vmax-vmin)*1d-10

    call qag(integrand_v,&
         vmin+(vmax-vmin)*1d-10, vmax-(vmax-vmin)*1d-10,&
         1d0, 1d-5, 6, flux_integral2, err, neval, ier)       

    ! Normalisation
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    dpsidreff = 2d0/a*sqrt(s)
    flux_integral2 = dpsidreff**(-2)*flux_integral2/Dp

    !call disp("flux_integral: ds/dr = ", dpsidreff)
    !call disp("flux_integral: Dp    = ", Dp)
    
    v = v0
    eta = eta0
  end function flux_integral2
end module driftorbit
