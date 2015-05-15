! Resonant regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, 2015

module driftorbit
  use do_magfie_mod
  use spline
  
  implicit none
  save

  integer, parameter :: nvar = 4

  ! Orbit parameters
  real(8) :: v, eta
  real(8) :: taub, bounceavg(nvar)

  ! Plasma parameters
  real(8) :: Om_tE, vth, M_t, n0

  ! For splining Om_tB
  integer, parameter :: netaspl = 200
  real(8) :: Om_spl_coeff(netaspl-1, 5)

  ! Harmonics TODO: make dynamic, multiple harmonics
  integer, parameter :: mph = 3, mth = 0

  ! Flux surface TODO: make a dynamic
  real(8) :: fsa, B0, a
contains

  subroutine init
    call do_magfie_init
    call init_misc
    call init_Om_tB_spl
    call init_fsa
  end subroutine init
 
  subroutine init_Om_tB_spl
    real(8) :: etarange(netaspl), Om_tB_v(netaspl)
    real(8) :: etamin, etamax
    integer :: k
    v = vth    
    etamin = etatp()*(1d0+1d-10)
    etamax = etadt()*(1d0-1d-10)
    
    do k = 0, netaspl-1
       eta = etamin * (1d0 + (k/(netaspl-1d0))**4*(etamax/etamin-1d0))
       etarange(k+1) = eta
       call bounce
       Om_tB_v(k+1) = bounceavg(3)
    end do
    Om_spl_coeff = spline_coeff(etarange, Om_tB_v)
  end subroutine init_Om_tB_spl

  subroutine init_misc
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    x(1) = s
    x(2) = 0d0
    x(3) = 0d0
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    B0 = bmod
    a = 4.6d1
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

  function etatp()
    ! returns eta at trapped passing boundary
    real(8) :: etatp
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    x(1) = s
    x(2) = 0d0
    x(3) = pi
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    etatp = 1d0/bmod
  end function etatp
  
  function etadt()
    ! returns maximum eta for deeply trapped orbits
    real(8) :: etadt
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    x(1) = s
    x(2) = 0d0
    x(3) = 0d0
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    etadt = 1d0/bmod
  end function etadt
  
  subroutine jac 
  end subroutine jac

  subroutine timestep(neq, t, y, ydot)
    !
    !  Timestep function for orbit integration.
    !  Includes poloidal angle theta and parallel velocity.
    !  More integrands may be added starting from y(3)
    !
    
    integer :: neq
    real(8) :: t, y(neq)
    real(8), intent(out) :: ydot(neq)
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(8) :: Om_tB_v

    x(1) = s
    x(2) = 0d0
    x(3) = y(1)
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )

    Om_tB_v = mi*c/(2d0*qi*sqrtg*hctrvr(3)*bmod**2)*(&      ! Om_tB/v**2
         -(2d0-eta*bmod)*bmod*hder(1)&
         +2d0*(1d0-eta*bmod)*hctrvr(3)*&
           (dBthcovds+q*dBphcovds)&
     !         +Bphcov*dqds)& TODO: add extra term!
    )

    ydot(1) = y(2)*hctrvr(3)                                ! theta
    ydot(2) = -v**2*eta/2d0*hctrvr(3)*hder(3)*bmod          ! v_par
    ydot(3) = Om_tB_v                                       ! v_ph
    ydot(4) = (2d0 - eta*bmod)*cos(mph*y(1)/iota)           ! Hmn/(mv**2/2)
  end subroutine timestep

  function findroot(y0, dt, tol)
    !
    !  Finds the root of an orbit after the first turn
    !
    ! TODO: might use Brent's method
    
    real(8) :: y0(nvar), dt, tol, findroot(nvar+1)
    integer :: n

    integer :: k, state, rootstate
    integer :: iwork(1000)
    real(8) :: ti, told
    real(8) :: rwork(1000), y(nvar), yold(nvar)

    n = 500
    rootstate = -1
    
    y = y0
    yold = y0
    ti = 0d0
    state = 1
    do k = 2,n
       yold = y
       told = ti
       call dlsode(timestep, nvar, y, ti, ti+dt, 1, 1d-13, 1d-13, 1, &
            state, 0, rwork, 1000, iwork, 1000, jac, 10)
       !if (k > 90) then
       !   print *, y(1)
       !end if
       if (abs(y(1)) < tol) then
          rootstate = 0
          exit
       end if
       if (yold(1)<0 .and. y(1)>0) then       
          dt = -yold(1)/(y(1)-yold(1))*dt ! guess for time interval to root
          y = yold
          ti = told
          state = 1
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

  subroutine bounce
    ! calculate all bounce averages
    real(8) :: findroot_res(nvar+1)
    real(8) :: taub_wesson
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(8) :: y0(nvar)

    x(1) = s
    x(2) = 0d0
    x(3) = 0d0
    call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
    y0 = 0d0
    y0(2) = vpar(bmod)
    taub_wesson = 2.0*pi/(vperp(bmod)*iota/R0*&
         sqrt(eps/2d0))

    findroot_res = findroot(y0, taub_wesson/10d0, 1d-13)
    taub = findroot_res(1)
    bounceavg = findroot_res(2:)/taub
  end subroutine bounce
    
  subroutine Om_ph(Omph, dOmphdv, dOmphdeta)
    ! returns canonical toroidal frequency
    ! and derivatives w.r.t. v and eta
    real(8), intent(out) :: Omph, dOmphdv, dOmphdeta
    real(8) :: Om_ph_can(3)
    Om_ph_can = spline_val_0(Om_spl_coeff, eta)*v**2
    Omph = Om_ph_can(1) + Om_tE
    dOmphdv = 2*Om_ph_can(1)/v
    dOmphdeta = Om_ph_can(2)
  end subroutine Om_ph

  function driftorbit_nroot(maxlevel)
    integer :: driftorbit_nroot, maxlevel
    real(8) :: etamin, etamax
    real(8) :: Omph_etamin, Omph_etamax, dOmphdv, dOmphdeta,&
         Omth_etamin, Omth_etamax, res_etamin, res_etamax
    
    etamin = etatp()*(1d0+1d-10)
    etamax = etadt()*(1d0-1d-10)
    
    driftorbit_nroot = 0
    ! TODO: return number of possible roots instead of 0 and 1
    eta = etamax
    call Om_ph(Omph_etamin, dOmphdv, dOmphdeta)
    Omth_etamin = 2d0*pi/taub
    eta = etamin
    call Om_ph(Omph_etamax, dOmphdv, dOmphdeta)
    Omth_etamin = 2d0*pi/taub

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

  function driftorbit_root(tol)
    real(8) :: driftorbit_root(2)
    real(8) :: tol, res, res_old, eta0, eta_old
    real(8) :: Omph, dOmphdv, dOmphdeta
    integer :: maxit, k, state
    real(8) :: etamin, etamax
    maxit = 100
    state = -2
    eta0 = eta
    eta_old = 0d0
    res = 0d0
    etamin = etatp()*(1d0+1d-10)
    etamax = etadt()*(1d0-1d-10)
    eta = (etamax+etamin)/2d0
    if(driftorbit_nroot(1) == 0) then
       print *, "ERROR: driftorbit_root couldn't bracket 0 for v/vth = ", v/vth
       return
    end if
    do k = 1,maxit
       res_old = res
       call Om_ph(Omph, dOmphdv, dOmphdeta)
       res = mph*Omph + mth*2d0*pi/taub

       if (k>90) then
          print *, "driftorbit_root: ", k, res, tol, eta
       end if
       
       driftorbit_root(1) = eta

       if (abs(res) < tol) then
          state = 1
          !driftorbit_root(2) = dOmphdeta
          eta = eta*(1+1d-10)
          call Om_ph(Omph, dOmphdv, dOmphdeta)
          driftorbit_root(2) = dOmphdeta
          driftorbit_root(2) = (mph*Omph + mth*2d0*pi/taub - res)/(eta*1d-10)
          exit       
       elseif (res > 0) then
          etamax = eta
          eta_old = eta
          eta = (eta+etamin)/2d0
       else
          etamin = eta
          eta_old = eta
          eta = (eta+etamax)/2d0
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
  
  function flux_integral(vrange)
    real(8) :: vrange(:), v0, eta0
    real(8) :: flux_integral(size(vrange))
    real(8) :: eta_res(2), Dp, Om_tB_ref, dpsidreff, dresdeta
    integer :: k

    flux_integral = 0d0

    v0 = v
    eta0 = eta
    
    ! Integral contributions
    do k=lbound(vrange,1), ubound(vrange,1)
       v = vrange(k)
       eta_res = driftorbit_root(1d-8*abs(Om_tE))
       eta = eta_res(1)
       dresdeta = eta_res(2)

       flux_integral(k) = integrand(vrange(k), eta_res(1))*&
            1d0/abs(dresdeta)       
    end do

    ! Normalisation
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    dpsidreff = 2d0/a*sqrt(s)
    flux_integral = dpsidreff**(-2)*flux_integral/Dp

    call disp("flux_integral: ds/dr = ", dpsidreff)
    call disp("flux_integral: Dp    = ", Dp)
    
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
    v = vx
    vr(1) = vx
    eta = etax
    maxres = maxwellian(vr)
    call bounce
    integrand =  pi/(2d0*n0)*(2*pi)**3/fsa*mph**2*c/(qi*iota*psi_pr)*&
         taub*mi**3*v**3*c/(qi*2d0*pi)*(mi*v**2/2d0)**2*abs(bounceavg(4))**2*&
         maxres(1)
!    integrand = (pi*mi*c/(2d0*qi))**2*pi/(iota*psi_pr*fsa)&
!         *mph**2*v**7*taub*abs(bounceavg(4))**2*maxres(1)*mi**3/n0
    ! TODO: factor of 2
  end function integrand
end module driftorbit
