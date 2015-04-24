! Resonant regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, 2015

module driftorbit
  use do_magfie_mod
  
  implicit none
  save

  integer, parameter :: nvar = 3
  
  real(8) :: v, eta
  real(8) :: taub, bounceavg(nvar)

contains

  subroutine init
    call do_magfie_init
  end subroutine init
  
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

  subroutine timestep(neq, t, y, ydot)
    !
    !  Timestep function for orbit integration.
    !  Includes poloidal angle theta and parallel velocity.
    !  More integrands may be added starting from y(3)
    !
    
    integer :: neq
    real(8) :: t, y(neq)
    real(8), intent(out) :: ydot(neq)
    real(8) :: bmod, sqrtg
    real(8) :: x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
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
    ydot(3) = Om_tB_v
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

    n = 100
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
       print *, "ERROR: findroot did not converge in 100 iterations"
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

  ! Shaing2009-035009 - (8) for bounce averaged toroidal drift 
  function Omph_shaing()
    real(8) :: Omph_shaing, kappa, Eell, Kell         
    
    ! Elliptic integrals for Shaing formula
    !kappa = sqrt((E()-qi*Phi_E(fsind)-J_perp()*om_c(0.0d0))/(2*J_perp()*&
    !     om_c(pi/2)*eps(fsind)))
    kappa = sqrt((mi*v**2/2d0-J_perp()*om_c(0d0))/(2*J_perp()*&
         om_c(pi/2)*eps))
    !kappa = sqrt(0.83)               ! test case so that 2*E/K = 1
    Kell = ellf(pi/2d0, kappa)
    Eell = elle(pi/2d0, kappa)

    ! Shaing2009-035009 - (8) for bounce averaged toroidal drift 
    Omph_shaing = Om_tE +&
         J_perp()*B0mn(fsind,1)/mi*depsdr(fsind)*(2*Eell/Kell-1)
  end function Omph_shaing
end module driftorbit
