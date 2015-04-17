! Resonant regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, 2015

module driftorbit
  use common
  use neo_exchange, only: nper,b_min,b_max, &
       theta_bmin,theta_bmax,phi_bmin,phi_bmax
  use magfie_mod, only: magfie, stevvo, magfie_deallocate
  use neo_magfie_mod, only: magfie_spline, magfie_sarray
  use nrtype, only: twopi
  
  implicit none
  save

  integer, parameter :: nvar = 3
  
  real(8) :: s
  real(8) :: v, eta
  

contains

  subroutine init
    
    real(8) :: thetab, phib, boozer_s
    real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)
    
    magfie_spline = 1
    
    allocate(magfie_sarray(1))
    magfie_sarray = 1d-1

    x(1) = 1d-1
    x(2:3) = 0.0d0
    call magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl )
  end subroutine init
  
  function Jperp()
    real(8) :: Jperp
    Jperp = mi**2*c/(2d0*qi)*eta*v**2
  end function Jperp
  
  function vpar(bmod)
    !   parallel velocity
    real(8) :: bmod, vpar
    vpar = v*sqrt(1-eta*bmod)
    if (isnan(vpar)) then
       vpar = 0d0
    end if
  end function vpar
  
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
    real(8) :: x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)

    x(1) = s
    x(2) = 0d0
    x(3) = y(1)
    call magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl )

    ydot(1) = y(2)*hctrvr(3)                           ! theta
    ydot(2) = -qi*Jperp()/(mi**2*c)*hctrvr(3)*bder(3)  ! v_par
    ydot(3) = 1d0                                      ! t

    print *, Jperp()
  end subroutine timestep

  

  function findroot(y0, t0, dt, tol)
    !
    !  Finds the root of an orbit after the first turn
    !
    ! TODO: might use Brent's method
    
    real(8) :: y0(nvar), t0, dt, tol, findroot(nvar+1)
    integer :: n

    integer :: k, state, rootstate
    integer :: iwork(1000)
    real(8) :: ti, told
    real(8) :: rwork(1000), y(nvar), yold(nvar)

    n = 100
    rootstate = -1
    
    y = y0
    yold = y0
    ti = t0
    state = 1
    do k = 2,n
       yold = y
       told = ti
       call dlsode(timestep, nvar, y, ti, ti+dt, 1, 1d-6, 1.0d-20, 1, &
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
end module driftorbit
