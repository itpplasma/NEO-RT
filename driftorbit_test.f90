PROGRAM main
  use driftorbit
  implicit none

  call init
  
  !call test_magfie
  call test_bounce
contains

  subroutine test_magfie
    integer, parameter :: nth = 50
    integer :: k
    real(8) :: thmin, thmax
    real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)

    print *, "test_magfie"
    open(unit=1, file='test_magfie.dat', recl=1024)

    thmin = -pi
    thmax = pi
    x(1) = 1d-1
    x(2) = 0d0
    x(3) = 0d0
    
    do k = 0, nth-1
       x(3) = thmin + k*(thmax-thmin)/(nth-1)
       call magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl )
       write(1,*) x(3), bmod
    end do
    
    close(unit=1)
  end subroutine test_magfie
  
  subroutine test_bounce
    real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)
    
    real(8) :: y0(nvar), t0, dt
    integer, parameter :: n = 100

    integer :: k, state
    integer :: iwork(1000)
    real(8) :: t
    real(8) :: rwork(1000), y(nvar)
    
    print *, "test_bounce"
    open(unit=1, file='test_bounce.dat', recl=1024)

    x(1) = 1d-1
    x(2) = 0d0
    x(3) = 0d0
    call magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl )

    v   = 1d8
    eta = (1d0-1d-2)/bmod
    
    y0(1) = 0d0
    y0(2) = vpar(bmod)
    y0(3) = 0d0

    print *, bmod
    
    y = y0    
    t = 0d0
    dt = 1d-6
    write(1, *) t, y(1)
    state = 1
    do k = 2,n
       call dlsode(timestep, nvar, y, t, t+dt, 1, 0.5d-5, 1.0d-20, 1, &
            state, 0, rwork, 1000, iwork, 1000, jac, 10)
       write(1, *) t, y(1)
    end do
    
    close(unit=1)
  end subroutine test_bounce
end PROGRAM main
