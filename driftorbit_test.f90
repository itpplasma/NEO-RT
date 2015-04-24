PROGRAM main
  use driftorbit
  use do_magfie_mod
  implicit none

  call init_test

  call test_magfie
  call test_bounce
contains

  subroutine init_test
    s = .16406d0
    v = 1d8
    call init
  end subroutine init_test

  subroutine test_magfie
    integer, parameter :: nth = 50
    integer :: k
    real(8) :: thmin, thmax
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

    print *, "test_magfie"
    open(unit=1, file='test_magfie.dat', recl=1024)

    thmin = -pi
    thmax = pi
    x(1) = s
    x(2) = 0d0
    x(3) = 0d0
    
    print *, "test_magfie: R0        = ", R0
    print *, "test_magfie: eps       = ", eps
    print *, "test_magfie: psi_pr    = ", psi_pr
    print *, "test_magfie: Bthcov    = ", Bthcov
    print *, "test_magfie: Bphcov    = ", Bphcov
    print *, "test_magfie: dBthcovds = ", dBthcovds
    print *, "test_magfie: dBphcovds = ", dBphcovds
    print *, "test_magfie: q         = ", q
    print *, "test_magfie: iota      = ", iota
    
    do k = 0, nth-1
       x(3) = thmin + k*(thmax-thmin)/(nth-1)
       call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
       write(1,*) x(3), bmod, sqrtg, hder(1), hder(2), hder(3), hcovar(1),&
            hcovar(2), hcovar(3), hctrvr(1), hctrvr(2), hctrvr(3),&
            boozer_curr_pol_hat_s
    end do
    
    close(unit=1)
  end subroutine test_magfie
  
  subroutine test_bounce
    real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)
    
    real(8) :: dt
    integer, parameter :: n = 100

    integer :: k, state
    integer :: iwork(1000)
    real(8) :: t
    real(8) :: rwork(1000), y(nvar)
    
    print *, "test_bounce"
    open(unit=1, file='test_bounce.dat', recl=1024)

    s = .16406d0
    x(1) = s
    x(2) = 0d0
    x(3) = 0d0
    call do_magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl ) 
    
    v   = 1d8
    eta = (1d0-1.86727d-1)/bmod
    call bounce

    print *, taub
    print *, bmod
    
    y = 0d0
    y(2) = vpar(bmod)
    t = 0d0
    dt = taub/(n-1)
    print *, "taub = ", taub
    write(1, *) t, y(1)
    state = 1
    do k = 2,n
       call dlsode(timestep, nvar, y, t, t+dt, 1, 1d-15, 1d-15, 1, &
            state, 0, rwork, 1000, iwork, 1000, jac, 10)
       write(1, *) t, y(1)
    end do
    
    close(unit=1)
  end subroutine test_bounce

  
  subroutine test_torfreq    
    print *, "test_torfreq"
    open(unit=1, file='test_torfreq.dat', recl=1024)
    ! TODO
    close(unit=1)
  end subroutine test_torfreq
end PROGRAM main
