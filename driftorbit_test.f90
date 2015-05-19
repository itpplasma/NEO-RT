PROGRAM main
  use driftorbit
  use do_magfie_mod
  use rkf45
  implicit none

  call init_test

  call test_magfie
  !call test_bounce
  call test_torfreq
  !call test_resline
  call test_flux
  !call test_driftorbit
contains

  subroutine init_test
    !s = .16406d0   ! flux surface no. 10
    !s = 5.4688d-2  ! flux surface no. 3
    !s = 1.58d-5     ! eps=1e-3
    !s = 1.5765d-3  ! eps=1e-2
    s = .1547      ! eps=1e-1
    !s = .3
    M_t = 1d-5   ! set Mach number M_t = Om_tE*R0/vth
    n0 = 1d22    ! particle density
    vth = 1d0

    call init
      
    ! set thermal velocity so that ExB = reference toroidal drift
    vth = abs(2*M_t*qi*psi_pr/(mi*c*R0)) ! thermal velocity
    Om_tE = vth*M_t/R0                   ! toroidal ExB drift frequency
  end subroutine init_test

  subroutine test_magfie
    integer, parameter :: nth = 50
    integer :: k
    real(8) :: thmin, thmax
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

    print *, "test_magfie"
    open(unit=9, file='test_magfie.dat', recl=1024)

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
    
    close(unit=9)
  end subroutine test_magfie
  
  subroutine test_bounce
    real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)
    
    real(8) :: dt
    integer, parameter :: n = 100

    integer :: k, state
    !integer :: iwork(1000)
    real(8) :: t, y(nvar)
    !real(8) :: rwork(1000)
    
    print *, "test_bounce"
    open(unit=9, file='test_bounce.dat', recl=1024)

    !s = .16406d0
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
    write(9, *) t, y(1)
    state = 1
    do k = 2,n
       !call dlsode(timestep, nvar, y, t, t+dt, 1, 1d-15, 1d-15, 1, &
       !     state, 0, rwork, 1000, iwork, 1000, jac, 10)
       !write(9, *) t, y(1)
    end do
    close(unit=9)
  end subroutine test_bounce
  
  subroutine test_torfreq
    integer, parameter :: n = 200
    integer :: k
    real(8) :: etamin, etamax
    real(8) :: Omph, dOmphdv, dOmphdeta
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: a, b, delta
    
    v = 1.75*vth
    
    call disp("test_torfreq: vth        = ", vth)
    call disp("test_torfreq: v/vth      = ", v/vth)
    call disp("test_torfreq: Om_tE      = ", Om_tE)
    call disp("test_torfreq: Om_tB_ref  = ", c*mi*vth**2/(2*qi*psi_pr))

    etamin = etatp()
    etamax = etadt()
    
    !etamin = etatp()*(1d0+1d-9)
    !etamax = etadt()*(1d0-1d-15)
    !etamax = etamin*(1d0+1d-8)

    delta = 1d-9   ! smallest relative distance to etamin
    b = log(delta)
    a = 1d0/(n-1d0)*(log(etamax/etamin - 1d0) - b)
    
    open(unit=9, file='test_torfreq.dat', recl=1024)
    do k = 0, n-1
       !eta = etamin + k/(n-1d0)*(etamax-etamin)
       eta = etamin*(1d0 + exp(a*k+b))
       !print *, etamin, etamax, eta
       call Om_ph(Omph, dOmphdv, dOmphdeta)
       call Om_th(Omth, dOmthdv, dOmthdeta)
       write(9, *) log((eta-etamin)/(etamax-etamin)), mph*(Omph-Om_tE),&
            mph*Om_tE, mth*Omth, mph*Omph + mth*Omth
       !write(9, *) log(eta-etamin), log(mph*dOmphdeta), log(mth*dOmthdeta),&
       !     log(mph*dOmphdeta + mth*dOmthdeta)
    end do
    close(unit=9)
  end subroutine test_torfreq
  
  subroutine test_Om_spline
    ! TODO: write test routine for splined canonical frequencies over eta
    ! Om_tB/v^2 and Omth/v are independent of v and splined in init_Om_tB_spl.
  end subroutine test_Om_spline
  
  subroutine test_resline
    integer, parameter :: n = 50
    integer :: k
    real(8) :: vmin, vmax
    real(8) :: etares(2)

    vmin = 1d-10*vth
    vmax = 3.5d1*vth

    vmin = find_vmin(vmin, vmax)
    
    call disp("test_resline: vmin/vth        = ", vmin/vth)
    call disp("test_resline: vmax/vth        = ", vmax/vth)

    
    open(unit=9, file='test_resline.dat', recl=1024)
    do k = 0, n-1
       v = vmin + k/(n-1d0)*(vmax-vmin)
       etares = driftorbit_root(1d-8*abs(Om_tE))
       write(9, *) v, log(etares(1)-etatp())
    end do
    close(unit=9)
  end subroutine test_resline
  
  subroutine test_flux
    integer, parameter :: n = 300
    integer :: k
    real(8) :: vrange(n), fluxint(n)
    real(8) :: vmin, vmax, dv

    vmin = 1d-10*vth
    vmax = 6d0*vth

    vmin = find_vmin(vmin, vmax)
    
    call disp("test_flux: vmin/vth        = ", vmin/vth)
    call disp("test_flux: vmax/vth        = ", vmax/vth)

    do k = 1, n
       vrange(k) = vmin + (k-1)*(vmax-vmin)/(n-1)
    end do

    dv = vrange(2) - vrange(1)
    fluxint = flux_integral(vrange)
    
    open(unit=9, file='test_flux.dat', recl=1024)
    do k = 1, n
       write(9, *) vrange(k)/vth, fluxint(k)*vth ! f(v)*vth
    end do
    close(unit=9)
    call disp("test_flux: D11/Dp          = ", sum(fluxint)*dv)
  end subroutine test_flux

  subroutine test_driftorbit
    integer, parameter :: n = 300
    integer :: k
    real(8) :: etamin, etamax
    real(8) :: OmtB, dOmtBdv, dOmtBdeta
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: pa, qa
    real(8) :: aa, ba, ca
    real(8) :: v1, v2
    
    v = 1*vth

    etamin = etatp()*(1d0+1d-10)
    etamax = etadt()*(1d0-1d-10)
    
    open(unit=9, file='test_driftorbit.dat', recl=1024)
    do k = 2, n-1
       eta = etamin + k/(n-1d0)*(etamax-etamin)
       call Om_tB(OmtB, dOmtBdv, dOmtBdeta)
       call Om_th(Omth, dOmthdv, dOmthdeta)

       aa = mph*OmtB/v**2
       ba = mth*Omth/v
       ca = mph*Om_tE
       
       pa = v*mth*Omth/(mph*OmtB)
       qa = v**2*Om_tE/OmtB

       v1 = (-ba + sqrt(ba**2 - 4*aa*ca))/(2*aa)
       v2 = (-ba - sqrt(ba**2 - 4*aa*ca))/(2*aa)
       
       write(9, *) eta, pa, qa, v1/vth, v2/vth, ba**2 - 4*aa*ca
       !write(9, *) eta, OmtB/v
    end do
    close(unit=9)
  end subroutine test_driftorbit
end PROGRAM main
