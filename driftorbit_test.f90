PROGRAM main
  use driftorbit
  use do_magfie_mod
  use rkf45
  implicit none

  call init_test

  !call test_magfie
  !call test_bounce
  !call test_torfreq
  !call test_resline
  !call test_flux
  !call test_driftorbit
  !call test_torfreq_pass
  !call test_machrange
  !call test_integral
  call test_machrange2
contains

  subroutine init_test
    !s = .16406d0   ! flux surface no. 10
    !s = 5.4688d-2  ! flux surface no. 3
    !s = 1.578d-5   ! eps=1e-3
    !s = 1.5765d-3  ! eps=1e-2
    s = .1547      ! eps=1e-1
    !s = .3
    !M_t = 1d-5   ! set Mach number M_t = Om_tE*R0/vth
    !M_t = -3.2d-2 ! set Mach number M_t = Om_tE*R0/vth
    !M_t = 0.2d-2 ! set Mach number M_t = Om_tE*R0/vth
    !M_t = 5.6d-5
    M_t = 2.8d-2
    n0 = 1d22     ! particle density
    vth = 1d8

    call init
    if (.not. nobdrift) then
       ! set thermal velocity so that ExB = reference toroidal drift
       vth = abs(2*M_t*qi*psi_pr/(mi*c*R0)) ! thermal velocity
    end if
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
    print *, "test_magfie: B0        = ", B0
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
       print *, sqrtg*bmod*hctrvr(3), psi_pr/q, psi_pr*s/q
    end do
    
    close(unit=9)
  end subroutine test_magfie
  
  subroutine test_bounce
    real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)
    
    real(8) :: dt
    integer, parameter :: n = 1000

    integer :: k, state
    !integer :: iwork(1000)
    real(8) :: t, y(nvar), yp(nvar), abserr, relerr
    !real(8) :: rwork(1000)

    mth = -115
    
    print *, "test_bounce"
    open(unit=9, file='test_bounce.dat', recl=1024)

    !s = .16406d0
    x(1) = s
    x(2) = 0d0
    x(3) = 0d0
    call do_magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl ) 
    
    v   = 1d8
    !eta = (1d0-1.86727d-1)/bmod
    !eta = etatp*(1+1d-2)
    eta = etatp*0.1
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
       relerr = 1d-9
       abserr = 1d-10
       call r8_rkf45 ( timestep, nvar, y, yp, t, t+dt, relerr, abserr, state )
       write(9, *) t, y(4), y(5)
    end do
    close(unit=9)
  end subroutine test_bounce
  
  subroutine test_torfreq
    integer, parameter :: n = 200
    integer :: k
    real(8) :: Omph, dOmphdv, dOmphdeta
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: a, b, delta

    mth = -3
    
    v = 3*vth
    
    call disp("test_torfreq: vth        = ", vth)
    call disp("test_torfreq: v/vth      = ", v/vth)
    call disp("test_torfreq: Om_tE      = ", Om_tE)
    call disp("test_torfreq: Om_tB_ref  = ", c*mi*vth**2/(2*qi*psi_pr))

    etamin = etatp*1d-7
    etamax = etatp*(1-1d-1)

    eta = etamax*(1d0-1d-7)
    call Om_th(Omth, dOmthdv, dOmthdeta)

    call disp("test_torfreq: Om_th_approx    = ", v/(q*R0*sqrt(2d0/eps)))
    call disp("test_torfreq: Om_th_deeptrap  = ", Omth)

    delta = 1d-9   ! smallest relative distance to etamin
    b = log(delta)
    a = 1d0/(n-1d0)*(log(etamax/etamin - 1d0) - b)
    
    open(unit=9, file='test_torfreq.dat', recl=1024)
    do k = 0, n-1
       eta = etamin + k/(n-1d0)*(etamax-etamin)
       !eta = etamin*(1d0 + exp(a*k+b))
       !print *, etamin, etamax, eta
       call Om_ph(Omph, dOmphdv, dOmphdeta)
       call Om_th(Omth, dOmthdv, dOmthdeta)
       write(9, *) eta, Omth, v*iota*sqrt(1-B0*eta)/R0, dOmthdeta,&
            -v*B0/(2*q*R0*sqrt(1-eta*B0))
       !write(9, *) log((eta-etamin)/(etamax-etamin)), mph*(Omph-Om_tE),&
       !     mph*Om_tE, mth*Omth, mph*Omph + mth*Omth
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
    integer, parameter :: n = 100
    integer :: k
    real(8) :: vmin, vmax
    real(8) :: etarest(2), etaresp(2)

    mth = -115
    M_t = 2.8d-2
    Om_tE = vth*M_t/R0

    vmin = 1d-6*vth
    vmax = 1d3*vth
    !vmax = 5d0*vth
    !vmin = 1d-10*vth
    !vmax = 3.5d1*vth
    !vmax = 4d-3*vth


    call find_vlim(vmin, vmax)
    call disp("test_resline: vmin/vth              = ", vmin/vth)
    call disp("test_resline: vmax/vth              = ", vmax/vth)
    
    !vmax = min(3.5d0*vth, vmax)
    !vmax = 2*vth
    
    open(unit=9, file='test_resline.dat', recl=1024)
    do k = 0, n-1
       v = vmin + k/(n-1d0)*(vmax-vmin)
       if (driftorbit_nroot(1, 1d-7*etatp, (1-1d-7)*etatp) > 0) then
          ! trapped resonance (passing)
          etaresp = driftorbit_root(1d-8*abs(Om_tE), 1d-7*etatp, (1-1d-7)*etatp)
       end if
       if (driftorbit_nroot(1, (1+1d-8)*etatp, (1-1d-8)*etadt) > 0) then
          ! passing resonance (trapped)
          etarest=driftorbit_root(1d-8*abs(Om_tE),(1+1d-8)*etatp,(1-1d-8)*etadt)
       end if
       write(9, *) v/vth, etaresp(1)/etatp, etarest(1)/etatp
    end do
    close(unit=9)
  end subroutine test_resline
  
  subroutine test_flux
    Integer, parameter :: n = 2000
    integer :: k
    real(8) :: vrange(n), fluxintp(n), fluxintt(n)
    real(8) :: vmin, vmax, dv
    real(8) :: rhol, Dp, Drp

    M_t = 0.0134
    Om_tE = vth*M_t/R0
    
    mth = -1

    fluxintp = 0d0
    fluxintt = 0d0
    
    vmin = 1d-6*vth
    vmax = 1d2*vth
    !vmax = 5d0*vth

    call find_vlim(vmin, vmax)

    !vmin = 2.8*vth
    vmax = min(vmax,5d0*vth)
    
    !vmax = .5*vth
    
    call disp("test_flux: vmin/vth        = ", vmin/vth)
    call disp("test_flux: vmax/vth        = ", vmax/vth)

    do k = 1, n
       vrange(k) = vmin + (k-1)*(vmax-vmin)/(n-1)
    end do

    dv = vrange(2) - vrange(1)
    
    open(unit=9, file='test_flux.dat', recl=1024)
    do k = 1, n
       v = vrange(k)
       if (driftorbit_nroot(1, 1d-7*etatp, (1-1d-7)*etatp) > 0) then
          ! passing resonance (transit)
          etamin = 1d-7*etatp
          etamax = (1-1d-7)*etatp
          fluxintp(k) = flux_integral(vrange(k))
       end if
       if (driftorbit_nroot(1, (1+1d-7)*etatp, (1-1d-7)*etadt) > 0) then
          ! trapped resonance (bounce)
          etamin = (1+1d-7)*etatp
          etamax = (1-1d-7)*etadt
          fluxintt(k) = flux_integral(vrange(k))
       end if
       write(9, *) vrange(k)/vth, fluxintp(k)*vth, fluxintt(k)*vth ! f(v)*vth
    end do
    close(unit=9)
    call disp("test_flux: D11/Dp passing  = ", sum(fluxintp)*dv)
    call disp("test_flux: D11/Dp trapped  = ", sum(fluxintt)*dv)
    call disp("test_flux: D11/Dp total    = ", (sum(fluxintt)+sum(fluxintp))*dv)


    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    rhol = vth*mi*c/(B0*qi)
    Drp = sqrt(pi)*mph*(1/iota)**2*(1/eps)**2*vth*rhol**2/(4*R0)
    
    call disp("test_flux: Drp/Dp           = ", Drp/Dp)
  end subroutine test_flux

  subroutine test_driftorbit
    integer, parameter :: n = 500
    integer :: k
    real(8) :: Omph, dOmphdv, dOmphdeta
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: vmin, vmax

    mth = -1
    
    vmin = 1d-6*vth
    vmax = 3.5d0*vth
    !vmax = 5d0*vth
    !vmin = 1d-10*vth
    !vmax = 3.5d1*vth
    !vmax = 4d-3*vth

    call find_vlim(vmin, vmax)
    call disp("test_driftorbit: vmin/vth        = ", vmin/vth)
    call disp("test_driftorbit: vmax/vth        = ", vmax/vth)
    v = 0.7*vth

    etamin = 1d-8*etatp
    !etamin = etatp*(1d0-1d-1)
    etamax = etadt*(1d0-1d-8)

    !etamin = etatp*(1d0-1d-7)
    !etamax = etatp*(1d0+1d-7)
    
    open(unit=9, file='test_driftorbit.dat', recl=1024)
    do k = 2, n-1
       eta = etamin + k/(n-1d0)*(etamax-etamin)
       call Om_ph(Omph, dOmphdv, dOmphdeta)
       call Om_th(Omth, dOmthdv, dOmthdeta)
       
       write(9, *) eta/etatp, mth*Omth, -mph*Om_tE
       !write(9, *) eta, OmtB/v
    end do
    close(unit=9)
  end subroutine test_driftorbit

   subroutine test_machrange     
    integer, parameter :: m=100, n=2000
    integer :: j, k
    real(8) :: vrange(n), fluxintp(n), fluxintt(n)
    real(8) :: mtrange(m), fluxresp(m), fluxrest(m)
    real(8) :: mtmin, mtmax, vmin, vmax, dv

    mth = -3

    mtmin = 0.1d-2
    mtmax = 5.6d-2
    do j = 1, n
      mtrange(j) = mtmin + (j-1)*(mtmax-mtmin)/(m-1)
    end do
    
    
    open(unit=9, file='test_machrange.dat', recl=1024)
    do j = 1, m
       M_t = mtrange(j)
       Om_tE = vth*M_t/R0
       
       fluxintp = 0d0
       fluxintt = 0d0

       vmin = 1d-6*vth
       vmax = 1d2*vth
       !vmax = 5d0*vth

       call find_vlim(vmin, vmax)

       !vmin = 2.8*vth
       vmax = min(vmax,5d0*vth)

       !vmax = .5*vth

       call disp("test_machrange: Mt       = ", M_t)
       call disp("test_machrange: vmin/vth = ", vmin/vth)
       call disp("test_machrange: vmax/vth = ", vmax/vth)

       do k = 1, n
          vrange(k) = vmin + (k-1)*(vmax-vmin)/(n-1)
       end do
       dv = vrange(2) - vrange(1)
       do k = 1, n
          v = vrange(k)
          if (driftorbit_nroot(1, 1d-7*etatp, (1-1d-7)*etatp) > 0) then
             ! passing resonance (transit)
             etamin = 1d-7*etatp
             etamax = (1-1d-7)*etatp
             fluxintp(k) = flux_integral(vrange(k))
          end if
          if (driftorbit_nroot(1, (1+1d-8)*etatp, (1-1d-8)*etadt) > 0) then
             ! trapped resonance (bounce)
             etamin = (1+1d-8)*etatp
             etamax = (1-1d-8)*etadt
             fluxintt(k) = flux_integral(vrange(k))
          end if
          fluxresp(j) = sum(fluxintp*dv)
          fluxrest(j) = sum(fluxintt*dv)
       end do
       call disp("test_machrange: Dp = ", fluxrest(j) + fluxresp(j))
       write(9, *) mtrange(j), fluxresp(j), fluxrest(j)
    end do
    close(unit=9)
  end subroutine test_machrange

  
   subroutine test_integral   
    integer, parameter :: m = 30
    integer :: j
    real(8) :: fluxresp(m), fluxrest(m)
    real(8) :: vminp, vmaxp, vmint, vmaxt
    real(8) :: rhol, Dp, Drp

    M_t = 0.0136
    Om_tE = vth*M_t/R0

    fluxrest = 0d0
    fluxresp = 0d0
    
    call disp("test_integral: Mach num Mt   = ", M_t)
    call disp("test_integral: Poloidal m0   = ", 1d0*m0)
    call disp("test_integral: Toroidal n0   = ", 1d0*mph)
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    rhol = vth*mi*c/(B0*qi)
    Drp = sqrt(pi)*mph*(1/iota)**2*(1/eps)**2*vth*rhol**2/(4*R0)
    
    call disp("test_integral: Drp/Dp           = ", Drp/Dp)
    
    open(unit=9, file='test_integral.dat', recl=1024)
    do j = 1, m
       mth = -j
       
       vminp = 1d-6*vth
       vmaxp = 1d2*vth
       vmint = 1d-6*vth
       vmaxt = 1d2*vth

       call find_vlim_p(vminp, vmaxp)
       call find_vlim_t(vmint, vmaxt)

       ! passing resonance (passing)
       if ((.not. nopassing) .and. (vmaxp > 0)) then
          etamin = 1d-7*etatp
          etamax = (1-1d-7)*etatp
          fluxresp(j) = flux_integral2(vminp, vmaxp)
       end if
          
       ! trapped resonance (trapped)
       etamin = (1+1d-8)*etatp
       etamax = (1-1d-8)*etadt
       fluxrest(j) = flux_integral2(vmint, vmaxt)
          
       !call disp("test_integral: D11/Dp pass  = ", fluxresp(j))
       !call disp("test_integral: D11/Dp trap  = ", fluxrest(j))
       !call disp("test_integral: D11/Dp total = ", fluxrest(j) + fluxresp(j))
       print *, "test_integral: mth = ", mth, " D11/Dp = ", fluxresp(j),&
            fluxrest(j), fluxresp(j) + fluxrest(j)
       write(9, *) j, fluxresp(j), fluxrest(j), fluxresp(j) + fluxrest(j)
    end do
    write(9, *) Drp/Dp, 0d0, 0d0, 0d0
    close(unit=9)
    call disp("test_integral: Dp total = ", sum(fluxrest + fluxresp))
  end subroutine test_integral

  subroutine test_machrange2  
    integer, parameter :: m = 30, n = 101
    integer :: j, k
    real(8) :: fluxresp(m), fluxrest(m)
    real(8) :: vminp, vmaxp, vmint, vmaxt
    real(8) :: Mtmin, Mtmax

    Mtmin = 0d0
    Mtmax = 1d-1
            
    open(unit=9, file='test_machrange2.dat', recl=1024)
    do k = 1, n
       if (k == 1) then
          M_t = 1d-5
       else
          M_t = mtmin + (k-1)*(mtmax-mtmin)/(n-1)
       end if
       Om_tE = vth*M_t/R0
       fluxrest = 0d0
       fluxresp = 0d0
       do j = 1, m
          mth = -j

          vminp = 1d-6*vth
          vmaxp = 1d2*vth
          vmint = 1d-6*vth
          vmaxt = 1d2*vth

          call find_vlim_p(vminp, vmaxp)
          call find_vlim_t(vmint, vmaxt)

          ! passing resonance (passing)
          if ((.not. nopassing) .and. (vmaxp > 0)) then
             etamin = 1d-7*etatp
             etamax = (1-1d-7)*etatp
             fluxresp(j) = flux_integral2(vminp, vmaxp)
          end if

          ! trapped resonance (trapped)
          etamin = (1+1d-8)*etatp
          etamax = (1-1d-8)*etadt
          fluxrest(j) = flux_integral2(vmint, vmaxt)

          print *, "test_flux: Mt = ", M_t, ", mth = ", mth, ", D11/Dp = ",&
               fluxresp(j), fluxrest(j), fluxresp(j) + fluxrest(j)
       end do
       write(9, *) M_t, sum(fluxresp), sum(fluxrest), sum(fluxresp+fluxrest)
       flush(9)
    end do
    close(unit=9)
  end subroutine test_machrange2
end PROGRAM main
