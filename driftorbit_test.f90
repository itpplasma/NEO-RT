program main
  use driftorbit
  use lineint
  use do_magfie_mod
  implicit none

  integer :: Mtnum, mthnum
  real(8) :: Mtmin, Mtmax
  logical :: odeint
  character(len=1024) :: controlfile
  
  call get_command_argument(1, controlfile)

  call read_control
  call init_test
  call test_magfie

  if (Mtnum < 1) return
  
  !call test_bounce
  !call test_boundaries
  !call test_torfreq
  !call test_resline
  !call test_flux
  !call test_driftorbit
  !call test_torfreq_pass
  !call test_machrange
  !call test_boundaries
  !call test_Hmn
  call test_integral
  !call test_machrange2
contains

  subroutine read_control
    character(1)          :: dummy
    real(8)               :: qs, ms

    open(unit=9,file=controlfile,status='old',form='formatted')
    read (9,*) dummy
    read (9,*) dummy
    read (9,*) dummy
    read (9,*) s
    read (9,*) M_t
    read (9,*) qs
    read (9,*) ms
    read (9,*) n0
    read (9,*) vth
    read (9,*) epsmn
    read (9,*) m0
    read (9,*) mph   
    read (9,*) mth   
    read (9,*) mthnum   
    read (9,*) Mtmin   
    read (9,*) Mtmax  
    read (9,*) Mtnum    
    read (9,*) supban
    read (9,*) magdrift
    read (9,*) nopassing  
    read (9,*) calcflux
    read (9,*) noshear
    read (9,*) pertfile
    read (9,*) odeint

    qi = qs*qe
    mi = ms*mu
  end subroutine read_control

  subroutine init_test

    call init
    !if ((abs(M_t) > 1d-12) .and. (.not. nobdrift)) then
    ! set thermal velocity so that ExB = reference toroidal drift
       !print *, "Set vth from: ", vth
       !vth = abs(2*M_t*qi*psi_pr/(mi*c*R0)) ! thermal velocity
       !print *, "Set vth to: ", vth
    !end if
    Om_tE = vth*M_t/R0                   ! toroidal ExB drift frequency

    etamin = (1+epst)*etatp
    etamax = (1-epst)*etadt
  end subroutine init_test

  subroutine test_magfie
    integer, parameter :: nth = 50
    integer :: k
    real(8) :: thmin, thmax
    real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
    real(8) :: Drp
    complex(8) :: bn

    Drp = 4*mph*q/(eps**2*sqrt(pi));

    open(unit=9, file='test_magfie_param.dat', recl=1024)

    thmin = -pi
    thmax = pi
    x(1) = s
    x(2) = 0d0
    x(3) = 0d0

    write(9,*) "-------------------------"
    write(9,*) "test_magfie: R0        = ", R0
    write(9,*) "test_magfie: a         = ", a
    write(9,*) "test_magfie: eps       = ", eps
    write(9,*) "test_magfie: A         = ", 1/eps
    write(9,*) "test_magfie: psi_pr    = ", psi_pr
    write(9,*) "test_magfie: B0        = ", B0
    !write(9,*) "test_magfie: B0h       = ", B0h
    !write(9,*) "test_magfie: B00       = ", B00
    write(9,*) "test_magfie: Bthcov    = ", Bthcov
    write(9,*) "test_magfie: Bphcov    = ", Bphcov
    write(9,*) "test_magfie: dBthcovds = ", dBthcovds
    write(9,*) "test_magfie: dBphcovds = ", dBphcovds
    write(9,*) "test_magfie: q         = ", q
    write(9,*) "test_magfie: iota      = ", iota
    write(9,*) "test_magfie: M_t       = ", M_t
    write(9,*) "test_magfie: Om_tE     = ", Om_tE
    write(9,*) "test_magfie: vth       = ", vth
    write(9,*) "test_magfie: T [eV]    = ", mi/2*vth**2/eV
    write(9,*) "test_magfie: m0        = ", 1d0*m0
    write(9,*) "test_magfie: n0        = ", 1d0*mph
    write(9,*) "test_magfie: Drp       = ", Drp
    write(9,*) "-------------------------"
    write(9,*) "test_magfie: pertfile  = ", pertfile
    write(9,*) "-------------------------"
    

    close(unit=9)
    
    open(unit=9, file='test_magfie.dat', recl=1024)
    do k = 0, nth-1
       x(3) = thmin + k*(thmax-thmin)/(nth-1)
       call do_magfie( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
       if (pertfile) then
          call neo_magfie_pert_amp( x, bn )
          bn = bn/bmod
       else
          bn = epsmn*exp(imun*m0*x(3))
       end if
       write(9,*) x(3), bmod, sqrtg, hder(1), hder(2), hder(3), hcovar(1),&
            hcovar(2), hcovar(3), hctrvr(1), hctrvr(2), hctrvr(3),&  ! 8
            real(bn), aimag(bn), real(epsmn*exp(imun*m0*x(3))),& !13
            aimag(epsmn*exp(imun*m0*x(3))) !16
    end do

    close(unit=9)
  end subroutine test_magfie

  subroutine test_bounce
    use dvode_f90_m2
    real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)

    real(8) :: dt
    integer, parameter :: n = 5000

    integer :: k, state
    !integer :: iwork(1000)
    real(8) :: t, y(nvar), y2(nvar), y2dot(nvar)
    !real(8) :: rwork(1000)

    real(8) :: atol(nvar), rtol, tout
    integer :: neq, itask, istate
    type (vode_opts) :: options

    neq = nvar
    t = 0.0d0
    rtol = 1d-9
    atol = 1d-10
    itask = 1
    istate = 1
    options = set_normal_opts(abserr_vector=atol, relerr=rtol)

    !mth = -1

    print *, "test_bounce"
    open(unit=9, file='test_bounce.dat', recl=1024)

    !s = .16406d0
    x(1) = s
    x(2) = 0d0
    x(3) = 0d0
    call do_magfie( x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl ) 

    v   = vth
    !eta = (1d0-1.86727d-1)/bmod
    !eta = etatp*(1+1d-2)
    !eta = etatp*(1+1d-7)
    eta = etadt*(1-1d-2)
    call bounce
    call disp("bounce: Om_tB = ", vth**2*bounceavg(3))
    call disp("bounce: taub  = ", taub)
    call disp("bounce: bmod  = ", bmod)


    taub = taub
    y = 0d0
    y(2) = vpar(bmod)
    y2 = 0d0
    y2(2) = vpar(bmod)
    t = 0d0
    dt = taub/(n-1)
    print *, "taub_est = ", 2.0*pi/(vperp(bmod)*iota/R0*sqrt(eps/2d0))
    print *, "taub     = ", taub
    write(9, *) t, y(1)
    state = 1
    do k = 2,n
       !relerr = 1d-9
       !abserr = 1d-10
       !call r8_rkf45 ( timestep, nvar, y, yp, t, t+dt, relerr, abserr, state)
       call timestep(t, y2, y2dot)
       y2 = y2 + y2dot*dt
       tout = t+dt
       call dvode_f90(timestep2, neq, y, t, tout, itask, istate, options)
       write(9, *) t/taub, y(1)/taub, y(3)*vth**2/taub,&
            y2(1)/taub, y2(3)*vth**2/taub
    end do
    close(unit=9)
  end subroutine test_bounce

  subroutine test_torfreq
    integer, parameter :: n = 1000
    integer :: k
    real(8) :: Omph, dOmphdv, dOmphdeta
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: OmtB, dOmtBdv, dOmtBdeta
    real(8) :: a, b, delta

    !v = vth*1.7571463448202738

    !mth = 2
    v = .25*vth

    call disp("test_torfreq: vth        = ", vth)
    call disp("test_torfreq: v/vth      = ", v/vth)
    call disp("test_torfreq: mph        = ", 1d0*mph)
    call disp("test_torfreq: mth        = ", 1d0*mth)
    call disp("test_torfreq: Om_tE      = ", Om_tE)
    call disp("test_torfreq: Om_tB_ref  = ", c*mi*vth**2/(2*qi*psi_pr))
    call bounce
    call disp("test_torfreq: Om_tB_ba   = ", vth**2*bounceavg(3))

    !etamin = etatp*(1-10*epst_spl)
    !etamax = etatp*(1+10*epst_spl)
    !etamin = etatp*(1+epst)
    !etamax = etatp*(1+10*epst_spl)
    !etamin = etatp*(1-10*epst_spl)
    !etamax = etatp*(1-epst)
    etamin = etatp*(1+epst)
    etamax = etadt*(1-epst)
    !etamin = etatp*(1-2*epst)
    !etamax = etatp*(1+2*epst)
    call disp("test_torfreq: etamin = ", etamin/etadt)
    call disp("test_torfreq: etamin = ", etamax/etadt)
    !etamax = etadt*(1-1d-7) 

    !eta = etamax*(1d0-1d-7)
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
       call Om_tB(OmtB, dOmtBdv, dOmtBdeta)
       write(9, *) (eta-etatp)/(etadt-etatp),&
            Om_tE, OmtB, dOmtbdv, dOmtbdeta,&
            Omth, dOmthdv, dOmthdeta,&
            Omph, dOmphdv, dOmphdeta, mph*Omph+mth*Omth
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
    real(8) :: roots(nlev, 3)
    integer :: nroots, kr

    vmin = 1d-6*vth
    vmax = 1d0*vth
    !vmax = 5d0*vth
    !vmin = 1d-10*vth
    !vmax = 3.5d1*vth
    !vmax = 4d-3*vth

    !call find_vlim(vmin, vmax)
    !vmax = (1+1d-3)*vmin
    call disp("test_resline: vmin/vth              = ", vmin/vth)
    call disp("test_resline: vmax/vth              = ", vmax/vth)

    !vmax = min(3.5d0*vth, vmax)
    !vmax = 2*vth

    etaresp = etatp
    etarest = etatp

    open(unit=9, file='test_resline.dat', recl=1024)
    do k = 0, n-1
       v = vmin + k/(n-1d0)*(vmax-vmin)
       
       if (.not. nopassing) then
          ! resonance (passing)

          call driftorbit_coarse(etatp*epsp, etatp*(1-epsp), roots, nroots)
          do kr = 1,nroots
             etaresp = driftorbit_root2(1d-8*abs(Om_tE), roots(kr,1), roots(kr,2))
             !write(9, *) v/vth, kr, roots(kr,1), roots(kr,2), etaresp(1)
          end do
       end if
       
       ! resonance (trapped)
       call driftorbit_coarse(etatp*(1+epst), etadt*(1-epst), roots, nroots)
       !print *, "trapped roots: ", v/vth, nroots
       do kr = 1,nroots
          etarest = driftorbit_root2(1d-8*abs(Om_tE), roots(kr,1), roots(kr,2))
          write(9, *) v/vth, kr, roots(kr,1), roots(kr,2), etarest(1)
       end do
       !etarest=driftorbit_root(1d-8*abs(Om_tE),(1+epst)*etatp,(1-epst)*etadt)
       
       !write(9, *) v/vth, (etaresp(1)-etatp)/(etadt-etatp),&
       !     (etarest(1)-etatp)/(etadt-etatp)
    end do
    close(unit=9)
  end subroutine test_resline

  subroutine test_driftorbit
    integer, parameter :: n = 500
    integer :: k
    real(8) :: Omph, dOmphdv, dOmphdeta
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: vmin, vmax

    !mth = -1

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

  subroutine test_machrange2  
    integer :: j, k
    real(8) :: fluxresp(2), fluxrest(2), fluxp(2), fluxt(2)
    real(8) :: vminp, vmaxp, vmint, vmaxt
    integer :: mthmin, mthmax
    character(1024) :: tmp, fn1, fn2

    fluxp = 0d0
    fluxt = 0d0

    write(tmp, *) s
    write(fn1, *) 'driftorbit_'//trim(adjustl(tmp))//'.out'
    write(fn2, *) 'driftorbit_integral_'//trim(adjustl(tmp))//'.out'

    open(unit=10, file=trim(adjustl(fn2)), recl=1024)
    do k = 1, Mtnum
       if (Mtnum > 1) M_t = Mtmin + (k-1)*(Mtmax-Mtmin)/(Mtnum-1)

       Om_tE = vth*M_t/R0
       fluxp = 0d0
       fluxt = 0d0

       if (M_t > 0) then
          mthmin = 1
          mthmax = ceiling(mph*q + mthnum)
       else
          mthmin = ceiling(-mph*q)
          mthmax = mthnum
       end if
       
       do j = mthmin, mthmax
          fluxrest = 0d0
          fluxresp = 0d0
          
          if (M_t < 0) then
             mth = j
          else
             mth = -j
          end if

          vminp = 1d-6*vth
          vmaxp = 5d0*vth
          vmint = 1d-6*vth
          vmaxt = 5d0*vth

          ! superbanana resonance
          if (supban) then
             mth = 0
             vmint = 0.1*vth
             vmaxt = 5*vth
             etamin = (1+100*epst)*etatp
             etamax = etatp + (1-100*epst)*(etadt-etatp)
             call find_vlim(vmint, vmaxt)
             print *, "VLIM: ", vmint/vth, vmaxt/vth
             vmint = 1.1111*vmint
             vmaxt = 0.9999*vmaxt
             etamin = (1+10*epst)*etatp
             etamax = (1-10*epst)*etadt
             if (odeint) then
                fluxrest = flux_integral_ode(vmint, vmaxt)
             else
                fluxrest = flux_integral(vmint, vmaxt)
             end if
             fluxt = fluxt + fluxrest

          else
             ! passing resonance (passing)
             if ((.not. nopassing) .and. &
                  ((M_t < 0 .and. mth > -mph*q) .or.&
                  (M_t > 0 .and. mth < -mph*q)) .and.&
                  (abs(mth+mph*q) <= mthnum)) then
                !print *, "BEFORE FIND_VLIM_P", mth
                !call find_vlim_p(vminp, vmaxp)
                !print *, "AFTER FIND_VLIM_P", vminp/vth, vmaxp/vth
                etamin = epsp*etatp
                etamax = (1-epsp)*etatp
                if (odeint) then
                   fluxresp = flux_integral_ode(vminp, vmaxp)
                else
                   fluxresp = flux_integral(vminp, vmaxp)
                end if
                fluxp = fluxp + fluxresp
             end if

             ! trapped resonance (trapped)
             if (M_t*mth < 0 .and. abs(mth) <= mthnum) then
                !call find_vlim_t(vmint, vmaxt)
                etamin = (1+epst)*etatp
                etamax = (1-epst)*etadt
                if (odeint) then
                   fluxrest = flux_integral_ode(vmint, vmaxt)
                else
                   fluxrest = flux_integral(vmint, vmaxt)
                end if
                fluxt = fluxt + fluxrest
             end if
          end if
                 
          print *, ''
          print *, "test_flux: Mt = ", M_t, ", mth = ", mth
          write(*,'(3ES12.2,2F12.2)') fluxresp(1), fluxrest(1),&
               fluxresp(1) + fluxrest(1), vminp/vth, vmint/vth
          write(*,'(3ES12.2,2F12.2)') fluxresp(2), fluxrest(2),&
               fluxresp(2) + fluxrest(2), vmaxp/vth, vmaxt/vth

          if (k == 1 .and. j == 1) then
             open(unit=10, file=trim(adjustl(fn2)), recl=1024)
          end if
          write(10, *) M_t, mth, fluxresp(1), fluxrest(1),&
               fluxresp(1) + fluxrest(1), fluxresp(2), fluxrest(2),&
               fluxresp(2) + fluxrest(2), vminp/vth, vmaxp/vth,&
               vmint/vth, vmaxt/vth
          flush(10)
          
          if (supban) then
             exit
          end if
       end do
       if (k == 1) then
          open(unit=9, file=trim(adjustl(fn1)), recl=1024)
       end if
       write(9, *) M_t, fluxp(1), fluxt(1),&
            fluxp(1) + fluxt(1),&
            fluxp(2), fluxt(2),&
            fluxp(2) + fluxt(2)
       flush(9)
    end do
    close(unit=9)
    close(unit=10)
  end subroutine test_machrange2

  subroutine test_integral
    real(8) :: eta_res(2)
    real(8) :: D11, D12, Dp, vold, dsdreff, D11sum, D12sum
    integer :: ku, nu
    real(8) :: xs, kappa2s
    real(8) :: roots(nlev, 3)
    integer :: nroots, kr

    call disp("test_integral: Mach num Mt         = ", M_t)
    call disp("test_integral: Poloidal mode mth   = ", 1d0*mth)

    nu = 300

    !mth = 2
    
    vmin2 = .1*vth
    vmax2 = 5*vth

    ! trapped
    etamin = (1+epst)*etatp
    etamax = (1-epst)*etadt
    !call find_vlim_t(vmin2, vmax2)
    !print *, 'vrange: ', vmin2/vth, vmax2/vth
    
    ! passing
    !etamin = epsp*etatp
    !etamax = (1-epsp)*etatp
    !call find_vlim_p(vmin2, vmax2)

    print *, 'vrange: ', vmin2/vth, vmax2/vth
    
    print *, 'D11/Dp (INT) = ', flux_integral(vmin2, vmax2)
    !print *, 'D11/Dp (ODE) = ', flux_integral_ode(vmin2, vmax2)
    vmin2 = vmin2 + 1d-10*(vmax2-vmin2)
    vmax2 = vmax2 - 1d-10*(vmax2-vmin2)

    D11  = 0d0
    D12  = 0d0
    D11sum = 0d0
    D12sum = 0d0
    dsdreff = 2d0/a*sqrt(s)
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
        
    open(unit=10, file='test_integral2.dat', recl=1024)
    v = vmax2
    vold = vmax2
    do ku = 0, nu-1
       vold = v
       v = vmax2 - ku*(vmax2-vmin2)/(nu-1)
       
       call driftorbit_coarse(etamin, etamax, roots, nroots)
       if (nroots == 0) cycle

       do kr = 1,nroots
          eta_res = driftorbit_root2(max(1d-9*abs(Om_tE),1d-12), etamin, etamax)
          eta = eta_res(1)
          D11 = (vold-v)/vth*D11int_u(v/vth)*dsdreff**(-2)
          D12 = (vold-v)/vth*D12int_u(v/vth)*dsdreff**(-2)
          if ((ku == 1) .OR. (ku == nu-1)) then
             D11 = D11/2d0
             D12 = D12/2d0
          end if
          D11sum = D11sum + D11
          D12sum = D12sum + D12
       end do

       xs = (v/vth)**2
       kappa2s = (1d0/eta_res(1)-B0*(1d0-eps))/(2d0*B0*eps)
       
       write(10, *) v/vth, eta_res(1), D11/Dp, D12/Dp,&
            (eta_res(1)-etamin)/(etamax-etamin),&
            vth, B0, xs, kappa2s
    end do
    print *, 'D11/Dp (SUM) = ', D11sum/Dp, D12sum/Dp
    close(unit=10)
  end subroutine test_integral

  subroutine test_supban
    real(8) :: vmin, vmax
    real(8) :: D1xDp(2)

    etamin = (1+1d-8)*etatp
    etamax = (1-1d-8)*etadt

    vmin = 1d-6*vth
    vmax = 10*vth

    call find_vlim(vmin, vmax)

    D1xdp = 0d0
    D1xdp = flux_integral(vmin,vmax)

    print *, s, eps, q, D1xdp(1), dqds, D1xdp(2)
  end subroutine test_supban

  
  subroutine test_supban2
    real(8) :: vmin, vmax
    real(8) :: D1xDp(2)

    etamin = (1+1d-8)*etatp
    etamax = (1-1d-8)*etadt

    vmin = 1d-6*vth
    vmax = 10*vth

    call find_vlim(vmin, vmax)

    D1xdp = 0d0
    D1xdp = flux_integral_ode(vmin,vmax)

    print *, s, eps, q, D1xdp(1), dqds, D1xdp(2)
  end subroutine test_supban2

  subroutine test_intstep2
    use dvode_f90_m2
    real(8) :: sk, ds, eta_res(2)
    real(8) :: D11, D12, Dp, D11u, vold, dsdreff
    integer :: ks, ns, ku, nu

    real(8) :: y(4), t, tout, rtol(4), atol(4), rwork(128)
    integer :: neq, itask, istate, iopt, itol, iwork(32), ng, jt
    TYPE(VODE_OPTS) :: OPTIONS

    call disp("test_intstep: Mach num Mt   = ", M_t)

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

    nu = 1

    vmin2 = 1d-2*vth
    vmax2 = 10*vth

    ! trapped
    etamin = (1+1d-8)*etatp
    etamax = (1-1d-8)*etadt
    call find_vlim(vmin2, vmax2)
    
    ! passing
    !etamin = 1d-8*etatp
    !etamax = (1-1d-8)*etatp
    !call find_vlim_p(vmin2, vmax2)

    !print *, 'D11/Dp (INT) = ', flux_integral(vmin2, vmax2)
    call disp("test_intstep: vmin/vt= ", vmin2/vth)
    call disp("test_intstep: vmax/vt= ", vmax2/vth)
    
    vmin2 = vmin2 + 1d-10*(vmax2-vmin2)
    vmax2 = vmax2 - 1d-8*(vmax2-vmin2)

    D11  = 0d0
    D12  = 0d0
    D11u = 0d0
    dsdreff = 2d0/a*sqrt(s)
    Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
    
    v = vmax2
    vold = vmax2
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
    !print *, y

    
    open(unit=10, file='test_intstep2.dat', recl=1024)
    v = vmin2
    vold = vmin2
    do ku = 0, nu-1
       vold = v
       v = vmax2 - ku*(vmax2-vmin2)/(nu-1)
       D11u = D11u + (vold-v)/vth*D11int_u(v/vth)*dsdreff**(-2)
       write(10, *) v/vth, D11u/Dp
    end do
    print *, 'D11/Dp (SUM) = ', D11u/Dp
    close(unit=10)

    OPTIONS = SET_OPTS(DENSE_J=.TRUE.,RELERR_VECTOR=RTOL, &
       ABSERR_VECTOR=ATOL,NEVENTS=NG)
    
    open(unit=9, file='test_intstep.dat', recl=1024)
    !write(9, *)                    't',&
    !     '                       vbar',&
    !     '                     etabar',&
    !     '                          I',&
    !     '                      v/vth',&
    !     '                      etark',&
    !     '                     etabar',&
    !     '                      error',&
    !     '                     D11/Dp'   
    do ks = 1, ns
       tout = t+ds
       call dvode_f90(intstep,neq,y,t,tout,itask,istate,options,g_fcn=bbox)
       if (istate == 3) then
          print *, 'D11/Dp (ODE) = ', D11/Dp, D12/Dp
          print *, 'reached boundary'
          print *, 'y = ', y
          exit
       end if
       vold = v
       v = vmin2 + y(1)*(vmax2-vmin2)
       eta_res = driftorbit_root(max(1d-9*abs(Om_tE),1d-12), etamin, etamax)
       D11 = dsdreff**(-2)*y(3) 
       D12 = dsdreff**(-2)*y(4) 
       D11u = D11u + (v-vold)/vth*D11int_u(v/vth)*dsdreff**(-2)
       write(9, *) t, y(1), y(2), y(3), &
            (vmin2 + (vmax2-vmin2)*y(1))/vth, etamin + (etamax-etamin)*y(2), &
            (eta_res(1)-etamin)/(etamax-etamin),&
            (etamin + (etamax-etamin)*y(2)-eta_res(1))/eta_res(1),&
            D11/Dp, D12/Dp
       flush(9)
    end do
    close(unit=9)
  end subroutine test_intstep2
  
  subroutine test_Hmn
    integer :: keta, neta
    real(8) :: vmin2, vmax2
    real(8) :: eta_res(2), H2mean
    !real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

    neta = 1000

    v = vth


    vmin2 = 1d-6*vth
    vmax2 = 1d1*vth
    call find_vlim_p(vmin2, vmax2)
    print *, vmin2/vth, vmax2/vth
    vmin2 = 1d-6*vth
    vmax2 = 1d1*vth
    call find_vlim_t(vmin2, vmax2)
    print *, vmin2/vth, vmax2/vth

    !etamin = 1d-6*etadt
    etamin = etatp*(1+1d-6)
    etamax = etadt*(1-1d-6)
    !etamax = etatp*(1-1d-6)
    eta_res=driftorbit_root(1d-8*abs(c*mi*vth**2/(2*qi*psi_pr)),&
         etamin, etamax)
    print *, 'eta_res_bar = ', (eta_res(1)-etatp)/(etadt-etatp)
    call bounce

    H2mean = 8*pi/mph*q*R0**2/(v**2*taub**2)*epsmn**2/&
         sqrt(eta_res(1)**2*eps**2*B0**2 - (eta_res(1)*B0-1)**2)

    print *, taub, (1/(eta_res(1)*B0) - 1d0)**2/eps**2

    !print *, 'eta res: Hmn = ', H2mean, bounceavg(4)**2+bounceavg(5)**2

    open(unit=9, file='test_Hmn.dat', recl=1024)
    do keta = 1, neta
       eta = etamin + (etamax-etamin)/(neta-1)*(keta-1)
       call bounce
       write(9, *) (eta-etatp)/(etadt-etatp), bounceavg(4)**2+bounceavg(5)**2,&
            H2mean
       flush(9)
    end do
    close(unit=9)
  end subroutine test_Hmn

  subroutine test_boundaries
    real(8) :: v1, v2, disc
    real(8) :: ae, be, ce, dtp
    real(8) :: Omph, dOmphdv, dOmphdeta
    real(8) :: Omth, dOmthdv, dOmthdeta
    real(8) :: Om_tB_bar, Om_tE_bar, Om_th_bar

    print *
    write(*,'(A)') 'test_boundaries:'
    
    dtp = 1d0
    mth = int(-sign(1d0, Om_tE)*ceiling(abs(mph*q)))

    write(*,'(A,I5)') 'mth = ', mth
    
    v1 = 1d-3*vth
    v2 = 1d3*vth
    call find_vlim_t(v1, v2)
    write(*,'(A,2ES10.2)') 'vlim_trap = ', v1/vth, v2/vth
  
    v1 = 1d-3*vth
    v2 = 1d3*vth
    call find_vlim_p(v1, v2)
    write(*,'(A,2ES10.2)') 'vlim_pass = ', v1/vth, v2/vth
    
    !eta = etatp*(1+1d-2)
    !eta = (etatp+etadt)/2d0
    eta = etatp*(1-1d-7)
    !eta = etatp*(1d-5)
    call Om_ph(Omph, dOmphdv, dOmphdeta)
    call Om_th(Omth, dOmthdv, dOmthdeta)
  
    Om_tB_bar = (Omph - Om_tE)/v**2
    Om_tE_bar = Om_tE
    Om_th_bar = Omth/v
    
    ae = Om_tB_bar
    be = (mth/mph+q*dtp)*Om_th_bar
    ce = Om_tE_bar

    disc = be**2 - 4*ae*ce
    print *, mph, mth
    print *, ae, be, ce, disc
    if (disc < 0) then
       print *, 'test_boundaries: discriminant < 0'
    else
       v1 = -2d0*ce/(be+sqrt(disc))
       v2 = -2d0*ce/(be-sqrt(disc))
       print *, 'vrange: u1 = ', v1/vth, ', u2 = ', v2/vth
       print *, 'test_boundaries: u1 = ', v1/vth, ', u2 = ', v2/vth
    end if
  end subroutine test_boundaries
end program main
