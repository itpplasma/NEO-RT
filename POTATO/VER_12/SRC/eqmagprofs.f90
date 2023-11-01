!
  module eqmagprofs_mod
  integer :: nsurf=300
  double precision :: sig_pol,sig_tor
  double precision, dimension(:,:), allocatable :: profs_mag
  end module eqmagprofs_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine load_eqmagprofs
!
! Preloads the profiles for eqmagprofs (see below) using field line integration
!
  use field_eq_mod, only   : psif
  use poicut_mod, only     : npc,rpc_arr,zpc_arr,rmagaxis,zmagaxis
  use eqmagprofs_mod, only : nsurf,sig_pol,sig_tor,profs_mag
!
  implicit none
!
  integer, parameter :: ndim=5, niter=20, iunit=3142
  double precision, parameter :: twopi=3.14159265358979d0*2.d0
  double precision, parameter :: relerr=1d-10
  logical :: firstpass, write_surf=.false., write_profs=.true.
  integer :: isurf,iter
  double precision :: rbeg,rend,zbeg,zend,w,phi0,dphi,phi, &
                      delr,delz,sign_delz,ddphi
  double precision, dimension(ndim) :: y,dy
!
  allocate(profs_mag(7,0:nsurf))
!
  rbeg=rmagaxis
  rend=rpc_arr(npc)
  zbeg=zmagaxis
  zend=zpc_arr(npc)
  delr=rend-rbeg
  delz=zend-zbeg
!
  phi0=0.d0
  dphi=0.1d0
  y(1)=rbeg
  y(2)=zbeg
  y(3:5)=0.d0
!
  call rhs(phi0,y,dy)
!
  profs_mag(:,0)=0.d0
  profs_mag(3,0)=psif ! poloidal flux on axis
!
  if(write_surf) open(iunit,file='magsurf.dat')
!
  do isurf=1,nsurf
    w=dble(isurf)/dble(nsurf)
    y(1)=rbeg+delr*w
    y(2)=zbeg+delz*w
    y(3:5)=0.d0
    if(write_surf) write(iunit,*) y(1:2)
!
    call odeint_allroutines(y,ndim,phi0,dphi,relerr,rhs)
!
    phi=dphi
    if(write_surf) write(iunit,*) y(1:2)
    sign_delz=sign(1.d0,y(2)-zbeg-delz*(y(1)-rbeg)/delr)
    firstpass=.true.
!
    do
!
      call odeint_allroutines(y,ndim,phi0,dphi,relerr,rhs)
!
      phi=phi+dphi
      if(sign_delz*(y(2)-zbeg-delz*(y(1)-rbeg)/delr).lt.0.d0) then
        if(firstpass) then
          firstpass=.false.
          sign_delz=-sign_delz
        else
          exit
        endif
      endif
      if(write_surf) write(iunit,*) y(1:2)
    enddo
!
    do iter=1,niter
!
      call rhs(phi0,y,dy)
!
      ddphi=(zbeg+delz*(y(1)-rbeg)/delr-y(2))/dy(2)
!
      call odeint_allroutines(y,ndim,phi0,ddphi,relerr,rhs)
!
      phi=phi+ddphi
      if(abs(dphi).lt.relerr) exit
    enddo
!
    if(write_surf) write(iunit,*) y(1:2)
    if(write_surf) write(iunit,*) 'NaN,NaN'
    if(write_surf) write(iunit,*) ' '
!
    profs_mag(3,isurf)=psif        ! poloidal flux
    profs_mag(5,isurf)=phi/twopi   ! safety factor
    profs_mag(6,isurf)=y(4)/y(3)   ! averaged nabla psi
    profs_mag(7,isurf)=y(5)        ! surface area
  enddo
!
! extrapolate safety factor to the axis:
  profs_mag(5,0)=3.d0*(profs_mag(5,1)-profs_mag(5,2))+profs_mag(5,3)
!
! toroidal flux (integrate $d \psi_{tor} = q d \psi_{pol}$):
  do isurf=1,nsurf
    profs_mag(4,isurf)=profs_mag(4,isurf-1)+0.5d0*(profs_mag(5,isurf)+profs_mag(5,isurf-1)) &
                      *(profs_mag(3,isurf)-profs_mag(3,isurf-1))
  enddo
!
! poloidal radius:
  profs_mag(1,:)=sqrt((profs_mag(3,:)-profs_mag(3,0)) &
                /(profs_mag(3,nsurf)-profs_mag(3,0)))
! toroidal radius:
  profs_mag(2,:)=sqrt((profs_mag(4,:)-profs_mag(4,0)) &
                /(profs_mag(4,nsurf)-profs_mag(4,0)))
!
  sig_pol=sign(1.d0,profs_mag(3,nsurf)-profs_mag(3,0))
  sig_tor=sign(1.d0,profs_mag(4,nsurf)-profs_mag(4,0))
!
  if(write_profs) then
    open(iunit,file='eqmagprofs.dat')
    do isurf=0,nsurf
      write (iunit,*) profs_mag(:,isurf)
    enddo
    close(iunit)
  endif
!
  profs_mag(3,:)=profs_mag(3,:)*sig_pol
  profs_mag(4,:)=profs_mag(4,:)*sig_tor
!
!------------  
!
  contains
!
!------------  
!
  subroutine rhs(phi,y,dy)
!
  use field_eq_mod, only : dpsidr,dpsidz
!
  implicit none
!
  double precision :: phi
  double precision, dimension(ndim) :: y,dy
  double precision :: bmod,sqrtg
  double precision, dimension(3)    :: x,bder,hcovar,hctrvr,hcurl,derphi
!
  x(1)=y(1)
  x(2)=0.d0
  x(3)=y(2)
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  dy(1)=hctrvr(1)/hctrvr(2)
  dy(2)=hctrvr(3)/hctrvr(2)
  dy(3)=1.d0/(hctrvr(2)*bmod)
  dy(4)=sqrt(dpsidr**2+dpsidz**2)/(hctrvr(2)*bmod)
  dy(5)=twopi*y(1)*sqrt(dy(1)**2+dy(2)**2)
!
  end subroutine rhs
!
!------------  
!
  end subroutine load_eqmagprofs
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine eqmagprofs(label_sw,rho_pol,rho_tor,psi_pol,psi_tor, &
                        qsaf,avnablapsi,surfarea,                 &
                        drho_pol,drho_tor,dpsi_pol,dpsi_tor,      &
                        dqsaf,davnablapsi,dsurfarea)
!
! Computes various flux surface labels, safety factor, average nabla psi and flux surface area
! as functions of a chosen flux surface label.
! Computes also first derivatives of those quantities over the chosen flux surface label.
! Arguments:
! label_sw (in)     - switch for the label used as independent variable, uses:
!                 1 - poloidal radius
!                 2 - toroidal radius
!                 3 - poloidal flux (with a sign - can be decreasing)
!                 4 - toroidal flux (with a sign - can be decreasing)
! rho_pol (in/out)  - poloidal radius
! rho_tor (in/out)  - toroidal radius
! psi_pol (in/out)  - poloidal flux
! psi_tor (in/out)  - toroidal flux
! avnablapsi (out)  - flux surface averaged nabla psi, $<|\nabla \psi|>$
! surfarea   (out)  - flux surface area
!
! drho_pol,drho_tor,dpsi_pol,dpsi_tor,dqsaf,davnablapsi,dsurfarea - derivatives of the above 
! quantities over flux surface label chosen for independent variable
!
  use eqmagprofs_mod, only : nsurf,sig_pol,sig_tor,profs_mag
!
  implicit none
!
  integer, parameter :: nplag = 4 !stencil for Largange polynomial interpolation (polynomial order + 1)
  integer, parameter :: nder  = 1 !number of derivatives from Largange polynomial interpolation
  integer            :: label_sw,ibeg,iend
  double precision   :: rho_pol,rho_tor,psi_pol,psi_tor,     &
                        qsaf,avnablapsi,surfarea,            &
                        drho_pol,drho_tor,dpsi_pol,dpsi_tor, &
                        dqsaf,davnablapsi,dsurfarea
  double precision   :: abscissa
  double precision, dimension(7) :: vec,dvec
  double precision, dimension(0:nder,nplag) :: coef
!
  select case(label_sw)
  case(1)
    abscissa=rho_pol
  case(2)
    abscissa=rho_tor
  case(3)
    abscissa=psi_pol*sig_pol
  case(4)
    abscissa=psi_tor*sig_tor
  case default
    print *,'eqmagprofs: wrong label switch'
    return
  end select
!
  call binsrc(profs_mag(label_sw,0:nsurf),0,nsurf,abscissa,ibeg)
!
  ibeg=max(1,ibeg-nplag/2)
  iend=ibeg+nplag-1
  if(iend.gt.nsurf) then
    iend=nsurf
    ibeg=iend+1-nplag
  endif
!
  call plag_coeff(nplag,nder,abscissa,profs_mag(label_sw,ibeg:iend),coef)
!
  vec=matmul(profs_mag(:,ibeg:iend),coef(0,:))
  dvec=matmul(profs_mag(:,ibeg:iend),coef(1,:))
!
  select case(label_sw)
  case(3)
    dvec=dvec*sig_pol
  case(4)
    dvec=dvec*sig_tor
  end select
!
  vec(3)=vec(3)*sig_pol
  dvec(3)=dvec(3)*sig_pol
  vec(4)=vec(4)*sig_tor
  dvec(4)=dvec(4)*sig_tor
!
  if(label_sw.ne.1) rho_pol=vec(1)
  if(label_sw.ne.2) rho_tor=vec(2)
  if(label_sw.ne.3) psi_pol=vec(3)
  if(label_sw.ne.4) psi_tor=vec(4)
  qsaf=vec(5)
  avnablapsi=vec(6)
  surfarea=vec(7)
!
  drho_pol=dvec(1)
  drho_tor=dvec(2)
  dpsi_pol=dvec(3)
  dpsi_tor=dvec(4)
  dqsaf=dvec(5)
  davnablapsi=dvec(6)
  dsurfarea=dvec(7)
!
  end subroutine eqmagprofs
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine test_eqmagprofs
!
  use eqmagprofs_mod, only : nsurf,sig_pol,sig_tor,profs_mag
!
  implicit none
!
  integer            :: label_sw,npoi=1000,i
  double precision   :: rho_pol,rho_tor,psi_pol,psi_tor,     &
                        qsaf,avnablapsi,surfarea,            &
                        drho_pol,drho_tor,dpsi_pol,dpsi_tor, &
                        dqsaf,davnablapsi,dsurfarea,         &
                        volume
!
  label_sw=1
  do i=0,npoi
    rho_pol=dble(i)/dble(npoi)
    call eqmagprofs(label_sw,rho_pol,rho_tor,psi_pol,psi_tor, &
                    qsaf,avnablapsi,surfarea,                 &
                    drho_pol,drho_tor,dpsi_pol,dpsi_tor,      &
                    dqsaf,davnablapsi,dsurfarea)
    write(2001,*) rho_pol,rho_tor,psi_pol,psi_tor,          &
                  qsaf,avnablapsi,surfarea,                 &
                  drho_pol,drho_tor,dpsi_pol,dpsi_tor,      &
                  dqsaf,davnablapsi,dsurfarea
  enddo
!
  label_sw=2
  do i=0,npoi
    rho_tor=dble(i)/dble(npoi)
    call eqmagprofs(label_sw,rho_pol,rho_tor,psi_pol,psi_tor, &
                    qsaf,avnablapsi,surfarea,                 &
                    drho_pol,drho_tor,dpsi_pol,dpsi_tor,      &
                    dqsaf,davnablapsi,dsurfarea)
    write(2002,*) rho_pol,rho_tor,psi_pol,psi_tor,          &
                  qsaf,avnablapsi,surfarea,                 &
                  drho_pol,drho_tor,dpsi_pol,dpsi_tor,      &
                  dqsaf,davnablapsi,dsurfarea
  enddo
!
  label_sw=3
  volume=0.d0
  do i=0,npoi
    psi_pol=(profs_mag(3,0)+dble(i)/dble(npoi)*(profs_mag(3,nsurf)-profs_mag(3,0)))*sig_pol
    call eqmagprofs(label_sw,rho_pol,rho_tor,psi_pol,psi_tor, &
                    qsaf,avnablapsi,surfarea,                 &
                    drho_pol,drho_tor,dpsi_pol,dpsi_tor,      &
                    dqsaf,davnablapsi,dsurfarea)
    write(2003,*) rho_pol,rho_tor,psi_pol,psi_tor,          &
                  qsaf,avnablapsi,surfarea,                 &
                  drho_pol,drho_tor,dpsi_pol,dpsi_tor,      &
                  dqsaf,davnablapsi,dsurfarea
    volume=volume+surfarea/avnablapsi
  enddo
  volume=volume*(profs_mag(3,nsurf)-profs_mag(3,0))/dble(npoi)
  print *,'plasma volume = ',sngl(volume),' cm^{-3},', &
          sngl(volume/500.d0),' half-liters,',sngl(volume*1d-6),' cubic meters'
!
  label_sw=4
  do i=0,npoi
    psi_tor=(profs_mag(4,0)+dble(i)/dble(npoi)*(profs_mag(4,nsurf)-profs_mag(4,0)))*sig_tor
    call eqmagprofs(label_sw,rho_pol,rho_tor,psi_pol,psi_tor, &
                    qsaf,avnablapsi,surfarea,                 &
                    drho_pol,drho_tor,dpsi_pol,dpsi_tor,      &
                    dqsaf,davnablapsi,dsurfarea)
    write(2004,*) rho_pol,rho_tor,psi_pol,psi_tor,          &
                  qsaf,avnablapsi,surfarea,                 &
                  drho_pol,drho_tor,dpsi_pol,dpsi_tor,      &
                  dqsaf,davnablapsi,dsurfarea
  enddo
!
  stop
  end subroutine test_eqmagprofs
