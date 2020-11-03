!
  module phielec_of_psi_mod
    integer, parameter :: npolyphi=1
    double precision, dimension(0:npolyphi) :: polyphi
  end module phielec_of_psi_mod
!
!cccccccccccccccccccccc
!
  program classify_orbits
!
  use parmot_mod, only : rmu,ro0
! 14.11.2011  use field_eq_mod, only : btf,rtf
! collisions
  use collis_alp, only : swcoll,iswmod
! end collisions
  use vparzero_line_mod, only : nvpline,s_line,vpar_line
  use orbit_dim_mod, only : write_orb,iunit1
!
  implicit none
!
  integer, parameter          :: next=3
!
  double precision, parameter :: pi=3.14159265358979d0
  double precision,parameter  :: c=2.9979d10
  double precision,parameter  :: e_charge=4.8032d-10
  double precision,parameter  :: e_mass=9.1094d-28
  double precision,parameter  :: p_mass=1.6726d-24
  double precision,parameter  :: ev=1.6022d-12
!
  logical :: closed_vpl,trapped,copass,ctrpass,triplet,wrorb
!
  integer          :: npoi,icase,ierr,L1i,nper,npoiper,i,j,k,ntimstep,ntestpart
  integer          :: ipart,notrace_passing,loopskip,iskip,ilost
  integer          :: nline,npoi_pot,iun1,iun2
  double precision :: dphi,rbeg,phibeg,zbeg,bmod00,rcham,rlarm,bmax,bmin
  double precision :: tau,dtau,xi,v0,bmod_ref,E_alpha,trace_time
  double precision :: toten,perpinv,rmin,rmax,zmin,zmax
  integer,          dimension(:),   allocatable :: ibinsrc_x
!
  integer, parameter          :: ndim=5
  double precision, parameter :: relerr=1d-10 !8
  double precision :: dL2_pol,dL2_pol_start,dtau_newt,taub,delphi
  double precision :: tau0,RNorm,ZNorm,vnorm,dnorm,vel_pol
  double precision :: sigma,dR_ds,dZ_ds,rend,zend
  double precision, dimension(ndim) :: z,z_start,vz
  double precision :: bmod,phi_elec,p2,sst,rst,zst,alam2,raxis,zaxis,hr,hz,eps
  double precision :: sst2,rst2,zst2,rxpoi2,zxpoi2,sigma_x2
  double precision :: raxis_co,zaxis_co,raxis_ctr,zaxis_ctr,p_phi
  double precision :: rxpoi,zxpoi,sigma_x,s_tp,s_tp2,potato,s,drds,dzds
  double precision, dimension(next) :: extraset
!
  external velo
!
  rmu=1.d30
  swcoll=.false.
!
  E_alpha=5d3
! alpha particle velocity, cm/s
  v0=sqrt(E_alpha*ev/(2.d0*p_mass))
! 14.04.2013 end
!
! Larmor radius:
  bmod_ref=1.d0
  rlarm=v0*2.d0*p_mass*c/(e_charge*bmod_ref)
! normalized slowing down time:
  tau=1d4 !trace_time*v0
! normalized time step:
  ro0=rlarm
  z=0.d0
  z(1)=171.d0
  z(3)=4.d0
  z(4)=1.d0
  z(5)=0.d0 !0.25d0
!
  call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
  toten=z(4)**2+phi_elec
!
  perpinv=z(4)**2*(1.d0-z(5)**2)/bmod
!perpinv=1.03d0*perpinv
!perpinv=0.9794d0*perpinv
!******
!perpinv=0.97d0*perpinv !uncomment this and one of 3 below
!perpinv=1.05d0*perpinv !problem !!! no convergence
!perpinv=1.07d0*perpinv !problem !!! jump
!perpinv=1.09d0*perpinv !problem !!! stagnation_point ierr =            1
!*****
!perpinv=0.97d0*perpinv
!perpinv=0.8152d0*perpinv
!perpinv=0.81d0*perpinv               !!!!
!perpinv=0.98989d0*perpinv

!perpinv=1.20d0*perpinv
!perpinv=1.40d0*perpinv
!perpinv=0.98232d0*perpinv !for ampl=0.25d0

perpinv=0.99d0*perpinv
!
  nline=200
!
  Rmin=140.d0
  Rmax=200.d0
  Zmin=-50.d0
  Zmax=50.d0
!
  ntimstep=100
  dtau=Rmax/dble(ntimstep)
!
  call find_vparzero_line(nline,toten,perpinv,rmin,rmax,zmin,zmax,icase)
!
  if(icase.eq.1) stop
  if(icase.ne.4) then
    do i=0,nvpline
      write(21,*) vpar_line(:,i),s_line(i)
    enddo
!
    do i=0,1000
      s=1.d-3*s_line(nvpline)*dble(i)
      call vparzero_line(s,rbeg,zbeg,drds,dzds)
      z(1)=rbeg
      z(3)=zbeg
      call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
      write(22,*) rbeg,zbeg,drds,dzds,toten-perpinv*bmod-phi_elec
    enddo
!
    call stagnation_point(icase,toten,perpinv,            &
                          sst,rst,zst,sst2,rst2,zst2,ierr)
!
    if(ierr.ne.0) then
      print *,'stagnation_point ierr = ',ierr
      stop
    endif
    print *,'stagnation point 1: ',rst,zst
    if(icase.eq.3) print *,'stagnation point 2: ',rst2,zst2
    trapped=.true.
  endif
!
! Invariant axes:
!
  call find_bstar_axes(icase,toten,perpinv,                   &
                       sst,rst,zst,sst2,rst2,zst2,dtau,       &
                       trapped,copass,ctrpass,triplet,        &
                       raxis_co,zaxis_co,raxis_ctr,zaxis_ctr, &
                       rxpoi,zxpoi,sigma_x,s_tp,              &
                       rxpoi2,zxpoi2,sigma_x2,s_tp2,ierr)
!
  if(ierr.ne.0) then
    print *,'find_bstar_axes ierr = ',ierr
    stop
  endif
!
  if(sigma_x.ne.0.d0) then
    print *,'right X-point: ',rxpoi,zxpoi, 'sigma: ',sigma_x
  else
    print *,'no right X-point'
  endif
!
  if(sigma_x2.ne.0.d0) then
    print *,'left X-point: ',rxpoi2,zxpoi2, 'sigma: ',sigma_x2
  else
    print *,'no left X-point'
  endif
!
  open(1,file='invar_axes.dat')
!
  if(copass) then
    print *,'co-passing axis: ',raxis_co,zaxis_co
    write(1,*) raxis_co,zaxis_co
  else
    print *,'no co-passing axis'
  endif
!
  if(ctrpass) then
    print *,'counter-passing axis: ',raxis_ctr,zaxis_ctr
    write(1,*) raxis_ctr,zaxis_ctr
  else
    print *,'no counter-passing axis'
  endif
!
  if(sigma_x.ne.0.d0) then
    write(1,*) rxpoi,zxpoi
!
    call vparzero_line(s_tp,z(1),z(3),dR_ds,dZ_ds)
!
    write(1,*) z(1),z(3)
  endif
!
  if(sigma_x2.ne.0.d0) then
    write(1,*) rxpoi2,zxpoi2
!
    call vparzero_line(s_tp2,z(1),z(3),dR_ds,dZ_ds)
!
    write(1,*) z(1),z(3)
  endif
!
  close(1)
!
! End invariant axes
!
! Plot orbits:
!
  eps=1.d-2
!
  npoi=15 !15 !0 !0 !2500
  npoi_pot=5
  iun1=1
  iun2=2
  wrorb=.true. !.false.
!
  call plot_orbit_set(wrorb,iun1,iun2,npoi,npoi_pot,dtau,icase,toten,perpinv,        &
                      trapped,copass,ctrpass,triplet,sst,rst,zst,sst2,rst2,zst2,     &
                      sigma_x,s_tp,rxpoi,zxpoi,sigma_x2,s_tp2,rxpoi2,zxpoi2,         &
                      raxis_ctr,zaxis_ctr,raxis_co,zaxis_co)
!
  end program classify_orbits
!
!ccccccccc
!
  subroutine phielec_of_psi(psi,phi_elec,dPhi_dpsi)
!
  use phielec_of_psi_mod, only : npolyphi,polyphi
  use field_eq_mod,       only : psi_sep
!
  implicit none
!
  logical :: prop=.true.
  integer :: i
  double precision :: psi,phi_elec,dPhi_dpsi
  double precision :: ampl
!
  ampl=1.12d0
!
  if(prop) then
    prop=.false.
    polyphi(0)=0.d0
    polyphi(1)=-ampl/psi_sep
  endif
!
  phi_elec=0.d0
  dPhi_dpsi=0.d0
!
  do i=npolyphi,0,-1
    phi_elec=polyphi(i)+phi_elec*psi
  enddo
!
  do i=npolyphi,1,-1
    dPhi_dpsi=polyphi(i)*dble(i)+dPhi_dpsi*psi
  enddo
!
  end subroutine phielec_of_psi
