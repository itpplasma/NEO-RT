!
  program classify_orbits
!
  use parmot_mod, only : rmu,ro0
! 14.11.2011  use field_eq_mod, only : btf,rtf
! collisions
  use collis_alp, only : swcoll,iswmod
! end collisions
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
  logical :: trapped,copass,ctrpass,triplet,wrorb
!
  integer          :: npoi,ierr,L1i,nper,npoiper,i,j,k,ntimstep,ntestpart
  integer          :: ipart,notrace_passing,loopskip,iskip,ilost
  integer          :: nline,npoi_pot,iun1,iun2
  double precision :: dphi,rbeg,phibeg,zbeg,bmod00,rcham,rlarm,bmax,bmin
  double precision :: tau,dtau,xi,v0,bmod_ref,E_alpha,trace_time
  double precision :: toten,perpinv,rmin,rmax,zmin,zmax
  integer,          dimension(:),   allocatable :: ibinsrc_x
  double precision, dimension(:),   allocatable :: R_line,Z_line
  double precision, dimension(:,:), allocatable :: xstart,zstart
  double precision, dimension(:),   allocatable :: confpart_trap,confpart_pass
!
  integer, parameter          :: ndim=5
  double precision, parameter :: relerr=1d-10 !8
  double precision :: dL2_pol,dL2_pol_start,dtau_newt,taub,delphi
  double precision :: tau0,RNorm,ZNorm,vnorm,dnorm,vel_pol
  double precision :: sigma,dR_dZ,rend,zend
  double precision, dimension(ndim) :: z,z_start,vz
  double precision :: bmod,phi_elec,p2,rst,zst,alam2,raxis,zaxis,hr,hz,eps
  double precision :: raxis_co,zaxis_co,raxis_ctr,zaxis_ctr,p_phi
  double precision :: rxpoi,zxpoi,sigma_x,z_tp,potato
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
perpinv=0.97d0*perpinv
!perpinv=0.8152d0*perpinv
!perpinv=0.81d0*perpinv

!
  nline=200
!
  allocate(R_line(-nline:nline),Z_line(-nline:nline))
!
  Rmin=140.d0
  Rmax=200.d0
  Zmin=-50.d0
  Zmax=50.d0
!
  ntimstep=100
  dtau=Rmax/dble(ntimstep)
!
  call find_vparzero_line(nline,toten,perpinv,rmin,rmax,zmin,zmax, &
                                R_line,Z_line,ierr)
print *,'find_vparzero_line ierr = ',ierr
!
  do i=-nline,nline
    write(21,*) R_line(i),Z_line(i)
  enddo
!
  if(ierr.eq.0) then
    do i=0,1000
      zbeg=zmin+1.d-3*(zmax-zmin)*dble(i)
      call vparzero_line(zbeg,rbeg,delphi)
      write(22,*) rbeg,zbeg,delphi,    toten-perpinv*bmod
    enddo
!
    call stagnation_point(toten,perpinv,rst,zst,ierr)
!
    if(ierr.ne.0) then
      print *,'stagnation_point ierr = ',ierr
      stop
    endif
    print *,'stagnation point: ',rst,zst
    trapped=.true.
  elseif(ierr.eq.1) then
    rst=R_line(0)
    zst=Z_line(0)
    print *,'stagnation point replaced with: ',rst,zst
    trapped=.false.
  else
    print *,'no orbits'
    stop
  endif
!
! Invariant axes:
!
  call find_bstar_axes(toten,perpinv,rst,zst,dtau,trapped,    &
                       copass,ctrpass,triplet,                &
                       raxis_co,zaxis_co,raxis_ctr,zaxis_ctr, &
                       rxpoi,zxpoi,sigma_x,z_tp,ierr)
!
  if(ierr.ne.0) then
    print *,'find_bstar_axes ierr = ',ierr
    stop
  endif
!
  if(sigma_x.ne.0.d0) then
    print *,'X-point: ',rxpoi,zxpoi
  else
    print *,'no X-point'
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
  if(triplet) then
    write(1,*) rxpoi,zxpoi
    z(3)=z_tp
!
    call vparzero_line(z(3),z(1),dR_dZ)
!
    write(1,*) z(1),z(3)
  else
    z_tp=zst
  endif
!
  close(1)
!
! End invariant axes
!
! Plot orbits:
!
  eps=1.d-3
!
  npoi=15 !0 !0 !2500
  npoi_pot=5
  iun1=1
  iun2=2
  wrorb=.true. !.false.
!
  call plot_orbit_set(wrorb,iun1,iun2,npoi,npoi_pot,dtau,eps,toten,perpinv, &
                      trapped,copass,ctrpass,triplet,rst,zst,sigma_x,Z_tp,  &
                      rxpoi,zxpoi,raxis_ctr,zaxis_ctr,raxis_co,zaxis_co)
!
  end program classify_orbits
