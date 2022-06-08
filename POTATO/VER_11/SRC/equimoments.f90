!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine fit_profiles
!
! Computes radial profiled of flux surface averaged guiding center density and toroidal flow
! density of the guiding centers using polynomial fit via Eq.(30) and Eq.(73).
!
  use field_eq_mod,      only : psif,nrad,nzet,rad,zet,psi_sep
  use poicut_mod,        only : rmagaxis,zmagaxis,psimagaxis,psi_bou,rhopol_bou
  use global_invariants, only : toten,perpinv
  use poicut_mod,        only : Rbou_lfs,Zbou_lfs
  use get_matrix_mod,    only : iclass
  use form_classes_doublecount_mod, only : nclasses
  use orbit_dim_mod,     only : numbasef

  implicit none
!
  integer,          parameter :: numbasef_fix=8 !6 !5
  double precision, parameter :: pi=3.14159265358979d0
  logical :: classes_talk
  integer :: nr,nz,ir,iz,i,k,iperp,nperp,ierr,nprof,ienerg,nenerg
  double precision :: rbeg,hr,zbeg,hz,weight,psi,psipow
  double precision :: bmod,phi_elec,phi_elec_min,phi_elec_max
  double precision :: toten_min,toten_max,thermen_max,toten_range
  double precision :: dens,omdens,trapez_fac,perpinv_max
  double precision :: xjperp,xenerg
  double precision, dimension(3) :: x
  integer,          dimension(numbasef_fix)              :: ipiv                 !<= change to allocatable
  double precision, dimension(numbasef_fix,numbasef_fix) :: amat,addmat          !<= change to allocatable
  double precision, dimension(numbasef_fix,2)            :: bvec,resints,ydummy  !<= change to allocatable
!
  numbasef=numbasef_fix
!
  nr=5000 !<= change to smaller number for debugging
  nz=5000 !<= change to smaller number for debugging
  rbeg=rad(1)
  hr=(rad(nrad)-rad(1))/dble(nr)
  zbeg=zet(1)
  hz=(zet(nzet)-zet(1))/dble(nz)
!
  x(2)=0.d0
!
  amat=0.d0
!
  do ir=0,nr
    x(1)=rbeg+hr*dble(ir)
    do iz=0,nz
      x(3)=zbeg+hz*dble(iz)
!
      call get_bmod_and_Phi(x,bmod,phi_elec)
      call thetafun(weight)
!
      psi=(psif-psimagaxis)/(psi_sep-psimagaxis)
      addmat(1,1)=weight*x(1)
      do i=2,numbasef
        addmat(i,1)=addmat(i-1,1)*psi
      enddo
      do i=2,numbasef
        addmat(:,i)=addmat(:,i-1)*psi
      enddo
      amat=amat+addmat
    enddo
  enddo
!
  amat=(2.d0*pi*hr*hz)*amat
!
  open(12345,file='amat.dat')
  do i=1,numbasef_fix
    write(12345,*) amat(i,:)
  enddo
  close(12345)
!
  call find_Phiminmax(phi_elec_min,phi_elec_max)
!
  thermen_max=6.d0
!
  toten_min=phi_elec_min
  toten_max=thermen_max+phi_elec_max
  toten_range=toten_max-toten_min
!
  print *,'miminum potential energy = ',phi_elec_min
  print *,'maximum potential energy = ',phi_elec_max
  print *,'miminum total energy = ',toten_min
  print *,'maximum total energy = ',toten_max
!
  nperp=50 !100
  nenerg=60 !40
!
  bvec=0.d0
!
  !$omp parallel default(firstprivate)
  do ienerg=1,nenerg
!  do ienerg=20,20 !<=fix energy for debugging
    xenerg=(dble(ienerg)-0.5d0)/dble(nenerg)
    toten=toten_min+toten_range*xenerg
print *,'toten = ',toten
!
    call find_Jperpmax(perpinv_max)
!
    do iperp=1,nperp
!    do iperp=25,25 !<=fix Jperp for debugging
      if(iperp.eq.nperp) then
        trapez_fac=0.5d0
      else
        trapez_fac=1.d0
      endif
      xjperp=dble(iperp)/dble(nperp)
      trapez_fac=trapez_fac*2.d0*xjperp
      perpinv=perpinv_max*(1.d0-xjperp**2)
!print *,' perpinv=',perpinv
!
      call find_bounds_fixpoints(ierr)
!
      if(ierr.ne.0) then
        print *,'find_bounds_fixpoints ierr = ',ierr
        cycle
      endif
!
!      classes_talk=.true.
      classes_talk=.false.
!
      call form_classes_doublecount(classes_talk,ierr)
!
      if(ierr.ne.0) then
        print *,'form_classes ierr = ',ierr
        cycle
      endif
!
      do iclass=1,nclasses
!
        call sample_class_doublecount(1000+iclass,ierr)
!
        if(ierr.eq.0) then
!
          call integrate_class_doublecount(2000+iclass,resints)
!
          bvec=bvec+perpinv_max*trapez_fac*resints
        else
          print *,'sample_class_doublecount error',ierr
        endif
      enddo
print *,'perpinv:',iperp,'/',nperp,' toten:',ienerg,'/',nenerg
    enddo
  enddo
  !$omp end parallel
!
  bvec=0.5d0*sqrt(pi)*toten_range*bvec/dble(nperp*nenerg)
!
  open(12345,file='bvec.dat')
  do i=1,numbasef_fix
    write(12345,*) bvec(i,:)
  enddo
  close(12345)
!
  call dgesv(numbasef,2,amat,numbasef,ipiv,bvec,numbasef,ierr)
!
  open(1900,file='compeqpar.dat')
!
  nprof=100
  do i=0,nprof
    psi=rhopol_bou**2*dble(i)/dble(nprof)
    psipow=1.d0
    dens=0.d0
    omdens=0.d0
    do k=1,numbasef
      dens=dens+bvec(k,1)*psipow
      omdens=omdens+bvec(k,2)*psipow
      psipow=psipow*psi
    enddo
    write (1900,*) psi,dens,omdens
  enddo
!
  close(1900)
!
  end subroutine fit_profiles
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine plotrot(Phi_eff,rho_pol)
!
! Computes and writes down density parameter of the Maxwellin
! and analytical result for the toroidal flow density, Eq.(82)
!
  use field_eq_mod,       only : psi_sep,psi_axis
!
  implicit none
!
  integer :: i,npoi
  double precision :: Phi_eff,rho_pol
  double precision :: psi,dens,temp,ddens,dtemp,phi_elec,dPhi_dpsi, &
                      psistep,epsdif,hpsi,torfl
!
  npoi=1000
  psistep=(psi_sep-psi_axis)*rho_pol**2/dble(npoi)
  epsdif=1.d-6
  hpsi=epsdif*abs(psi_sep-psi_axis)
!
  open(337,file='testeqpar.dat')
!
  do i=0,npoi
    psi=psi_axis+dble(i)*psistep
!
    call phielec_of_psi(psi,phi_elec,dPhi_dpsi)
!
    call denstemp_of_psi(psi,dens,temp,ddens,dtemp)
!
    torfl=Phi_eff*(dens*dPhi_dpsi+ddens)
    write(337,*) rho_pol**2*dble(i)/dble(npoi),dens,torfl
  enddo
!
  close(337)
!
  end subroutine plotrot
