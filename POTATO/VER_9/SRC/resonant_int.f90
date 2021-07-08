!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! modules
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module resint_mod
    integer          :: nmodes
    double precision :: twopim2,rm3,taub_new,delphi_new
    integer, dimension(:), allocatable :: marr,narr
    double precision, dimension(:), allocatable :: delint_mode
  end module resint_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! rotines
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine velo_res(dtau,z,vz)
!
! Computes the extended RHS for equations of motion and Fourier amplitude
! of the perturbed Hamiltonian $\hat H_\bm$, Eq.(93)
!
  use orbit_dim_mod,     only : neqm,next,numbasef
  use global_invariants, only : toten,perpinv
  use resint_mod,        only : twopim2,rm3,taub_new,delphi_new
!
  implicit none
!
  double complex, parameter :: imun=(0.d0,1.d0)
!
  double precision :: dtau,bmod,phi_elec
  double precision, dimension(neqm+next) :: z,vz
  double complex :: bmod_n,comfac
!
  call velo(dtau,z(1:neqm),vz(1:neqm))
  call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
  call bmod_pert(z(1),z(3),bmod_n)
!
  comfac=(2.d0*(toten-phi_elec)/bmod-perpinv)*bmod_n &
        *exp(imun*(rm3*z(2)-(twopim2+delphi_new*rm3)*z(6)/taub_new))
!
  vz(6)=1.d0
  vz(7)=dble(comfac)
  vz(8)=dimag(comfac)
!
  end subroutine velo_res
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine pertham(z,absHn2)
!
! Computes modulus squared of the Fourier amplitude of the normalized perturbed 
! Hamiltoninan, $|\hat H_\bm|^2$, with $\hat H_\bm$ defined by Eq.(93).
!
  use orbit_dim_mod,     only : neqm,next     ,write_orb,iunit1
  use global_invariants, only : toten,perpinv,dtau
  use resint_mod,        only : taub_new,delphi_new
!
!
  implicit none
!
  double complex, parameter :: imun=(0.d0,1.d0)
!
  double precision :: absHn2,bmod,phi_elec,taub,delphi
  double precision, dimension(neqm) :: z
  double precision, dimension(:), allocatable :: extraset
!
  external :: velo,velo_res
!
  call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
  toten=z(4)**2+phi_elec
  perpinv=z(4)**2*(1.d0-z(5)**2)/bmod
!
  next=0
  allocate(extraset(next))
!
  call find_bounce(next,velo,dtau,z,taub,delphi,extraset)
!
  taub_new=taub
  delphi_new=delphi
  deallocate(extraset)
!
  next=3
  allocate(extraset(next))
  extraset=0.d0
!
  call find_bounce(next,velo_res,dtau,z,taub,delphi,extraset)
!
  absHn2=(extraset(2)/taub)**2+(extraset(3)/taub)**2
!
  deallocate(extraset)
!
  end subroutine pertham
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine integrate_class_resonances(xint_using_delta)
!
! Computes sum over resonances $x=x^{res}_{(\bm,k)}$ in Eq.(94) for a given class $k$
!
  use find_all_roots_mod, only : customgrid,ncustom,xcustom,nroots,roots
  use get_matrix_mod,     only : relmargin,iclass
  use form_classes_doublecount_mod, only : ifuntype,R_class_beg,R_class_end,sigma_class
  use resint_mod,         only : nmodes,marr,narr,twopim2,rm3,delint_mode
  use orbit_dim_mod,      only : neqm
  use global_invariants,  only : toten,perpinv,Phi_eff
  use sample_matrix_mod,  only : npoi,xarr
!
  implicit none
!
  double precision, parameter :: twopi=3.14159265358979d0*2.d0
!
  integer          :: mode,iroot,ierr
  double precision :: xint_using_delta,relmargin_loc,widthclass,xbeg,xend
  double precision :: rescond,dresconddx,psiast,dpsiastdx,taub,delphi
  double precision :: one_res,sigma,delta_R,Rst,xi,dxi_dx,dpsiast_dRst,absHn2
  double precision :: toten_loc,perpinv_loc,fmaxw,A1ast,A2ast
  double precision :: dens,temp,ddens,dtemp,phi_elec,dPhi_dpsi
  double precision, dimension(neqm) :: z
!
  toten_loc=toten
  perpinv_loc=perpinv
!
  sigma=sigma_class(iclass)
  delta_R=R_class_end(iclass)-R_class_beg(iclass)
  relmargin_loc=1.d-8
  widthclass=1.d0
!
  call classbounds(ifuntype(iclass),relmargin_loc,widthclass,xbeg,xend)
!
  customgrid=.true.
  ncustom=npoi
  allocate(xcustom(ncustom))
  xcustom=xarr
!
  xint_using_delta=0.d0
!
  do mode=1,nmodes
    twopim2=twopi*dble(marr(mode))
    rm3=dble(narr(mode))
    delint_mode(mode)=0.d0
!
    call find_all_roots(get_rescond,xbeg,xend,ierr)
!
    if(ierr.ne.0) then
      print *,'integrate_class_resonances: error in find_all_roots'
      return
    endif
!
    do iroot=1,nroots
!
      call get_rescond(roots(iroot),rescond,dresconddx)
      call xi_func(ifuntype(iclass),roots(iroot),xi,dxi_dx)
!
      Rst=R_class_beg(iclass)+delta_R*xi
!
      call starter_doublecount(toten,perpinv,sigma,Rst,   &
                               psiast,dpsiast_dRst,z,ierr)
!
      if(ierr.ne.0) then
        print *,'integrate_class_resonances: error in starter_doublecount'
        cycle
      endif
!
      dpsiastdx=dpsiast_dRst*delta_R*dxi_dx     !$\difp{\psi^\ast}{x}$
!
      call pertham(z,absHn2)
      call equilmaxw(psiast,fmaxw)
      call denstemp_of_psi(psiast,dens,temp,ddens,dtemp)
      call phielec_of_psi(psiast,phi_elec,dPhi_dpsi)
!
      toten=toten_loc
      perpinv=perpinv_loc
!
! Non-local thermodynamic forces Eq.(87):
      A2ast=dtemp/temp
      A1ast=ddens/dens+dPhi_dpsi/temp-1.5d0*A2ast
!
! Expression under summation sign except the last line in Eq.(94):
      one_res=abs(dpsiastdx/dresconddx)*absHn2*fmaxw            &
             *(Phi_eff*taub*(A1ast+A2ast*(toten-phi_elec)/temp)+delphi/temp)
!
!
! emulator of box average Eq.(89) (or Eq.(88)):
      one_res=one_res*Phi_eff*taub
! end emulator of box average
!
!
      xint_using_delta=xint_using_delta+one_res
      delint_mode(mode)=delint_mode(mode)+one_res
    enddo
  enddo
!
  customgrid=.false.
  deallocate(xcustom)
!
!------------
  contains
!------------
!
  subroutine get_rescond(x,rescond,dresconddx)
!
! Computes resonance condition $F(x)=\Delta\varphi_b+2\pi m_2/m_3$ 
! and its derivative $F^\prime(x)$ for $F(x)=0$ root finding.
! Computes as by-products normalized toroidal momentum $\psi^\ast$, 
! bounce time $\tau_b$ and toroidal displacement $\Delta\varphi_b$.
!
  use sample_matrix_mod, only : n1
!
  implicit none
!
  double precision :: x,rescond,dresconddx
  double precision, dimension(n1) :: vec,dvec
!
  call interpolate_class_doublecount(x,vec,dvec)
!
  psiast=vec(1)               ! $\psi^\ast$
  taub=vec(2)                 ! $\tau_b$
  delphi=vec(3)               ! $\Delta\varphi_b$
  rescond=delphi+twopim2/rm3
  dresconddx=dvec(3)
!
  end subroutine get_rescond
!
!------------
!
  end subroutine integrate_class_resonances
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine resonant_flux
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
  use resint_mod,        only : nmodes,delint_mode

  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
  logical :: classes_talk
  integer :: nr,nz,ir,iz,i,k,iperp,nperp,ierr,nprof,ienerg,nenerg
  double precision :: rbeg,hr,zbeg,hz,weight,psi,psipow
  double precision :: bmod,phi_elec,phi_elec_min,phi_elec_max
  double precision :: toten_min,toten_max,thermen_max,toten_range
  double precision :: dens,omdens,trapez_fac,perpinv_max
  double precision :: totresflux,xint_using_delta
  double precision :: xjperp,xenerg,totxint
  double precision, dimension(:), allocatable :: totdelint
!
  numbasef=0 !no extra integrals sampled, pure orbit integration
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
  nperp=5000 !100
  nenerg=60 !40
!
  totresflux=0.d0
  allocate(totdelint(nmodes))
!
!  do ienerg=1,nenerg
  do ienerg= 10,10 !20,20 !<=fix energy for debugging
    xenerg=(dble(ienerg)-0.5d0)/dble(nenerg)
    toten=toten_min+toten_range*xenerg
print *,'toten = ',toten
!
    call find_Jperpmax(perpinv_max)
!
open(1901,file='subintegrand94_vsJperp.dat')
    do iperp=1,nperp
!    do iperp=1,1 !<=fix Jperp for debugging
      if(iperp.eq.nperp) then
        trapez_fac=0.5d0
      else
        trapez_fac=1.d0
      endif
      xjperp=dble(iperp)/dble(nperp)
      trapez_fac=trapez_fac*2.d0*xjperp
      perpinv=perpinv_max*(1.d0-xjperp**2)
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
      totxint=0.d0
      totdelint=0.d0
!
      do iclass=1,nclasses
!
        call sample_class_doublecount(1000+iclass,ierr)
!
        if(ierr.eq.0) then
!
          call integrate_class_resonances(xint_using_delta)
!
          totxint=totxint+xint_using_delta
          totdelint=totdelint+delint_mode
        else
          print *,'sample_class_doublecount error',ierr
        endif
      enddo
!
      totresflux=totresflux+perpinv_max*trapez_fac*totxint
write(1901,*) perpinv,totxint,nclasses,totdelint
print *,'perpinv:',iperp,'/',nperp,' toten:',ienerg,'/',nenerg
    enddo
close(1901)
  enddo
!
  totresflux=0.25d0*sqrt(pi)*pi*toten_range*totresflux/dble(nperp*nenerg)
!
  print *,'resonant flux = ',totresflux
!
  end subroutine resonant_flux
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
