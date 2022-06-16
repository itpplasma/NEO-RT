!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! modules
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module resint_mod
!
    type respoints_fix_jperp_mode_class
      integer :: nrespoi
      double precision :: toten_res,perpinv_res
      double precision, dimension(:),   allocatable :: w_res
      double precision, dimension(:,:), allocatable :: z_res
      double precision, dimension(:),   allocatable :: taub
    end type
!
    type respoint_single
      double precision :: toten_res,perpinv_res,w_res
      double precision, dimension(5) :: z_res
      double precision :: taub
    end type
!
    type respoints_fix_jperp
      integer :: nclasses
      type(respoints_fix_jperp_mode_class), dimension(:,:), allocatable :: respoints_jp
    end type
!
    integer          :: nmodes,nperp_max=0
    double precision :: twopim2,rm3,taub_new,delphi_new
    integer, dimension(:), allocatable :: marr,narr
    double precision, dimension(:), allocatable :: delint_mode
!
    type(respoints_fix_jperp_mode_class), dimension(:,:), allocatable :: respoints_jp
    type(respoints_fix_jperp),            dimension(:),   allocatable :: respoints_all, &
                                                                         respoints_all_tmp
    type(respoint_single),                dimension(:),   allocatable :: respoint
  end module resint_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! rotines
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine velo_res(dtau,z,vz)
!
! Computes the extended RHS for equations of motion and Fourier amplitude
! of the perturbed Hamiltonian $\hat H_\bm$, Eq.(102) (former Eq.(93))
!
  use orbit_dim_mod,     only : neqm,next,numbasef
  use global_invariants, only : toten,perpinv
  use resint_mod,        only : twopim2,rm3,taub_new,delphi_new
!
  implicit none
!
  complex(8), parameter :: imun=(0.d0,1.d0)
!
  double precision :: dtau,bmod,phi_elec
  double precision, dimension(neqm+next) :: z,vz
  complex(8) :: bmod_n,comfac
!
  call velo(dtau,z(1:neqm),vz(1:neqm))
  call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
  call bmod_pert(z(1),z(3),bmod_n)
!
  comfac=(2.d0*(toten-phi_elec)/bmod-perpinv)*bmod_n &
        *exp(imun*(rm3*z(2)-(twopim2+delphi_new*rm3)*z(6)/taub_new))
!
  vz(6)=1.d0
  vz(7)=real(comfac)
  vz(8)=aimag(comfac)
!
  end subroutine velo_res
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine pertham(z,absHn2)
!
! Computes modulus squared of the Fourier amplitude of the normalized perturbed
! Hamiltoninan, $|\hat H_\bm|^2$, with $\hat H_\bm$ defined by Eq.(102) (former Eq.(93)).
!
  use orbit_dim_mod,     only : neqm,next     ,write_orb,iunit1
  use global_invariants, only : toten,perpinv,dtau
  use resint_mod,        only : taub_new,delphi_new
!
!
  implicit none
!
  complex(8), parameter :: imun=(0.d0,1.d0)
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
  subroutine integrate_class_resonances
!
! Computes sum over resonances $x=x^{res}_{(\bm,k)}$ in Eq.(104) for a given class $k$
!
  use find_all_roots_mod, only : customgrid,ncustom,xcustom,nroots,roots
  use get_matrix_mod,     only : relmargin,iclass
  use form_classes_doublecount_mod, only : ifuntype,R_class_beg,R_class_end,sigma_class
  use resint_mod,         only : nmodes,marr,narr,twopim2,rm3,delint_mode,respoints_jp
  use orbit_dim_mod,      only : neqm
  use global_invariants,  only : toten,perpinv,cE_ref,Phi_eff
  use sample_matrix_mod,  only : npoi,xarr
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0, twopi=2.d0*pi, &
                                 pi32_over4m=-0.25d0*pi**1.5d0
!
  integer          :: mode,iroot,ierr
  double precision :: relmargin_loc,widthclass,xbeg,xend
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
!old=>  relmargin_loc=1.d-8
!old=>  widthclass=1.d0
  relmargin_loc=relmargin                                      !<=new
  widthclass=abs(R_class_end(iclass)/R_class_beg(iclass)-1.d0) !<=new
!
  call classbounds(ifuntype(iclass),relmargin_loc,widthclass,xbeg,xend)
!
  customgrid=.true.
  ncustom=npoi
  allocate(xcustom(ncustom))
  xcustom=xarr
!
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
      customgrid=.false.
      deallocate(xcustom)
      return
    endif
!
    respoints_jp(mode,iclass)%nrespoi=nroots
    respoints_jp(mode,iclass)%toten_res=toten
    respoints_jp(mode,iclass)%perpinv_res=perpinv
    if(nroots.eq.0) cycle
    allocate(respoints_jp(mode,iclass)%w_res(nroots), &
             respoints_jp(mode,iclass)%z_res(5,nroots), &
             respoints_jp(mode,iclass)%taub(nroots))
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
      if(.true.) then
        write(31415,*) toten,perpinv,psiast,marr(mode),narr(mode)   !<=resonant line for plotting
      endif
!
      respoints_jp(mode,iclass)%z_res(:,iroot)=z
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
! Non-local thermodynamic forces Eq.(95) (former Eq.(87)):
      A2ast=dtemp/temp
      A1ast=ddens/dens+dPhi_dpsi/temp-1.5d0*A2ast
!
! Expression under summation signs except the last line in Eq.(104) (former Eq.(94)):
      one_res=abs(dpsiastdx/dresconddx)*absHn2*fmaxw            &
!ERROR, SEE in RED=>             *(Phi_eff*taub*(A1ast+A2ast*(toten-phi_elec)/temp)+delphi/temp)
             *Phi_eff*taub*(A1ast+A2ast*(toten-phi_elec)/temp)  !<=ERROR CORRECTED
!
! emulator of box average (Heaviside function replaced with one in (104) - result is integral torque in
! the whole volume normalized by the reference energy $\cE_{ref}$):
      one_res=one_res*taub
! end emulator of box average
!
! multiply expression under summation over modes with toroidal mode number, with factor $-\pi^{3/2}/4$
! and with reference energy:
      one_res=one_res*rm3*pi32_over4m*cE_ref
!
      respoints_jp(mode,iclass)%w_res(iroot)=one_res
      respoints_jp(mode,iclass)%taub(iroot)=taub
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
  subroutine get_matrix_res
!
  use sample_matrix_out_mod,        only : n1,n2,x,amat,icount
  use global_invariants,            only : toten,perpinv
  use form_classes_doublecount_mod, only : nclasses
  use get_matrix_mod,               only : iclass
  use resint_mod,                   only : nmodes,nperp_max,respoints_jp, &
                                           respoints_all,respoints_all_tmp
  use cc_mod,                       only : wrbounds,dowrite
  use orbit_dim_mod,                only : write_orb
!
  logical :: classes_talk
!
  integer :: ierr,mode
!
  wrbounds=.false.
  dowrite=.false.
  write_orb=.false.
  classes_talk=.false.
!
  perpinv=x
!
  call find_bounds_fixpoints(ierr)
!
  if(ierr.ne.0) then
    print *,'get_matrix_res: find_bounds_fixpoints ierr = ',ierr
    stop
  endif
!
  call form_classes_doublecount(classes_talk,ierr)
!
  if(ierr.ne.0) then
    print *,'get_matrix_res: form_classes ierr = ',ierr
    stop
  endif
!
  allocate(respoints_jp(nmodes,nclasses))
  amat=0.d0
!
  do iclass=1,nclasses
!
    call sample_class_doublecount(1,ierr)
!
    if(ierr.eq.0) then
!
      call integrate_class_resonances
!
      do mode=1,nmodes
! amat(:,1) - sum over classes and resonances within each class:
        if(respoints_jp(mode,iclass)%nrespoi.gt.0) then
          amat(mode,1)=amat(mode,1)+sum(respoints_jp(mode,iclass)%w_res)
        endif
      enddo
    else
      print *,'get_matrix_res: sample_class_doublecount error',ierr
      do mode=1,nmodes
        respoints_jp(mode,iclass)%nrespoi=0
        respoints_jp(mode,iclass)%toten_res=toten
        respoints_jp(mode,iclass)%perpinv_res=perpinv
      enddo
    endif
  enddo
!
  icount=icount+1
  if(icount.gt.nperp_max) then
    if(nperp_max.le.0) then
      nperp_max=1
      if(allocated(respoints_all)) deallocate(respoints_all)
      allocate(respoints_all(nperp_max))
    else
      allocate(respoints_all_tmp(nperp_max))
      respoints_all_tmp=respoints_all
      deallocate(respoints_all)
      allocate(respoints_all(2*nperp_max))
      respoints_all(1:nperp_max)=respoints_all_tmp
      deallocate(respoints_all_tmp)
      nperp_max=2*nperp_max
    endif
  endif
  respoints_all(icount)%nclasses=nclasses
  allocate(respoints_all(icount)%respoints_jp(nmodes,nclasses))
  respoints_all(icount)%respoints_jp=respoints_jp
  deallocate(respoints_jp)
!print *,icount  !uncomment to see some life
!
  end subroutine get_matrix_res
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine resonant_torque
!
!
  use field_eq_mod,      only : psif,nrad,nzet,rad,zet,psi_sep
  use poicut_mod,        only : rmagaxis,zmagaxis,psimagaxis,psi_bou,rhopol_bou
  use global_invariants, only : toten,perpinv
  use poicut_mod,        only : Rbou_lfs,Zbou_lfs
  use get_matrix_mod,    only : iclass
  use form_classes_doublecount_mod, only : nclasses
  use orbit_dim_mod,     only : numbasef
  use resint_mod,        only : nmodes,delint_mode,respoints_jp,respoints_all,nperp_max, &
                                respoints_all_tmp,respoint
  use sample_matrix_out_mod, only : nlagr,n1,n2,npoi,itermax,x,amat,icount,xbeg,xend,eps, &
                                    ind_hist,xarr,amat_arr

  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
  logical :: adaptive_jperp
  integer :: nr,nz,ir,iz,i,k,iperp,nperp,ierr,nprof,ienerg,nenerg
  integer :: nrespoints
  double precision :: rbeg,hr,zbeg,hz,weight,psi,psipow
  double precision :: bmod,phi_elec,phi_elec_min,phi_elec_max
  double precision :: toten_min,toten_max,thermen_max,toten_range
  double precision :: omdens,trapez_fac,perpinv_max
  double precision :: torque_int,torque_int_loc
  double precision :: xjperp,xenerg,totxint,step_energ
  double precision :: time_beg,time_end
  double precision :: dens, temp, ddens, dtemp
  double precision, dimension(:), allocatable :: torque_int_modes

  ! Number of boxes in radius to localize torque contributions
  integer, parameter :: nbox = 100
  double precision :: sbox(nbox)
  double precision :: taubox(nbox)
!
!Distribution of integral torque over boxes:
  double precision, dimension(nbox) :: torquebox
!
  external :: get_matrix_res
!
  adaptive_jperp=.true.    !use adaptive integration over J_perp
!  adaptive_jperp=.false.   !use non-adaptive integration over J_perp
!
! size of result matrix in get_matrix_res:
  n1=nmodes
  n2=1
! Lagrange polynomial order for sampling adaptive grid:
  nlagr=3
! Error of sampling:
  eps=1d-2
! Maximum number of refinements:
  itermax=5
!
  numbasef=0 !no extra integrals sampled, pure orbit integration
  call linspace(1d0/nbox, 1d0, nbox, sbox)
!
! Find minimum and maximum values of electrostatic potential in the computation domain:
!
  call find_Phiminmax(phi_elec_min,phi_elec_max)
!
  thermen_max=6.d0   !maximum kinetic energy in units of temperature on the axis (used for upper energy integration limit)
!
  call denstemp_of_psi(psimagaxis, dens, temp, ddens, dtemp)
!
  thermen_max=thermen_max*temp  !maximum kinetic energy in units of reference energy
!
! Energy integration limits:
  toten_min=phi_elec_min
  toten_max=thermen_max+phi_elec_max
  toten_range=toten_max-toten_min
!
  print *,'maximum kinetic energy = ',thermen_max
  print *,'miminum potential energy = ',phi_elec_min
  print *,'maximum potential energy = ',phi_elec_max
  print *,'miminum total energy = ',toten_min
  print *,'maximum total energy = ',toten_max
!
  nenerg=60 !40 !grid size for energy integration (number of energy levels)
!
  torque_int=0.d0
  allocate(torque_int_modes(nmodes))
  torque_int_modes=0.d0
!
  torquebox=0.d0
!
  step_energ=toten_range/dble(nenerg) !integration step over total energy
!step_energ=0.22521463755624047d0
!
  if(adaptive_jperp) then
    open(1901,file='subint_ofH0int_104_vsJperp_fromresp.dat')
    open(1902,file='subint_ofH0int_104_vsJperp_adapt.dat')
  else
    open(1902,file='subint_ofH0int_104_vsJperp_equi.dat')
  endif
!
!  do ienerg=1,nenerg
  do ienerg=2,nenerg
!  do ienerg= 20,20 !10,10 !20,20 !<=fix energy for debugging
    xenerg=(dble(ienerg)-0.5d0)/dble(nenerg)
    toten=toten_min+toten_range*xenerg
!toten =   -3.1211921097605737d0
    print *,'toten = ',toten
!
    call cpu_time(time_beg)
!
! find maximum possible value of J_perp for given total energy:
!
    call find_Jperpmax(perpinv_max)
!
    if(adaptive_jperp) then
!
! New, adaptive integration:
!
      nperp_max=0   !initial size of respoints_all(nperp_max) - will be increased by sample_matrix_out
      xbeg=0.d0                 !lower integration limit over normalized J_perp
      xend=perpinv_max*0.9999d0 !upper integration limit over normalized J_perp
      npoi=51                   !initial equidistant grid size over J_perp
      icount=0                  !initialize the counter of resonant points
!
! Generate J_perp grid with function values:
!
      call sample_matrix_out(get_matrix_res,ierr)
!
      allocate(respoints_all_tmp(npoi))
!
! reorder J_perp points to the increasing sequence:
      do i=1,npoi
        respoints_all_tmp(i)=respoints_all(ind_hist(i))
      enddo
!
      deallocate(respoints_all)
      allocate(respoints_all(npoi))
      respoints_all=respoints_all_tmp
      deallocate(respoints_all_tmp)
!
! count the number of resonant points and update their weights in accordance
! with trapezoidal integration rule over J_perp:
      nrespoints=0
      do iperp=1,npoi
        if(iperp.eq.1) then
          trapez_fac=0.5d0*(xarr(2)-xarr(1))
        elseif(iperp.eq.npoi) then
          trapez_fac=0.5d0*(xarr(npoi)-xarr(npoi-1))
        else
          trapez_fac=0.5d0*(xarr(iperp+1)-xarr(iperp-1))
        endif
        do iclass=1,respoints_all(iperp)%nclasses
          do i=1,nmodes
            if(respoints_all(iperp)%respoints_jp(i,iclass)%nrespoi.gt.0) then
              respoints_all(iperp)%respoints_jp(i,iclass)%w_res = &
              respoints_all(iperp)%respoints_jp(i,iclass)%w_res*trapez_fac
              nrespoints=nrespoints &
                        +respoints_all(iperp)%respoints_jp(i,iclass)%nrespoi
            endif
          enddo
        enddo
      enddo
!
      allocate(respoint(nrespoints))
!
! extract separate resonant points and update their weight in accorance with
! integration over total energy:
      nrespoints=0
      do iperp=1,npoi
        do iclass=1,respoints_all(iperp)%nclasses
          do i=1,nmodes
            do k=1,respoints_all(iperp)%respoints_jp(i,iclass)%nrespoi
              nrespoints=nrespoints+1
              respoint(nrespoints)%toten_res                             &
                =respoints_all(iperp)%respoints_jp(i,iclass)%toten_res
              respoint(nrespoints)%perpinv_res                           &
                =respoints_all(iperp)%respoints_jp(i,iclass)%perpinv_res
              respoint(nrespoints)%w_res                                 &
                =respoints_all(iperp)%respoints_jp(i,iclass)%w_res(k)    &
                *step_energ
              respoint(nrespoints)%z_res &
                =respoints_all(iperp)%respoints_jp(i,iclass)%z_res(:,k)
              respoint(nrespoints)%taub &
                  =respoints_all(iperp)%respoints_jp(i,iclass)%taub(k)
            enddo
          enddo
        enddo
      enddo
!
! Computation of the integral torque:
      torque_int_loc=0.d0
!
      do i=1,nrespoints
! Here box counter should be used. Below we compute an integral torque as a sum of orbit weights:
        torque_int_loc=torque_int_loc+respoint(i)%w_res
! Quantity torque_int_loc is a contribution of one energy level to the integral torque. Width of energy
! interval is already included in the weight.
! For box counter use the starting coordinates of the orbit at the Poincare cut. They are stored in
! respoint(i)%z_res(1:5)
! Weight of the orbit is stored in
! respoint(i)%w_res
! If needed, there are also total energy and perpendicular invariant available in
! respoint(i)%toten_res
! and
! respoint(i)%perpinv_res
! respectively
        call time_in_box(respoint(i)%z_res(1:5), nbox, sbox, &
          respoint(i)%taub, taubox)
        write(1901,*) respoint(i)%toten_res, &
          respoint(i)%perpinv_res, &
          ! tormom_of_RZ(respoint(i)%toten_res, respoint(i)%perpinv_res, TODO ), &
          torque_int_loc    !,taubox/respoint(i)%taub
!
        torquebox=torquebox+respoint(i)%w_res*taubox/respoint(i)%taub   !<=sum up resonances in boxes
!
      enddo
      write(1901,*) ' '
!
      deallocate(respoint)
      print *,'number of J_perp points = ',npoi
!
! Sum up contributions of energy levels:
      torque_int=torque_int+torque_int_loc
      print *,'method 1, torque_int_loc = ',torque_int_loc
!
! Alternative computation via sampling matrix. Here distribution over resonances is also computed:
      do i=1,npoi
!        write(10000,*) xarr(i),amat_arr(:,1,i)
        if(i.eq.1) then
          torque_int_loc=0.d0
        else
          torque_int_loc=torque_int_loc &
                    +0.5d0*(sum(amat_arr(:,1,i))+sum(amat_arr(:,1,i-1))) &
                    *(xarr(i)-xarr(i-1))*step_energ
          torque_int_modes=torque_int_modes &
                    +0.5d0*(amat_arr(:,1,i)+amat_arr(:,1,i-1)) &
                    *(xarr(i)-xarr(i-1))*step_energ
        endif
        write(1902,*) xarr(i),torque_int_loc,torque_int_modes
      enddo
      write(1902,*) ' '
      print *,'method 2, torque_int_loc = ',torque_int_loc
!
    else
!
! Old, non-adaptive integration:
! Not prepared for box counting, use for testing only.
      nperp=2500 !5000 !100    !size of the integration grid over normalized J_perp
      if(.not.allocated(amat)) allocate(amat(n1,n2))
      icount=0
      nperp_max=0
      do iperp=1,nperp
        if(iperp.eq.nperp) then
          trapez_fac=0.5d0
        else
          trapez_fac=1.d0
        endif
        xjperp=dble(iperp)/dble(nperp)
        trapez_fac=trapez_fac*2.d0*xjperp/dble(nperp)
        x=perpinv_max*(1.d0-xjperp**2)
!
        call get_matrix_res
!
        torque_int=torque_int+perpinv_max*trapez_fac*sum(amat(:,1))*step_energ
        torque_int_modes=torque_int_modes+perpinv_max*trapez_fac*amat(:,1)*step_energ
! Subintegrand of dimensional integral over energy for a total torque as function of J_perp:
        write(1902,*) perpinv,torque_int,torque_int_modes
        print *,'perpinv:',iperp,'/',nperp,' toten:',ienerg,'/',nenerg
      enddo
    endif
!
    call cpu_time(time_end)
    print *,' toten:',ienerg,'/',nenerg,' cpu time = ',time_end-time_beg,' sec'
!
  enddo
!
  if(adaptive_jperp) close(1901)
  close(1902)
!
! Integral torque:
  print *,'resonant torque  = ',torque_int
  open(1,file='integral_torque.dat')
  write(1,*) torque_int
  close(1)
!
!Write box-counted integral torque:
  if(adaptive_jperp) then
    open(1,file='boxcounted_torque.dat')
    write(1,*) torquebox(1), 0d0, sbox(1)
    do i=2,nbox
      write(1,*) torquebox(i), sbox(i), sbox(i+1)
    enddo
    close(1)
  endif
!
  end subroutine resonant_torque
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
