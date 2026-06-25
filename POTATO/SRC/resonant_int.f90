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
    !$omp threadprivate(twopim2,rm3,taub_new,delphi_new)
! Energy and perpendicular invariant of the resonant orbit, set by pertham and
! read by velo_res for the perturbed-Hamiltonian Fourier integral.  Held separate
! from the global class invariants toten,perpinv so the per-mode resonance loop
! never writes those -- they stay the shared, read-only class values, which lets
! the loop be parallelized without making toten,perpinv threadprivate.
    double precision :: toten_orb,perpinv_orb
    !$omp threadprivate(toten_orb,perpinv_orb)
! By-products of the resonance condition get_rescond ($\psi^\ast$, bounce time,
! toroidal shift), recovered for the orbit weight.  get_rescond is an internal
! procedure of integrate_class_resonances that is also passed as the dummy root
! function to find_all_roots; with gfortran's trampoline for that dummy-procedure
! call, a host-local PRIVATE variable written through it is not seen as private,
! so these must be threadprivate module state instead of privatized host locals.
    double precision :: psiast_res,taub_res,delphi_res
    !$omp threadprivate(psiast_res,taub_res,delphi_res)
    integer, dimension(:), allocatable :: marr,narr
    double precision, dimension(:), allocatable :: delint_mode
!
    type(respoints_fix_jperp_mode_class), dimension(:,:), allocatable :: respoints_jp
    type(respoints_fix_jperp),            dimension(:),   allocatable :: respoints_all, &
                                                                         respoints_all_tmp
    type(respoint_single),                dimension(:),   allocatable :: respoint
! Per-energy-slice resonance bookkeeping: the energy loop (resonant_torque) runs
! slices in parallel, and each slice fills its own J_perp grid of resonant points
! (respoints_all/_tmp, respoints_jp, respoint) and per-mode integral (delint_mode),
! and sizes them with its own nperp_max.  Threadprivate so slices do not collide.
! nmodes,marr,narr are set once before the loop and read only, so they stay shared.
    !$omp threadprivate(nperp_max,delint_mode,respoints_jp,respoints_all, &
    !$omp               respoints_all_tmp,respoint)
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
  use resint_mod,        only : twopim2,rm3,taub_new,delphi_new,toten_orb,perpinv_orb
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
  comfac=(2.d0*(toten_orb-phi_elec)/bmod-perpinv_orb)*bmod_n &
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
  use global_invariants, only : dtau
  use resint_mod,        only : taub_new,delphi_new,toten_orb,perpinv_orb
!
!
  implicit none
!
  complex(8), parameter :: imun=(0.d0,1.d0)
!
  integer :: ierr
  double precision :: absHn2,bmod,phi_elec,taub,delphi
  double precision, dimension(neqm) :: z
  double precision, dimension(:), allocatable :: extraset
!
  external :: velo,velo_res
!
  call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
  toten_orb=z(4)**2+phi_elec
  perpinv_orb=z(4)**2*(1.d0-z(5)**2)/bmod
!
  next=0
  allocate(extraset(next))
!
  call find_bounce(next,velo,dtau,z,taub,delphi,extraset,ierr)
  if(ierr.ne.0) then
    absHn2 = 0.d0
    deallocate(extraset)
    return
  endif
!
  taub_new=taub
  delphi_new=delphi
  deallocate(extraset)
!
  next=3
  allocate(extraset(next))
  extraset=0.d0
!
  call find_bounce(next,velo_res,dtau,z,taub,delphi,extraset,ierr)
  if(ierr.ne.0) then
    absHn2 = 0.d0
    deallocate(extraset)
    return
  endif
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
  use resint_mod,         only : nmodes,marr,narr,twopim2,rm3,delint_mode,respoints_jp, &
                                 psiast_res,taub_res,delphi_res
  use orbit_dim_mod,      only : neqm
  use global_invariants,  only : toten,perpinv,cE_ref,Phi_eff
  use sample_matrix_mod,  only : npoi,xarr
  use logging_mod,        only : tee_message
  use interp_cache_mod,   only : interp_cache_reset
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0, twopi=2.d0*pi, &
                                 pi32_over4m=-0.25d0*pi**1.5d0
!
  integer          :: mode,iroot,ierr
  double precision :: relmargin_loc,widthclass,xbeg,xend
  double precision :: rescond,dresconddx,dpsiastdx
  double precision :: one_res,sigma,delta_R,Rst,xi,dxi_dx,dpsiast_dRst,absHn2
  double precision :: fmaxw,A1ast,A2ast
  double precision :: dens,temp,ddens,dtemp,phi_elec,dPhi_dpsi
  double precision, dimension(neqm) :: z
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
! Parallel over modes: each mode is independent.  It runs its own find_all_roots
! (nroots,roots threadprivate) and per-root starter+pertham (twopim2,rm3,taub_new,
! delphi_new,toten_orb,perpinv_orb,next, and the interp cache threadprivate), and
! writes only its disjoint column respoints_jp(mode,iclass) and delint_mode(mode).
! toten,perpinv stay shared read-only (the class invariants); customgrid,ncustom,
! xcustom and the grid arrays are shared read-only too.  The threadprivate
! form-class bounds it reads (ifuntype,R_class_beg,R_class_end,sigma_class) are
! copyin'd.  The get_rescond by-products psiast_res,taub_res,delphi_res are
! threadprivate module state (resint_mod), NOT host-local privates: get_rescond
! is passed as the dummy root function to find_all_roots, and gfortran's
! trampoline for that call does not honor a host-local private.  reslines
! write(31415) and tee_message go in a critical section.  Each thread resets its
! interp cache at entry so it drops the previous class's grid before reusing it.
  !$omp parallel default(shared) &
  !$omp   private(mode,iroot,ierr,rescond,dresconddx,dpsiastdx) &
  !$omp   private(one_res,Rst,xi,dxi_dx,dpsiast_dRst,absHn2,fmaxw,A1ast,A2ast) &
  !$omp   private(dens,temp,ddens,dtemp,phi_elec,dPhi_dpsi,z) &
  !$omp   copyin(ifuntype,R_class_beg,R_class_end,sigma_class)
  call interp_cache_reset
  !$omp do schedule(dynamic)
  do mode=1,nmodes
    twopim2=twopi*dble(marr(mode))
    rm3=dble(narr(mode))
    delint_mode(mode)=0.d0
!
    call find_all_roots(get_rescond,xbeg,xend,ierr)
!
    if(ierr.ne.0) then
      !$omp critical (reslines_log)
      call tee_message( &
        'integrate_class_resonances: error in find_all_roots')
      !$omp end critical (reslines_log)
      respoints_jp(mode,iclass)%nrespoi=0
      respoints_jp(mode,iclass)%toten_res=toten
      respoints_jp(mode,iclass)%perpinv_res=perpinv
      cycle
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
                               psiast_res,dpsiast_dRst,z,ierr)
!
      if(ierr.ne.0) then
        !$omp critical (reslines_log)
        call tee_message( &
          'integrate_class_resonances: error in starter_doublecount')
        !$omp end critical (reslines_log)
        cycle
      endif
!
      if(.true.) then
        !$omp critical (reslines_log)
        write(31415,*) toten,perpinv,psiast_res,marr(mode),narr(mode)   !<=resonant line for plotting
        !$omp end critical (reslines_log)
      endif
!
      respoints_jp(mode,iclass)%z_res(:,iroot)=z
!
      dpsiastdx=dpsiast_dRst*delta_R*dxi_dx     !$\difp{\psi^\ast}{x}$
!
      call pertham(z,absHn2)
      call equilmaxw(psiast_res,fmaxw)
      call denstemp_of_psi(psiast_res,dens,temp,ddens,dtemp)
      call phielec_of_psi(psiast_res,phi_elec,dPhi_dpsi)
!
! toten,perpinv are never written in this loop (pertham writes toten_orb,
! perpinv_orb instead), so they still hold the class invariants set on entry --
! no save/restore needed.
!
! Non-local thermodynamic forces Eq.(95) (former Eq.(87)):
      A2ast=dtemp/temp
      A1ast=ddens/dens+dPhi_dpsi/temp-1.5d0*A2ast
!
! Expression under summation signs except the last line in Eq.(104) (former Eq.(94)):
      one_res=abs(dpsiastdx/dresconddx)*absHn2*fmaxw            &
!ERROR, SEE in RED=>             *(Phi_eff*taub_res*(A1ast+A2ast*(toten-phi_elec)/temp)+delphi_res/temp)
             *Phi_eff*taub_res*(A1ast+A2ast*(toten-phi_elec)/temp)  !<=ERROR CORRECTED
!
! emulator of box average (Heaviside function replaced with one in (104) - result is integral torque in
! the whole volume normalized by the reference energy $\cE_{ref}$):
      one_res=one_res*taub_res
! end emulator of box average
!
! multiply expression under summation over modes with toroidal mode number, with factor $-\pi^{3/2}/4$
! and with reference energy:
      one_res=one_res*rm3*pi32_over4m*cE_ref
!
      respoints_jp(mode,iclass)%w_res(iroot)=one_res
      respoints_jp(mode,iclass)%taub(iroot)=taub_res
      delint_mode(mode)=delint_mode(mode)+one_res
    enddo
  enddo
  !$omp end do
  !$omp end parallel
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
  psiast_res=vec(1)           ! $\psi^\ast$
  taub_res=vec(2)             ! $\tau_b$
  delphi_res=vec(3)           ! $\Delta\varphi_b$
  rescond=delphi_res+twopim2/rm3
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
  use logging_mod,                  only : tee_message, close_logging
  use bounds_fixpoints_mod,         only : region_set_t
!
  logical :: classes_talk
!
  integer :: ierr,mode
  character(len=256) :: msg
! Region set owned per J_perp node (hence per energy slice, hence per thread).
  type(region_set_t) :: regions
!
  wrbounds=.false.
  dowrite=.false.
  write_orb=.false.
  classes_talk=.false.
!
  perpinv=x
!
  call find_bounds_fixpoints(regions,ierr)
!
  if(ierr.ne.0) then
    write(msg, '(A,I0)') &
      'get_matrix_res: find_bounds_fixpoints ierr = ', ierr
    call tee_message(trim(msg))
    call close_logging()
    stop
  endif
!
  call form_classes_doublecount(regions,classes_talk,ierr)
!
  if(ierr.ne.0) then
    write(msg, '(A,I0)') &
      'get_matrix_res: form_classes ierr = ', ierr
    call tee_message(trim(msg))
    call close_logging()
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
      write(msg, '(A,I0)') &
        'get_matrix_res: sample_class_doublecount error ', ierr
      call tee_message(trim(msg))
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
  use field_sub, only : psif
  use field_eq_mod, only : nrad,nzet,rad,zet,psi_sep
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
  use potato_input_mod,  only : nbox, nenerg_input => nenerg, &
                                thermen_max_input => thermen_max, &
                                enkin_min_over_temp, &
                                adaptive_jperp, npoi_init, nlagr_sampling, &
                                eps_sampling, itermax_sampling
  use logging_mod,       only : tee_message
  !$ use omp_lib, only : omp_set_max_active_levels

  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
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
  character(len=256) :: msg
  double precision, dimension(:), allocatable :: torque_int_modes
  double precision, dimension(:), allocatable :: sbox
  double precision, dimension(:), allocatable :: taubox
  double precision, dimension(:), allocatable :: torquebox
! Per-resonance-point box times, filled by the parallel time_in_box pass and
! consumed by the serial accumulate/write pass below.
  double precision, dimension(:,:), allocatable :: taubox_all
!
  external :: get_matrix_res
!
  allocate(sbox(nbox), taubox(nbox), torquebox(nbox))
!
! n1,n2,nlagr,eps,itermax (sample_matrix_out_mod) are threadprivate per-slice grid
! config; they are set inside the energy loop so each thread seeds its own.
!
  numbasef=0 !no extra integrals sampled, pure orbit integration
  call linspace(1d0/nbox, 1d0, nbox, sbox)
!
! Find minimum and maximum values of electrostatic potential in the computation domain:
!
  call find_Phiminmax(phi_elec_min,phi_elec_max)
!
  thermen_max=thermen_max_input
!
  call denstemp_of_psi(psimagaxis, dens, temp, ddens, dtemp)
!
  thermen_max=thermen_max*temp  !maximum kinetic energy in units of reference energy
!
! Energy integration limits.  Start the toten grid so the lowest slice keeps at
! least enkin_min_over_temp*T of kinetic energy at the deepest potential well
! (phi_elec_max), cutting off the cold near-boundary band that otherwise breeds
! spurious low-energy resonances.
  toten_min=phi_elec_max+enkin_min_over_temp*temp
  toten_max=thermen_max+phi_elec_max
  toten_range=toten_max-toten_min
!
  write(msg, '(A,ES14.6)') &
    'maximum kinetic energy = ', thermen_max
  call tee_message(trim(msg))
  write(msg, '(A,ES14.6)') &
    'minimum potential energy = ', phi_elec_min
  call tee_message(trim(msg))
  write(msg, '(A,ES14.6)') &
    'maximum potential energy = ', phi_elec_max
  call tee_message(trim(msg))
  write(msg, '(A,ES14.6)') &
    'minimum total energy = ', toten_min
  call tee_message(trim(msg))
  write(msg, '(A,ES14.6)') &
    'maximum total energy = ', toten_max
  call tee_message(trim(msg))
!
  nenerg=nenerg_input
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
! Coarse-grain parallelism: the total-energy slices are independent, so run them
! in parallel.  All per-slice scratch (the class and J_perp grids, the resonant-
! point bookkeeping, the orbit invariants toten,perpinv,sigma,dtau) is
! threadprivate, so each slice keeps its own; the torque sums reduce.  Cap the
! active level at one so the inner per-class/per-mode regions stay serial inside
! each slice instead of oversubscribing.  File writes go to libgfortran's
! per-unit lock; resonance lines are compared sorted, so their interleaving is
! immaterial.
  !$ call omp_set_max_active_levels(1)
  !$omp parallel do default(shared) schedule(dynamic) &
  !$omp   private(xenerg,perpinv_max,trapez_fac,nrespoints,i,k,iperp,ierr,nperp) &
  !$omp   private(time_beg,time_end,xjperp,torque_int_loc,msg,taubox_all) &
  !$omp   reduction(+:torque_int,torque_int_modes,torquebox)
  do ienerg=1,nenerg
    xenerg=(dble(ienerg)-0.5d0)/dble(nenerg)
    toten=toten_min+toten_range*xenerg
! per-slice J_perp grid config (sample_matrix_out_mod is threadprivate):
    n1=nmodes
    n2=1
    nlagr=nlagr_sampling
    eps=eps_sampling
    itermax=itermax_sampling
!toten =   -3.1211921097605737d0
    write(msg, '(A,ES22.14)') 'toten = ', toten
    call tee_message(trim(msg))
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
      npoi=npoi_init            !initial equidistant grid size over J_perp
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
! Box counter (time_in_box) integrates each resonant orbit through the radial
! boxes via dvode -- the dominant per-point cost.  Points are independent and
! dvode_f90_m is threadprivate/thread-safe, so run them in parallel into the
! per-point taubox_all.  z_res(1:5) is the orbit start at the Poincare cut.
      allocate(taubox_all(nbox,nrespoints))
      !$omp parallel do default(shared) private(i) schedule(dynamic)
      do i=1,nrespoints
        call time_in_box(respoint(i)%z_res(1:5), nbox, sbox, &
          respoint(i)%taub, taubox_all(:,i))
      enddo
      !$omp end parallel do
!
! Serial accumulate/write pass: torque_int_loc is a running prefix sum written
! per row to fort.1901, so it stays serial and in order; torquebox sums the
! per-point box times.  Width of the energy interval is already in the weight.
      do i=1,nrespoints
        torque_int_loc=torque_int_loc+respoint(i)%w_res
        write(1901,*) respoint(i)%toten_res, &
          respoint(i)%perpinv_res, &
          torque_int_loc
        torquebox=torquebox+respoint(i)%w_res*taubox_all(:,i)/respoint(i)%taub
      enddo
      deallocate(taubox_all)
      write(1901,*) ' '
!
      deallocate(respoint)
      write(msg, '(A,I0)') &
        'number of J_perp points = ', npoi
      call tee_message(trim(msg))
!
! Sum up contributions of energy levels:
      torque_int=torque_int+torque_int_loc
      write(msg, '(A,ES22.14)') &
        'method 1, torque_int_loc = ', torque_int_loc
      call tee_message(trim(msg))
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
      write(msg, '(A,ES22.14)') &
        'method 2, torque_int_loc = ', torque_int_loc
      call tee_message(trim(msg))
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
        write(msg, '(A,I0,A,I0,A,I0,A,I0)') &
          'perpinv:', iperp, '/', nperp, &
          ' toten:', ienerg, '/', nenerg
        call tee_message(trim(msg))
      enddo
    endif
!
    call cpu_time(time_end)
    write(msg, '(A,I0,A,I0,A,F10.2,A)') &
      ' toten:', ienerg, '/', nenerg, &
      ' cpu time = ', time_end-time_beg, ' sec'
    call tee_message(trim(msg))
!
  enddo
  !$omp end parallel do
!
  if(adaptive_jperp) close(1901)
  close(1902)
!
! Integral torque:
  write(msg, '(A,ES22.14)') &
    'resonant torque  = ', torque_int
  call tee_message(trim(msg))
  open(1,file='integral_torque.dat')
  write(1,*) torque_int
  close(1)
!
!Write box-counted integral torque:
  if(adaptive_jperp) then
    open(1,file='boxcounted_torque.dat')
    write(1,*) torquebox(1), 0d0, sbox(1)
    do i=2,nbox
      write(1,*) torquebox(i), sbox(i-1), sbox(i)
    enddo
    close(1)
  endif
!
  end subroutine resonant_torque
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
