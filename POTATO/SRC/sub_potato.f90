!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Modules:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module phielec_of_psi_mod
    integer, parameter :: npolyphi=9
    double precision, dimension(0:npolyphi) :: polyphi
    double precision, dimension(0:npolyphi) :: polytemp
    double precision, dimension(0:npolyphi) :: polydens
  end module phielec_of_psi_mod
!
  module cc_mod
    logical :: dowrite=.false. !write canonical frequencies and bounce integrals for sampling grid and interpolation
    logical :: wrbounds=.false. !write vpar^2 and psi^* curves vs cut parameter R_c, extremum and boundary points
! These controls are set by get_matrix_res for each parallel energy slice.
! Keeping them per-thread prevents a diagnostic flag write in one slice from
! racing with a read in another slice.
    !$omp threadprivate(dowrite,wrbounds)
  end module cc_mod
!
  module global_invariants
! integration step, normalized total energy, normalized perpendicular invariant, velocity sign:
    double precision :: dtau,toten,perpinv,sigma
! reference energy $\cE_{ref}$, effective potential $\Phi_{eff}$:
    double precision :: cE_ref,Phi_eff
! Total-energy slices run independently in resonant_torque.  These invariants
! are the mutable per-slice/per-orbit state read by the orbit and class builders.
! Keep the normalization constants shared; they are set once at startup.
    !$omp threadprivate(dtau,toten,perpinv,sigma)
  end module global_invariants
!
  module poicut_mod
! number of points for Poincare cut:
    integer :: npc
! R value at the LFS Poincare cut boundary, grid step size over R (equidistant):
    double precision :: rpc_beg,h_rpc
! coordinates (R,Z) and poloidal flux value at the magnetic axis:
    double precision :: rmagaxis,zmagaxis,psimagaxis
! normalzied poloidal radius and dimensional poloidal flux of the last closed surface limiting computation domain:
    double precision :: rhopol_bou,psi_bou
! coordinates (R,Z) at the LFS and HFS ends of Poincare cut segment located in the computation domain:
    double precision :: Rbou_lfs,Zbou_lfs,Rbou_hfs,Zbou_hfs
! arrays of R and Z coordinates at the Poincare cut:
    double precision, dimension(:), allocatable :: rpc_arr,zpc_arr
  end module poicut_mod
!
  module bounds_fixpoints_mod
    type allowed_region
! flag whether the region contains the orbits partly or fully within rho_pol boundary
      logical :: within_rhopol
! boudary type, .true. - inner boundary where vpar2=0, .false. - outer boundary where vpar2>0:
      logical :: inner_b,inner_e
! number of O- and X-points:
      integer :: n_o,n_x,n_sep
! cylindrical coordinates of the boundaries of the allowed region vpar2>0:
      double precision :: R_b,Z_b,R_e,Z_e
! values of normalized canonical momentum $\psi^\ast$ at the boundaries:
      double precision :: psiast_b,psiast_e
! cylindrical coordinates and normalized canonical momentum at O-points:
      double precision, dimension(:), allocatable :: R_o,Z_o,psiast_o
! cylindrical coordinates and normalized canonical momentum at X-points:
      double precision, dimension(:), allocatable :: R_x,Z_x,psiast_x
! Radii of separatrix crossings with Poincare cut (includes X-point radii):
      double precision, dimension(:), allocatable :: R_sep
! label of X-point among separatrix crossings:
      logical,          dimension(:), allocatable :: xpoint
    end type
!
! The region set is returned to the caller instead of held as module state, so
! concurrent energy slices do not share derived-type allocatable scratch.
    type region_set_t
      integer :: nregions = 0
      type(allowed_region), dimension(:,:), allocatable :: all_regions
    end type
  end module bounds_fixpoints_mod
!
  module form_classes_doublecount_mod
! number of classes for fixed H_0 and J_perp:
    integer :: nclasses
! type of class parameterization function:
    integer,          dimension(:), allocatable :: ifuntype
! class interval limits over cut parameter R_c and parallel velocity sign:
    double precision, dimension(:), allocatable :: R_class_beg,R_class_end,sigma_class
    !$omp threadprivate(nclasses, ifuntype, R_class_beg,R_class_end,sigma_class)
  end module form_classes_doublecount_mod
!
!------------------------------------------------------
!
  module orbit_dim_mod
    integer, parameter  :: neqm=5
    logical             :: write_orb=.false.
    logical             :: clip_resonance_classes=.true.
    integer             :: iunit1=100,next,numbasef
    double precision    :: Rorb_max
    logical             :: orbit_wall_loss=.false.
    integer             :: orbit_failure_stage=0
! Rorb_max is per-orbit scratch: find_bounce sets it to the orbit's maximum R as
! it integrates.  next is the extra-integral count: in the parallel resonance
! mode loop pertham writes it (next=0 then next=3 for its two find_bounce calls)
! and velo_res reads it for array dims, so each thread needs its own.  The
! find_bounce-heavy node-fill loops run in parallel, so both must be threadprivate.
! The grid-build regions never write next; they copyin the master value to keep
! the prior shared semantics.  numbasef is set before any parallel region and
! read only, so it stays shared.
    !$omp threadprivate(write_orb,Rorb_max,next,orbit_wall_loss,orbit_failure_stage)
  end module orbit_dim_mod
!
!------------------------------------------------------
!
  module get_matrix_mod
    double precision :: relerror=1.d-4, relmargin=1.d-4
! Class domain is trimmed to where |delphi_b| <= delphi_max: no resonance lives
! past 2*pi*m_max/n, and the divergent X-point endpoint beyond it defeats the
! adaptive sampler (sample_matrix ierr=2 -> class dropped).  delphi_max <= 0
! disables the trim (default, preserves non-torque paths).
    double precision :: delphi_max=0.d0
    integer          :: iclass
! Relerror and relmargin are configured by sample_class_doublecount for each
! energy slice; they are not immutable global input once the torque loop starts.
    !$omp threadprivate(relerror,relmargin,delphi_max,iclass)
  end module get_matrix_mod
!
!------------------------------------------------------
!
  module potato_boundary_scan_mod
! State used only by the bounded one-dimensional callbacks below.  A topology
! certificate is built independently in each energy-slice thread; these
! controls must therefore not be shared between concurrent scans.
    integer :: fixedpoint_scan_sigma=1,fixedpoint_scan_branch=1
    double precision :: fixedpoint_scan_left=0.d0,fixedpoint_scan_right=0.d0
    integer, parameter :: jperpmax_success=0
    integer, parameter :: jperpmax_unresolved=1
    integer, parameter :: jperpmax_invalid_domain=2
    integer :: jperpmax_status=jperpmax_unresolved
    logical :: jperpmax_certified=.false.
    double precision :: jperpmax_witness=0.d0,jperpmax_upper_bound=huge(1.d0)
    !$omp threadprivate(fixedpoint_scan_sigma,fixedpoint_scan_branch, &
    !$omp&                fixedpoint_scan_left,fixedpoint_scan_right)
    !$omp threadprivate(jperpmax_status,jperpmax_certified, &
    !$omp&                jperpmax_witness,jperpmax_upper_bound)
  end module potato_boundary_scan_mod
!
  module interp_cache_mod
! Exact memoization of interpolate_class_doublecount, keyed on the class
! coordinate x.  Within one integrate_class_resonances call the grid and iclass
! are fixed, but the modes scan the same x values (only the mode-dependent m,n
! combination changes), so the interpolation -- whose cost is dominated by
! plag_coeff -- is otherwise recomputed identically across modes.  A hit returns
! the stored result bit-for-bit, so resonance points are unchanged.  The cache is
! per-thread scratch: the parallel mode loop reads and fills it from every thread,
! so each thread keeps its own buffer and calls interp_cache_reset inside the
! parallel region to drop the previous class's grid before reusing it.
    integer :: ncache=0, n1cache=0
    double precision, dimension(:),   allocatable :: xcache
    double precision, dimension(:,:), allocatable :: veccache, dveccache
    !$omp threadprivate(ncache,n1cache,xcache,veccache,dveccache)
  contains
    subroutine interp_cache_reset
      ncache=0
    end subroutine interp_cache_reset
  end module interp_cache_mod
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Routines:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
  subroutine find_bounce(next,velo_ext,dtau_in,z_eqm,taub,delphi,extraset,ierr)
!
! Integrates the orbit over one bounce time (finds this time). If needed
! (write_orb=.true.) writes it to the file with unit number "iunit1".
! Besides orbit equations integrates along the orbit and extra set of
! functions of phase space coordinates.
! Agruments:
! velo_ext       - external routine to integrate (input)
! dtau_in        - maximum value of time step (input)
! z_eqm(5)       - phase space variables (input)
! taub           - bounce time (output)
! delphi         - toroidal shift per bounce time (output)
! extraset(next) - extra integrals along the orbit (inout)
! ierr           - error flag, 0 = success, 1 = orbit left domain (output)
!
  use orbit_dim_mod, only : neqm,write_orb,iunit1,Rorb_max,orbit_wall_loss, &
                            orbit_failure_stage
  use field_eq_mod, only : ierrfield
!
  implicit none
!
! nousecut - way to close the orbit during primary search:
! .true. - without using Poincare cut, can be used for general Phi distributions
! .false. - with using the cut, valid for Phi=Phi(psi) only
!  logical, parameter :: nousecut=.true.
  logical, parameter :: nousecut=.false.
!
! maximum number of Newton iterations for closing the orbit:
  integer, parameter :: niter=20
  integer, parameter :: max_primary_steps=100000
!
! relative error of orbit integrator:
  double precision, parameter :: relerr=1d-10 !8

  integer, intent(in) :: next
  double precision, intent(in) :: dtau_in
  double precision, dimension(neqm), intent(in) :: z_eqm
  double precision, dimension(next), intent(inout) :: extraset
  double precision, intent(out) :: taub, delphi
  integer, intent(out) :: ierr

  logical :: firstpass,primary_closed,newton_converged
  integer :: ndim, iter,primary_steps
  double precision :: dtau
  double precision :: dL2_pol,dL2_pol_start,dtau_newt,r_prev,z_prev
  double precision :: tau0,RNorm,ZNorm,vnorm,dnorm,vel_pol,dL2_pol_min
  double precision :: dZ_dR,sign_delZ,Z_tmp
  double precision, dimension(neqm+next) :: z,z_start,vz
!
  external velo_ext
!
  ndim = neqm+next
  ierr = 0
  primary_closed=.false.
  newton_converged=.false.
  primary_steps=0
  orbit_wall_loss=.false.
  orbit_failure_stage=0
  ierrfield=0
!
  z(1:neqm)=z_eqm
  if(next.gt.0) then
    z(neqm+1:ndim)=extraset
  endif
  Rorb_max=z(1)
!
!  dtau=dtau_in/max(z(4),1d-3)
  dtau=dtau_in
!
! Primary search:
!
  z_start=z
!
  orbit_failure_stage=1
  call velo_ext(dtau,z,vz)
  if(ierrfield.ne.0) then
    ierr = 1
    return
  endif
!
! unit 2D vector in the direction of the guiding center velocity in RZ-plane:
!
  vel_pol=sqrt(vz(1)**2+vz(3)**2)
  if(vel_pol.ne.vel_pol) then
    ierr=2
    return
  endif
  if(vel_pol.le.0.d0) then
    ierr=2
    return
  endif
  RNorm=vz(1)/vel_pol
  ZNorm=vz(3)/vel_pol
!
  tau0=0.d0
!
  if(write_orb) write (iunit1,*) z(1:neqm),vz(5)
!
! first step:
!
  orbit_failure_stage=2
  call odeint_allroutines(z,ndim,tau0,dtau,relerr,velo_ext)
  if(ierrfield.ne.0) then
    ierr = 1
    return
  endif
!
  taub=dtau
  if(nousecut) then
! initialize sqrt(2) of the poloidal length of the step
    dL2_pol=2.d0*(z(1)-z_start(1))**2+(z(3)-z_start(3))**2
  else
! initialize Poincare cut crossing check
!
    call get_poicut(z(1),Z_tmp,dZ_dR)
!
    sign_delZ=sign(1.d0,z(3)-Z_tmp)
    firstpass=.true.
    dL2_pol=2.d0*(z(1)-z_start(1))**2+(z(3)-z_start(3))**2
  endif
!
  if(write_orb) then
!
    call velo_ext(dtau,z,vz)
!
    write (iunit1,*) z(1:neqm),vz(5)
  endif
!
  do
    primary_steps=primary_steps+1
    if(primary_steps.gt.max_primary_steps) then
      ierr=2
      return
    endif
    r_prev=z(1)
    z_prev=z(3)
!
    orbit_failure_stage=3
    call odeint_allroutines(z,ndim,tau0,dtau,relerr,velo_ext)
    if(ierrfield.ne.0) then
      ierr = 1
      return
    endif
!
    taub=taub+dtau
    if(nousecut) then
! check if poloidal distance to the starting point is smaller than
! sqrt(2) of the poloidal length of the step:
      dL2_pol=2.d0*((z(1)-r_prev)**2+(z(3)-z_prev)**2)
      dL2_pol_start=(z(1)-z_start(1))**2+(z(3)-z_start(3))**2
      if(dL2_pol_start.lt.dL2_pol) then
        primary_closed=.true.
        exit
      endif
    else
! check if Poincare cut has been crossed
!
      call get_poicut(z(1),Z_tmp,dZ_dR)
!
      if(sign_delZ*(z(3)-Z_tmp).lt.0.d0) then
        if(firstpass) then
! first crossing (continue integration)
          firstpass=.false.
          sign_delZ=-sign_delZ
      else
! second crossing (stop integration)
          primary_closed=.true.
          exit
        endif
      endif
    endif
!
    if(write_orb) then
!
      call velo_ext(dtau,z,vz)
!
      write (iunit1,*) z(1:neqm),vz(5)
    endif
    Rorb_max=max(Rorb_max,z(1))
  enddo
!
! End primary search
!
  if(.not.primary_closed) then
    ierr=2
    return
  endif
! Newton adjustment
!
  dtau_newt=0.d0
  do iter=1,niter
!
    orbit_failure_stage=4
    call velo_ext(dtau,z,vz)
!
    vnorm=vz(1)*RNorm+vz(3)*ZNorm
    dnorm=(z_start(1)-z(1))*RNorm+(z_start(3)-z(3))*ZNorm
    if(dnorm**2.lt.dL2_pol*relerr) then
      newton_converged=.true.
      exit
    endif
    if(vnorm.ne.vnorm .or. abs(vnorm).le.64.d0*epsilon(1.d0)) then
      ierr=2
      return
    endif
    dtau_newt=dnorm/vnorm
    if(dtau_newt.ne.dtau_newt .or. abs(dtau_newt).gt.huge(dtau_newt)) then
      ierr=2
      return
    endif
!
    orbit_failure_stage=5
    call odeint_allroutines(z,ndim,tau0,dtau_newt,relerr,velo_ext)
    if(ierrfield.ne.0) then
      ierr = 1
      return
    endif
!
    taub=taub+dtau_newt
  enddo

  if(.not.newton_converged) then
    ierr=2
    return
  endif
  if(primary_steps.gt.10000 .or. abs(dtau_newt).gt.dtau .or. taub.gt.1.d6) then
    print *, 'find_bounce: large return correction, steps,dtau,dtau_newt,taub=', &
        primary_steps,dtau,dtau_newt,taub
  endif
!
  if(next.gt.0) then
    extraset=z(neqm+1:ndim)
  endif
!
  if(write_orb) then
!
    call velo_ext(dtau,z,vz)
!
    write (iunit1,*) z(1:neqm),vz(5)
    write (iunit1,*) 'NaN NaN NaN NaN NaN NaN'
  endif
!
! End Newton adjustment
!
  delphi=z(2)-z_start(2)
!
  end subroutine find_bounce
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine first_return_map(sigma_beg,R_beg,sigma_end,R_end, &
                              tau_fr,dphi_fr,ierr)
!
! Computes first return point to the Poincare cut.
! Input arguments:
! sigma_beg  - starting v_parallel sign
! R_beg      - starting R on the Poncare cut
! Output arguments:
! sigma_end  - v_parallel sign after first return
! R_end      - R on the Poncare cut after first return
!
  use orbit_dim_mod,     only : neqm
  use global_invariants, only : toten,perpinv,sigma,dtau
!
  implicit none
!
  integer, parameter :: niter=20
  double precision, parameter :: relerr=1d-10 !8
!
  integer :: iter, ierr
  double precision :: sigma_beg,R_beg,sigma_end,R_end,tau_fr,dphi_fr
  double precision :: tau0,dZ_dR,sign_delZ,Z_tmp
  double precision :: dtau_newt,r_prev,z_prev
  double precision, dimension(neqm) :: z,vz
!
  external velo
!
!
  z(1)=R_beg
  z(2)=0.d0
!
  call get_poicut(z(1),z(3),dZ_dR)
!
  call get_z45(toten,perpinv,sigma_beg,z,ierr)
!
  if(ierr.ne.0) return
!
! Primary search:
!
  tau0=0.d0
  tau_fr=0.d0
  Z_tmp=z(3)
!
  call odeint_allroutines(z,neqm,tau0,dtau,relerr,velo)
!
  tau_fr=tau_fr+dtau
!
  sign_delZ=sign(1.d0,z(3)-Z_tmp)
!
  do
!
    call odeint_allroutines(z,neqm,tau0,dtau,relerr,velo)
!
    tau_fr=tau_fr+dtau
!
    call get_poicut(z(1),Z_tmp,dZ_dR)
!
    if((z(3)-Z_tmp)*sign_delZ.lt.0.d0) exit
  enddo
!
! End primary search
!
! Newton adjustment
!
  do iter=1,niter
!
    call velo(dtau,z,vz)
!
    call get_poicut(z(1),Z_tmp,dZ_dR)
!
    dtau_newt=(Z_tmp-z(3))/vz(3)
!
    call odeint_allroutines(z,neqm,tau0,dtau_newt,relerr,velo)
!
    tau_fr=tau_fr+dtau_newt
!
    if(abs(Z_tmp-z(3)).lt.z(1)*relerr) exit
  enddo
!
  R_end=z(1)
  sigma_end=sign(1.d0,z(5))
  dphi_fr=z(2)
!
  end subroutine first_return_map
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_bmod_and_Phi(x,bmod,phi_elec)
!
! Wrapper for comuting module of B and electrostatic potential
!
  implicit none
!
  double precision :: bmod,sqrtg,phi_elec
  double precision, dimension(3)    :: x,bder,hcovar,hctrvr,hcurl,derphi
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  call elefie(x,phi_elec,derphi)
!
  end subroutine get_bmod_and_Phi
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_tormom(z,p_phi)
!
! Computes the normalized toroidal moment p_phi=psi^* for given phase space
! coordinates of alpha_lifetime
!
  use field_sub, only : psif,dpsidr,dpsidz
  use parmot_mod,   only : ro0
!
  implicit none
!
  integer, parameter :: neqm=5
  double precision :: p_phi
  double precision :: bmod,sqrtg,phi_elec
  double precision, dimension(3)    :: bder,hcovar,hctrvr,hcurl,derphi
  double precision, dimension(neqm) :: z
!
  call magfie(z(1:3),bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  p_phi=ro0*z(4)*z(5)*hcovar(2)+psif
!
  end subroutine get_tormom
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_z45(toten,perpinv,sigma,z,ierr)
!
! Computes momentum module z(4) and pitch parameter z(5)
! for given coordinates z(1:3), total energy "toten",
! perpendicular invariant "perpinv" and velocity sign "sigma".
! Error code "ierr": 0 - OK, 1 - negative kinetic energy,
! 2 - invalid perpendicular action or negative parallel kinetic energy
!
  use pitch_boundary_mod, only : resolve_pitch_squared
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
!
  implicit none
!
  integer, parameter :: neqm=5
  integer :: ierr
  double precision   :: toten,perpinv,sigma,bmod,phi_elec,p2,alam2
  double precision, dimension(neqm) :: z
!
  ierr=0
  z(4:5)=0.d0

! J_perp is a positive action.  The generated positive-part expression is
! reserved for the outer domain envelope; it must never turn a negative
! physical sample into a class/orbit by clamping it to zero.
  if(.not.ieee_is_finite(perpinv) .or. perpinv.lt.0.d0) then
    ierr=2
    return
  endif
!
  call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
  p2=toten-phi_elec
  if(p2.le.0.d0) then
    ierr=1
    return
  endif
  z(4)=sqrt(p2)
!
  call resolve_pitch_squared(perpinv*bmod/p2,alam2,ierr)
  if(ierr.ne.0) return
  z(5)=sign(sqrt(alam2),sigma)
!
  end subroutine get_z45
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine tormom_of_RZ(toten,perpinv,sigma,R_in,Z_in,p_phi,ierr)
!
! Computes the normalized toroidal moment p_phi=psi^* for given cylindrical
! coordinates (R_in,Z_in), invariants of motion (toten,perpinv) and parallel
! velocity sign sigma
!
  implicit none
!
  integer, parameter :: neqm=5
  integer :: ierr
  double precision :: toten,perpinv,sigma,R_in,Z_in,p_phi
  double precision, dimension(neqm) :: z
!
  ierr=0
  z(1)=R_in
  z(2)=0.d0
  z(3)=Z_in
!
  call get_z45(toten,perpinv,sigma,z,ierr)
!
  if(ierr.ne.0) then
    print *,'tormom_of_RZ ierr = ',ierr
    return
  endif
!
  call get_tormom(z,p_phi)
!
  end subroutine tormom_of_RZ
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine starter_doublecount(toten,perpinv,sigma,Rst,   &
                                 psiast,dpsiast_dRst,z,ierr)
!
! computes alpha_lifetime phase space coordinates z as functions of
! the position on the Poincare cut Rst, invariants of motion toten,perpinv
! and parallel velocity sign sigma. Computes also the normalized
! toroidal moment psiast=psi^* and its derivative along the cut with
! fixed invariants of motion dpsiast_dRst.
!
  use parmot_mod, only : gradpsiast,dpsiast_dR,dpsiast_dZ
  use poicut_mod, only : npc,rpc_arr
!
  implicit none
!
  integer, parameter :: neqm=5
  double precision, parameter :: dtau=0.d0
  integer :: ierr
  double precision :: toten,perpinv,sigma,Rst,psiast,dpsiast_dRst,dZ_dR
  double precision, dimension(neqm) :: z,vz
!
  if(Rst.lt.rpc_arr(0)) then
    ierr=11
    return
  elseif(Rst.gt.rpc_arr(npc)) then
    ierr=12
    return
  else
    ierr=0
  endif
!
  z(1)=Rst
  z(2)=0.d0
!
  call get_poicut(z(1),z(3),dZ_dR)
!
  call get_z45(toten,perpinv,sigma,z,ierr)
!
  if(ierr.ne.0) return
!
  call get_tormom(z,psiast)
!
  gradpsiast=.true.
!
  call velo(dtau,z,vz)
!
  gradpsiast=.false.
!
  dpsiast_dRst=dpsiast_dR+dpsiast_dZ*dZ_dR
!
  end subroutine starter_doublecount
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine velo_pphint(dtau_ext,z_ext,vz_ext)
!
! Computes the extended RHS for equations of motion and weighted
! integrals along the orbit of the powers of poloidal flux which are
! base functions phi_k for minimization of functional Eq.(29)
!
  use field_sub, only : psif
  use field_eq_mod, only : psi_sep
  use poicut_mod,    only : psimagaxis
  use orbit_dim_mod, only : neqm,next,numbasef
!
  implicit none
!
  integer :: i,k,kk
  double precision :: dtau_ext,psipow,psihat,weight
  double precision, dimension(neqm+next) :: z_ext,vz_ext
!
  call velo(dtau_ext,z_ext(1:neqm),vz_ext(1:neqm))
  call thetafun(weight)
!
  psihat=(psif-psimagaxis)/(psi_sep-psimagaxis)
  psipow=1.d0
  kk=neqm+numbasef
!
  do i=1,numbasef
    k=neqm+i
    vz_ext(k)=psipow*weight
    vz_ext(kk+i)=vz_ext(2)*vz_ext(k)
    psipow=psipow*psihat
  enddo
!
  end subroutine velo_pphint
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine thetafun(weight)
!
! generalized Theta-function Theta_w(x), Eq.(32)
! where x=1-(psi-psi_axis)(psi_b-psi_axis) and psi_b is poloidal flux
! at rho_pol boundary
!
  use field_sub, only : psif
  use poicut_mod,   only : psimagaxis,psi_bou
!
  implicit none
!
  integer, parameter :: kbou=3
  double precision, parameter :: boumarg=0.8d0
  double precision, parameter :: addbou=(1.-boumarg)**kbou
  double precision :: weight,dweight,ddweight,x1,x2,x,xx
!
  x=(psif-psimagaxis)/(psi_bou-psimagaxis)
!
  if(x.lt.1.d0) then
    xx=1.d0-x
    weight=xx**kbou/(xx**kbou+addbou)
  else
    weight=0.d0
  endif
!
  end subroutine thetafun
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine gpsigb_and_ders(R,Z,gpgb,dgpgb_dr,dgpgb_dz)
!
! Computes phi component of the vector product of gradient psi with gradient
! of B and the numerical derivatives of this product over R and Z.
! Used for finding the Poincare cut which is a line where this product is zero.
!
  implicit none
!
  double precision :: epsdif=1.d-6
  double precision :: R,Z,gpgb,dgpgb_dr,dgpgb_dz,RR,ZZ,hdif
!
  hdif=R*epsdif
!
  RR=R+hdif
  ZZ=Z
!
  call gradpsi_times_gradb
!
  dgpgb_dr=gpgb
  RR=R-hdif
!
  call gradpsi_times_gradb
!
  dgpgb_dr=(dgpgb_dr-gpgb)/(2.d0*hdif)
  RR=R
  ZZ=Z+hdif
!
  call gradpsi_times_gradb
!
  dgpgb_dz=gpgb
  ZZ=Z-hdif
!
  call gradpsi_times_gradb
!
  dgpgb_dz=(dgpgb_dz-gpgb)/(2.d0*hdif)
  ZZ=Z
!
  call gradpsi_times_gradb
!
!------------
  contains
!------------
!
  subroutine gradpsi_times_gradb
!
  use field_sub, only : psif,dpsidr,dpsidz
!
  implicit none
!
  double precision :: bmod,sqrtg
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
!
  x(1)=RR
  x(2)=0.d0
  x(3)=ZZ
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  gpgb=dpsidr*bder(3)-dpsidz*bder(1)
!
  end subroutine gradpsi_times_gradb
!-------------
  end subroutine gpsigb_and_ders
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_poicut(rho_pol,npline)
!
! Finds a line where $\nabla \psi \times \nabla B = 0$. Line starts and ends
! at the separatrix. In up-down symmeric configuration it is reduced to
! a midplane poloidal cross-section $Z=0$. This line contains all fixed
! (O and X) points of the flow as well as the magnetic axis and can be
! used as a Poincare cut.
!
  use field_sub, only : psif,dpsidr,dpsidz
  use field_eq_mod, only : psi_axis,psi_sep,rtf
  use poicut_mod, only   : npc,rpc_beg,h_rpc,rpc_arr,zpc_arr
!
  implicit none
!
  integer, parameter :: nsplit=100,niter=100
  double precision, parameter :: relerr=1d-12
  integer :: i,j,npline
  double precision :: rho_pol,psils,sigpsi
  double precision :: h_R,det,delR,delZ,stepfac,err_dist,stepmod
  double precision :: R,Z,gpgb,dgpgb_dr,dgpgb_dz
  double precision :: Rb,Zb,Re,Ze
!
  npc=npline
  allocate(rpc_arr(0:npc),zpc_arr(0:npc))
!
!
  R=1.d0
  Z=0.d0
!
  call gpsigb_and_ders(R,Z,gpgb,dgpgb_dr,dgpgb_dz)
!
  err_dist=rtf*relerr
  h_R=rtf/dble(nsplit)
  sigpsi=sign(1.d0,psi_sep-psi_axis)
  psils=psi_axis+rho_pol**2*(psi_sep-psi_axis)
  R=rtf
!
  do i=2,nsplit
    R=R-h_R
!
    call gpsigb_and_ders(R,Z,gpgb,dgpgb_dr,dgpgb_dz)
!
    if((psif-psils)*(psi_sep-psi_axis).gt.0.d0) exit
  enddo
!
  do i=1,niter
!
    call gpsigb_and_ders(R,Z,gpgb,dgpgb_dr,dgpgb_dz)
!
    det=dpsidr*dgpgb_dz-dpsidz*dgpgb_dr
    delR=((psif-psils)*dgpgb_dz-dpsidz*gpgb)/det
    delZ=(dpsidr*gpgb-(psif-psils)*dgpgb_dr)/det
    stepmod=sqrt(delR**2+delZ**2)
    stepfac=h_R/max(h_R,stepmod)
    R=R-delR*stepfac
    Z=Z-delZ*stepfac
    if(stepmod.lt.err_dist) exit
  enddo
!
  Rb=R
  Zb=Z
!
  R=rtf
!
  do i=2,nsplit
    R=R+h_R
!
    call gpsigb_and_ders(R,Z,gpgb,dgpgb_dr,dgpgb_dz)
!
    if((psif-psils)*(psi_sep-psi_axis).gt.0.d0) exit
  enddo
!
  do i=1,niter
!
    call gpsigb_and_ders(R,Z,gpgb,dgpgb_dr,dgpgb_dz)
!
    det=dpsidr*dgpgb_dz-dpsidz*dgpgb_dr
    delR=((psif-psils)*dgpgb_dz-dpsidz*gpgb)/det
    delZ=(dpsidr*gpgb-(psif-psils)*dgpgb_dr)/det
    stepmod=sqrt(delR**2+delZ**2)
    stepfac=h_R/max(h_R,stepmod)
    R=R-delR*stepfac
    Z=Z-delZ*stepfac
    if(stepmod.lt.err_dist) exit
  enddo
!
  Re=R
  Ze=Z
!
  rpc_beg=Rb
  h_rpc=(Re-Rb)/dble(npc)
  rpc_arr(0)=Rb
  zpc_arr(0)=Zb
  rpc_arr(npc)=Re
  zpc_arr(npc)=Ze
!
  R=Rb
  Z=Zb
!
  do j=1,npline-1
    R=R+h_rpc
!
    do i=1,niter
!
      call gpsigb_and_ders(R,Z,gpgb,dgpgb_dr,dgpgb_dz)
!
      delZ=gpgb/dgpgb_dz
      Z=Z-delZ
      if(abs(delZ).lt.err_dist) exit
    enddo
!
    rpc_arr(j)=R
    zpc_arr(j)=Z
  enddo
!
  call find_magaxis
!
  if(.true.) then
    print *,'cut range:',rpc_arr(0),rpc_arr(npc)
  endif
!
  end subroutine find_poicut
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_poicut(R,Z,dZ_dR)
!
! Interpolates a pre-calculated line where $\nabla \psi \times \nabla B = 0$,
! i.e. computes function $Z(R)$ and its derivative dZ_dR.
!
  use poicut_mod, only : npc,rpc_beg,h_rpc,rpc_arr,zpc_arr

  implicit none
!
  integer, parameter :: nplag=4, nder=1
  integer :: ibeg,iend
  double precision :: s,Z,R,dZ_dR
  double precision, dimension(0:nder,nplag) :: coef
!
  ibeg=int((R-rpc_beg)/h_rpc)
  ibeg=max(0,ibeg-nplag/2)
  iend=ibeg+nplag-1
  if(iend.gt.npc) then
    iend=npc
    ibeg=iend+1-nplag
  endif
!
  call plag_coeff(nplag,nder,R,rpc_arr(ibeg:iend),coef)
!
  Z=sum(coef(0,:)*zpc_arr(ibeg:iend))
  dZ_dR=sum(coef(1,:)*zpc_arr(ibeg:iend))
!
  end subroutine get_poicut
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_magaxis
!
! Find cylindrical coordinates of the magnetic axis
!
  use field_sub, only : psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  use poicut_mod, only   : npc,rpc_arr,rmagaxis,zmagaxis,psimagaxis
!
  implicit none
!
  integer, parameter :: nsplit=100,niter=100
  double precision, parameter :: relerr=1d-12
  integer :: i,iter
  double precision :: R,Z,dZ_dR,bmod,phi_elec
  double precision :: errdist2,h_R,dpsi_dl,dpsi_dl_prev,del_R,del_Z,det
  double precision,dimension(3) :: x
!
  errdist2=(relerr*rpc_arr(npc))**2
  h_R=(rpc_arr(npc)-rpc_arr(0))/dble(nsplit)
  R=rpc_arr(0)
!
  call get_poicut(R,Z,dZ_dR)
!
  x(1)=R
  x(2)=0.d0
  x(3)=Z
!
  call get_bmod_and_Phi(x,bmod,phi_elec)
!
  dpsi_dl=dpsidr+dpsidz*dZ_dR
!
  do i=1,nsplit
    dpsi_dl_prev=dpsi_dl
    R=rpc_arr(0)+h_R*dble(i)
!
    call get_poicut(R,Z,dZ_dR)
!
    x(1)=R
    x(3)=Z
!
    call get_bmod_and_Phi(x,bmod,phi_elec)
!
    dpsi_dl=dpsidr+dpsidz*dZ_dR
    if(dpsi_dl_prev*dpsi_dl.lt.0.d0) exit
  enddo
!
  do iter=1,niter
    x(1)=R
    x(3)=Z
!
    call get_bmod_and_Phi(x,bmod,phi_elec)
!
    det=d2psidr2*d2psidz2-d2psidrdz**2
    del_R=(d2psidrdz*dpsidz-d2psidz2*dpsidr)/det
    del_Z=(d2psidrdz*dpsidr-d2psidr2*dpsidz)/det
    R=R+del_R
    Z=Z+del_Z
    if(del_R**2+del_Z**2.lt.errdist2) exit
  enddo
!
  rmagaxis=R
  zmagaxis=Z
  psimagaxis=psif
!
  if(.true.) then
    print *,'magnetic axis:'
    print *,rmagaxis,zmagaxis
  endif
!
  end subroutine find_magaxis
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine rhopol_boundary(rho_pol)
!
! Finds poloidal flux psi_bou for a given poloidal radius rho_pol and cylindrical
! coordinates of LFS and HFS intersections of the Poincare cut with the flux surface
! with given rho_pol. Results are stored in the module poicut_mod
!
  use field_sub, only : psif,dpsidr,dpsidz
  use field_eq_mod, only : psi_axis,psi_sep
  use poicut_mod, only   : npc,rpc_arr,rmagaxis,zmagaxis,psimagaxis, &
                           rhopol_bou,psi_bou,Rbou_lfs,Zbou_lfs,Rbou_hfs,Zbou_hfs
!
  implicit none
!
  integer, parameter :: niter=100
  double precision, parameter :: relerr=1d-12
  integer :: i,iter
  double precision :: rho_pol,psipol,R_rhopol,Z_rhopol
  double precision :: R,Z,dZ_dR,bmod,phi_elec
  double precision :: errdist,dpsi_dl,del_R
  double precision,dimension(3) :: x
!
  rhopol_bou=rho_pol
  psi_bou=psimagaxis+(psi_sep-psimagaxis)*rho_pol**2
!
  errdist=relerr*rpc_arr(npc)
!
  Rbou_lfs=rmagaxis+(rpc_arr(npc)-rmagaxis)*rho_pol
  x(2)=0.d0
!
  do iter=1,niter
!
    call get_poicut(Rbou_lfs,Zbou_lfs,dZ_dR)
!
    x(1)=Rbou_lfs
    x(3)=Zbou_lfs
!
    call get_bmod_and_Phi(x,bmod,phi_elec)
!
    dpsi_dl=dpsidr+dpsidz*dZ_dR
    del_R=(psi_bou-psif)/dpsi_dl
    Rbou_lfs=Rbou_lfs+del_R
    if(abs(del_R).lt.errdist) exit
  enddo
!
  call get_poicut(Rbou_lfs,Zbou_lfs,dZ_dR)
!
  Rbou_hfs=rmagaxis+(rpc_arr(0)-rmagaxis)*rho_pol
!
  do iter=1,niter
!
    call get_poicut(Rbou_hfs,Zbou_hfs,dZ_dR)
!
    x(1)=Rbou_hfs
    x(3)=Zbou_hfs
!
    call get_bmod_and_Phi(x,bmod,phi_elec)
!
    dpsi_dl=dpsidr+dpsidz*dZ_dR
    del_R=(psi_bou-psif)/dpsi_dl
    Rbou_hfs=Rbou_hfs+del_R
    if(abs(del_R).lt.errdist) exit
  enddo
!
  call get_poicut(Rbou_hfs,Zbou_hfs,dZ_dR)
!
  if(.true.) then
    print *,'outer boundary LFS:'
    print *,Rbou_lfs,Zbou_lfs
    print *,'outer boundary HFS:'
    print *,Rbou_hfs,Zbou_hfs
  endif
!
  end subroutine rhopol_boundary
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine poltor_field_dir(ifdir_type)
!
! Determines mutual direction of poloidal and toroidal fields -
! If $B^\theta*B^\varphi > 0$ (typical AUG case) idir_type=1, otherwise idir_type=2
! In case 1 the non-vanishing axis (located at the LFS from the magnetic axis)
! corresponds to the counter-passing orbits (vpar<0), in case 2 - for the co-passing (vpar>0)
!
  use poicut_mod, only : npc,rpc_arr,zpc_arr
!
  implicit none
!
  integer :: ifdir_type
  double precision :: bmod,sqrtg
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
!
  x(1)=rpc_arr(npc)
  x(2)=0.d0
  x(3)=zpc_arr(npc)
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  if(hctrvr(2)*hctrvr(3).gt.0.d0) then
    ifdir_type=1
  else
    ifdir_type=2
  endif
!
  end subroutine poltor_field_dir
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_bounds_fixpoints(regions,ierr)
!
  use global_invariants,    only : toten,perpinv,sigma
  use find_all_roots_mod,   only : nroots,roots,relerr_allroots,          &
                                   root_eval_valid,root_eval_error,        &
                                   root_invalid_domain,root_no_intersection
  use poicut_mod,           only : npc,rpc_arr,zpc_arr,rmagaxis, &
                                   Rbou_lfs,Zbou_lfs,Rbou_hfs,Zbou_hfs
  use field_sub, only : psif
  use bounds_fixpoints_mod, only : allowed_region,region_set_t
  use cc_mod, only : wrbounds
  use orbit_dim_mod, only : clip_resonance_classes
!
  implicit none
!
  type(region_set_t), intent(inout) :: regions
  integer :: nbounds,nregions
  double precision,     dimension(:),   allocatable :: R_bo,Z_bo,psiast_bo
  type(allowed_region), dimension(:,:), allocatable :: all_regions
  logical, parameter :: doublecount=.true.
  integer, parameter :: niter=100
  double precision, parameter :: relerr=1d-12
  logical :: separin_b,separin_e,within_rhopol,start_minmax
  integer :: ierr,ireg0
  integer :: i,k,iter,isig,ireg,n_o,n_x,nxp_tot,ixp_tot,nsc,nsc_raw
  double precision :: R,Z,vpar2,dvpar2_dR,dZ_dR,errdist2,R_b,R_e,dummy
  double precision :: R_b_in,R_e_in,p_phi,sigma_new,R_new,R_old,tau_fr,dphi_fr
  double precision :: gpgb,dgpgb_dr,dgpgb_dz,det,del_R,del_Z,delR_marg
  double precision :: pphi_b,pphi_e,pphi_min_reg,pphi_max_reg,pphi_min,pphi_max
  double precision :: separatrix_root_tolerance
  double precision :: pphi_minmax
  double precision, dimension(2) :: gradvpar2,Rblfs,Rbhfs,Zblfs,Zbhfs
  logical,          dimension(:),   allocatable :: opoint_extr,logdummy
  integer,          dimension(:,:), allocatable :: ipoi_x_tot
  integer,          dimension(:),   allocatable :: ipoi_tmp
  double precision, dimension(:),   allocatable :: R_extr,Z_extr,psiast_extr, &
                                                   psiast_x_tot,R_x_tot,rsc_tmp
!------------
!
  ierr=0
  regions%nregions=0
  if(allocated(regions%all_regions)) deallocate(regions%all_regions)

  if(wrbounds) then
! write vpar2 on the cut as function of R
    open(2001,file='vpar2_vs_R.dat')
    do i=1,10000
      R=rpc_arr(0)+(rpc_arr(npc)-rpc_arr(0))*dble(i)/10001.d0
!
      call vparzero1D(R,vpar2,Z)
!
      write(2001,*) R,vpar2,Z
    enddo
    close(2001)
  endif
!
! Find boundaries of allowed regions
!
! check if allowed region extends to the separatrix at the HFS:
!
  call vparzero1D(rpc_arr(0),vpar2,dummy)
!
  if(vpar2.gt.0.d0) then
    separin_b=.true.
  else
    separin_b=.false.
  endif
!
!
! check if allowed region extends to the separatrix at the LFS:
!
  call vparzero1D(rpc_arr(npc),vpar2,dummy)
!
  if(vpar2.gt.0.d0) then
    separin_e=.true.
  else
    separin_e=.false.
  endif
!
! find inner boundaries of allowed regions:
  relerr_allroots=1.d-12
!
  call find_all_roots(vparzero1D,rpc_arr(0),rpc_arr(npc),ierr)
!
  if(ierr.ne.0) then
    print *,'find_bounds_fixpoints: error in find_all_roots, inner boundaries'
    return
  endif
!
! Keep the exact turning root.  p_parallel=0 is an allowed phase-space
! boundary; moving it by an unconstrained relative offset can leave the
! Poincare interval and fabricate an evaluation in forbidden space.
  if(nroots.gt.0) then
    if(allocated(R_bo)) deallocate(R_bo,Z_bo,psiast_bo)
    allocate(R_bo(nroots),Z_bo(nroots),psiast_bo(nroots))
    R_bo=roots
  endif
!
! store R, Z and psi of inner boundaries:
!
  do i=1,nroots
    R=R_bo(i)
!
    call get_poicut(R,Z,dZ_dR)
!
! compute poloidal flux psif at the boundary point
!
    call vparzero_vec(R,Z,vpar2,gradvpar2)
!
    Z_bo(i)=Z
    psiast_bo(i)=psif
  enddo
!
  nbounds=nroots
  if(wrbounds) then
!write the boundary data:
    open(3001,file='vpar2_zeropoints.dat')
    do i=1,nroots
      write(3001,*) R_bo(i),0.d0,Z_bo(i),psiast_bo(i)
    enddo
    close(3001)
  endif
!
! End find boundaries of allowed regions
!
! Determine the number of allowed (vpar2>0) regions:
!
! general case, vpar2 may be negative at the LFS separatrix crossing,
! regions are numbered from 1
  if(nbounds.eq.0) then
!no inner boundaries
    if(separin_e) then
! vpar2>0 in the whole domain, one region over the whole domain
      nregions=1
    else
! vpar2<0 in the whole domain, no regions in the domain
      nregions=0
      regions%nregions=0
      ! No allowed orbit at this J_perp is an empty physical class set, not
      ! a topology failure.  The matrix callback records zero contribution.
      ierr=0
      return
    endif
  elseif(modulo(nbounds,2).eq.0) then
! even number of inner boundaries
    if(separin_e) then
! vpar2>0 at domain boundaries
      nregions=nbounds/2+1
    else
! vpar2<0 at domain boundaries
      nregions=nbounds/2
    endif
  else
! odd number of inner boundaries
    nregions=(nbounds+1)/2
  endif
!
! End determine the number of vpar2>0 regions
!
!...........................
!
! Form allowed regions for the whole Poincare cut:
! general case, vpar2 may be negative at the LFS separatrix crossing,
! regions are numbered from 1
!
  if(allocated(all_regions)) deallocate(all_regions)
  allocate(all_regions(2,nregions))
!
! start point of the first region:
  if(separin_b) then
    all_regions(:,1)%R_b=rpc_arr(0)
    all_regions(:,1)%Z_b=zpc_arr(0)
    R=rpc_arr(0)
    Z=zpc_arr(0)
!
    do isig=1,2
      if(isig.eq.1) then
        sigma=1.d0
      else
        sigma=-1.d0
      endif
!
      call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
!
      all_regions(isig,1)%psiast_b=p_phi
      all_regions(isig,1)%inner_b=.false.
    enddo
!
  else
    all_regions(:,1)%R_b=R_bo(1)
    all_regions(:,1)%Z_b=Z_bo(1)
    all_regions(:,1)%psiast_b=psiast_bo(1)
    all_regions(:,1)%inner_b=.true.
  endif
!
! end point of the last region:
  if(separin_e) then
    all_regions(:,nregions)%R_e=rpc_arr(npc)
    all_regions(:,nregions)%Z_e=zpc_arr(npc)
    R=rpc_arr(npc)
    Z=zpc_arr(npc)
!
    do isig=1,2
      if(isig.eq.1) then
        sigma=1.d0
      else
        sigma=-1.d0
      endif
!
      call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
!
      all_regions(isig,nregions)%psiast_e=p_phi
      all_regions(isig,nregions)%inner_e=.false.
    enddo
!
  else
    all_regions(:,nregions)%R_e=R_bo(nbounds)
    all_regions(:,nregions)%Z_e=Z_bo(nbounds)
    all_regions(:,nregions)%psiast_e=psiast_bo(nbounds)
    all_regions(:,nregions)%inner_e=.true.
  endif
!
! store all data for inner boundaries:
  do ireg=1,nregions-1
    if(separin_b) then
      k=ireg*2-1
    else
      k=ireg*2
    endif
    all_regions(:,ireg)%R_e=R_bo(k)
    all_regions(:,ireg)%Z_e=Z_bo(k)
    all_regions(:,ireg)%psiast_e=psiast_bo(k)
    all_regions(:,ireg)%inner_e=.true.
    k=k+1
    all_regions(:,ireg+1)%R_b=R_bo(k)
    all_regions(:,ireg+1)%Z_b=Z_bo(k)
    all_regions(:,ireg+1)%psiast_b=psiast_bo(k)
    all_regions(:,ireg+1)%inner_b=.true.
  enddo
!
! End form allowed regions for the whole Poincare cut
!
!...........................
!
! Find minimum and maximum values of p_phi in the domain limited
! by the boundary flux surface with given rho_pol:
!
  if(clip_resonance_classes) then
  start_minmax=.true.
!
  do isig=1,2
!
    if(isig.eq.1) then
      sigma=1.d0
    else
      sigma=-1.d0
    endif
!
    do ireg=1,nregions
      within_rhopol=.false.
      R_b_in=all_regions(isig,ireg)%R_b
      R_e_in=all_regions(isig,ireg)%R_e
!
      if(R_b_in.gt.Rbou_hfs.and.R_e_in.lt.Rbou_lfs) then
! allowed region is fully within rho_pol domain
        within_rhopol=.true.
        pphi_b=all_regions(isig,ireg)%psiast_b
        pphi_e=all_regions(isig,ireg)%psiast_e
      elseif(R_b_in.lt.Rbou_hfs.and.R_e_in.gt.Rbou_lfs) then
! rho_pol domain is fully within allowed region
        within_rhopol=.true.
        R_b_in=Rbou_hfs
!
        call tormom_of_RZ(toten,perpinv,sigma,Rbou_hfs,Zbou_hfs,pphi_b,ierr)
!
        R_e_in=Rbou_lfs
!
        call tormom_of_RZ(toten,perpinv,sigma,Rbou_lfs,Zbou_lfs,pphi_e,ierr)
!
      elseif(R_b_in.lt.Rbou_hfs.and.R_e_in.gt.Rbou_hfs) then
! HFS boundary of the rho_pol domain is within allowed region
        within_rhopol=.true.
        R_b_in=Rbou_hfs
!
        call tormom_of_RZ(toten,perpinv,sigma,Rbou_hfs,Zbou_hfs,pphi_b,ierr)
!
        pphi_e=all_regions(isig,ireg)%psiast_e
      elseif(R_b_in.lt.Rbou_lfs.and.R_e_in.gt.Rbou_lfs) then
! LFS boundary of the rho_pol domain is within allowed region
        within_rhopol=.true.
        pphi_b=all_regions(isig,ireg)%psiast_b
        R_e_in=Rbou_lfs
!
        call tormom_of_RZ(toten,perpinv,sigma,Rbou_lfs,Zbou_lfs,pphi_e,ierr)
!
      endif
!
      if(within_rhopol) then
!
        call classify_extrema(R_b_in,R_e_in,ierr)
        if(ierr.ne.0) then
          print *,'find_bounds_fixpoints: unresolved extrema for range, ierr = ',ierr
          return
        endif
!
        if(nroots.gt.0) then
          pphi_min_reg=min(pphi_b,pphi_e,minval(psiast_extr))
          pphi_max_reg=max(pphi_b,pphi_e,maxval(psiast_extr))
        else
          pphi_min_reg=min(pphi_b,pphi_e)
          pphi_max_reg=max(pphi_b,pphi_e)
        endif
!
        if(start_minmax) then
          start_minmax=.false.
          pphi_min=pphi_min_reg
          pphi_max=pphi_max_reg
        else
          pphi_min=min(pphi_min,pphi_min_reg)
          pphi_max=max(pphi_max,pphi_max_reg)
        endif
!
      endif
!
    enddo
  enddo
!
  if(start_minmax) then
    print *,'minimum and maximum values of p_phi are not determined H,J = ', &
        toten,perpinv
    ! Every allowed region lies outside the clipped rho_pol domain.  There is
    ! therefore no resonance class to construct at this J_perp; return the
    ! certified empty class set instead of continuing with uninitialized
    ! p_phi limits.
    regions%nregions=0
    if(allocated(all_regions)) deallocate(all_regions)
    ierr=0
    return
  endif

!
! End find minimum and maximum values of p_phi in the domain limited
! by the boundary flux surface with given rho_pol
!
!...........................
!
! Cut out from the allowed regions the segments occupied by the orbits
! never visiting rho_pol domain. Segments are cut only at the edges:
!
  do isig=1,2
!
    if(isig.eq.1) then
      sigma=1.d0
    else
      sigma=-1.d0
    endif
!
    do ireg=1,nregions
      R_b_in=all_regions(isig,ireg)%R_b
      R_e_in=all_regions(isig,ireg)%R_e

!
! Left boundary:
!
      if(all_regions(isig,ireg)%psiast_b.lt.pphi_min) then
        pphi_minmax=pphi_min
        relerr_allroots=1.d-11
!
        call find_all_roots(boucross_with_endpoint_limit,R_b_in,R_e_in,ierr)
!
          if(ierr.ne.0) then
            print *,'find_bounds_fixpoints: error in find_all_roots, cut left boundary 1'
            return
          endif
!
        if(nroots.gt.0) then
          within_rhopol=.true.
          R=roots(1)
!
          call get_poicut(R,Z,dZ_dR)
          call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
          if(ierr.ne.0) return
!
          all_regions(isig,ireg)%inner_b=.false.
          all_regions(isig,ireg)%R_b=R
          all_regions(isig,ireg)%Z_b=Z
          all_regions(isig,ireg)%psiast_b=p_phi
        else
! the whole region should be cut out
          within_rhopol=.false.
        endif
      elseif(all_regions(isig,ireg)%psiast_b.gt.pphi_max) then
        pphi_minmax=pphi_max
        relerr_allroots=1.d-11
!
        call find_all_roots(boucross_with_endpoint_limit,R_b_in,R_e_in,ierr)
!
          if(ierr.ne.0) then
            print *,'find_bounds_fixpoints: error in find_all_roots, cut left boundary 2'
            return
          endif
!
        if(nroots.gt.0) then
          within_rhopol=.true.
          R=roots(1)
!
          call get_poicut(R,Z,dZ_dR)
          call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
          if(ierr.ne.0) return
!
          all_regions(isig,ireg)%inner_b=.false.
          all_regions(isig,ireg)%R_b=R
          all_regions(isig,ireg)%Z_b=Z
          all_regions(isig,ireg)%psiast_b=p_phi
        else
! the whole region should be cut out
          within_rhopol=.false.
        endif
      else
        within_rhopol=.true.
      endif
!
! Right boundary:
!
      if(within_rhopol) then
        if(all_regions(isig,ireg)%psiast_e.lt.pphi_min) then
          pphi_minmax=pphi_min
          relerr_allroots=1.d-11
!
          call find_all_roots(boucross_with_endpoint_limit,R_b_in,R_e_in,ierr)
!
          if(ierr.ne.0) then
            print *,'find_bounds_fixpoints: error in find_all_roots, cut right boundary 1'
            return
          endif
!
          if(nroots.gt.0) then
            R=roots(nroots)
!
            call get_poicut(R,Z,dZ_dR)
            call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
            if(ierr.ne.0) return
!
            all_regions(isig,ireg)%inner_e=.false.
            all_regions(isig,ireg)%R_e=R
            all_regions(isig,ireg)%Z_e=Z
            all_regions(isig,ireg)%psiast_e=p_phi
          else
! A right endpoint outside the global range must have a certified crossing.
! A successful scan with no crossing is an unresolved topology condition,
! not permission to omit torque from this region.
            print *,'cutter: no certified right boundary'
            ierr=root_no_intersection
            return
          endif
        elseif(all_regions(isig,ireg)%psiast_e.gt.pphi_max) then
          pphi_minmax=pphi_max
          relerr_allroots=1.d-11
!
          call find_all_roots(boucross_with_endpoint_limit,R_b_in,R_e_in,ierr)
!
          if(ierr.ne.0) then
            print *,'find_bounds_fixpoints: error in find_all_roots, cut right boundary 2'
            return
          endif
!
          if(nroots.gt.0) then
            R=roots(nroots)
!
            call get_poicut(R,Z,dZ_dR)
            call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
            if(ierr.ne.0) return
!
            all_regions(isig,ireg)%inner_e=.false.
            all_regions(isig,ireg)%R_e=R
            all_regions(isig,ireg)%Z_e=Z
            all_regions(isig,ireg)%psiast_e=p_phi
          else
! A right endpoint outside the global range must have a certified crossing.
! A successful scan with no crossing is an unresolved topology condition.
            print *,'cutter: no certified right boundary'
            ierr=root_no_intersection
            return
          endif
        endif
      endif
!
      all_regions(isig,ireg)%within_rhopol=within_rhopol
    enddo
  enddo

!
! End cut out from the allowed regions the segments occupied by the orbits
! never visiting rho_pol domain.
!
  else
    do isig=1,2
      do ireg=1,nregions
        all_regions(isig,ireg)%within_rhopol=.true.
      enddo
    enddo
  endif
!
!...........................
!
! Find and classify extremum points of psi^* within each region:
!
  if(wrbounds) then
! write psi^* dependencies along the cut (vs R):
    open(5001,file='psiast_vsR_p.dat')
    open(5002,file='psiast_vsR_m.dat')
  endif
!
  do ireg=1,nregions
    do isig=1,2
      if(isig.eq.1) then
        sigma=1.d0
      else
        sigma=-1.d0
      endif
      R_b_in=all_regions(isig,ireg)%R_b
      R_e_in=all_regions(isig,ireg)%R_e
!
      call classify_extrema(R_b_in,R_e_in,ierr)
      if(ierr.ne.0) then
        print *,'find_bounds_fixpoints: unresolved extremum topology, ierr = ',ierr
        return
      endif
!
      n_o=0
      n_x=0
!
      do i=1,nroots
        if(opoint_extr(i)) then
          n_o=n_o+1
        else
          n_x=n_x+1
        endif
      enddo
!
      all_regions(isig,ireg)%n_o=n_o
      all_regions(isig,ireg)%n_x=n_x
!
      if(n_o.gt.0) then
        allocate(all_regions(isig,ireg)%R_o(n_o))
        allocate(all_regions(isig,ireg)%Z_o(n_o))
        allocate(all_regions(isig,ireg)%psiast_o(n_o))
      endif
!
      if(n_x.gt.0) then
        allocate(all_regions(isig,ireg)%R_x(n_x))
        allocate(all_regions(isig,ireg)%Z_x(n_x))
        allocate(all_regions(isig,ireg)%psiast_x(n_x))
      endif
!
      n_o=0
      n_x=0
!
      do i=1,nroots
        if(opoint_extr(i)) then
          n_o=n_o+1
          all_regions(isig,ireg)%R_o(n_o)=R_extr(i)
          all_regions(isig,ireg)%Z_o(n_o)=Z_extr(i)
          all_regions(isig,ireg)%psiast_o(n_o)=psiast_extr(i)
        else
          n_x=n_x+1
          all_regions(isig,ireg)%R_x(n_x)=R_extr(i)
          all_regions(isig,ireg)%Z_x(n_x)=Z_extr(i)
          all_regions(isig,ireg)%psiast_x(n_x)=psiast_extr(i)
        endif
      enddo
    enddo
  enddo
!
  if(wrbounds) then
    close(5001)
    close(5002)
  endif
!
! End find and classify extremum points of psi^* within each region
!
!
!...........................
!
! Find separatrix crossings of the Poincare cut:
!
  nxp_tot=0
!
  do isig=1,2
    do ireg=1,nregions
      if(all_regions(isig,ireg)%within_rhopol) then
        nxp_tot=nxp_tot+all_regions(isig,ireg)%n_x
      endif
    enddo
  enddo
!
  allocate(psiast_x_tot(nxp_tot),R_x_tot(nxp_tot))
  allocate(ipoi_x_tot(2,nxp_tot)) !isig and ireg for X-point in the total list
!
  nxp_tot=0
!
  do isig=1,2
    do ireg=1,nregions
      if(all_regions(isig,ireg)%within_rhopol) then
        do i=1,all_regions(isig,ireg)%n_x
          nxp_tot=nxp_tot+1
          ipoi_x_tot(1,nxp_tot)=isig
          ipoi_x_tot(2,nxp_tot)=ireg
          psiast_x_tot(nxp_tot)=all_regions(isig,ireg)%psiast_x(i)
          R_x_tot(nxp_tot)=all_regions(isig,ireg)%R_x(i)
        enddo
      endif
    enddo
  enddo
!
  do isig=1,2
!
    if(isig.eq.1) then
      sigma=1.d0
    else
      sigma=-1.d0
    endif
!
! protective margin around separatrix crossings (set by hands after experimenting)
! needed to exclude the X-point (singular separarix crossing) from the numerical
! search of regular separatrix crossings
    delR_marg=1.d-6*rpc_arr(npc)
!
    do ireg=1,nregions
      if(all_regions(isig,ireg)%within_rhopol) then
        nsc=0
!
        do ixp_tot=1,nxp_tot
          R_b=all_regions(isig,ireg)%R_b
          R_e=all_regions(isig,ireg)%R_e
!
          if(ipoi_x_tot(1,ixp_tot).eq.isig .and. &
             ipoi_x_tot(2,ixp_tot).eq.ireg) then
! X-point belongs to this region, exclude it from the search:
            if(R_b.lt.R_x_tot(ixp_tot)-delR_marg) then
              relerr_allroots=1.d-12
!
              call find_all_roots(sepcross,R_b,R_x_tot(ixp_tot)-delR_marg,ierr)
!
              if(ierr.ne.0) then
                print *,'find_bounds_fixpoints: error in find_all_roots, sepcross 1, H,J,sigma,ireg,Rlo,Rhi,ierr = ', &
                    toten,perpinv,sigma,ireg,R_b_in,R_x_tot(ixp_tot)-delR_marg,ierr
                return
              endif
!
              nsc=nsc+nroots
            endif
            if(R_x_tot(ixp_tot)+delR_marg.lt.R_e) then
              relerr_allroots=1.d-12
!
              call find_all_roots(sepcross,R_x_tot(ixp_tot)+delR_marg,R_e,ierr)
!
              if(ierr.ne.0) then
                print *,'find_bounds_fixpoints: error in find_all_roots, sepcross 2'
                return
              endif
!
              nsc=nsc+nroots
            endif
          else
! X-point belongs to a different region:
            relerr_allroots=1.d-12
!
            call find_all_roots(sepcross,R_b,R_e,ierr)
!
            if(ierr.ne.0) then
              print *,'find_bounds_fixpoints: error in find_all_roots, sepcross 3'
              return
            endif
!
            nsc=nsc+nroots
          endif
        enddo
!
        n_x=all_regions(isig,ireg)%n_x
        nsc=nsc+n_x
        allocate(rsc_tmp(nsc),ipoi_tmp(nsc))
        nsc=0
!
        do ixp_tot=1,nxp_tot
          R_b=all_regions(isig,ireg)%R_b
          R_e=all_regions(isig,ireg)%R_e
!
          if(ipoi_x_tot(1,ixp_tot).eq.isig .and. &
             ipoi_x_tot(2,ixp_tot).eq.ireg) then
! X-point belongs to this region, exclude it from the search:
            if(R_b.lt.R_x_tot(ixp_tot)-delR_marg) then
              relerr_allroots=1.d-12
!
              call find_all_roots(sepcross,R_b,R_x_tot(ixp_tot)-delR_marg,ierr)
!
              if(ierr.ne.0) then
                print *,'find_bounds_fixpoints: error in find_all_roots, sepcross 4'
                return
              endif
!
              if(nroots.gt.0) then
                rsc_tmp(nsc+1:nsc+nroots)=roots
                nsc=nsc+nroots
              endif
            endif
            if(R_x_tot(ixp_tot)+delR_marg.lt.R_e) then
              relerr_allroots=1.d-12
!
              call find_all_roots(sepcross,R_x_tot(ixp_tot)+delR_marg,R_e,ierr)
!
              if(ierr.ne.0) then
                print *,'find_bounds_fixpoints: error in find_all_roots, sepcross 5'
                return
              endif
!
              if(nroots.gt.0) then
                rsc_tmp(nsc+1:nsc+nroots)=roots
                nsc=nsc+nroots
              endif
            endif
          else
! X-point belongs to a different region:
            relerr_allroots=1.d-12
!
            call find_all_roots(sepcross,R_b,R_e,ierr)
!
            if(ierr.ne.0) then
              print *,'find_bounds_fixpoints: error in find_all_roots, sepcross 6'
              return
            endif
!
            if(nroots.gt.0) then
              rsc_tmp(nsc+1:nsc+nroots)=roots
              nsc=nsc+nroots
            endif
          endif
!
        enddo
!
        allocate(logdummy(nsc+n_x))
        logdummy=.false.
!
! add X-points:
        if(n_x.gt.0) then
          rsc_tmp(nsc+1:nsc+n_x)=all_regions(isig,ireg)%R_x
          logdummy(nsc+1:nsc+n_x)=.true.
        endif
        nsc=nsc+n_x
!
        call sortin(rsc_tmp,ipoi_tmp,nsc)

! The same regular crossing can be returned by adjacent root brackets when it
! is a double root.  Keep one copy, but never merge an X-point with a regular
! crossing: those are different boundary types at the exact event.
        separatrix_root_tolerance=max(relerr_allroots* &
            (rpc_arr(npc)-rpc_arr(0)),256.d0*epsilon(1.d0)* &
            max(1.d0,abs(rpc_arr(0)),abs(rpc_arr(npc))))
        nsc_raw=nsc
        nsc=0
        do i=1,nsc_raw
          ixp_tot=ipoi_tmp(i)
          if(nsc.gt.0) then
            if(logdummy(ixp_tot).eqv.logdummy(ipoi_tmp(nsc)) .and. &
               abs(rsc_tmp(ixp_tot)-rsc_tmp(ipoi_tmp(nsc))).le. &
               separatrix_root_tolerance) cycle
          endif
          nsc=nsc+1
          ipoi_tmp(nsc)=ixp_tot
        enddo
!
        all_regions(isig,ireg)%n_sep=nsc
        allocate(all_regions(isig,ireg)%R_sep(nsc),all_regions(isig,ireg)%xpoint(nsc))
        all_regions(isig,ireg)%R_sep=rsc_tmp(ipoi_tmp)
        all_regions(isig,ireg)%xpoint=logdummy(ipoi_tmp)
        deallocate(rsc_tmp,ipoi_tmp,logdummy)
      endif
    enddo
  enddo
!
! End find separatrix crossings of the Poincare cut
!
!...........................
!
  if(wrbounds) then
!
    ireg0=1
!
    open(5003,file='opoints.dat')
    open(5004,file='xpoints.dat')
    open(5005,file='boundaries_beg.dat')
    open(5006,file='boundaries_end.dat')
    do ireg=ireg0,nregions
      do isig=1,2
        do i=1,all_regions(isig,ireg)%n_o
          write(5003,*) all_regions(isig,ireg)%R_o(i),all_regions(isig,ireg)%psiast_o(i)
        enddo
        do i=1,all_regions(isig,ireg)%n_x
          write(5004,*) all_regions(isig,ireg)%R_x(i),all_regions(isig,ireg)%psiast_x(i)
        enddo
        write(5005,*) all_regions(isig,ireg)%R_b,all_regions(isig,ireg)%psiast_b
        write(5006,*) all_regions(isig,ireg)%R_e,all_regions(isig,ireg)%psiast_e
      enddo
    enddo
    close(5003)
    close(5004)
    close(5005)
    close(5006)
  endif
!
!------------
!
  regions%nregions=nregions
  call move_alloc(all_regions,regions%all_regions)
!
!------------
!
  contains
!
!------------
!
  subroutine classify_extrema(R_b_in,R_e_in,ierr_classify)
!
! Find all extremum points of psi^* on the Poincare cut in the interval R_b_in < R < R_e_in,
! determines their types (X- or O-point). Results - cylindrical coordinates (R,Z), values
! of psi^* and types of fixed points - are stored in the arrays R_extr, Z_extr, psiast_extr
! and opoint_extr (logical array with .true. for O-points and .false. for X-points), respectively.
! Dimension of the arrays is nroots.
!
  use potato_topology_mod, only : root_has_two_sided_neighborhood,root_is_open_interval
  implicit none
!
  integer, intent(out) :: ierr_classify
  integer :: ierr, ierr_type, nroots1,nroots2,nvalid,i
  double precision :: R_b_in,R_e_in
  double precision :: vvrt,dvvrt_dx,dvrdz,p_phi
  double precision, dimension(:), allocatable :: R_tmp1,R_tmp2,R_tmp
  double precision, dimension(:), allocatable :: R_valid,Z_valid,psi_valid
  logical, dimension(:), allocatable :: opoint_valid
  logical :: opoint
!
  ierr_classify=0
  nroots=0
  if(allocated(opoint_extr)) then
    deallocate(opoint_extr,R_extr,Z_extr,psiast_extr)
  endif
  R_e=0.5d0*(R_b_in+R_e_in)
!
  R_b=R_b_in
  relerr_allroots=1.d-8
!
  call find_all_roots(get_vvert,0.d0,1.d0,ierr)
!
  if(ierr.ne.0) then
    ierr_classify=ierr
    return
  endif
!
  nroots1=nroots
  if(nroots1.gt.0) allocate(R_tmp1(nroots1))
!
  do i=1,nroots
    R_tmp1(i)=R_b+(R_e-R_b)*roots(i)**2
  enddo
!
  if(wrbounds) then
    do i=1,1000
      R=R_b+(R_e-R_b)*(1e-3*dble(i))**2
      call get_poicut(R,Z,dZ_dR)
      call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
      write(5000+isig,*) R,p_phi
    enddo
    write(5000+isig,*) 'NaN NaN'
  endif
!
  R_b=R_e_in
  relerr_allroots=1.d-8
!
  call find_all_roots(get_vvert,0.d0,1.d0,ierr)
!
  if(ierr.ne.0) then
    if(allocated(R_tmp1)) deallocate(R_tmp1)
    ierr_classify=ierr
    return
  endif
!
  nroots2=nroots
  if(nroots2.gt.0) allocate(R_tmp2(nroots2))
!
  do i=1,nroots
    R_tmp2(i)=R_b+(R_e-R_b)*roots(i)**2
  enddo
!
  if(wrbounds) then
    do i=1,1000
      R=R_b+(R_e-R_b)*(1e-3*dble(i))**2
      call get_poicut(R,Z,dZ_dR)
      call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
      write(5000+isig,*) R,p_phi
    enddo
    write(5000+isig,*) 'NaN NaN'
  endif
!
  nroots=nroots1+nroots2
  if(nroots.gt.0) allocate(R_tmp(nroots))
!
  if(nroots1.gt.0) then
    R_tmp(1:nroots1)=R_tmp1
    deallocate(R_tmp1)
  endif
!
  if(nroots2.gt.0) then
    R_tmp(nroots1+1:nroots)=R_tmp2(nroots2:1:-1)
    deallocate(R_tmp2)
  endif
!
  nvalid=0
  if(nroots.gt.0) then
    allocate(R_valid(nroots),Z_valid(nroots),psi_valid(nroots),opoint_valid(nroots))
    do i=1,nroots
      R=R_tmp(i)

! A zero of the cut velocity at an allowed-region endpoint is the physical
! turning/separatrix boundary, not an O- or X-point in the open interval.
! Keep the exact boundary in all_regions; do not pass it to the two-sided
! fixed-point classifier, whose finite differences would necessarily probe
! the forbidden side and report a false topology failure.
      if(.not.root_is_open_interval(R,R_b_in,R_e_in)) then
        call get_poicut(R,Z,dZ_dR)
        print *,'classify_extrema: endpoint stationary root H,J,sigma,Rlo,Rhi,R,Z = ', &
                toten,perpinv,sigma,R_b_in,R_e_in,R,Z
        cycle
      endif

      call get_poicut(R,Z,dZ_dR)
      call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
      if(ierr.ne.0) then
        deallocate(R_tmp,R_valid,Z_valid,psi_valid,opoint_valid)
        ierr_classify=ierr
        return
      endif
      call determine_fixpoint_type(R,Z,R_b_in,R_e_in,opoint,ierr_type)
      if(ierr_type.ne.0) then
        deallocate(R_tmp,R_valid,Z_valid,psi_valid,opoint_valid)
        ierr_classify=ierr_type
        return
      endif
      nvalid=nvalid+1
      R_valid(nvalid)=R
      Z_valid(nvalid)=Z
      psi_valid(nvalid)=p_phi
      opoint_valid(nvalid)=opoint
    enddo
    deallocate(R_tmp)
  endif

  nroots=nvalid
  if(nvalid.gt.0) then
    allocate(opoint_extr(nvalid),R_extr(nvalid),Z_extr(nvalid),psiast_extr(nvalid))
    opoint_extr=opoint_valid(1:nvalid)
    R_extr=R_valid(1:nvalid)
    Z_extr=Z_valid(1:nvalid)
    psiast_extr=psi_valid(1:nvalid)
  endif
  if(allocated(R_valid)) deallocate(R_valid,Z_valid,psi_valid,opoint_valid)
!
  end subroutine classify_extrema
!
!------------
!
  subroutine vparzero_vec(R,Z,vpar2,gradvpar2)
!
! Computes the square of parallel velocity and
! its gradient
!
  implicit none
!
  double precision :: bmod,sqrtg,phi_elec,R,Z,vpar2
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl,derphi
  double precision, dimension(2) :: gradvpar2
!
  x(1)=R
  x(2)=0.d0
  x(3)=Z
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
  call elefie(x,phi_elec,derphi)
!
  vpar2=toten-phi_elec-perpinv*bmod
  gradvpar2(1)=-derphi(1)-perpinv*bmod*bder(1)
  gradvpar2(2)=-derphi(3)-perpinv*bmod*bder(3)
!
  end subroutine vparzero_vec
!
!------------
!
  subroutine vparzero1D(R,vpar2,dvpar2_dR)
!
! Computes the square of parallel velocity at the Poincare cut
! and its derivative along the cut.
! To be used as formal argument in subroutine find_all_roots
! for finding forbidden boundaries
!
  implicit none
!
  double precision :: R,Z,vpar2,dvpar2_dR,dZ_dR
  double precision, dimension(2) :: gradvpar2
!
  call get_poicut(R,Z,dZ_dR)
!
  call vparzero_vec(R,Z,vpar2,gradvpar2)
!
  dvpar2_dR=gradvpar2(1)+gradvpar2(2)*dZ_dR
!
  end subroutine vparzero1D
!
!------------
!
  subroutine vvert(x,vvrt)
!
! Computes the normal guiding center velocity
! with respect to the Poincare cut
!
  use find_all_roots_mod, only : root_eval_valid,root_eval_error, &
                                  root_invalid_domain

  implicit none
!
  integer :: ierr
  double precision, parameter :: dtau=0.d0
  double precision :: x,vvrt,dZ_dR,dummy
  double precision, dimension(5) :: z,vz
  double precision, dimension(2) :: gradvpar2
!
  z=0.d0
  z(1)=R_b+(R_e-R_b)*x**2
!
  call get_poicut(z(1),z(3),dZ_dR)
  call get_z45(toten,perpinv,sigma,z,ierr)
  if(ierr.ne.0) then
    vvrt=0.d0
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
    return
  endif
  call velo(dtau,z,vz)
!
  vvrt=vz(3)-dZ_dR*vz(1)
!
  end subroutine vvert
!
!------------
!
  subroutine get_vvert(x,vvrt,dvvrt_dx)
!
! Computes numerically the derivatives of normal velocity
! (see above) and the velocity itself for root finiding.
! To be used as formal argument in subroutine find_all_roots.
!
  implicit none
!
  double precision, parameter :: hdiff=1.d-6
  double precision :: x,vvrt,dvvrt_dx,xm,xp,xx,vvrtm,vvrtp
!
  xm=max(0.d0,x-hdiff)
!
  call vvert(xm,vvrtm)
!
  xp=min(1.d0,x+hdiff)
!
  call vvert(xp,vvrtp)
!
  dvvrt_dx=(vvrtp-vvrtm)/(xp-xm)
!
  call vvert(x,vvrt)
!
  end subroutine get_vvert
!
!------------
!
  subroutine determine_fixpoint_type(R_in,Z_in,R_lo,R_hi,opoint,ierr_out)
!
! Determines type of fixed point with coordinates (R_in,Z_in)
! opoint = .true.   - O-point
! opoint = .false.  - X-point
!
  use potato_topology_mod, only : choose_two_sided_step
  implicit none
!
  double precision, parameter :: dtau=0.d0, hdiff=1.d-6, boundary_safety=0.5d0
  logical :: opoint
  integer :: ierr,ierr_out
  logical :: resolved_step
  double precision :: R_in,Z_in,R_lo,R_hi,dvrdr,dvrdz,dvzdr,dvzdz,hstep
  double precision, dimension(5) :: z,vz
!
  ierr_out=0
  opoint=.false.
  z(1)=R_in
  z(2)=0.d0
  z(3)=Z_in
!
  call get_z45(toten,perpinv,sigma,z,ierr)
  if(ierr.ne.0) then
    ierr_out=ierr
    return
  endif
  call choose_two_sided_step(R_in,R_lo,R_hi,hdiff*R_in,boundary_safety,hstep,resolved_step)
  if(.not.resolved_step) then
    print *,'determine_fixpoint_type: root lacks strict two-sided neighborhood H,J,sigma,Rlo,Rhi,R,Z,h,ierr = ', &
            toten,perpinv,sigma,R_lo,R_hi,R_in,Z_in,hstep,2
    ierr_out=2
    return
  endif

! A fixed-point type requires certified two-sided finite differences.  A
! forbidden perturbation is unresolved topology, not a one-sided derivative.
  z(1)=R_in-hstep
  call get_z45(toten,perpinv,sigma,z,ierr)
  if(ierr.ne.0) then
    print *,'determine_fixpoint_type: invalid R- neighborhood H,J,sigma,Rlo,Rhi,R,Z,h,ierr = ', &
            toten,perpinv,sigma,R_lo,R_hi,R_in,Z_in,hstep,ierr
    ierr_out=ierr
    return
  endif
  call velo(dtau,z,vz)
  dvrdr=vz(1)
  dvzdr=vz(3)

  z(1)=R_in+hstep
  call get_z45(toten,perpinv,sigma,z,ierr)
  if(ierr.ne.0) then
    print *,'determine_fixpoint_type: invalid R+ neighborhood H,J,sigma,Rlo,Rhi,R,Z,h,ierr = ', &
            toten,perpinv,sigma,R_lo,R_hi,R_in,Z_in,hstep,ierr
    ierr_out=ierr
    return
  endif
  call velo(dtau,z,vz)
  dvrdr=(vz(1)-dvrdr)/(2.d0*hstep)
  dvzdr=(vz(3)-dvzdr)/(2.d0*hstep)

  z(1)=R_in
  z(3)=Z_in-hstep
  call get_z45(toten,perpinv,sigma,z,ierr)
  if(ierr.ne.0) then
    print *,'determine_fixpoint_type: invalid Z- neighborhood H,J,sigma,Rlo,Rhi,R,Z,h,ierr = ', &
            toten,perpinv,sigma,R_lo,R_hi,R_in,Z_in,hstep,ierr
    ierr_out=ierr
    return
  endif
  call velo(dtau,z,vz)
  dvrdz=vz(1)
  dvzdz=vz(3)

  z(3)=Z_in+hstep
  call get_z45(toten,perpinv,sigma,z,ierr)
  if(ierr.ne.0) then
    print *,'determine_fixpoint_type: invalid Z+ neighborhood H,J,sigma,Rlo,Rhi,R,Z,h,ierr = ', &
            toten,perpinv,sigma,R_lo,R_hi,R_in,Z_in,hstep,ierr
    ierr_out=ierr
    return
  endif
  call velo(dtau,z,vz)
  dvrdz=(vz(1)-dvrdz)/(2.d0*hstep)
  dvzdz=(vz(3)-dvzdz)/(2.d0*hstep)
!
  if(4.d0*dvrdz*dvzdr+(dvrdr-dvzdz)**2.lt.0.d0) then
    opoint=.true.
  else
    opoint=.false.
  endif
!
  end subroutine determine_fixpoint_type
!
!------------
!
  subroutine sepcross(Rst,delpphi,ddelpphi_dRst)
!
! Roots of this function (delpphi(Rst)) are the values of cut parameter
! at separatrix crossings with the Poincare cut. To be used as formal
! argument in subroutine find_all_roots.
!
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use find_all_roots_mod, only : root_eval_valid,root_eval_error, &
                                 root_invalid_domain
  implicit none
!
  integer :: ierr
  double precision :: Rst,delpphi,ddelpphi_dRst
  double precision, dimension(5) :: z
!
  call starter_doublecount(toten,perpinv,sigma,Rst,   &
                           delpphi,ddelpphi_dRst,z,ierr)
!
  if(ierr.ne.0) then
    delpphi=0.d0
    ddelpphi_dRst=0.d0
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
    return
  endif
!
  delpphi=delpphi-psiast_x_tot(ixp_tot)
  if(.not.ieee_is_finite(delpphi) .or. &
     .not.ieee_is_finite(ddelpphi_dRst)) then
    delpphi=0.d0
    ddelpphi_dRst=0.d0
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
  endif
!
  end subroutine sepcross
!
!------------
!
  subroutine boucross(Rst,delpphi,ddelpphi_dRst)
!
! Roots of this function (delpphi(Rst)) are the values of cut parameter
! at the "most external with respect to rho_pol" orbit crossings with the
! Poincare cut. To be used as formal argument in subroutine find_all_roots.
!
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use find_all_roots_mod, only : root_eval_valid,root_eval_error, &
                                 root_invalid_domain
  implicit none
  integer :: ierr
  double precision :: Rst,delpphi,ddelpphi_dRst
  double precision, dimension(5) :: z
!
  call starter_doublecount(toten,perpinv,sigma,Rst,   &
                           delpphi,ddelpphi_dRst,z,ierr)
!
  if(ierr.ne.0) then
    delpphi=0.d0
    ddelpphi_dRst=0.d0
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
    return
  endif
!
  delpphi=delpphi-pphi_minmax
  if(.not.ieee_is_finite(delpphi) .or. &
     .not.ieee_is_finite(ddelpphi_dRst)) then
    delpphi=0.d0
    ddelpphi_dRst=0.d0
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
  endif
!
  end subroutine boucross
!
!------------
!
  subroutine boucross_with_endpoint_limit(Rst,delpphi,ddelpphi_dRst)
!
! At a turning endpoint starter_doublecount is undefined although psi^*
! has a finite one-sided limit.  Supplying that limit keeps the clipping
! root search continuous at the endpoint.
!
  implicit none
  double precision :: Rst,delpphi,ddelpphi_dRst
!
  if(Rst.eq.R_b_in) then
    delpphi=all_regions(isig,ireg)%psiast_b-pphi_minmax
    ddelpphi_dRst=0.d0
    root_eval_valid=.true.
    root_eval_error=0
    return
  endif
  if(Rst.eq.R_e_in) then
    delpphi=all_regions(isig,ireg)%psiast_e-pphi_minmax
    ddelpphi_dRst=0.d0
    root_eval_valid=.true.
    root_eval_error=0
    return
  endif
!
  call boucross(Rst,delpphi,ddelpphi_dRst)
!
  end subroutine boucross_with_endpoint_limit
!
!------------
!
  end subroutine find_bounds_fixpoints
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine form_classes_doublecount(regions,classes_talk,ierr)
!
! Here "class" means one connected allowed interval on the Poincare cut, with
! its two boundary/homoclinic-pair records and the parallel-velocity sign.  It
! is not a trapped/passing/winding label.  The routine determines these
! connected segments and their boundary types (usual boundary, vpar2=0
! boundary or separatrix crossing); types only select the parameterization.
!
  use form_classes_doublecount_mod, only : nclasses,ifuntype,sigma_class,  &
                                           R_class_beg,R_class_end
  use bounds_fixpoints_mod,         only : region_set_t
  use poicut_mod,                   only : Rbou_lfs,Zbou_lfs,Rbou_hfs,Zbou_hfs
  use global_invariants,            only : toten,perpinv,sigma
!
  implicit none
!
  type(region_set_t), intent(in) :: regions
  logical :: classes_talk
  integer          :: ierr,isig,ireg,nxp,ixp,icl,i
  integer, dimension(:), allocatable :: ibt_b,ibt_e
!
  associate(nregions => regions%nregions, all_regions => regions%all_regions)
!
  ierr=0
!
! Count the number of classes:
!
  nclasses=0
!
  do ireg=1,nregions
    do isig=1,2
      if(all_regions(isig,ireg)%within_rhopol) then
! Number of classes per allowed motion region is number of separatrix crossings
! (including X-points) plus one:
        nclasses=nclasses+all_regions(isig,ireg)%n_sep+1
      endif
    enddo
  enddo
!
! End count the number of classes
!
  if(allocated(ifuntype)) deallocate(ifuntype,sigma_class,R_class_beg,R_class_end)
  allocate(ifuntype(nclasses),sigma_class(nclasses),  &
           R_class_beg(nclasses),R_class_end(nclasses))
  allocate(ibt_b(nclasses),ibt_e(nclasses))
!
! Set the boundaries of class domains (cut parameter values [R_b,R_e]
! and boundary types [ibt_b,ibt_e]:
! Boundary types: 1 - rho_pol bondary with vpar2 > 0
!                 2 - inner bondary with vpar2 = 0
!                 3 - regular separatrix crossing (trace of X-point)
!                 4 - X-point (singular separatrix crossing)
!
!
  nclasses=0
!
  do ireg=1,nregions
    do isig=1,2
!
      if(isig.eq.1) then
        sigma=1.d0
      else
        sigma=-1.d0
      endif
!
      if(all_regions(isig,ireg)%within_rhopol) then
        nclasses=nclasses+1
        sigma_class(nclasses)=sigma
        R_class_beg(nclasses)=all_regions(isig,ireg)%R_b
        if(all_regions(isig,ireg)%inner_b) then
          ibt_b(nclasses)=2
        else
          ibt_b(nclasses)=1
        endif
        nxp=all_regions(isig,ireg)%n_sep
        do ixp=1,nxp
          R_class_end(nclasses)=all_regions(isig,ireg)%R_sep(ixp)
          if(all_regions(isig,ireg)%xpoint(ixp)) then
            ibt_e(nclasses)=4
            ibt_b(nclasses+1)=4
          else
            ibt_e(nclasses)=3
            ibt_b(nclasses+1)=3
          endif
          nclasses=nclasses+1
          sigma_class(nclasses)=sigma
          R_class_beg(nclasses)=all_regions(isig,ireg)%R_sep(ixp)
        enddo
        R_class_end(nclasses)=all_regions(isig,ireg)%R_e
        if(all_regions(isig,ireg)%inner_e) then
          ibt_e(nclasses)=2
        else
          ibt_e(nclasses)=1
        endif
      endif
!
    enddo
  enddo
!
! End set the boundaries of class domains
!
! Set the index of class parameterization function $\xi_{ij}$, Eq.(57)
! tens - left boundary type, units - right boundary type:
  do icl=1,nclasses
    ifuntype(icl)=ibt_b(icl)*10+ibt_e(icl)
  enddo
!
  if(classes_talk) then
    print *,'number of classes = ',nclasses
    do i=1,nclasses
      print *,R_class_beg(i),R_class_end(i),ifuntype(i)
    enddo
  endif
!
  deallocate(ibt_b,ibt_e)
!
  end associate
!
  end subroutine form_classes_doublecount
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine xi_func(ifuntype,x,xi,dxi_dx)
!
! Class parameterization functions $\xi_{ij}$, Eq.(57)
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
!
  integer :: ifuntype
  double precision :: x,xi,dxi_dx
!
  select case(ifuntype)
  case(11)
! two rho_pol boundaries, 0<x<1
    xi=x
    dxi_dx=1.d0
  case(12)
! left- rho_pol boundary, right - inner boundary, 0<x<1
    xi=1.d0-(1.d0-x)**2
    dxi_dx=2.d0*(1.d0-x)
  case(13,14)
! left- rho_pol boundary, right - X-point, 0<x<inf
    xi=tanh(x)
    dxi_dx=1.d0/cosh(x)**2
  case(21)
! left- inner boundary, right - rho_pol boundary, 0<x<1
    xi=x**2
    dxi_dx=2.d0*x
  case(22)
! two inner boundaries, 0<x<1
    xi=0.5d0*(1.d0-cos(pi*x))
    dxi_dx=0.5d0*pi*sin(pi*x)
  case(23,24)
! left- inner boundary, right - X-point, 0<x<inf
    xi=tanh(x)**2
    dxi_dx=2.d0*tanh(x)/cosh(x)**2
  case(31,41)
! left- X-point, right - rho_pol boundary, -inf<x<0
    xi=tanh(x)+1.d0
    dxi_dx=1.d0/cosh(x)**2
  case(32,42)
! left- X-point, right - inner boundary, -inf<x<0
    xi=1.d0-tanh(x)**2
    dxi_dx=-2.d0*tanh(x)/cosh(x)**2
  case(33,34,43,44)
! two X-points, -inf<x<inf
    xi=0.5d0*(1.d0+tanh(x))
    dxi_dx=0.5d0/cosh(x)**2
  case default
    print *,'xi_func: wrong function type'
    stop
  end select
!
  return
!
  end subroutine xi_func
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine classbounds(ifuntype,relmargin,widthclass,xbeg,xend)
!
! Sets the boundaries over class parameter x for
!
  implicit none
!
  integer :: ifuntype
  double precision :: relmargin,widthclass,xbeg,xend,turn_x,turn_margin

  ! relmargin is a physical cut-coordinate margin.  The inner-boundary
  ! charts use xi~x**2, so using relmargin directly in x would place the
  ! orbit at relmargin**2 and make the ODE start indistinguishable from the
  ! turning point.  These are the exact inverse chart coordinates.
  turn_margin=sqrt(relmargin)
  turn_x=atanh(turn_margin)
!
  select case(ifuntype)
  case(11)
! two rho_pol boundaries, 0<x<1
    xbeg=0.d0
    xend=1.d0
  case(12)
! left- rho_pol boundary, right - inner boundary, 0<x<1
    xbeg=0.d0
    xend=1.d0-turn_margin
  case(13)
! left- rho_pol boundary, right - X-point, 0<x<inf
    xbeg=0.d0
    xend=-0.5d0*log(relmargin/widthclass)
  case(14)
! left- rho_pol boundary, right - X-point, 0<x<inf
    xbeg=0.d0
    xend=-0.5d0*log(relmargin/widthclass)*0.5d0
  case(21)
! left- inner boundary, right - rho_pol boundary, 0<x<1
    xbeg=turn_margin
    xend=1.d0
  case(22)
! two inner boundaries, 0<x<1
    turn_x=acos(1.d0-2.d0*relmargin)/acos(-1.d0)
    xbeg=turn_x
    xend=1.d0-turn_x
  case(23)
! left- inner boundary, right - X-point, 0<x<inf
    xbeg=turn_x
    xend=-0.5d0*log(relmargin/widthclass)
  case(24)
! left- inner boundary, right - X-point, 0<x<inf
    xbeg=turn_x
    xend=-0.25d0*log(relmargin/widthclass)
  case(31)
! left- X-point, right - rho_pol boundary, -inf<x<0
    xbeg=0.5d0*log(relmargin/widthclass)
    xend=0.d0
  case(32)
! left- X-point, right - inner boundary, -inf<x<0
    xbeg=0.5d0*log(relmargin/widthclass)
    xend=-turn_x
  case(33)
! two X-points, -inf<x<inf
    xbeg=0.5d0*log(relmargin/widthclass)
    xend=-xbeg
  case(34)
! two X-points, -inf<x<inf
    xbeg=0.5d0*log(relmargin/widthclass)
    xend=-xbeg*0.5d0
  case(41)
! left- X-point, right - rho_pol boundary, -inf<x<0
    xbeg=0.25d0*log(relmargin/widthclass)
    xend=0.d0
  case(42)
! left- X-point, right - inner boundary, -inf<x<0
    xbeg=0.25d0*log(relmargin/widthclass)
    xend=-turn_x
  case(43)
! two X-points, -inf<x<inf
    xend=-0.5d0*log(relmargin/widthclass)
    xbeg=-xend*0.5d0
  case(44)
! two X-points, -inf<x<inf
    xbeg=0.25d0*log(relmargin/widthclass)
    xend=-xbeg
  case default
    xbeg=0.d0
    xend=0.d0
  end select
!
  end subroutine classbounds
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine get_matrix_doublecount
!
! Computation of normalized toroidal moment, normalizaed bounce time,
! toroidal displacement per bounce time and bounce integrals as functions
! of class parameter x for adaptive refinement of interpolation grid over x
!
  use sample_matrix_mod, only : n1,n2,x,amat,matrix_eval_valid, &
                                matrix_eval_error,matrix_eval_success, &
                                matrix_eval_starter_failure, &
                                matrix_eval_orbit_failure,matrix_eval_wall_loss, &
                                matrix_eval_nonfinite, &
                                matrix_boundary_error
  use orbit_dim_mod,     only : neqm,next,orbit_wall_loss,orbit_failure_stage
  use get_matrix_mod,    only : iclass
  use global_invariants, only : dtau,toten,perpinv
  use form_classes_doublecount_mod, only : ifuntype,R_class_beg,R_class_end,sigma_class
!
  implicit none
!
  logical :: fullbounce
  integer :: ierr,k
  double precision :: psiast,dpsiast_dRst,taub,delphi,xi,dxi_dx,Rst,sigma,delta_R, &
                      tau_fr,dphi_fr
  double precision, dimension(neqm) :: z
  double precision, dimension(next) :: extraset
!
  external :: velo,velo_pphint

  matrix_eval_valid=.true.
  matrix_eval_error=matrix_eval_success
!
  sigma=sigma_class(iclass)
  delta_R=R_class_end(iclass)-R_class_beg(iclass)
!
  call xi_func(ifuntype(iclass),x,xi,dxi_dx)
!
!
  Rst=R_class_beg(iclass)+delta_R*xi
!
  call starter_doublecount(toten,perpinv,sigma,Rst,   &
                           psiast,dpsiast_dRst,z,ierr)
!
  if(ierr.ne.0) then
    matrix_eval_valid=.false.
    if(orbit_wall_loss) then
      matrix_eval_error=matrix_eval_wall_loss
    else
      matrix_eval_error=matrix_eval_starter_failure
    endif
    print *,'get_matrix_doublecount: starter failure ierr = ',ierr
    return
  endif
!
! Temporary A/B experiment: use the existing two-pass first-return map.
! This commit is diagnostic only and will be reverted after the probe run.
  fullbounce=.false.
!
  if(next.eq.0) then
    if(fullbounce) then
!
      call find_bounce(next,velo,dtau,z,taub,delphi,extraset,ierr)
      if(ierr.ne.0) then
        matrix_eval_valid=.false.
        if(orbit_wall_loss) then
          matrix_eval_error=matrix_eval_wall_loss
        else
          matrix_eval_error=matrix_eval_orbit_failure
        endif
        return
      endif
!
    else
!
      call first_return_map(sigma,Rst,sigma,Rst,taub,delphi,ierr)
      if(ierr.ne.0) then
        matrix_eval_valid=.false.
        matrix_eval_error=matrix_eval_orbit_failure
        return
      endif
      call first_return_map(sigma,Rst,sigma,Rst,tau_fr,dphi_fr,ierr)
      if(ierr.ne.0) then
        matrix_eval_valid=.false.
        matrix_eval_error=matrix_eval_orbit_failure
        return
      endif
!
      taub=taub+tau_fr
      delphi=delphi+dphi_fr
    endif
  else
    extraset=0.d0
!
    call find_bounce(next,velo_pphint,dtau,z,taub,delphi,extraset,ierr)
    if(ierr.ne.0) then
      matrix_eval_valid=.false.
      if(orbit_wall_loss) then
        matrix_eval_error=matrix_eval_wall_loss
      else
        matrix_eval_error=matrix_eval_orbit_failure
      endif
      return
    endif
!
  endif
!
  amat(1,1)=psiast
  amat(2,1)=taub
  amat(3,1)=delphi
!
  if(next.gt.0) then
    amat(4:3+next,1)=extraset
  endif

  if(any(amat.ne.amat)) then
    matrix_eval_valid=.false.
    matrix_eval_error=matrix_eval_nonfinite
    return
  endif
  if(any(abs(amat).gt.huge(amat))) then
    matrix_eval_valid=.false.
    matrix_eval_error=matrix_eval_nonfinite
    return
  endif
!
  end subroutine get_matrix_doublecount
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine sample_class_doublecount(iunit,ierr)
!
! Formation of adaptive grid over class parameter x together with data
! for adaptive interpolation of normalized toroidal moment, normalizaed bounce time,
! toroidal displacement per bounce time and bounce integrals
!
  use sample_matrix_mod, only : nlagr,n1,n2,itermax,eps,xbeg,xend,  &
                                npoi,x,xarr,amat_arr,matrix_boundary_error, &
                                matrix_eval_error, &
                                matrix_eval_success
  use sample_class_status_mod, only : sample_class_success, &
                                      sample_class_no_resonance
  use orbit_dim_mod,     only : next,numbasef,orbit_failure_stage
  use get_matrix_mod,    only : relerror,relmargin,iclass,delphi_max
  use global_invariants, only : dtau,toten,perpinv,sigma
  use form_classes_doublecount_mod, only : ifuntype,sigma_class, &
                                           R_class_beg,R_class_end
  use cc_mod, only : dowrite
  use interp_cache_mod,  only : interp_cache_reset
!
  implicit none
!
  integer :: iunit,ierr,i
  double precision :: psiastbeg,psiastend,widthclass
  logical :: empty_class
!
  external :: get_matrix_doublecount
!
  nlagr=7         !<= temporary place, should be moved out for centralized input
  relerror=1.d-3  !<= temporary place, should be moved out for centralized input
  relmargin=1.d-7 !<= temporary place, should be moved out for centralized input
  itermax=20      !<= temporary place, should be moved out for centralized input
!
  next=2*numbasef
  n1=3+next
  n2=1
!
  ierr=sample_class_success
  matrix_boundary_error=matrix_eval_success
!
  widthclass=abs(R_class_end(iclass)/R_class_beg(iclass)-1.d0)
!
  if(widthclass/relmargin.lt.1.d0) then
    print *,'ignore class'
    ierr=sample_class_no_resonance
    return
  endif
!
  call classbounds(ifuntype(iclass),relmargin,widthclass,xbeg,xend)
!
!
  empty_class=.false.
  if(delphi_max.gt.0.d0) then
    call bound_class_delphi(ifuntype(iclass),xbeg,xend,empty_class)
  endif
  if(matrix_boundary_error.ne.matrix_eval_success) then
    ierr=matrix_boundary_error
    return
  endif
  if(empty_class) then
    ! The finite harmonic guard certifies that this class has no roots.  It is
    ! a zero contribution, not an adaptive-sampler failure.
    ierr=sample_class_no_resonance
    call interp_cache_reset
    return
  endif
!
  call bound_class_return(ifuntype(iclass),xbeg,xend)
  if(matrix_boundary_error.ne.matrix_eval_success) then
    ierr=matrix_boundary_error
    return
  endif
!
  call bound_class_wall(xbeg,xend)
  if(matrix_boundary_error.ne.matrix_eval_success) then
    ierr=matrix_boundary_error
    return
  endif
!
  eps=relerror
!
  call sample_matrix(get_matrix_doublecount,ierr)
  if(ierr.ne.0) then
    print *, 'sample_class_doublecount: class=',iclass, &
      ' iftype=',ifuntype(iclass), &
      ' sigma=',sigma_class(iclass), &
      ' Rbeg,Rend=',R_class_beg(iclass),R_class_end(iclass), &
      ' xbeg,xend=',xbeg,xend, &
      ' delphi_max=',delphi_max, &
      ' sample_x,error,ierr,stage=',x,matrix_eval_error,ierr,orbit_failure_stage
    print *, 'sample_nodes npoi=',npoi
    do i=max(1,npoi-8),npoi
      print *, 'sample_node=',i,xarr(i),amat_arr(:,1,i)
    enddo
  endif
!
! The grid (xarr,amat_arr,npoi) and iclass just changed; drop memoized entries
! so interpolate_class_doublecount never returns a value from the old grid.
  call interp_cache_reset
!
  if(dowrite) then
    print *,'npoi = ',npoi
    do i=1,npoi
      write(iunit,*) xarr(i),amat_arr(:,1,i)
    enddo
  endif
!
  end subroutine sample_class_doublecount
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine bound_class_delphi(iftype,xb,xe,empty_class)
!
! Trim the X-point side(s) of a class so the sampled domain stays where
! |delphi_b| <= delphi_max.  delphi_b diverges logarithmically toward the
! X-point; past 2*pi*m_max/n no resonance exists, yet that endpoint is what
! makes sample_matrix exceed itermax (ierr=2 -> the class is dropped) and floods
! the root search with non-physical dense roots.  The endpoint is located by
! bisection on |delphi_b(x)| using get_matrix_doublecount as the orbit evaluator.
!
  use sample_matrix_mod, only : n1,n2,x,amat,matrix_eval_valid, &
                                matrix_eval_error,matrix_eval_success, &
                                matrix_eval_wall_loss,matrix_eval_nonfinite, &
                                matrix_boundary_error, &
                                matrix_failure_is_open_boundary
  use get_matrix_mod,    only : delphi_max
  use resonance_mode_bounds_mod, only : resonance_no_root_for_any, &
                                        harmonic_guard_success
!
  implicit none
!
  integer :: iftype,ib_beg,ib_end
  double precision :: xb,xe,xmid
  logical :: amat_was_alloc,empty_class,empty_beg,empty_end
!
  external :: get_matrix_doublecount
!
  empty_class=.false.
  amat_was_alloc=allocated(amat)
  if(.not.amat_was_alloc) allocate(amat(n1,n2))
!
! Boundary types 3 (regular separatrix crossing) and 4 (X-point) are the
! divergent ends; ifuntype packs them as ib_beg*10+ib_end.
  ib_beg=iftype/10
  ib_end=mod(iftype,10)
!
  if(ib_beg.ge.3 .and. ib_end.ge.3) then
    xmid=0.5d0*(xb+xe)
    empty_end=.false.
    empty_beg=.false.
    call trim_endpoint(xmid,xe,empty_end)
    call trim_endpoint(xmid,xb,empty_beg)
    empty_class=empty_beg .and. empty_end
  elseif(ib_end.ge.3) then
    call trim_endpoint(xb,xe,empty_class)
  elseif(ib_beg.ge.3) then
    call trim_endpoint(xe,xb,empty_class)
  endif
!
  if(.not.amat_was_alloc) deallocate(amat)
!
  contains
!
  subroutine trim_endpoint(xsafe,xdiv,empty_interval)
  double precision :: xsafe,xdiv,xroot,xnoroot,xm
  double precision :: delphi_safe,delphi_div,delphi_mid
  integer :: it,guard_ierr
  logical :: no_root,empty_interval
  !
  empty_interval=.false.
  delphi_div=eval_delphi(xdiv)
  if(matrix_boundary_error.ne.matrix_eval_success) return
  call resonance_no_root_for_any(delphi_div,no_root,guard_ierr)
  if(guard_ierr.ne.harmonic_guard_success) then
    matrix_boundary_error=matrix_eval_nonfinite
    return
  endif
  if(.not.no_root) return

  delphi_safe=eval_delphi(xsafe)
  if(matrix_boundary_error.ne.matrix_eval_success) return
  call resonance_no_root_for_any(delphi_safe,no_root,guard_ierr)
  if(guard_ierr.ne.harmonic_guard_success) then
    matrix_boundary_error=matrix_eval_nonfinite
    return
  endif
  if(no_root) then
    empty_interval=.true.
    return
  endif
!
  ! Keep the two sides named by their physical meaning.  Class coordinates
  ! need not increase toward the divergent endpoint.
  xroot=xsafe
  xnoroot=xdiv
  do it=1,40
    xm=0.5d0*(xroot+ xnoroot)
    delphi_mid=eval_delphi(xm)
    if(matrix_boundary_error.ne.matrix_eval_success) return
    call resonance_no_root_for_any(delphi_mid,no_root,guard_ierr)
    if(guard_ierr.ne.harmonic_guard_success) then
      matrix_boundary_error=matrix_eval_nonfinite
      return
    endif
    if(no_root) then
      xnoroot=xm
    else
      xroot=xm
    endif
    if(abs(xnoroot-xroot).lt.1.d-3*max(1.d0,abs(xdiv))) exit
  enddo
  ! The sampled interval ends on the resonant side of the bracket.  The
  ! no-root side is the X-point-facing numerical singularity we are excluding.
  xdiv=xroot
  end subroutine trim_endpoint
!
  double precision function eval_delphi(xval)
  double precision :: xval
  matrix_eval_valid=.true.
  matrix_eval_error=matrix_eval_success
  matrix_boundary_error=matrix_eval_success
  x=xval
  call get_matrix_doublecount
  if(matrix_failure_is_open_boundary(matrix_eval_error)) then
    ! An orbit-return or wall failure is the expected invalid side of a
    ! separatrix-facing class endpoint.  Other failures remain fatal.
    eval_delphi=sign(huge(1.d0),1.d0)
  elseif(matrix_eval_error.ne.matrix_eval_success) then
    matrix_boundary_error=matrix_eval_error
    eval_delphi=0.d0
  else
    eval_delphi=amat(3,1)
  endif
  end function eval_delphi
!
  end subroutine bound_class_delphi
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine bound_class_return(iftype,xb,xe)
!
! A type-1 cut endpoint is a computational rho boundary, not a guarantee
! that a finite-width orbit closes.  Keep all endpoints whose full orbit
! returns, and trim only an endpoint interval with an orbit-closure failure.
! This also preserves SOL excursions when edge_extension is enabled: a
! leaving-and-returning orbit evaluates successfully and is never trimmed.
!
  use sample_matrix_mod, only : n1,n2,x,amat,matrix_eval_valid, &
                                matrix_eval_error,matrix_eval_success, &
                                matrix_eval_orbit_failure, &
                                matrix_boundary_error
!
  implicit none
!
  integer :: iftype,ib_beg,ib_end
  double precision :: xb,xe,xmid
  logical :: amat_was_alloc
!
  external :: get_matrix_doublecount
!
  ib_beg=iftype/10
  ib_end=mod(iftype,10)
  if(ib_beg.ne.1 .and. ib_end.ne.1) return
!
  amat_was_alloc=allocated(amat)
  if(.not.amat_was_alloc) allocate(amat(n1,n2))
  xmid=0.5d0*(xb+xe)
  if(orbit_return_invalid(xmid)) then
    matrix_boundary_error=matrix_eval_orbit_failure
  elseif(ib_end.eq.1) then
    call trim_return(xmid,xe)
    if(matrix_boundary_error.eq.matrix_eval_success .and. ib_beg.eq.1) then
      call trim_return(xmid,xb)
    endif
  elseif(ib_beg.eq.1) then
    call trim_return(xmid,xb)
  endif
!
  if(.not.amat_was_alloc) deallocate(amat)
!
  contains
!
  subroutine trim_return(xsafe,xdiv)
  double precision :: xsafe,xdiv,xvalid,xinvalid,xm
  integer :: it
!
  if(.not.orbit_return_invalid(xdiv)) return
  if(orbit_return_invalid(xsafe)) then
    matrix_boundary_error=matrix_eval_orbit_failure
    return
  endif
  xvalid=xsafe
  xinvalid=xdiv
  do it=1,40
    xm=0.5d0*(xvalid+xinvalid)
    if(orbit_return_invalid(xm)) then
      xinvalid=xm
    else
      xvalid=xm
    endif
    if(abs(xinvalid-xvalid).lt.1.d-3*max(1.d0,abs(xdiv))) exit
  enddo
  xdiv=xvalid
  end subroutine trim_return
!
  logical function orbit_return_invalid(xval)
  double precision :: xval
!
  matrix_eval_valid=.true.
  matrix_eval_error=matrix_eval_success
  matrix_boundary_error=matrix_eval_success
  x=xval
  call get_matrix_doublecount
  orbit_return_invalid=(matrix_eval_error.eq.matrix_eval_orbit_failure)
  if(matrix_eval_error.ne.matrix_eval_success .and. &
     matrix_eval_error.ne.matrix_eval_orbit_failure) then
    matrix_boundary_error=matrix_eval_error
    orbit_return_invalid=.false.
  endif
  end function orbit_return_invalid
!
  end subroutine bound_class_return
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine bound_class_wall(xb,xe)
!
! Trim the class interval to where the orbit closes inside the limiter wall.
! A wall-crossing orbit fails in get_matrix_doublecount (amat(3,1) left at the
! huge sentinel); without this trim a single lost sample makes sample_matrix
! exceed itermax and the whole class is dropped (ierr=2) instead of only the
! lost orbits.  Trims both endpoints by bisection on orbit validity, so the
! inside-wall part of the class is still sampled and its resonances kept.
!
  use sample_matrix_mod, only : n1,n2,x,amat,matrix_eval_valid, &
                                matrix_eval_error,matrix_eval_success, &
                                matrix_eval_wall_loss,matrix_boundary_error
  use wall_loss_mod,     only : wall_loaded
!
  implicit none
!
  double precision :: xb,xe,xmid
  logical :: amat_was_alloc
!
  external :: get_matrix_doublecount
!
  if(.not.wall_loaded) return
!
  amat_was_alloc=allocated(amat)
  if(.not.amat_was_alloc) allocate(amat(n1,n2))
!
  xmid=0.5d0*(xb+xe)
  if(.not.orbit_lost(xmid)) then
    call trim_lost(xmid,xe)
    call trim_lost(xmid,xb)
  endif
!
  if(.not.amat_was_alloc) deallocate(amat)
!
  contains
!
  subroutine trim_lost(xsafe,xdiv)
  double precision :: xsafe,xdiv,xlo,xhi,xm
  integer :: it
!
  if(.not.orbit_lost(xdiv)) return
!
  xlo=xsafe
  xhi=xdiv
  do it=1,40
    xm=0.5d0*(xlo+xhi)
    if(orbit_lost(xm)) then
      xhi=xm
    else
      xlo=xm
    endif
    if(abs(xhi-xlo).lt.1.d-3*max(1.d0,abs(xdiv))) exit
  enddo
  xdiv=xlo
  end subroutine trim_lost
!
  logical function orbit_lost(xval)
  double precision :: xval
  matrix_eval_valid=.true.
  matrix_eval_error=matrix_eval_success
  x=xval
  call get_matrix_doublecount
  orbit_lost=(matrix_eval_error.eq.matrix_eval_wall_loss)
  if(matrix_eval_error.ne.matrix_eval_success .and. &
     matrix_eval_error.ne.matrix_eval_wall_loss) then
    matrix_boundary_error=matrix_eval_error
    orbit_lost=.false.
  endif
  end function orbit_lost
!
  end subroutine bound_class_wall
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine interpolate_class_doublecount(x,vec,dvec)
!
! Interpolation on the adaptive grid over class parameter x of normalized toroidal
! moment, normalizaed bounce time, toroidal displacement per bounce time and bounce
! integrals (vec) and their derivatives (dvec)
!
  use sample_matrix_mod, only : nlagr,n1,npoi,xarr,amat_arr
  use get_matrix_mod,    only : iclass
  use form_classes_doublecount_mod, only : ifuntype
  use interp_cache_mod,  only : ncache,xcache,veccache,dveccache
!
  implicit none
!
  integer, parameter :: nder=1
  integer            :: k,npoilag,nshift,ibeg,iend,ic
  double precision   :: x,xin,xi,xib,dxi_dx,dxi_dxb,dpphi_dxib,xi_inf,a,b
!
  double precision, dimension(n1)             :: vec,dvec
  double precision, dimension(0:nder,nlagr+1) :: coef
!
! Exact memoization (see interp_cache_mod): a hit at this x returns the stored
! vec,dvec bit-for-bit; the scan is far cheaper than the plag_coeff/matmul work.
  do ic=1,ncache
    if(xcache(ic).eq.x) then
      vec=veccache(1:n1,ic)
      dvec=dveccache(1:n1,ic)
      return
    endif
  enddo
!
  npoilag=nlagr+1
  nshift=nlagr/2
!
  call binsrc(xarr,1,npoi,x,k)
!
  ibeg=max(1,k-nshift)
  iend=min(npoi,ibeg+nlagr)
  ibeg=iend-nlagr
!
  xin=max(xarr(1),min(xarr(npoi),x))
!
  call plag_coeff(npoilag,nder,xin,xarr(ibeg:iend),coef)
!
  vec=matmul(amat_arr(:,1,ibeg:iend),coef(0,:))
  dvec=matmul(amat_arr(:,1,ibeg:iend),coef(1,:))
!
  if(x.lt.xarr(1)) then
    vec=amat_arr(:,1,1)+(x-xarr(1))*dvec
!
    select case(ifuntype(iclass))
    case(31)
! left- separatrix, right - rho_pol boundary, -inf<x<0
      xi=tanh(x)+1.d0
      xib=tanh(xarr(1))+1.d0
      dxi_dx=1.d0/cosh(x)**2
      dxi_dxb=1.d0/cosh(xarr(1))**2
      dpphi_dxib=dvec(1)/dxi_dxb
      vec(1)=amat_arr(1,1,1)+dpphi_dxib*(xi-xib)
      dvec(1)=dpphi_dxib*dxi_dx
    case(41)
! left- X-point, right - rho_pol boundary, -inf<x<0
      xi=tanh(x)+1.d0
      xib=tanh(xarr(1))+1.d0
      dxi_dx=1.d0/cosh(x)**2
      dxi_dxb=1.d0/cosh(xarr(1))**2
      xi_inf=0.d0
      b=dvec(1)/(2.d0*(xib-xi_inf)*dxi_dxb)
      a=amat_arr(1,1,1)-b*(xib-xi_inf)**2
      vec(1)=a+b*(xi-xi_inf)**2
      dvec(1)=2.d0*b*(xi-xi_inf)*dxi_dx
    case(32)
! left- separatrix, right - inner boundary, -inf<x<0
      xi=1.d0-tanh(x)**2
      xib=1.d0-tanh(xarr(1))**2
      dxi_dx=-2.d0*tanh(x)/cosh(x)**2
      dxi_dxb=-2.d0*tanh(xarr(1))/cosh(xarr(1))**2
      dpphi_dxib=dvec(1)/dxi_dxb
      vec(1)=amat_arr(1,1,1)+dpphi_dxib*(xi-xib)
      dvec(1)=dpphi_dxib*dxi_dx
    case(42)
! left- X-point, right - inner boundary, -inf<x<0
      xi=1.d0-tanh(x)**2
      xib=1.d0-tanh(xarr(1))**2
      dxi_dx=-2.d0*tanh(x)/cosh(x)**2
      dxi_dxb=-2.d0*tanh(xarr(1))/cosh(xarr(1))**2
      xi_inf=0.d0
      b=dvec(1)/(2.d0*(xib-xi_inf)*dxi_dxb)
      a=amat_arr(1,1,1)-b*(xib-xi_inf)**2
      vec(1)=a+b*(xi-xi_inf)**2
      dvec(1)=2.d0*b*(xi-xi_inf)*dxi_dx
    case(33,34)
! left - separatrix, right - separatrix or X-point, -inf<x<inf
      xi=0.5d0*(1.d0+tanh(x))
      xib=0.5d0*(1.d0+tanh(xarr(1)))
      dxi_dx=0.5d0/cosh(x)**2
      dxi_dxb=0.5d0/cosh(xarr(1))**2
      dpphi_dxib=dvec(1)/dxi_dxb
      vec(1)=amat_arr(1,1,1)+dpphi_dxib*(xi-xib)
      dvec(1)=dpphi_dxib*dxi_dx
    case(43,44)
! left - X-point, right - separatrix or X-point, -inf<x<inf
      xi=0.5d0*(1.d0+tanh(x))
      xib=0.5d0*(1.d0+tanh(xarr(1)))
      dxi_dx=0.5d0/cosh(x)**2
      dxi_dxb=0.5d0/cosh(xarr(1))**2
      xi_inf=0.d0
      b=dvec(1)/(2.d0*(xib-xi_inf)*dxi_dxb)
      a=amat_arr(1,1,1)-b*(xib-xi_inf)**2
      vec(1)=a+b*(xi-xi_inf)**2
      dvec(1)=2.d0*b*(xi-xi_inf)*dxi_dx
    end select
  elseif(x.gt.xarr(npoi)) then
    vec=amat_arr(:,1,npoi)+(x-xarr(npoi))*dvec
!
    select case(ifuntype(iclass))
    case(13)
! left- rho_pol boundary, right - separatrix, 0<x<inf
      xi=tanh(x)
      xib=tanh(xarr(npoi))
      dxi_dx=1.d0/cosh(x)**2
      dxi_dxb=1.d0/cosh(xarr(npoi))**2
      dpphi_dxib=dvec(1)/dxi_dxb
      vec(1)=amat_arr(1,1,npoi)+dpphi_dxib*(xi-xib)
      dvec(1)=dpphi_dxib*dxi_dx
    case(14)
! left- rho_pol boundary, right - X-point, 0<x<inf
      xi=tanh(x)
      xib=tanh(xarr(npoi))
      dxi_dx=1.d0/cosh(x)**2
      dxi_dxb=1.d0/cosh(xarr(npoi))**2
      xi_inf=1.d0
      b=dvec(1)/(2.d0*(xib-xi_inf)*dxi_dxb)
      a=amat_arr(1,1,npoi)-b*(xib-xi_inf)**2
      vec(1)=a+b*(xi-xi_inf)**2
      dvec(1)=2.d0*b*(xi-xi_inf)*dxi_dx
    case(23)
! left- inner boundary, right - separatrix, 0<x<inf
      xi=tanh(x)**2
      xib=tanh(xarr(npoi))**2
      dxi_dx=2.d0*tanh(x)/cosh(x)**2
      dxi_dxb=2.d0*tanh(xarr(npoi))/cosh(xarr(npoi))**2
      dpphi_dxib=dvec(1)/dxi_dxb
      vec(1)=amat_arr(1,1,npoi)+dpphi_dxib*(xi-xib)
      dvec(1)=dpphi_dxib*dxi_dx
    case(24)
! left- inner boundary, right - X-point, 0<x<inf
      xi=tanh(x)**2
      xib=tanh(xarr(npoi))**2
      dxi_dx=2.d0*tanh(x)/cosh(x)**2
      dxi_dxb=2.d0*tanh(xarr(npoi))/cosh(xarr(npoi))**2
      xi_inf=1.d0
      b=dvec(1)/(2.d0*(xib-xi_inf)*dxi_dxb)
      a=amat_arr(1,1,npoi)-b*(xib-xi_inf)**2
      vec(1)=a+b*(xi-xi_inf)**2
      dvec(1)=2.d0*b*(xi-xi_inf)*dxi_dx
    case(33,43)
! left - separatrix or X-point, right -separatrix, -inf<x<inf
      xi=0.5d0*(1.d0+tanh(x))
      xib=0.5d0*(1.d0+tanh(xarr(npoi)))
      dxi_dx=0.5d0/cosh(x)**2
      dxi_dxb=0.5d0/cosh(xarr(npoi))**2
      dpphi_dxib=dvec(1)/dxi_dxb
      vec(1)=amat_arr(1,1,npoi)+dpphi_dxib*(xi-xib)
      dvec(1)=dpphi_dxib*dxi_dx
    case(34,44)
! left - separatrix or X-point, right -X-point, -inf<x<inf
      xi=0.5d0*(1.d0+tanh(x))
      xib=0.5d0*(1.d0+tanh(xarr(npoi)))
      dxi_dx=0.5d0/cosh(x)**2
      dxi_dxb=0.5d0/cosh(xarr(npoi))**2
      xi_inf=1.d0
      b=dvec(1)/(2.d0*(xib-xi_inf)*dxi_dxb)
      a=amat_arr(1,1,npoi)-b*(xib-xi_inf)**2
      vec(1)=a+b*(xi-xi_inf)**2
      dvec(1)=2.d0*b*(xi-xi_inf)*dxi_dx
    end select
  endif
!
! Store this (x,vec,dvec) for reuse by later modes at the same x.
  call interp_cache_store(x,vec,dvec)
!
  end subroutine interpolate_class_doublecount
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine interp_cache_store(x,vec,dvec)
!
  use sample_matrix_mod, only : n1
  use interp_cache_mod,  only : ncache,n1cache,xcache,veccache,dveccache
!
  implicit none
!
  double precision, intent(in) :: x
  double precision, dimension(n1), intent(in) :: vec,dvec
!
  integer :: ncap,newcap
  double precision, dimension(:),   allocatable :: xtmp
  double precision, dimension(:,:), allocatable :: vectmp,dvectmp
!
! Row count change forces a fresh buffer (and would mean a new grid anyway).
  if(allocated(xcache)) then
    if(n1cache.ne.n1) then
      deallocate(xcache,veccache,dveccache)
    endif
  endif
!
  if(.not.allocated(xcache)) then
    newcap=64
    allocate(xcache(newcap),veccache(n1,newcap),dveccache(n1,newcap))
    n1cache=n1
    ncache=0
  endif
!
  ncap=size(xcache)
  if(ncache.ge.ncap) then
    newcap=2*ncap
    allocate(xtmp(newcap),vectmp(n1,newcap),dvectmp(n1,newcap))
    xtmp(1:ncap)=xcache
    vectmp(:,1:ncap)=veccache
    dvectmp(:,1:ncap)=dveccache
    call move_alloc(xtmp,xcache)
    call move_alloc(vectmp,veccache)
    call move_alloc(dvectmp,dveccache)
  endif
!
  ncache=ncache+1
  xcache(ncache)=x
  veccache(1:n1,ncache)=vec
  dveccache(1:n1,ncache)=dvec
!
  end subroutine interp_cache_store
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine integrate_class_doublecount(iunit,resint)
!
! Computes integrals over x in Eq.(74) using adaptive polynomial interpolation
! of bounce integrals.
!
  use sample_matrix_mod, only : n1
  use orbit_dim_mod,     only : next,numbasef
  use get_matrix_mod,    only : relmargin,iclass
  use form_classes_doublecount_mod, only : ifuntype
  use cc_mod, only : dowrite
!
  implicit none
!
  integer :: iunit,ierr,i,npoiint,ib1,ie1,ib2,ie2
  double precision                :: x,xbeg,xend,hint,fmaxw
  double precision, dimension(n1) :: vec,dvec
  double precision, dimension(numbasef,2) :: resint
!
  npoiint=500
!
  ierr=0
!
  select case(ifuntype(iclass))
  case(11)
! two rho_pol boundaries, 0<x<1
    xbeg=0.d0
    xend=1.d0
  case(12)
! left- rho_pol boundary, right - inner boundary, 0<x<1
    xbeg=0.d0
    xend=1.d0
  case(13,14)
! left- rho_pol boundary, right - X-point, 0<x<inf
    xbeg=0.d0
    xend=-log(relmargin)
  case(21)
! left- inner boundary, right - rho_pol boundary, 0<x<1
    xbeg=0.d0
    xend=1.d0
  case(22)
! two inner boundaries, 0<x<1
    xbeg=0.d0
    xend=1.d0
  case(23,24)
! left- inner boundary, right - X-point, 0<x<inf
    xbeg=0.d0
    xend=-log(relmargin)
  case(31,41)
! left- X-point, right - rho_pol boundary, -inf<x<0
    xbeg=log(relmargin)
    xend=0.d0
  case(32,42)
! left- X-point, right - inner boundary, -inf<x<0
    xbeg=log(relmargin)
    xend=0.d0
  case(33,34,43,44)
! two X-points, -inf<x<inf
    xbeg=log(relmargin)
    xend=-xbeg
  end select
!
  ib1=4
  ie1=3+numbasef
  ib2=4+numbasef
  ie2=3+2*numbasef
!
  hint=(xend-xbeg)/npoiint
!
  x=xbeg
!
  call interpolate_class_doublecount(x,vec,dvec)
!
  call equilmaxw(vec(1),fmaxw)
!
  resint(:,1)=abs(dvec(1))*vec(ib1:ie1)*fmaxw
  resint(:,2)=abs(dvec(1))*vec(ib2:ie2)*fmaxw
!
  x=xend
!
  call interpolate_class_doublecount(x,vec,dvec)
!
  call equilmaxw(vec(1),fmaxw)
!
  resint(:,1)=resint(:,1)+abs(dvec(1))*vec(ib1:ie1)*fmaxw
  resint(:,2)=resint(:,2)+abs(dvec(1))*vec(ib2:ie2)*fmaxw
!
  resint=0.5d0*resint
!
  do i=1,npoiint-1
    x=xbeg+hint*dble(i)
!
    call interpolate_class_doublecount(x,vec,dvec)
!
    call equilmaxw(vec(1),fmaxw)
!
    resint(:,1)=resint(:,1)+abs(dvec(1))*vec(ib1:ie1)*fmaxw
    resint(:,2)=resint(:,2)+abs(dvec(1))*vec(ib2:ie2)*fmaxw
!
    if(dowrite) then
      write(iunit,*) x,vec
      write(iunit+1000,*) x,dvec
    endif
  enddo
!
  resint=hint*resint
!
  end subroutine integrate_class_doublecount
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine equilmaxw(psiast,fmaxw)
!
! Computes normalized equilibrium Maxwellian $f_0$ - factor in front of bounce
! time integral in Eq.(73)
!
  use global_invariants, only : toten
!
  implicit none
!
  double precision :: psiast,fmaxw,dens,temp,ddens,dtemp,phi_elec,dPhi_dpsi
!
!
  call denstemp_of_psi(psiast,dens,temp,ddens,dtemp)
  call phielec_of_psi(psiast,phi_elec,dPhi_dpsi)
!
!  fmaxw=dens*exp(phi_elec-toten)/temp**1.5d0
  fmaxw=dens*exp((phi_elec-toten)/temp)/temp**1.5d0
!
  end subroutine equilmaxw
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_minmax_bsc(sw_min_max,fun,xb,xe,f_minmax)
!
! Finds global minimum (sw_min_max=.true.) or global maximum (sw_min_max=.false.)
! of the function specified by "fun" in the inerval [xb,xe] using bisection
! Result - f_minmax
!
  implicit none
!
  integer, parameter :: nprimsrc=1000 !initial interval splitting for crude search
  double precision, parameter :: epserr=1.d-12
!
  logical :: sw_min_max
  integer :: i
  double precision :: xb,xe,f_minmax,hx,sigma,x,f,x1,x2,f1,f2,errdist
  integer, dimension(1) :: iextr
  double precision, dimension(0:nprimsrc) :: xarr,farr
  external :: fun
!
  errdist=abs(xe-xb)*epserr
!
  if(sw_min_max) then
    sigma=1.d0
  else
    sigma=-1.d0
  endif
!
  hx=(xe-xb)/dble(nprimsrc)
!
  do i=0,nprimsrc
    xarr(i)=xb+hx*dble(i)
!
    call fun(xarr(i),farr(i))
!
  enddo
!
  farr=sigma*farr
!
  iextr=minloc(farr)
!
  i=iextr(1)-1
!
  if(i.eq.0) then
    x1=xarr(0)
    x2=xarr(1)
    f1=farr(0)
    f2=farr(1)
  elseif(i.eq.nprimsrc) then
    x1=xarr(nprimsrc-1)
    x2=xarr(nprimsrc)
    f1=farr(nprimsrc-1)
    f2=farr(nprimsrc)
  else
    x1=xarr(i-1)
    x2=xarr(i+1)
    f1=farr(i-1)
    f2=farr(i+1)
  endif
!
  do
!
    x=0.5d0*(x2+x1)
!
    call fun(x,f)
!
    f=f*sigma
    if(f1.lt.f2) then
      x2=x
      f2=f
    else
      x1=x
      f1=f
    endif
!
    if((x2-x1).lt.errdist) exit
  enddo
!
  f_minmax=sigma*min(f1,f2)
!
  end subroutine find_minmax_bsc
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_Phiminmax(phi_elec_min,phi_elec_max)
!
! Finds minimum and maximum values of the normalized electrostatic potential
! on the Poincare cut. They are needed for setting the limits of integration
! over total energy
!
  use poicut_mod,        only : Rbou_hfs,Rbou_lfs
!
  implicit none
!
  double precision :: phi_elec_min,phi_elec_max
!
  call find_minmax_bsc(.true.,Phioncut,Rbou_hfs,Rbou_lfs,phi_elec_min)
  call find_minmax_bsc(.false.,Phioncut,Rbou_hfs,Rbou_lfs,phi_elec_max)
!
!------------
  contains
!------------
!
  subroutine Phioncut(R,phi_elec)
!
! Computes the normalized electrostatic potential as fnction of cut parameter R_c
!
  implicit none
!
  double precision :: R,phi_elec,Z,dZ_dR,bmod
  double precision,dimension(3) :: x
!
  call get_poicut(R,Z,dZ_dR)
!
  x(1)=R
  x(2)=0.d0
  x(3)=Z
!
  call get_bmod_and_Phi(x,bmod,phi_elec)
!
  end subroutine Phioncut
!
!------------
!
  end subroutine find_Phiminmax
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_Jperpmax(perpinv_max)
!
! Find the operational outer J_perp envelope on the Poincare cut.  The
! established bounded 1-D maximizer is used here; the certified topology
! contribution bound is checked separately by sample_matrix_out.  This
! routine must not turn a numerical envelope into an EQDSK interval proof.
!
  use potato_boundary_scan_mod, only : jperpmax_success, &
      jperpmax_unresolved,jperpmax_invalid_domain,jperpmax_status, &
      jperpmax_certified,jperpmax_witness,jperpmax_upper_bound
  use poicut_mod, only : npc,rpc_arr
  use global_invariants, only : toten
!
  implicit none
!
  integer :: i
  double precision :: perpinv_max,value,range,tolerance
  logical :: ok,jperp_eval_failed
!
  perpinv_max=0.d0
  jperpmax_status=jperpmax_unresolved
  jperpmax_certified=.false.
  jperpmax_witness=0.d0
  jperpmax_upper_bound=huge(1.d0)
  if(npc.lt.1) then
    jperpmax_status=jperpmax_invalid_domain
    return
  endif
  range=rpc_arr(npc)-rpc_arr(0)
  if(range.lt.0.d0) then
    jperpmax_status=jperpmax_invalid_domain
    return
  endif
  tolerance=128.d0*epsilon(1.d0)*max(1.d0,abs(rpc_arr(0)),abs(rpc_arr(npc)))
!
! Check all supplied cut nodes for a finite, strictly positive B.  This is a
! domain sanity check, not an interior positivity certificate.
  do i=0,npc
    call evaluate_jperp(rpc_arr(i),value,ok)
    if(.not.ok) then
      jperpmax_status=jperpmax_invalid_domain
      return
    endif
    jperpmax_witness=max(jperpmax_witness,value)
  enddo
  if(range.le.tolerance) then
    perpinv_max=jperpmax_witness
    jperpmax_upper_bound=jperpmax_witness
    jperpmax_status=jperpmax_success
    jperpmax_certified=.true.
    return
  endif
!
! The maximizer's 1000-point bracketing scan is independent of the cut-node
! witness above.  A callback failure remains a hard domain error.
  jperp_eval_failed=.false.
  call find_minmax_bsc(.false.,jperp_oncut,rpc_arr(0),rpc_arr(npc),perpinv_max)
  if(jperp_eval_failed .or. perpinv_max.ne.perpinv_max .or. &
     abs(perpinv_max).ge.huge(perpinv_max)) then
    jperpmax_status=jperpmax_invalid_domain
    return
  endif
  tolerance=1024.d0*epsilon(1.d0)*max(1.d0,abs(perpinv_max), &
                                      abs(jperpmax_witness))
  if(perpinv_max.lt.jperpmax_witness-tolerance) then
    jperpmax_status=jperpmax_unresolved
    return
  endif
  perpinv_max=max(perpinv_max,jperpmax_witness)
  jperpmax_upper_bound=perpinv_max
  jperpmax_status=jperpmax_success
  jperpmax_certified=.true.
!
!------------
  contains
!------------
!
  subroutine evaluate_jperp(R,perpinv,ok_value)
!
  use global_invariants, only : toten
  use potato_symbolic_kernel_mod, only : potato_jperp_kernel
!
  implicit none
!
  double precision, intent(in) :: R
  double precision, intent(out) :: perpinv
  logical, intent(out) :: ok_value
  double precision :: phi_elec,Z,dZ_dR,bmod,sqrtg
  double precision :: qPhi_prime,magnetic_field_B_prime
  double precision :: candidate,envelope_bound,derivative
  double precision,dimension(3) :: x,bder,hcovar,hctrvr,hcurl,derphi
!
  call get_poicut(R,Z,dZ_dR)
  x(1)=R
  x(2)=0.d0
  x(3)=Z
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  call elefie(x,phi_elec,derphi)
  ok_value=(bmod.gt.0.d0 .and. bmod.eq.bmod .and. abs(bmod).lt.huge(bmod))
  if(.not.ok_value) then
    perpinv=0.d0
    return
  endif
  qPhi_prime=derphi(1)+derphi(3)*dZ_dR
  magnetic_field_B_prime=bmod*(bder(1)+bder(3)*dZ_dR)
  call potato_jperp_kernel(toten,phi_elec,bmod,qPhi_prime, &
                           magnetic_field_B_prime, &
                           candidate,envelope_bound,derivative)
! Only the generated positive part is used here, as an upper-bound helper for
! max_D[(H-qPhi)/B]_+.  No physical class is built from this clamped value.
  ok_value=(envelope_bound.eq.envelope_bound .and. &
            abs(envelope_bound).lt.huge(envelope_bound))
  if(ok_value) then
    perpinv=envelope_bound
  else
    perpinv=0.d0
  endif
  end subroutine evaluate_jperp
!
  subroutine jperp_oncut(R,perpinv)
    double precision, intent(in) :: R
    double precision, intent(out) :: perpinv
    logical :: ok_value

    call evaluate_jperp(R,perpinv,ok_value)
    if(.not.ok_value) then
      jperp_eval_failed=.true.
      perpinv=0.d0
    endif
  end subroutine jperp_oncut
!
!------------
!
  end subroutine find_Jperpmax

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_jperp_topology_boundaries(candidates,nmax,ncandidates,ierr)
!
! Return the certified outer-J_perp topology candidates for one total energy.
! On the Poincare cut
!
!   v_parallel^2 = H - Phi(R) - J_perp B(R)
!                 = B(R) [ q(R;H) - J_perp ],
!   q(R;H) = [H-Phi(R)]/B(R).
!
! The allowed-region root set can change only at a cut endpoint or at a
! stationary root of q.  The class/fixpoint root set has two additional
! one-dimensional discriminants: a fixed point can be born or die when the
! fixed-point branch J_perp(R_c) is stationary, and a fixed point can enter or
! leave the rho_pol-clipped cut at a clipping endpoint.  All of these roots
! are found from the equations used by find_bounds_fixpoints; the outer
! sampler must not infer them from finite midpoint signatures.
!
  use find_all_roots_mod, only : nroots,roots,nsearch_min,relerr_allroots, &
                                 root_success,root_unresolved_separation
  use global_invariants, only : toten
  use poicut_mod, only : npc,rpc_arr,Rbou_hfs,Rbou_lfs
  use potato_boundary_scan_mod, only : fixedpoint_scan_sigma, &
                                       fixedpoint_scan_branch, &
                                       fixedpoint_scan_left, fixedpoint_scan_right
  implicit none
!
  integer, intent(in) :: nmax
  integer, intent(out) :: ncandidates,ierr
  double precision, intent(out) :: candidates(nmax)
  integer, parameter :: max_scan_roots=4096
  integer :: ierr_roots,nsearch_save,i,j,k,nq,ndisc,npart,nsigma,branch
  integer :: npart_initial,nfp_segments,nbd_segments,ierr_local
  integer :: collision_segment_left,collision_segment_right,collision_boundary_segment
  integer :: collision_boundary_left,collision_boundary_right
  integer :: boundary_stage
  double precision :: relerr_save,q,range,midpoint,jvalue,pstar,refined_root
  double precision :: collision_jlo,collision_jhi
  double precision :: boundary_stage_r,boundary_stage_j
  double precision, allocatable :: qroots(:),disc_roots(:),rpartition(:)
  integer, allocatable :: fp_sigma(:),fp_branch(:)
  double precision, allocatable :: fp_rlo(:),fp_rhi(:),fp_jlo(:),fp_jhi(:), &
                                   fp_plo(:),fp_phi(:)
  integer, allocatable :: bd_type(:),bd_sigma(:)
  double precision, allocatable :: bd_rlo(:),bd_rhi(:),bd_jlo(:),bd_jhi(:)
  logical :: ok
  external :: find_all_roots_bracketed,fixedpoint_discriminant, &
              fixedpoint_branch_stationary,fixedpoint_roots_at_R, &
              fixedpoint_branch_value,fixedpoint_turning_intersection
!
  ncandidates=0
  ierr=root_success
  candidates=0.d0
  boundary_stage=0
  boundary_stage_r=0.d0
  boundary_stage_j=0.d0
  collision_segment_left=0
  collision_segment_right=0
  collision_boundary_segment=0
  collision_boundary_left=0
  collision_boundary_right=0
  nfp_segments=0
  nbd_segments=0
  npart=0
  if(nmax.le.0 .or. npc.le.0) then
    ierr=3
    return
  endif
  range=rpc_arr(npc)-rpc_arr(0)
  if(range.le.0.d0) then
    ierr=3
    return
  endif
!
! Endpoint roots are boundary candidates even when they are not stationary.
  call jperp_value(rpc_arr(0),q,ok)
  if(.not.ok) then
    ierr=2
    return
  endif
  call add_candidate(q)
  if(ierr.ne.0) return
  call jperp_value(rpc_arr(npc),q,ok)
  if(.not.ok) then
    ierr=2
    return
  endif
  call add_candidate(q)
  if(ierr.ne.0) return
  if(Rbou_hfs.gt.rpc_arr(0) .and. Rbou_hfs.lt.rpc_arr(npc)) then
    call jperp_value(Rbou_hfs,q,ok)
    if(.not.ok) then
      ierr=2
      return
    endif
    call add_candidate(q)
    if(ierr.ne.0) return
  endif
  if(Rbou_lfs.gt.rpc_arr(0) .and. Rbou_lfs.lt.rpc_arr(npc)) then
    call jperp_value(Rbou_lfs,q,ok)
    if(.not.ok) then
      ierr=2
      return
    endif
    call add_candidate(q)
    if(ierr.ne.0) return
  endif
!
! Find every stationary q root on a deliberately fine bounded scan.  The
! callback has a two-sided finite-difference derivative only in the open cut;
! endpoint values were already added explicitly above.
  nsearch_save=nsearch_min
  relerr_save=relerr_allroots
  nsearch_min=max(nsearch_min,4096)
  relerr_allroots=1.d-11
  boundary_stage=10
  call find_all_roots_bracketed(jperp_stationary,rpc_arr(0),rpc_arr(npc),ierr_roots)
  if(ierr_roots.ne.root_success) then
    call fail_boundary(ierr_roots)
    return
  endif
  nq=nroots
  allocate(qroots(max(1,nq)))
  if(nq.gt.0) qroots(1:nq)=roots(1:nq)
!
  do i=1,nq
    call jperp_value(qroots(i),q,ok)
    if(.not.ok) then
      call fail_boundary(2)
      return
    endif
    call add_candidate(q)
    if(ierr.ne.0) then
      call fail_boundary(ierr)
      return
    endif
  enddo

! A fixed-point branch can end on the turning boundary u=0.  Substituting
! u=0 into the exact fixed-point equation a*u**2+b*u+c=0 gives c(R)=0.
! These intersections are not fixed-point discriminants and need their own
! R roots; their J_perp values are topology boundaries of the outer scan.
  boundary_stage=15
  call find_all_roots_bracketed(fixedpoint_turning_intersection, &
                                rpc_arr(0),rpc_arr(npc),ierr_roots)
  if(ierr_roots.ne.root_success) then
    call fail_boundary(ierr_roots)
    return
  endif
  do i=1,nroots
    call jperp_value(roots(i),q,ok)
    if(.not.ok) then
      call fail_boundary(2)
      return
    endif
    call add_candidate(q)
    if(ierr.ne.0) then
      call fail_boundary(ierr)
      return
    endif
    call add_rpartition(roots(i))
    if(ierr.ne.0) then
      call fail_boundary(ierr)
      return
    endif
  enddo
!
! The discriminant of the exact fixed-point equation is a second physical
! boundary set.  With u=p_parallel (signed), psi*=psi+rho0*h_phi*u and
! u^2=A-J_perp*B, the stationary condition d(psi*)/dR_c=0 is
!
!   a(R) u^2 + b(R) u + c(R) = 0 .
!
! Its discriminant partitions the cut into intervals on which each fixed-point
! branch is a single-valued function J_perp(R_c).  We then root-find
! dJ_perp/dR_c on every such interval, rather than sampling J_perp midpoints.
  boundary_stage=20
  call find_all_roots_bracketed(fixedpoint_discriminant, &
                                rpc_arr(0),rpc_arr(npc),ierr_roots)
  if(ierr_roots.ne.root_success) then
    call fail_boundary(ierr_roots)
    return
  endif
  ndisc=nroots
  allocate(disc_roots(max(1,ndisc)))
  if(ndisc.gt.0) disc_roots(1:ndisc)=roots(1:ndisc)
  do i=1,ndisc
    call add_fixedpoint_candidates(disc_roots(i))
    if(ierr.ne.0) then
      call fail_boundary(ierr)
      return
    endif
  enddo
  boundary_stage=25
  call add_fixedpoint_candidates(rpc_arr(0))
  if(ierr.ne.0) then
    call fail_boundary(ierr)
    return
  endif
  call add_fixedpoint_candidates(rpc_arr(npc))
  if(ierr.ne.0) then
    call fail_boundary(ierr)
    return
  endif
!
! Split the R_c scan at every allowed-region and fixed-point branch endpoint.
  allocate(rpartition(max_scan_roots))
  npart=0
  call add_rpartition(rpc_arr(0))
  call add_rpartition(rpc_arr(npc))
  do i=1,nq
    call add_rpartition(qroots(i))
  enddo
  do i=1,ndisc
    call add_rpartition(disc_roots(i))
  enddo
  if(Rbou_hfs.gt.rpc_arr(0) .and. Rbou_hfs.lt.rpc_arr(npc)) &
      call add_rpartition(Rbou_hfs)
  if(Rbou_lfs.gt.rpc_arr(0) .and. Rbou_lfs.lt.rpc_arr(npc)) &
      call add_rpartition(Rbou_lfs)
  if(ierr.ne.0) then
    call fail_boundary(ierr)
    return
  endif
  call sort_rpartition
!
! Add the projection discriminants of both signed fixed-point branches.
  npart_initial=npart
  do nsigma=-1,1,2
    fixedpoint_scan_sigma=nsigma
    do branch=1,2
      fixedpoint_scan_branch=branch
      do i=1,npart_initial-1
        boundary_stage=30
        boundary_stage_r=0.5d0*(rpartition(i)+rpartition(i+1))
        midpoint=0.5d0*(rpartition(i)+rpartition(i+1))
        call fixedpoint_branch_value(midpoint,nsigma,branch,jvalue,pstar,ok)
        if(.not.ok) cycle
        fixedpoint_scan_left=rpartition(i)
        fixedpoint_scan_right=rpartition(i+1)
        call find_all_roots_bracketed(fixedpoint_branch_stationary, &
                                      rpartition(i),rpartition(i+1),ierr_local)
        if(ierr_local.ne.root_success) then
          call fail_boundary(ierr_local)
          return
        endif
        do j=1,nroots
          call fixedpoint_branch_value(roots(j),nsigma,branch,jvalue,pstar,ok)
          if(.not.ok) then
            call fail_boundary(2)
            return
          endif
          call add_candidate(jvalue)
          if(ierr.ne.0) then
            call fail_boundary(ierr)
            return
          endif
          call add_rpartition(roots(j))
          if(ierr.ne.0) then
            call fail_boundary(ierr)
            return
          endif
        enddo
      enddo
    enddo
  enddo
  call sort_rpartition
! The separatrix root set also changes when two distinct fixed points have the
! same (J_perp,psi^*) value.  This critical-value collision is not a
! fixed-point birth/death: it is where a moving cut curve becomes tangent to
! another fixed-point level.  The branch projection is monotone on every
! R interval below, so each collision is reduced to a bounded root problem in
! J_perp; no midpoint topology inference is used.
  allocate(fp_sigma(max_scan_roots),fp_branch(max_scan_roots), &
      fp_rlo(max_scan_roots),fp_rhi(max_scan_roots), &
      fp_jlo(max_scan_roots),fp_jhi(max_scan_roots), &
      fp_plo(max_scan_roots),fp_phi(max_scan_roots))
  boundary_stage=40
  call collect_fixedpoint_segments
  if(ierr.ne.0) then
    call fail_boundary(ierr)
    return
  endif
  do i=1,nfp_segments-1
    do j=i+1,nfp_segments
      boundary_stage=50
      boundary_stage_r=0.5d0*(fp_rlo(i)+fp_rhi(i))
      collision_segment_left=i
      collision_segment_right=j
      collision_boundary_segment=0
      collision_boundary_left=0
      collision_boundary_right=0
      collision_jlo=max(0.d0,max(min(fp_jlo(i),fp_jhi(i)), &
                                 min(fp_jlo(j),fp_jhi(j))))
      collision_jhi=min(max(fp_jlo(i),fp_jhi(i)), &
                        max(fp_jlo(j),fp_jhi(j)))
      if(collision_jhi.lt.collision_jlo) cycle
      if(collision_jhi-collision_jlo.le.256.d0*epsilon(1.d0)* &
         max(1.d0,abs(collision_jlo),abs(collision_jhi))) then
        midpoint=0.5d0*(collision_jlo+collision_jhi)
        call fixedpoint_collision_value(midpoint,jvalue,ok)
        if(.not.ok) then
          call fail_boundary(2)
          return
        endif
        if(abs(jvalue).le.1.d-10*max(1.d0,abs(jvalue))) then
          call add_candidate(midpoint)
          if(ierr.ne.0) then
            call fail_boundary(ierr)
            return
          endif
        endif
      else
        call find_all_roots_bracketed(fixedpoint_collision, &
                                      collision_jlo,collision_jhi,ierr_local)
        if(ierr_local.ne.root_success) then
          call fail_boundary(ierr_local)
          return
        endif
        do k=1,nroots
          call refine_collision_candidate(roots(k),refined_root,ok)
          if(ok) then
            call add_candidate(refined_root)
          else
            call add_candidate(roots(k))
          endif
          if(ierr.ne.0) then
            call fail_boundary(ierr)
            return
          endif
        enddo
      endif
  enddo
  enddo
  allocate(bd_type(max_scan_roots),bd_sigma(max_scan_roots), &
      bd_rlo(max_scan_roots),bd_rhi(max_scan_roots), &
      bd_jlo(max_scan_roots),bd_jhi(max_scan_roots))
  boundary_stage=60
  call collect_boundary_segments
  if(ierr.ne.0) then
    call fail_boundary(ierr)
    return
  endif
! Two turning-boundary branches can meet at the same canonical-momentum
! level.  This is the limiting case of a separatrix entering or leaving the
! clipped cut; it is a boundary-boundary collision, not a fixed-point root.
  do i=1,nbd_segments-1
    do j=i+1,nbd_segments
      boundary_stage=65
      collision_segment_left=0
      collision_segment_right=0
      collision_boundary_segment=0
      collision_boundary_left=i
      collision_boundary_right=j
      collision_jlo=max(0.d0,max(min(bd_jlo(i),bd_jhi(i)), &
                                 min(bd_jlo(j),bd_jhi(j))))
      collision_jhi=min(max(bd_jlo(i),bd_jhi(i)), &
                        max(bd_jlo(j),bd_jhi(j)))
      if(collision_jhi.lt.collision_jlo) cycle
      if(collision_jhi-collision_jlo.le.256.d0*epsilon(1.d0)* &
         max(1.d0,abs(collision_jlo),abs(collision_jhi))) then
        midpoint=0.5d0*(collision_jlo+collision_jhi)
        call fixedpoint_collision_value(midpoint,jvalue,ok)
        if(.not.ok) then
          call fail_boundary(2)
          return
        endif
        if(abs(jvalue).le.1.d-10*max(1.d0,abs(jvalue))) then
          call add_candidate(midpoint)
          if(ierr.ne.0) then
            call fail_boundary(ierr)
            return
          endif
        endif
      else
        call find_all_roots_bracketed(fixedpoint_collision, &
                                      collision_jlo,collision_jhi,ierr_local)
        if(ierr_local.ne.root_success) then
          call fail_boundary(ierr_local)
          return
        endif
        do k=1,nroots
          call refine_collision_candidate(roots(k),refined_root,ok)
          if(ok) then
            call add_candidate(refined_root)
          else
            call add_candidate(roots(k))
          endif
          if(ierr.ne.0) then
            call fail_boundary(ierr)
            return
          endif
        enddo
      endif
    enddo
  enddo

  do i=1,nfp_segments
    do j=1,nbd_segments
      boundary_stage=70
      boundary_stage_r=0.5d0*(fp_rlo(i)+fp_rhi(i))
      collision_segment_left=i
      collision_segment_right=0
      collision_boundary_segment=j
      collision_boundary_left=0
      collision_boundary_right=j
      collision_jlo=max(0.d0,max(min(fp_jlo(i),fp_jhi(i)), &
                                 min(bd_jlo(j),bd_jhi(j))))
      collision_jhi=min(max(fp_jlo(i),fp_jhi(i)), &
                        max(bd_jlo(j),bd_jhi(j)))
      if(collision_jhi.lt.collision_jlo) cycle
      if(collision_jhi-collision_jlo.le.256.d0*epsilon(1.d0)* &
         max(1.d0,abs(collision_jlo),abs(collision_jhi))) then
        midpoint=0.5d0*(collision_jlo+collision_jhi)
        call fixedpoint_collision_value(midpoint,jvalue,ok)
        if(.not.ok) then
          call fail_boundary(2)
          return
        endif
        if(abs(jvalue).le.1.d-10*max(1.d0,abs(jvalue))) then
          call add_candidate(midpoint)
          if(ierr.ne.0) then
            call fail_boundary(ierr)
            return
          endif
        endif
      else
        call find_all_roots_bracketed(fixedpoint_collision, &
                                      collision_jlo,collision_jhi,ierr_local)
        if(ierr_local.ne.root_success) then
          call fail_boundary(ierr_local)
          return
        endif
        do k=1,nroots
          call refine_collision_candidate(roots(k),refined_root,ok)
          if(ok) then
            call add_candidate(refined_root)
          else
            call add_candidate(roots(k))
          endif
          if(ierr.ne.0) then
            call fail_boundary(ierr)
            return
          endif
        enddo
      endif
    enddo
  enddo
!
! A fixed point crossing the rho_pol clipping boundary changes the class set
! even when v_parallel^2 remains positive there.  The quadratic above gives
! those candidates directly at each clipping endpoint.
  if(Rbou_hfs.gt.rpc_arr(0) .and. Rbou_hfs.lt.rpc_arr(npc)) then
    boundary_stage=80
    call add_fixedpoint_candidates(Rbou_hfs)
    if(ierr.ne.0) then
      call fail_boundary(ierr)
      return
    endif
  endif
  if(Rbou_lfs.gt.rpc_arr(0) .and. Rbou_lfs.lt.rpc_arr(npc)) then
    boundary_stage=81
    call add_fixedpoint_candidates(Rbou_lfs)
    if(ierr.ne.0) then
      call fail_boundary(ierr)
      return
    endif
  endif
!
  call restore_scan_state
  deallocate(qroots,disc_roots,rpartition,fp_sigma,fp_branch,fp_rlo,fp_rhi, &
      fp_jlo,fp_jhi,fp_plo,fp_phi,bd_type,bd_sigma,bd_rlo,bd_rhi,bd_jlo,bd_jhi)
  return
!
contains

  subroutine fail_boundary(status)
    integer, intent(in) :: status

    print *,'find_jperp_topology_boundaries: failure stage,H,R,J,npart,nfp,nbd,ierr = ', &
        boundary_stage,toten,boundary_stage_r,boundary_stage_j,npart, &
        nfp_segments,nbd_segments,status
    ierr=status
    call restore_scan_state
    if(allocated(qroots)) deallocate(qroots)
    if(allocated(disc_roots)) deallocate(disc_roots)
    if(allocated(rpartition)) deallocate(rpartition)
    if(allocated(fp_sigma)) deallocate(fp_sigma,fp_branch,fp_rlo,fp_rhi, &
        fp_jlo,fp_jhi,fp_plo,fp_phi)
    if(allocated(bd_type)) deallocate(bd_type,bd_sigma,bd_rlo,bd_rhi,bd_jlo,bd_jhi)
  end subroutine fail_boundary

  subroutine restore_scan_state
    nsearch_min=nsearch_save
    relerr_allroots=relerr_save
    fixedpoint_scan_sigma=1
    fixedpoint_scan_branch=1
    fixedpoint_scan_left=0.d0
    fixedpoint_scan_right=0.d0
    collision_segment_left=0
    collision_segment_right=0
    collision_boundary_segment=0
    collision_boundary_left=0
    collision_boundary_right=0
  end subroutine restore_scan_state

  subroutine add_candidate(value)
    double precision, intent(in) :: value
    integer :: local_i
    double precision :: candidate_tolerance

    if(value.ne.value .or. abs(value).gt.huge(value)) then
      ierr=2
      return
    endif
    ! The outer integration variable is J_perp >= 0.  Fixed-point and
    ! boundary level equations are evaluated on a wider algebraic domain;
    ! roots below zero are valid equation roots but cannot change the
    ! physical topology certificate.
    if(value.lt.0.d0) return
    do local_i=1,ncandidates
      candidate_tolerance=max(32.d0*epsilon(1.d0)* &
          max(1.d-300,abs(value),abs(candidates(local_i))), &
          relerr_allroots*max(abs(value),abs(candidates(local_i))))
      if(abs(value-candidates(local_i)).le.candidate_tolerance) return
    enddo
    if(ncandidates.ge.nmax) then
      ierr=3
      return
    endif
    ncandidates=ncandidates+1
    candidates(ncandidates)=value
  end subroutine add_candidate

  subroutine add_fixedpoint_candidates(R)
    double precision, intent(in) :: R
    integer :: local_sigma,local_branch,nfixed
    double precision :: local_u(2),local_j(2),local_p(2)
    logical :: local_ok

    do local_sigma=-1,1,2
      call fixedpoint_roots_at_R(R,local_sigma,local_u,local_j,local_p, &
                                 nfixed,local_ok)
      if(.not.local_ok) then
        ierr=2
        return
      endif
      do local_branch=1,nfixed
        call add_candidate(local_j(local_branch))
        if(ierr.ne.0) return
      enddo
    enddo
  end subroutine add_fixedpoint_candidates

  subroutine add_rpartition(value)
    double precision, intent(in) :: value
    integer :: local_i
    double precision :: tolerance,scale

    if(value.ne.value .or. abs(value).gt.huge(value)) then
      ierr=2
      return
    endif
    if(value.lt.rpc_arr(0) .or. value.gt.rpc_arr(npc)) return
    scale=max(1.d0,abs(value))
    tolerance=128.d0*epsilon(1.d0)*scale
    do local_i=1,npart
      if(rpartition(local_i).eq.value) return
      if(abs(rpartition(local_i)-value).le.tolerance) then
        ! Distinct boundary equations that collapse at floating-point scale
        ! cannot be certified as one transition.  Merging them would omit a
        ! possible narrow topology interval.
        ierr=root_unresolved_separation
        return
      endif
    enddo
    if(npart.ge.max_scan_roots) then
      ierr=3
      return
    endif
    npart=npart+1
    rpartition(npart)=value
  end subroutine add_rpartition

  subroutine sort_rpartition
    integer :: local_i,local_j
    double precision :: local_key

    do local_i=2,npart
      local_key=rpartition(local_i)
      local_j=local_i-1
      do while(local_j.ge.1 .and. rpartition(local_j).gt.local_key)
        rpartition(local_j+1)=rpartition(local_j)
        local_j=local_j-1
      enddo
      rpartition(local_j+1)=local_key
    enddo
  end subroutine sort_rpartition

  subroutine collect_fixedpoint_segments
! On each final R interval the quadratic branch is single-valued and its
! projection J_perp(R) has no stationary root.  Store its bounded endpoint
! values for the pairwise critical-value collision scan.
    integer :: local_sigma,local_branch,local_i
    double precision :: local_mid,local_j,local_p,local_jl,local_jr, &
                        local_pl,local_pr,local_jq1,local_jq2, &
                        local_pq1,local_pq2,local_rlo,local_rhi,tolerance,width
    logical :: local_ok,local_okl,local_okr,local_okq1,local_okq2

    nfp_segments=0
    do local_sigma=-1,1,2
      do local_branch=1,2
        do local_i=1,npart-1
          width=rpartition(local_i+1)-rpartition(local_i)
          if(width.le.0.d0) then
            ierr=3
            return
          endif
          local_mid=0.5d0*(rpartition(local_i)+rpartition(local_i+1))
          boundary_stage_r=local_mid
          call fixedpoint_branch_value(local_mid,local_sigma,local_branch, &
                                       local_j,local_p,local_ok)
          if(.not.local_ok) cycle
          boundary_stage=41
          boundary_stage_j=local_j
          call fixedpoint_branch_endpoint(rpartition(local_i), &
              rpartition(local_i),rpartition(local_i+1),1.d0,local_sigma, &
              local_branch,local_jl,local_pl,local_rlo,local_okl)
          call fixedpoint_branch_endpoint(rpartition(local_i+1), &
              rpartition(local_i),rpartition(local_i+1),-1.d0,local_sigma, &
              local_branch,local_jr,local_pr,local_rhi,local_okr)
          if(.not.local_okl .or. .not.local_okr) then
            boundary_stage=42
            ierr=2
            return
          endif
          if(local_rhi-local_rlo.le.256.d0*epsilon(1.d0)* &
             max(1.d0,abs(local_rlo),abs(local_rhi))) then
            ! The endpoint contraction can leave a formally positive branch
            ! interval below the representable R resolution.  It has no open
            ! phase-space interval on which a collision can occur; retaining
            ! it would make its inverse callback fail at every midpoint.
            cycle
          endif
          call fixedpoint_branch_value(rpartition(local_i)+0.25d0*width, &
              local_sigma,local_branch,local_jq1,local_pq1,local_okq1)
          call fixedpoint_branch_value(rpartition(local_i)+0.75d0*width, &
              local_sigma,local_branch,local_jq2,local_pq2,local_okq2)
          if(.not.local_okq1 .or. .not.local_okq2) then
            boundary_stage=43
            ierr=2
            return
          endif
! Branches with no non-negative physical J_perp value do not belong to
! the outer integration domain.  Discard them before applying the strict
! monotonicity check; otherwise an entirely negative branch can be rejected
! for a harmless projection shape outside the certified physical domain.
          if(max(local_jl,local_jr,local_jq1,local_jq2).le.0.d0) cycle
          tolerance=1.d-10*max(1.d0,abs(local_jl),abs(local_jr), &
                                abs(local_jq1),abs(local_jq2))
          if(local_jq1.lt.min(local_jl,local_jr)-tolerance .or. &
             local_jq1.gt.max(local_jl,local_jr)+tolerance .or. &
             local_jq2.lt.min(local_jl,local_jr)-tolerance .or. &
             local_jq2.gt.max(local_jl,local_jr)+tolerance) then
            boundary_stage=44
            ierr=2
            return
          endif
          if(abs(local_jr-local_jl).le.256.d0*epsilon(1.d0)* &
             max(1.d0,abs(local_jl),abs(local_jr))) then
            ! A branch whose projection has zero resolvable J_perp width
            ! contributes no open interval to the outer partition.  It is a
            ! degenerate endpoint, not an unresolved topology transition.
            cycle
          endif
          if(nfp_segments.ge.max_scan_roots) then
            ierr=3
            return
          endif
          nfp_segments=nfp_segments+1
          fp_sigma(nfp_segments)=local_sigma
          fp_branch(nfp_segments)=local_branch
! The stored R interval is the certified open branch interval, not the
! physical endpoint at which a quadratic root may coalesce.  The exact
! physical endpoint remains in rpartition/candidate data above.  Using the
! contracted positions here prevents the inverse collision map from probing
! through an invalid endpoint during its bounded bisection.
          fp_rlo(nfp_segments)=local_rlo
          fp_rhi(nfp_segments)=local_rhi
          fp_jlo(nfp_segments)=local_jl
          fp_jhi(nfp_segments)=local_jr
          fp_plo(nfp_segments)=local_pl
          fp_phi(nfp_segments)=local_pr
        enddo
      enddo
    enddo
  end subroutine collect_fixedpoint_segments

  subroutine collect_boundary_segments
! Separatrix intersections can also be born at a physical boundary when an
! X-point level reaches that boundary value.  Include both kinds of boundary
! level curves used by find_bounds_fixpoints: v_parallel=0 turning roots and
! fixed cut/rho_pol endpoints.  The q(R) projection is monotone on the same
! final R partition, so its inverse is bounded just like a fixed-point branch.
    integer :: local_i
    double precision :: local_width,local_q0,local_q1,local_q25,local_q75, &
                        local_p,local_tolerance
    logical :: local_ok0,local_ok1,local_ok25,local_ok75

    nbd_segments=0
    do local_i=1,npart-1
      local_width=rpartition(local_i+1)-rpartition(local_i)
      if(local_width.le.0.d0) then
        ierr=3
        return
      endif
      call jperp_value(rpartition(local_i),local_q0,local_ok0)
      call jperp_value(rpartition(local_i+1),local_q1,local_ok1)
      call jperp_value(rpartition(local_i)+0.25d0*local_width, &
                       local_q25,local_ok25)
      call jperp_value(rpartition(local_i)+0.75d0*local_width, &
                       local_q75,local_ok75)
      if(.not.local_ok0 .or. .not.local_ok1 .or. .not.local_ok25 .or. &
         .not.local_ok75) then
        ierr=2
        return
      endif
      local_tolerance=1.d-10*max(1.d0,abs(local_q0),abs(local_q1), &
                                  abs(local_q25),abs(local_q75))
      if(local_q25.lt.min(local_q0,local_q1)-local_tolerance .or. &
         local_q25.gt.max(local_q0,local_q1)+local_tolerance .or. &
         local_q75.lt.min(local_q0,local_q1)-local_tolerance .or. &
         local_q75.gt.max(local_q0,local_q1)+local_tolerance) then
        ierr=2
        return
      endif
      if(max(local_q0,local_q1).gt.0.d0) then
        call add_boundary_segment(1,0,rpartition(local_i), &
            rpartition(local_i+1),local_q0,local_q1)
        if(ierr.ne.0) return
      endif
    enddo
    call add_fixed_boundary_segment(rpc_arr(0))
    call add_fixed_boundary_segment(rpc_arr(npc))
    if(Rbou_hfs.gt.rpc_arr(0) .and. Rbou_hfs.lt.rpc_arr(npc)) &
        call add_fixed_boundary_segment(Rbou_hfs)
    if(Rbou_lfs.gt.rpc_arr(0) .and. Rbou_lfs.lt.rpc_arr(npc)) &
        call add_fixed_boundary_segment(Rbou_lfs)
  end subroutine collect_boundary_segments

  subroutine add_boundary_segment(local_type,local_sigma,local_rlo,local_rhi, &
                                  local_jlo,local_jhi)
    integer, intent(in) :: local_type,local_sigma
    double precision, intent(in) :: local_rlo,local_rhi,local_jlo,local_jhi

    if(nbd_segments.ge.max_scan_roots) then
      ierr=3
      return
    endif
    nbd_segments=nbd_segments+1
    bd_type(nbd_segments)=local_type
    bd_sigma(nbd_segments)=local_sigma
    bd_rlo(nbd_segments)=local_rlo
    bd_rhi(nbd_segments)=local_rhi
    bd_jlo(nbd_segments)=local_jlo
    bd_jhi(nbd_segments)=local_jhi
  end subroutine add_boundary_segment

  subroutine add_fixed_boundary_segment(local_r)
    double precision, intent(in) :: local_r
    double precision :: local_j
    integer :: local_sigma
    logical :: local_ok

    call jperp_value(local_r,local_j,local_ok)
    if(.not.local_ok) then
      ierr=2
      return
    endif
    if(local_j.le.0.d0) return
    do local_sigma=-1,1,2
      call add_boundary_segment(2,local_sigma,local_r,local_r,0.d0,local_j)
      if(ierr.ne.0) return
    enddo
  end subroutine add_fixed_boundary_segment

  subroutine fixedpoint_branch_endpoint(endpoint,rlo,rhi,direction, &
                                         local_sigma,local_branch, &
                                         local_j,local_p,local_r,local_ok)
! A branch can coalesce at a discriminant endpoint, where the branch index
! deliberately has no exact two-root representation.  Contract only this
! certified endpoint inward; failure to find the branch is a certificate
! failure, not permission to drop a collision interval.
    double precision, intent(in) :: endpoint,rlo,rhi,direction
    integer, intent(in) :: local_sigma,local_branch
    double precision, intent(out) :: local_j,local_p,local_r
    logical, intent(out) :: local_ok
    integer :: local_k
    double precision :: distance,max_distance,probe,scale

    call fixedpoint_branch_value(endpoint,local_sigma,local_branch, &
                                 local_j,local_p,local_ok)
    local_r=endpoint
    if(local_ok) return
    local_ok=.false.
    local_j=0.d0
    local_p=0.d0
    scale=max(1.d0,abs(endpoint),abs(rlo),abs(rhi))
    max_distance=0.49d0*(rhi-rlo)
    distance=max(1.d-8*(rhi-rlo),256.d0*epsilon(1.d0)*scale)
    do local_k=1,32
      if(distance.gt.max_distance) distance=max_distance
      probe=endpoint+direction*distance
      if(probe.le.rlo .or. probe.ge.rhi) exit
      call fixedpoint_branch_value(probe,local_sigma,local_branch, &
                                   local_j,local_p,local_ok)
      if(local_ok) then
        local_r=probe
        return
      endif
      if(distance.ge.max_distance) exit
      distance=min(2.d0*distance,max_distance)
    enddo
  end subroutine fixedpoint_branch_endpoint

  subroutine fixedpoint_collision(jtest,difference,ddifference_dj)
! Difference of two critical values after bounded inversion of their monotone
! J_perp(R) projections.  The derivative is only a local diagnostic for the
! all-roots scanner; all function evaluations remain inside certified branch
! intervals.
    use find_all_roots_mod, only : root_eval_valid,root_eval_error, &
                                   root_invalid_domain
    double precision, intent(in) :: jtest
    double precision, intent(out) :: difference,ddifference_dj
    double precision :: h,dm,dp,scale
    logical :: ok0,okm,okp

    root_eval_valid=.true.
    root_eval_error=0
    call fixedpoint_collision_value(jtest,difference,ok0)
    boundary_stage_j=jtest
    if(.not.ok0) then
      root_eval_valid=.false.
      root_eval_error=root_invalid_domain
      difference=0.d0
      ddifference_dj=0.d0
      return
    endif
    scale=max(1.d0,abs(jtest),abs(collision_jlo),abs(collision_jhi))
! The overlap can be much narrower than the absolute J scale.  A nominal
! relative step then collapses to the roundoff floor even though a bounded
! two-sided quarter-interval is perfectly resolved.  Raise the provisional
! step above roundoff before clipping it to the certified overlap.
    h=max(1.d-7*(collision_jhi-collision_jlo),4096.d0*epsilon(1.d0)*scale)
    if(jtest.le.collision_jlo+ h) then
      h=min(h,0.25d0*(collision_jhi-collision_jlo))
      if(h.le.256.d0*epsilon(1.d0)*scale) then
        root_eval_valid=.false.
        root_eval_error=root_invalid_domain
        ddifference_dj=0.d0
        return
      endif
      call fixedpoint_collision_value(jtest+h,dp,okp)
      if(.not.okp) then
        root_eval_valid=.false.
        root_eval_error=root_invalid_domain
        difference=0.d0
        ddifference_dj=0.d0
        return
      endif
      ddifference_dj=(dp-difference)/h
      root_eval_valid=.true.
      root_eval_error=0
      return
    elseif(jtest.ge.collision_jhi-h) then
      h=min(h,0.25d0*(collision_jhi-collision_jlo))
      if(h.le.256.d0*epsilon(1.d0)*scale) then
        root_eval_valid=.false.
        root_eval_error=root_invalid_domain
        ddifference_dj=0.d0
        return
      endif
      call fixedpoint_collision_value(jtest-h,dm,okm)
      if(.not.okm) then
        root_eval_valid=.false.
        root_eval_error=root_invalid_domain
        difference=0.d0
        ddifference_dj=0.d0
        return
      endif
      ddifference_dj=(difference-dm)/h
      root_eval_valid=.true.
      root_eval_error=0
      return
    endif
    h=min(h,0.25d0*(jtest-collision_jlo),0.25d0*(collision_jhi-jtest))
    if(h.le.256.d0*epsilon(1.d0)*scale) then
      root_eval_valid=.false.
      root_eval_error=root_invalid_domain
      ddifference_dj=0.d0
      return
    endif
    call fixedpoint_collision_value(jtest-h,dm,okm)
    call fixedpoint_collision_value(jtest+h,dp,okp)
    if(.not.okm .or. .not.okp) then
        root_eval_valid=.false.
      root_eval_error=root_invalid_domain
      difference=0.d0
      ddifference_dj=0.d0
      return
    endif
    ddifference_dj=(dp-dm)/(2.d0*h)
    root_eval_valid=.true.
    root_eval_error=0
  end subroutine fixedpoint_collision

  subroutine fixedpoint_collision_value(jtest,difference,ok_value)
    double precision, intent(in) :: jtest
    double precision, intent(out) :: difference
    logical, intent(out) :: ok_value
    double precision :: rleft,rright,pleft,pright
    logical :: okleft,okright

    if(collision_boundary_left.gt.0) then
      call invert_boundary_segment(collision_boundary_left,jtest,rleft,pleft,okleft)
    else
      call invert_fixedpoint_branch(collision_segment_left,jtest,rleft,pleft,okleft)
    endif
    if(collision_boundary_right.eq.0 .and. collision_boundary_segment.eq.0) then
      call invert_fixedpoint_branch(collision_segment_right,jtest,rright,pright,okright)
    else
      if(collision_boundary_right.gt.0) then
        call invert_boundary_segment(collision_boundary_right,jtest,rright,pright,okright)
      else
        call invert_boundary_segment(collision_boundary_segment,jtest,rright,pright,okright)
      endif
    endif
    ok_value=okleft .and. okright
    if(ok_value) then
      difference=pleft-pright
    else
      difference=0.d0
    endif
  end subroutine fixedpoint_collision_value

  subroutine refine_collision_candidate(seed,refined,ok_refined)
! Refine a collision root only after the broad all-roots scan has found it.
! The refinement is bounded by the same certified overlap and is accepted
! only when a valid sign bracket exists.  This removes the scan-scale root
! error without turning an unrelated invalid branch endpoint into a guessed
! topology boundary.
    double precision, intent(in) :: seed
    double precision, intent(out) :: refined
    logical, intent(out) :: ok_refined
    integer :: local_it,local_expand
    double precision :: a,b,mid,fa,fb,fmid,step,scale
    logical :: oka,okb,okmid

    refined=seed
    ok_refined=.false.
    if(seed.lt.collision_jlo .or. seed.gt.collision_jhi) return
    scale=max(1.d-300,abs(seed),abs(collision_jlo),abs(collision_jhi))
    step=max(1024.d0*epsilon(1.d0)*scale, &
             1.d-12*max(1.d-300,collision_jhi-collision_jlo))
    do local_expand=1,32
      a=max(collision_jlo,seed-step)
      b=min(collision_jhi,seed+step)
      if(a.ge.b) exit
      call fixedpoint_collision_value(a,fa,oka)
      call fixedpoint_collision_value(b,fb,okb)
      if(oka .and. okb .and. (fa.eq.0.d0 .or. fb.eq.0.d0 .or. fa*fb.lt.0.d0)) exit
      step=2.d0*step
    enddo
    if(.not.oka .or. .not.okb) return
    if(fa.ne.0.d0 .and. fb.ne.0.d0 .and. fa*fb.gt.0.d0) return
    do local_it=1,200
      mid=0.5d0*(a+b)
      if(mid.eq.a .or. mid.eq.b) exit
      call fixedpoint_collision_value(mid,fmid,okmid)
      if(.not.okmid) then
        ok_refined=.false.
        return
      endif
      if(fmid.eq.0.d0) then
        a=mid
        b=mid
        exit
      elseif(fa.eq.0.d0 .or. fa*fmid.le.0.d0) then
        b=mid
        fb=fmid
      else
        a=mid
        fa=fmid
      endif
    enddo
    ! Use the lower representable side of the bracket.  The sampler then
    ! probes the open interval below the event instead of rounding the event
    ! upward into the newly created class.
    refined=min(a,b)
    ok_refined=.true.
  end subroutine refine_collision_candidate

  subroutine invert_fixedpoint_branch(segment,jtarget,rvalue,pvalue,ok_value)
    integer, intent(in) :: segment
    double precision, intent(in) :: jtarget
    double precision, intent(out) :: rvalue,pvalue
    logical, intent(out) :: ok_value
    integer :: local_it
    double precision :: rlo_local,rhi_local,jlo_local,jhi_local,rmid, &
                        jmid,pmid,tolerance,jlo_value,jhi_value, &
                        plo_value,phi_value,radial_resolution,jwidth,weight,jscale
    logical :: okmid,increasing,bracketed

    rvalue=0.d0
    pvalue=0.d0
    ok_value=.false.
    jlo_local=fp_jlo(segment)
    jhi_local=fp_jhi(segment)
    ! The inverse is evaluated on J_perp itself.  An absolute tolerance
    ! scaled by one is many orders too large for the small-J outer scan and
    ! can move a fixed-point/boundary collision across its true topology
    ! transition.  Keep the residual relative to this local J interval.
    jscale=max(1.d-300,abs(jtarget),abs(jlo_local),abs(jhi_local))
    tolerance=max(1024.d0*epsilon(1.d0)*jscale,1.d-12*jscale)
    if(jtarget.lt.min(jlo_local,jhi_local)-tolerance .or. &
       jtarget.gt.max(jlo_local,jhi_local)+tolerance) return
    if(abs(jtarget-jlo_local).le.tolerance) then
      rvalue=fp_rlo(segment)
      pvalue=fp_plo(segment)
      ok_value=.true.
      return
    endif
    if(abs(jtarget-jhi_local).le.tolerance) then
      rvalue=fp_rhi(segment)
      pvalue=fp_phi(segment)
      ok_value=.true.
      return
    endif
    rlo_local=fp_rlo(segment)
    rhi_local=fp_rhi(segment)
    jlo_value=jlo_local
    jhi_value=jhi_local
    plo_value=fp_plo(segment)
    phi_value=fp_phi(segment)
    increasing=jhi_local.gt.jlo_local
    do local_it=1,100
      rmid=0.5d0*(rlo_local+rhi_local)
      call fixedpoint_branch_value(rmid,fp_sigma(segment),fp_branch(segment), &
                                   jmid,pmid,okmid)
      if(.not.okmid) then
        return
      endif
      if(abs(jmid-jtarget).le.tolerance) then
        rvalue=rmid
        pvalue=pmid
        ok_value=.true.
        return
      endif
      if(increasing) then
        if(jmid.lt.jtarget) then
          rlo_local=rmid
          jlo_value=jmid
          plo_value=pmid
        else
          rhi_local=rmid
          jhi_value=jmid
          phi_value=pmid
        endif
      else
        if(jmid.gt.jtarget) then
          rlo_local=rmid
          jlo_value=jmid
          plo_value=pmid
        else
          rhi_local=rmid
          jhi_value=jmid
          phi_value=pmid
        endif
      endif
    enddo
    rvalue=rmid
    pvalue=pmid
    ok_value=abs(jmid-jtarget).le.10.d0*tolerance
    if(.not.ok_value) then
      radial_resolution=256.d0*epsilon(1.d0)*max(1.d0,abs(rlo_local), &
                                                   abs(rhi_local))
      bracketed=.false.
      if(abs(rhi_local-rlo_local).le.radial_resolution) then
        if(jtarget.ge.min(jlo_value,jhi_value)-tolerance) then
          if(jtarget.le.max(jlo_value,jhi_value)+tolerance) bracketed=.true.
        endif
      endif
      if(bracketed) then
        ! The continuous inverse lies between adjacent representable R values.
        ! Interpolate the level value across that machine-resolution bracket;
        ! a wider or unbracketed gap remains an invalid topology certificate.
        jwidth=jhi_value-jlo_value
        if(abs(jwidth).gt.tolerance) then
          weight=(jtarget-jlo_value)/jwidth
        else
          weight=0.5d0
        endif
        weight=max(0.d0,min(1.d0,weight))
        rvalue=rlo_local+weight*(rhi_local-rlo_local)
        pvalue=plo_value+weight*(phi_value-plo_value)
        ok_value=.true.
        return
      endif
    endif
  end subroutine invert_fixedpoint_branch

  subroutine invert_boundary_segment(segment,jtarget,rvalue,pvalue,ok_value)
    integer, intent(in) :: segment
    double precision, intent(in) :: jtarget
    double precision, intent(out) :: rvalue,pvalue
    logical, intent(out) :: ok_value
    integer :: local_it
    double precision :: rlo_local,rhi_local,jlo_local,jhi_local,rmid, &
                        jmid,pmid,tolerance
    logical :: okmid,increasing

    rvalue=0.d0
    pvalue=0.d0
    ok_value=.false.
    jlo_local=bd_jlo(segment)
    jhi_local=bd_jhi(segment)
    ! Match the fixed-point inverse: the boundary level must be resolved at
    ! the scale of the local J_perp interval, not at an absolute unit scale.
    tolerance=max(1024.d0*epsilon(1.d0)*max(1.d-300,abs(jtarget), &
        abs(jlo_local),abs(jhi_local)),1.d-12*max(1.d-300,abs(jtarget), &
        abs(jlo_local),abs(jhi_local)))
    if(jtarget.lt.min(jlo_local,jhi_local)-tolerance .or. &
       jtarget.gt.max(jlo_local,jhi_local)+tolerance) return
    if(bd_type(segment).eq.2) then
      call fixed_boundary_level(bd_rlo(segment),jtarget,bd_sigma(segment), &
                                pvalue,ok_value)
      rvalue=bd_rlo(segment)
      return
    endif
    if(abs(jtarget-jlo_local).le.tolerance) then
      rvalue=bd_rlo(segment)
      call turning_boundary_level(rvalue,pvalue,ok_value)
      return
    endif
    if(abs(jtarget-jhi_local).le.tolerance) then
      rvalue=bd_rhi(segment)
      call turning_boundary_level(rvalue,pvalue,ok_value)
      return
    endif
    rlo_local=bd_rlo(segment)
    rhi_local=bd_rhi(segment)
    increasing=jhi_local.gt.jlo_local
    do local_it=1,100
      rmid=0.5d0*(rlo_local+rhi_local)
      call jperp_value(rmid,jmid,okmid)
      if(.not.okmid) return
      if(abs(jmid-jtarget).le.tolerance) then
        rvalue=rmid
        call turning_boundary_level(rvalue,pvalue,ok_value)
        return
      endif
      if(increasing) then
        if(jmid.lt.jtarget) then
          rlo_local=rmid
        else
          rhi_local=rmid
        endif
      else
        if(jmid.gt.jtarget) then
          rlo_local=rmid
        else
          rhi_local=rmid
        endif
      endif
    enddo
    rvalue=rmid
    call turning_boundary_level(rvalue,pvalue,ok_value)
    ok_value=ok_value .and. abs(jmid-jtarget).le.10.d0*tolerance
  end subroutine invert_boundary_segment

  subroutine turning_boundary_level(local_r,local_p,local_ok)
    use field_sub, only : psif
    double precision, intent(in) :: local_r
    double precision, intent(out) :: local_p
    logical, intent(out) :: local_ok
    double precision :: local_z,local_dz,local_bmod,local_phi
    double precision, dimension(3) :: local_x

    call get_poicut(local_r,local_z,local_dz)
    local_x(1)=local_r
    local_x(2)=0.d0
    local_x(3)=local_z
    call get_bmod_and_Phi(local_x,local_bmod,local_phi)
    local_p=psif
    local_ok=(local_p.eq.local_p .and. abs(local_p).lt.huge(local_p))
  end subroutine turning_boundary_level

  subroutine fixed_boundary_level(local_r,local_j,local_sigma,local_p,local_ok)
    use parmot_mod, only : ro0
    double precision, intent(in) :: local_r,local_j
    integer, intent(in) :: local_sigma
    double precision, intent(out) :: local_p
    logical, intent(out) :: local_ok
    double precision :: local_a,local_b,local_ap,local_bp,local_ps,local_psp, &
                        local_h,local_hp,local_u2,tolerance

    call fixedpoint_geometry(local_r,local_a,local_b,local_ap,local_bp, &
                             local_ps,local_psp,local_h,local_hp,local_ok)
    if(.not.local_ok) then
      local_p=0.d0
      return
    endif
    local_u2=local_a-local_j*local_b
    tolerance=1.d-11*max(1.d0,abs(local_a),abs(local_j*local_b))
    if(local_u2.lt.-tolerance) then
      local_ok=.false.
      local_p=0.d0
      return
    endif
    local_p=local_ps+ro0*local_h*dble(local_sigma)*sqrt(max(0.d0,local_u2))
  end subroutine fixed_boundary_level

  subroutine jperp_value(R,value,ok_value)
    use potato_symbolic_kernel_mod, only : potato_jperp_kernel
    double precision, intent(in) :: R
    double precision, intent(out) :: value
    logical, intent(out) :: ok_value
    double precision :: Z,dZ_dR,bmod,sqrtg,phi_elec
    double precision :: qPhi_prime,magnetic_field_B_prime
    double precision :: candidate,positive,derivative
    double precision, dimension(3) :: xx,bder,hcovar,hctrvr,hcurl,derphi

    call get_poicut(R,Z,dZ_dR)
    xx(1)=R
    xx(2)=0.d0
    xx(3)=Z
    call magfie(xx,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    call elefie(xx,phi_elec,derphi)
    ok_value=(bmod.gt.0.d0 .and. bmod.eq.bmod .and. abs(bmod).lt.huge(bmod))
    if(.not.ok_value) then
      value=0.d0
      return
    endif
    qPhi_prime=derphi(1)+derphi(3)*dZ_dR
    magnetic_field_B_prime=bmod*(bder(1)+bder(3)*dZ_dR)
    call potato_jperp_kernel(toten,phi_elec,bmod,qPhi_prime, &
                             magnetic_field_B_prime, &
                             candidate,positive,derivative)
    ok_value=(candidate.eq.candidate .and. abs(candidate).lt.huge(candidate))
    if(ok_value) then
      value=candidate
    else
      value=0.d0
    endif
  end subroutine jperp_value

  subroutine jperp_value_and_derivative(R,value,dvalue,ok_value)
    use global_invariants, only : toten
    use potato_symbolic_kernel_mod, only : potato_jperp_kernel
    double precision, intent(in) :: R
    double precision, intent(out) :: value,dvalue
    logical, intent(out) :: ok_value
    double precision :: Z,dZ_dR,bmod,sqrtg,phi_elec
    double precision :: qPhi_prime,magnetic_field_B_prime
    double precision :: candidate,positive,derivative
    double precision, dimension(3) :: xx,bder,hcovar,hctrvr,hcurl,derphi

    call get_poicut(R,Z,dZ_dR)
    xx(1)=R
    xx(2)=0.d0
    xx(3)=Z
    call magfie(xx,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
    call elefie(xx,phi_elec,derphi)
    ok_value=(bmod.gt.0.d0 .and. bmod.eq.bmod .and. abs(bmod).lt.huge(bmod))
    if(.not.ok_value) then
      value=0.d0
      dvalue=0.d0
      return
    endif
    qPhi_prime=derphi(1)+derphi(3)*dZ_dR
    magnetic_field_B_prime=bmod*(bder(1)+bder(3)*dZ_dR)
    call potato_jperp_kernel(toten,phi_elec,bmod,qPhi_prime, &
                             magnetic_field_B_prime, &
                             candidate,positive,derivative)
    if(candidate.ne.candidate .or. abs(candidate).gt.huge(candidate) .or. &
       derivative.ne.derivative .or. abs(derivative).gt.huge(derivative)) then
      ok_value=.false.
      value=0.d0
      dvalue=0.d0
      return
    endif
    value=candidate
    dvalue=derivative
  end subroutine jperp_value_and_derivative

  subroutine jperp_stationary(R,stationary,dstationary)
    use find_all_roots_mod, only : root_eval_valid,root_eval_error, &
                                   root_invalid_domain
    double precision, intent(in) :: R
    double precision, intent(out) :: stationary,dstationary
    double precision :: value,dvalue,value_m,dvalue_m,value_p,dvalue_p,h
    logical :: ok_m,ok_p,ok

    root_eval_valid=.true.
    root_eval_error=0
    call jperp_value_and_derivative(R,value,stationary,ok)
    if(.not.ok) then
      root_eval_valid=.false.
      root_eval_error=root_invalid_domain
      stationary=0.d0
      dstationary=0.d0
      return
    endif
    h=max(1.d-8*range,256.d0*epsilon(1.d0)*max(1.d0,abs(R)))
    h=min(h,0.25d0*(R-rpc_arr(0)),0.25d0*(rpc_arr(npc)-R))
    if(h.le.0.d0) then
      dstationary=0.d0
      return
    endif
    call jperp_value_and_derivative(R-h,value_m,dvalue_m,ok_m)
    call jperp_value_and_derivative(R+h,value_p,dvalue_p,ok_p)
    if(.not.ok_m .or. .not.ok_p) then
      root_eval_valid=.false.
      root_eval_error=root_invalid_domain
      stationary=0.d0
      dstationary=0.d0
      return
    endif
    dstationary=(dvalue_p-dvalue_m)/(2.d0*h)
  end subroutine jperp_stationary

  end subroutine find_jperp_topology_boundaries
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine fixedpoint_geometry(R,A,B,Aprime,Bprime,psistar,psiprime, &
                                 hphi,hphiprime,ok)
!
! Geometry of the exact cut representation
!
!   psi^*(R_c;H,J,sigma) = psi(R_c) + rho0*h_phi(R_c)*u,
!   u^2 = A(R_c)-J*B(R_c),  A=H-Phi.
!
! The derivative of h_phi is evaluated on the bounded cut.  At a cut endpoint
! the one-sided geometric derivative is used only to evaluate a boundary
! equation; fixed-point classification itself still requires a strict
! two-sided neighbourhood in determine_fixpoint_type.
!
  use field_sub, only : psif,dpsidr,dpsidz
  use global_invariants, only : toten
  use poicut_mod, only : npc,rpc_arr
  implicit none
!
  double precision, intent(in) :: R
  double precision, intent(out) :: A,B,Aprime,Bprime,psistar,psiprime, &
                                   hphi,hphiprime
  logical, intent(out) :: ok
  double precision :: Z,dZ_dR,bmod,sqrtg,phi_elec,hminus,hplus,hstep
  double precision :: scale,range
  double precision, dimension(3) :: xx,bder,hcovar,hctrvr,hcurl,derphi
  logical :: ok_minus,ok_plus

  ok=.false.
  A=0.d0
  B=0.d0
  Aprime=0.d0
  Bprime=0.d0
  psistar=0.d0
  psiprime=0.d0
  hphi=0.d0
  hphiprime=0.d0
  if(npc.le.0 .or. R.lt.rpc_arr(0) .or. R.gt.rpc_arr(npc)) return
  range=rpc_arr(npc)-rpc_arr(0)
  if(range.le.0.d0) return

  call get_poicut(R,Z,dZ_dR)
  xx(1)=R
  xx(2)=0.d0
  xx(3)=Z
  call magfie(xx,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
  call elefie(xx,phi_elec,derphi)
  if(bmod.le.0.d0) return

  A=toten-phi_elec
  B=bmod
  Aprime=-derphi(1)-derphi(3)*dZ_dR
  Bprime=bmod*(bder(1)+bder(3)*dZ_dR)
  psistar=psif
  psiprime=dpsidr+dpsidz*dZ_dR
  hphi=hcovar(2)

  scale=max(1.d0,abs(R),abs(rpc_arr(0)),abs(rpc_arr(npc)))
  hstep=max(1.d-7*range,256.d0*epsilon(1.d0)*scale)
  if(R.gt.rpc_arr(0) .and. R.lt.rpc_arr(npc)) then
    hstep=min(hstep,0.25d0*(R-rpc_arr(0)), &
                    0.25d0*(rpc_arr(npc)-R))
  else
    hstep=min(hstep,0.25d0*range)
  endif
  if(hstep.le.0.d0) return
  ok_minus=.false.
  ok_plus=.false.
  if(R.gt.rpc_arr(0)) then
    call fixedpoint_hphi(R-hstep,hminus,ok_minus)
  endif
  if(R.lt.rpc_arr(npc)) then
    call fixedpoint_hphi(R+hstep,hplus,ok_plus)
  endif
  if(ok_minus .and. ok_plus) then
    hphiprime=(hplus-hminus)/(2.d0*hstep)
  elseif(ok_plus) then
    hphiprime=(hplus-hphi)/hstep
  elseif(ok_minus) then
    hphiprime=(hphi-hminus)/hstep
  else
    return
  endif
  ok=.true.

contains

  subroutine fixedpoint_hphi(Rin,hout,okout)
    double precision, intent(in) :: Rin
    double precision, intent(out) :: hout
    logical, intent(out) :: okout
    double precision :: Zin,dZin,bmodin,sqrtgin,phiin
    double precision, dimension(3) :: xin,bderin,hcovarin,hctrvrin,hcurlin,derphiin

    okout=.false.
    hout=0.d0
    if(Rin.lt.rpc_arr(0) .or. Rin.gt.rpc_arr(npc)) return
    call get_poicut(Rin,Zin,dZin)
    xin(1)=Rin
    xin(2)=0.d0
    xin(3)=Zin
    call magfie(xin,bmodin,sqrtgin,bderin,hcovarin,hctrvrin,hcurlin)
    call elefie(xin,phiin,derphiin)
    if(bmodin.le.0.d0) return
    hout=hcovarin(2)
    okout=.true.
  end subroutine fixedpoint_hphi

  end subroutine fixedpoint_geometry

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine fixedpoint_coefficients(R,acoef,bcoef,ccoef,energy_a,bfield, &
                                    psistar,hphi,ok)
!
! For u=sigma*sqrt(A-J*B), d psi^*/d R_c=0 reduces exactly to
!
!   a u^2+b u+c=0,
!
! with no square-root extrapolation across a reflecting boundary.
!
  use parmot_mod, only : ro0
  implicit none
  double precision, intent(in) :: R
  double precision, intent(out) :: acoef,bcoef,ccoef,energy_a,bfield, &
                                   psistar,hphi
  logical, intent(out) :: ok
  double precision :: Aprime,Bprime,psiprime,hprime

  call fixedpoint_geometry(R,energy_a,bfield,Aprime,Bprime,psistar,psiprime, &
                           hphi,hprime,ok)
  if(.not.ok) then
    acoef=0.d0
    bcoef=0.d0
    ccoef=0.d0
    hphi=0.d0
    return
  endif
  acoef=ro0*(2.d0*hprime*bfield+hphi*Bprime)
  bcoef=2.d0*bfield*psiprime
  ccoef=ro0*hphi*(Aprime*bfield-energy_a*Bprime)
  end subroutine fixedpoint_coefficients

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine fixedpoint_roots_at_R(R,sigma_value,uroots,jroots,proots, &
                                   nfixed,ok)
  use parmot_mod, only : ro0
  use potato_topology_mod, only : solve_fixedpoint_quadratic
  implicit none
  double precision, intent(in) :: R
  integer, intent(in) :: sigma_value
  double precision, intent(out) :: uroots(2),jroots(2),proots(2)
  integer, intent(out) :: nfixed
  logical, intent(out) :: ok
  double precision :: acoef,bcoef,ccoef,energy_a,bfield,psistar,hphi
  double precision :: u_tol
  double precision, dimension(2) :: all_u
  integer :: nall,i
  logical :: coeff_ok,quadratic_ok

  nfixed=0
  uroots=0.d0
  jroots=0.d0
  proots=0.d0
  call fixedpoint_coefficients(R,acoef,bcoef,ccoef,energy_a,bfield, &
                               psistar,hphi,coeff_ok)
  ok=coeff_ok
  if(.not.ok .or. bfield.le.0.d0) return
  call solve_fixedpoint_quadratic(acoef,bcoef,ccoef,all_u,nall,quadratic_ok)
  if(.not.quadratic_ok) then
    ok=.false.
    return
  endif
  u_tol=128.d0*epsilon(1.d0)*max(1.d0,sqrt(abs(energy_a)))
  do i=1,nall
    if(dble(sigma_value)*all_u(i).le.u_tol) cycle
    if(nfixed.gt.0) then
      if(abs(all_u(i)-uroots(nfixed)).le.128.d0*epsilon(1.d0)* &
         max(1.d0,abs(all_u(i)))) cycle
    endif
    nfixed=nfixed+1
    uroots(nfixed)=all_u(i)
    jroots(nfixed)=(energy_a-all_u(i)**2)/bfield
    proots(nfixed)=psistar+ro0*hphi*all_u(i)
  enddo
  end subroutine fixedpoint_roots_at_R

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine fixedpoint_branch_value(R,sigma_value,branch,jvalue,pstar,ok)
  implicit none
  double precision, intent(in) :: R
  integer, intent(in) :: sigma_value
  integer, intent(in) :: branch
  double precision, intent(out) :: jvalue,pstar
  logical, intent(out) :: ok
  double precision :: u(2),j(2),p(2)
  integer :: nfixed

  call fixedpoint_roots_at_R(R,sigma_value,u,j,p,nfixed,ok)
  if(.not.ok .or. branch.lt.1 .or. branch.gt.nfixed) then
    ok=.false.
    jvalue=0.d0
    pstar=0.d0
    return
  endif
  jvalue=j(branch)
  pstar=p(branch)
  ok=.true.
  end subroutine fixedpoint_branch_value

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine fixedpoint_discriminant(R,discriminant,ddiscriminant_dR)
  use find_all_roots_mod, only : root_eval_valid,root_eval_error, &
                                 root_invalid_domain
  use poicut_mod, only : npc,rpc_arr
  implicit none
  double precision, intent(in) :: R
  double precision, intent(out) :: discriminant,ddiscriminant_dR
  double precision :: dcenter,dminus,dplus,h,range
  logical :: ok,ok_m,ok_p

  root_eval_valid=.true.
  root_eval_error=0
  call fixedpoint_discriminant_value(R,dcenter,ok)
  if(.not.ok) then
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
    discriminant=0.d0
    ddiscriminant_dR=0.d0
    return
  endif
  range=rpc_arr(npc)-rpc_arr(0)
  h=max(1.d-7*range,256.d0*epsilon(1.d0)*max(1.d0,abs(R)))
  if(R.gt.rpc_arr(0) .and. R.lt.rpc_arr(npc)) then
    h=min(h,0.25d0*(R-rpc_arr(0)),0.25d0*(rpc_arr(npc)-R))
  else
    h=min(h,0.25d0*range)
  endif
  ok_m=.false.
  ok_p=.false.
  if(R.gt.rpc_arr(0)) then
    call fixedpoint_discriminant_value(R-h,dminus,ok_m)
  endif
  if(R.lt.rpc_arr(npc)) then
    call fixedpoint_discriminant_value(R+h,dplus,ok_p)
  endif
  if(ok_m .and. ok_p) then
    ddiscriminant_dR=(dplus-dminus)/(2.d0*h)
  elseif(ok_p) then
    ddiscriminant_dR=(dplus-dcenter)/h
  elseif(ok_m) then
    ddiscriminant_dR=(dcenter-dminus)/h
  else
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
    discriminant=0.d0
    ddiscriminant_dR=0.d0
  endif
  discriminant=dcenter
  end subroutine fixedpoint_discriminant

  subroutine fixedpoint_turning_intersection(R,intersection, &
                                             dintersection_dR)
  use find_all_roots_mod, only : root_eval_valid,root_eval_error, &
                                 root_invalid_domain
  use poicut_mod, only : npc,rpc_arr
  implicit none
  double precision, intent(in) :: R
  double precision, intent(out) :: intersection,dintersection_dR
  double precision :: acoef,bcoef,ccoef,energy_a,bfield,psistar,hphi
  double precision :: cminus,cplus,unused_a,unused_b,unused_e,unused_bfield, &
                      unused_psistar,unused_hphi,h,range
  logical :: ok,okminus,okplus

  root_eval_valid=.true.
  root_eval_error=0
  call fixedpoint_coefficients(R,acoef,bcoef,ccoef,energy_a,bfield, &
                               psistar,hphi,ok)
  if(.not.ok) then
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
    intersection=0.d0
    dintersection_dR=0.d0
    return
  endif
  intersection=ccoef
  range=rpc_arr(npc)-rpc_arr(0)
  h=max(1.d-7*range,256.d0*epsilon(1.d0)*max(1.d0,abs(R)))
  if(R.gt.rpc_arr(0) .and. R.lt.rpc_arr(npc)) then
    h=min(h,0.25d0*(R-rpc_arr(0)),0.25d0*(rpc_arr(npc)-R))
  else
    h=min(h,0.25d0*range)
  endif
  okminus=.false.
  okplus=.false.
  if(R.gt.rpc_arr(0)) then
    call fixedpoint_coefficients(R-h,unused_a,unused_b,cminus,unused_e, &
                                 unused_bfield,unused_psistar,unused_hphi, &
                                 okminus)
  endif
  if(R.lt.rpc_arr(npc)) then
    call fixedpoint_coefficients(R+h,unused_a,unused_b,cplus,unused_e, &
                                 unused_bfield,unused_psistar,unused_hphi, &
                                 okplus)
  endif
  if(okminus .and. okplus) then
    dintersection_dR=(cplus-cminus)/(2.d0*h)
  elseif(okplus) then
    dintersection_dR=(cplus-ccoef)/h
  elseif(okminus) then
    dintersection_dR=(ccoef-cminus)/h
  else
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
    intersection=0.d0
    dintersection_dR=0.d0
  endif
  end subroutine fixedpoint_turning_intersection

  subroutine fixedpoint_discriminant_value(R,discriminant,ok)
  implicit none
  double precision, intent(in) :: R
  double precision, intent(out) :: discriminant
  logical, intent(out) :: ok
  double precision :: acoef,bcoef,ccoef,energy_a,bfield,psistar,hphi

  call fixedpoint_coefficients(R,acoef,bcoef,ccoef,energy_a,bfield, &
                               psistar,hphi,ok)
  if(ok) then
    discriminant=bcoef*bcoef-4.d0*acoef*ccoef
  else
    discriminant=0.d0
  endif
  end subroutine fixedpoint_discriminant_value

!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  subroutine fixedpoint_branch_stationary(R,stationary,dstationary)
    use find_all_roots_mod, only : root_eval_valid,root_eval_error, &
                                 root_invalid_domain,root_left_endpoint_contracted, &
                                 root_right_endpoint_contracted,root_search_left, &
                                 root_search_right
  use potato_boundary_scan_mod, only : fixedpoint_scan_sigma, &
                                       fixedpoint_scan_branch, &
                                       fixedpoint_scan_left, fixedpoint_scan_right
  use poicut_mod, only : npc,rpc_arr
  implicit none
  double precision, intent(in) :: R
  double precision, intent(out) :: stationary,dstationary
  double precision :: j0,jm,jp,j2,h,range,left,right
  logical :: ok0,okm,okp,ok2

  root_eval_valid=.true.
  root_eval_error=0
  left=fixedpoint_scan_left
  right=fixedpoint_scan_right
  if(right.le.left) then
    left=rpc_arr(0)
    right=rpc_arr(npc)
  endif
! find_all_roots may have contracted an invalid physical endpoint.  The
! finite-difference classifier must use that certified open endpoint during
! subsequent refinement; probing against the original forbidden boundary
! would manufacture a second, spurious invalid-domain failure.
  if(root_left_endpoint_contracted) left=max(left,root_search_left)
  if(root_right_endpoint_contracted) right=min(right,root_search_right)
  if(right.le.left) then
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
    stationary=0.d0
    dstationary=0.d0
    return
  endif
  range=right-left
  call fixedpoint_branch_value(R,fixedpoint_scan_sigma,fixedpoint_scan_branch, &
                               j0,stationary,ok0)
  if(.not.ok0) then
    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
    stationary=0.d0
    dstationary=0.d0
    return
  endif
  h=max(1.d-7*range,256.d0*epsilon(1.d0)*max(1.d0,abs(R)))
  h=min(h,0.25d0*range)
  okm=.false.
  okp=.false.
  if(R.gt.left) then
    h=min(h,0.25d0*(R-left))
    call fixedpoint_branch_value(R-h,fixedpoint_scan_sigma, &
                                 fixedpoint_scan_branch,jm,stationary,okm)
  endif
  if(R.lt.right) then
    h=min(h,0.25d0*(right-R))
    call fixedpoint_branch_value(R+h,fixedpoint_scan_sigma, &
                                 fixedpoint_scan_branch,jp,stationary,okp)
  endif
  if(okm .and. okp) then
    stationary=(jp-jm)/(2.d0*h)
    dstationary=(jp-2.d0*j0+jm)/(h*h)
    return
  endif
  if(okp) then
    call fixedpoint_branch_value(R+2.d0*h,fixedpoint_scan_sigma, &
                                 fixedpoint_scan_branch,j2,stationary,ok2)
    if(ok2) then
      stationary=(jp-j0)/h
      dstationary=(j2-2.d0*jp+j0)/(h*h)
      return
    endif
  endif
  if(okm) then
    call fixedpoint_branch_value(R-2.d0*h,fixedpoint_scan_sigma, &
                                 fixedpoint_scan_branch,j2,stationary,ok2)
    if(ok2) then
      stationary=(j0-jm)/h
      dstationary=(j0-2.d0*jm+j2)/(h*h)
      return
    endif
  endif
  root_eval_valid=.false.
  root_eval_error=root_invalid_domain
  stationary=0.d0
  dstationary=0.d0
  end subroutine fixedpoint_branch_stationary

!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
