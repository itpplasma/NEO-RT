!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Modules:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module phielec_of_psi_mod
    integer, parameter :: npolyphi=3 !2 !1
    double precision, dimension(0:npolyphi) :: polyphi
  end module phielec_of_psi_mod
!
  module cc_mod
    logical :: dowrite=.false. !write canonical frequencies and bounce integrals for sampling grid and interpolation
    logical :: wrbounds=.false. !write vpar^2 and psi^* curves vs cut parameter R_c, extremum and boundary points
  end module cc_mod
!
  module global_invariants
! integration step, normalized total energy, normalized perpendicular invariant, velocity sign:
    double precision :: dtau,toten,perpinv,sigma
! reference energy $\cE_{ref}$, effective potential $\Phi_{eff}$:
    double precision :: cE_ref,Phi_eff
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
    integer :: nbounds,nregions
    double precision,     dimension(:),   allocatable :: R_bo,Z_bo,psiast_bo
    type(allowed_region), dimension(:,:), allocatable :: all_regions
  end module bounds_fixpoints_mod
!
  module form_classes_doublecount_mod
! number of classes for fixed H_0 and J_perp:
    integer :: nclasses
! type of class parameterization function:
    integer,          dimension(:), allocatable :: ifuntype
! class interval limits over cut parameter R_c and parallel velocity sign:
    double precision, dimension(:), allocatable :: R_class_beg,R_class_end,sigma_class
  end module form_classes_doublecount_mod
!
!------------------------------------------------------
!
  module orbit_dim_mod
    integer, parameter  :: neqm=5
    logical             :: write_orb=.true.
    integer             :: iunit1=100,next,numbasef
    double precision    :: Rorb_max
  end module orbit_dim_mod
!
!------------------------------------------------------
!
  module get_matrix_mod
    double precision :: relerror=1.d-4, relmargin=1.d-4
    integer          :: iclass
  end module get_matrix_mod
!
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Routines:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
!
  subroutine find_bounce(next,velo_ext,dtau_in,z_eqm,taub,delphi,extraset)
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
! extraset(next) - extra integrals along the orbit
!
  use orbit_dim_mod, only : neqm,write_orb,iunit1,Rorb_max
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
!
! relative error of orbit integrator:
  double precision, parameter :: relerr=1d-10 !8
!
  logical :: firstpass
  integer :: next, ndim, iter
  double precision :: dtau_in,dtau,taub,delphi
  double precision :: dL2_pol,dL2_pol_start,dtau_newt,r_prev,z_prev
  double precision :: tau0,RNorm,ZNorm,vnorm,dnorm,vel_pol,dL2_pol_min
  double precision :: dZ_dR,sign_delZ,Z_tmp
  double precision, dimension(neqm) :: z_eqm
  double precision, dimension(next) :: extraset
  double precision, dimension(neqm+next) :: z,z_start,vz
!
  external velo_ext
!
  ndim = neqm+next
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
  call velo_ext(dtau,z,vz)
!
! unit 2D vector in the direction of the guiding center velocity in RZ-plane:
!
  vel_pol=sqrt(vz(1)**2+vz(3)**2)
  RNorm=vz(1)/vel_pol
  ZNorm=vz(3)/vel_pol
!
  tau0=0.d0
!
  if(write_orb) write (iunit1,*) z(1:neqm),vz(5)
!
! first step:
!
  call odeint_allroutines(z,ndim,tau0,dtau,relerr,velo_ext)
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
    r_prev=z(1)
    z_prev=z(3)
!
    call odeint_allroutines(z,ndim,tau0,dtau,relerr,velo_ext)
!
    taub=taub+dtau
    if(nousecut) then
! check if poloidal distance to the starting point is smaller than
! sqrt(2) of the poloidal length of the step:
      dL2_pol=2.d0*((z(1)-r_prev)**2+(z(3)-z_prev)**2)
      dL2_pol_start=(z(1)-z_start(1))**2+(z(3)-z_start(3))**2
      if(dL2_pol_start.lt.dL2_pol) exit
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
! Newton adjustment
!
  do iter=1,niter
!
    call velo_ext(dtau,z,vz)
!
    vnorm=vz(1)*RNorm+vz(3)*ZNorm
    dnorm=(z_start(1)-z(1))*RNorm+(z_start(3)-z(3))*ZNorm
    if(dnorm**2.lt.dL2_pol*relerr) exit
    dtau_newt=dnorm/vnorm
!
    call odeint_allroutines(z,ndim,tau0,dtau_newt,relerr,velo_ext)
!
    taub=taub+dtau_newt
  enddo
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
  use field_eq_mod, only : psif     ,dpsidr,dpsidz
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
! 2 - negative parallel kinetic energy
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
  alam2=1.d0-perpinv*bmod/p2
  if(alam2.lt.0.d0) then
    ierr=2
    return
  endif
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
  use field_eq_mod,  only : psif,psi_sep
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
  use field_eq_mod, only : psif
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
  use field_eq_mod, only : psif,dpsidr,dpsidz
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
  use field_eq_mod, only : psif,dpsidr,dpsidz,psi_axis,psi_sep,rtf
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
  use field_eq_mod, only : psif,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
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
  use field_eq_mod, only : psif,dpsidr,dpsidz,psi_axis,psi_sep
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
  subroutine find_bounds_fixpoints(ierr)
!
  use global_invariants,    only : toten,perpinv,sigma
  use find_all_roots_mod,   only : nroots,roots,relerr_allroots
  use poicut_mod,           only : npc,rpc_arr,zpc_arr,rmagaxis, &
                                   Rbou_lfs,Zbou_lfs,Rbou_hfs,Zbou_hfs
  use field_eq_mod,         only : psif
  use bounds_fixpoints_mod, only : nbounds,R_bo,Z_bo,psiast_bo, &
                                   nregions,all_regions
  use cc_mod, only : wrbounds
!
  implicit none
!
  logical, parameter :: doublecount=.true.
  integer, parameter :: niter=100
  double precision, parameter :: relerr=1d-12
  logical :: separin_b,separin_e,within_rhopol,start_minmax
  integer :: ierr,ireg0
  integer :: i,k,iter,isig,ireg,n_o,n_x,nxp_tot,ixp_tot,nsc
  double precision :: R,Z,vpar2,dvpar2_dR,dZ_dR,errdist2,R_b,R_e,dummy
  double precision :: R_b_in,R_e_in,p_phi,sigma_new,R_new,R_old,tau_fr,dphi_fr
  double precision :: gpgb,dgpgb_dr,dgpgb_dz,det,del_R,del_Z,delR_marg
  double precision :: pphi_b,pphi_e,pphi_min_reg,pphi_max_reg,pphi_min,pphi_max
  double precision :: pphi_minmax
  double precision, dimension(2) :: gradvpar2,Rblfs,Rbhfs,Zblfs,Zbhfs
  logical,          dimension(:),   allocatable :: opoint_extr,logdummy
  integer,          dimension(:,:), allocatable :: ipoi_x_tot
  integer,          dimension(:),   allocatable :: ipoi_tmp
  double precision, dimension(:),   allocatable :: R_extr,Z_extr,psiast_extr, &
                                                   psiast_x_tot,R_x_tot,rsc_tmp
!------------
!
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
  ierr=0
!
! find inner boundaries of allowed regions:
  relerr_allroots=1.d-12
!
  call find_all_roots(vparzero1D,rpc_arr(0),rpc_arr(npc),ierr)
!
  if(ierr.ne.0) then
    print *,'find_bounds_fixpoints: error in find_all_roots, inner boundaries'
  endif
!
! displace the roots by root accuracy distance so that vpar2 > 0 always:
  do i=1,nroots
!
    call vparzero1D(roots(i),vpar2,dvpar2_dR)
!
    roots(i)=roots(i)*(1.d0+sign(relerr,dvpar2_dR))
  enddo
!
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
      ierr=2
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
        call classify_extrema(R_b_in,R_e_in)
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
    print *,'minimum and maximum values of p_phi are not determined'
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
        call find_all_roots(boucross,R_b_in,R_e_in,ierr)
!
        if(ierr.ne.0) then
          print *,'find_bounds_fixpoints: error in find_all_roots, cut left boundary 1'
        endif
!
        if(nroots.gt.0) then
          within_rhopol=.true.
          R=roots(1)
!
          call get_poicut(R,Z,dZ_dR)
          call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
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
        call find_all_roots(boucross,R_b_in,R_e_in,ierr)
!
        if(ierr.ne.0) then
          print *,'find_bounds_fixpoints: error in find_all_roots, cut left boundary 2'
        endif
!
        if(nroots.gt.0) then
          within_rhopol=.true.
          R=roots(1)
!
          call get_poicut(R,Z,dZ_dR)
          call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
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
          call find_all_roots(boucross,R_b_in,R_e_in,ierr)
!
          if(ierr.ne.0) then
            print *,'find_bounds_fixpoints: error in find_all_roots, cut right boundary 1'
          endif
!
          if(nroots.gt.0) then
            R=roots(nroots)
!
            call get_poicut(R,Z,dZ_dR)
            call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
!
            all_regions(isig,ireg)%inner_e=.false.
            all_regions(isig,ireg)%R_e=R
            all_regions(isig,ireg)%Z_e=Z
            all_regions(isig,ireg)%psiast_e=p_phi
          else
! this shold never happen
            print *,'cutter error'
            stop !no further processing, debug the code
          endif
        elseif(all_regions(isig,ireg)%psiast_e.gt.pphi_max) then
          pphi_minmax=pphi_max
          relerr_allroots=1.d-11
!
          call find_all_roots(boucross,R_b_in,R_e_in,ierr)
!
          if(ierr.ne.0) then
            print *,'find_bounds_fixpoints: error in find_all_roots, cut right boundary 2'
          endif
!
          if(nroots.gt.0) then
            R=roots(nroots)
!
            call get_poicut(R,Z,dZ_dR)
            call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
!
            all_regions(isig,ireg)%inner_e=.false.
            all_regions(isig,ireg)%R_e=R
            all_regions(isig,ireg)%Z_e=Z
            all_regions(isig,ireg)%psiast_e=p_phi
          else
! this shold never happen
            print *,'cutter error'
            stop !no further processing, debug the code
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
      call classify_extrema(R_b_in,R_e_in)
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
                print *,'find_bounds_fixpoints: error in find_all_roots, sepcross 1'
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
  contains
!
!------------
!
  subroutine classify_extrema(R_b_in,R_e_in)
!
! Find all extremum points of psi^* on the Poincare cut in the interval R_b_in < R < R_e_in,
! determines their types (X- or O-point). Results - cylindrical coordinates (R,Z), values
! of psi^* and types of fixed points - are stored in the arrays R_extr, Z_extr, psiast_extr
! and opoint_extr (logical array with .true. for O-points and .false. for X-points), respectively.
! Dimension of the arrays is nroots.
!
  implicit none
!
  integer :: ierr,nroots1,nroots2,i
  double precision :: R_b_in,R_e_in
  double precision :: vvrt,dvvrt_dx,dvrdz,p_phi
  double precision, dimension(:), allocatable :: R_tmp1,R_tmp2,R_tmp
!
  R_e=0.5d0*(R_b_in+R_e_in)
!
  R_b=R_b_in
  relerr_allroots=1.d-8
!
  call find_all_roots(get_vvert,0.d0,1.d0,ierr)
!
  if(ierr.ne.0) then
    print *,'classify_extrema: error in find_all_roots 1'
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
    print *,'classify_extrema: error in find_all_roots 2'
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
  if(nroots.gt.0) then
    if(allocated(opoint_extr)) then
      deallocate(opoint_extr,R_extr,Z_extr,psiast_extr)
    endif
    allocate(opoint_extr(nroots))
    allocate(R_extr(nroots),Z_extr(nroots),psiast_extr(nroots))
  endif
!
  do i=1,nroots
!
    R=R_tmp(i)
!
    call get_poicut(R,Z,dZ_dR)
    call tormom_of_RZ(toten,perpinv,sigma,R,Z,p_phi,ierr)
    call determine_fixpoint_type(R,Z,opoint_extr(i))
!
    R_extr(i)=R
    Z_extr(i)=Z
    psiast_extr(i)=p_phi
!
  enddo
!
  if(nroots.gt.0) deallocate(R_tmp)
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
  subroutine determine_fixpoint_type(R_in,Z_in,opoint)
!
! Determines type of fixed point with coordinates (R_in,Z_in)
! opoint = .true.   - O-point
! opoint = .false.  - X-point
!
  implicit none
!
  double precision, parameter :: dtau=0.d0, hdiff=1.d-6
  logical :: opoint
  integer :: ierr,iden_R,iden_Z
  double precision :: R_in,Z_in,dvrdr,dvrdz,dvzdr,dvzdz
  double precision, dimension(5) :: z,vz,vz_base
!
  z(1)=R_in
  z(2)=0.d0
  z(3)=Z_in
!
  iden_R=0
  iden_Z=0
!
  call get_z45(toten,perpinv,sigma,z,ierr)
  call velo(dtau,z,vz_base)
!
  z(1)=R_in-hdiff*R_in
!
  call get_z45(toten,perpinv,sigma,z,ierr)
!
  if(ierr.eq.0) then
!
    call velo(dtau,z,vz)
!
    iden_R=iden_R+1
    dvrdr=vz(1)
    dvzdr=vz(3)
  else
    dvrdr=vz_base(1)
    dvzdr=vz_base(3)
  endif
!
  z(1)=R_in+hdiff*R_in
!
  call get_z45(toten,perpinv,sigma,z,ierr)
!
  if(ierr.eq.0) then
!
    call velo(dtau,z,vz)
!
    iden_R=iden_R+1
    dvrdr=(vz(1)-dvrdr)/dble(iden_R)
    dvzdr=(vz(3)-dvzdr)/dble(iden_R)
  else
    if(iden_R.eq.0) then
      print *,'determine_fixpoint_type: iden_R = 0'
      stop
    endif
    dvrdr=(vz_base(1)-dvrdr)/dble(iden_R)
    dvzdr=(vz_base(3)-dvzdr)/dble(iden_R)
  endif
!
  z(1)=R_in
  z(3)=Z_in-hdiff*R_in
!
  call get_z45(toten,perpinv,sigma,z,ierr)
!
  if(ierr.eq.0) then
!
    call velo(dtau,z,vz)
!
    iden_Z=iden_Z+1
    dvrdz=vz(1)
    dvzdz=vz(3)
  else
    dvrdz=vz_base(1)
    dvzdz=vz_base(3)
  endif
!
  z(3)=Z_in+hdiff*R_in
!
  call get_z45(toten,perpinv,sigma,z,ierr)
!
  if(ierr.eq.0) then
!
    call velo(dtau,z,vz)
!
    iden_Z=iden_Z+1
    dvrdz=(vz(1)-dvrdz)/dble(iden_Z)
    dvzdz=(vz(3)-dvzdz)/dble(iden_Z)
  else
    if(iden_Z.eq.0) then
      print *,'determine_fixpoint_type: iden_Z = 0'
      stop
    endif
    dvrdz=(vz_base(1)-dvrdz)/dble(iden_Z)
    dvzdz=(vz_base(3)-dvzdz)/dble(iden_Z)
  endif
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
    print *,'sepcross: error in starter_doublecount = ',ierr
    stop !no further processing, debug the code
  endif
!
  delpphi=delpphi-psiast_x_tot(ixp_tot)
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
    print *,'boucross: error in starter_doublecount = ',ierr
    stop !no further processing, debug the code
  endif
!
  delpphi=delpphi-pphi_minmax
!
  end subroutine boucross
!
!------------
!
  end subroutine find_bounds_fixpoints
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine form_classes_doublecount(classes_talk,ierr)
!
! Here "class" means segments coming from splitting by X-points of the domains of allowed motion
! for particles starting at the Poincare with a given parallel velocity sign. The routine
! determines these segments and their boundary types (usual boundary, vpar2=0 boundary
! or separatrix crossing). Types are used for determination of class parameterization.
!
  use form_classes_doublecount_mod, only : nclasses,ifuntype,sigma_class,  &
                                           R_class_beg,R_class_end
  use bounds_fixpoints_mod,         only : nregions,all_regions
  use poicut_mod,                   only : Rbou_lfs,Zbou_lfs,Rbou_hfs,Zbou_hfs
  use global_invariants,            only : toten,perpinv,sigma
!
  implicit none
!
  logical :: classes_talk
  integer          :: ierr,isig,ireg,nxp,ixp,icl,i
  integer, dimension(:), allocatable :: ibt_b,ibt_e
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
  double precision :: relmargin,widthclass,xbeg,xend
!
  select case(ifuntype)
  case(11)
! two rho_pol boundaries, 0<x<1
    xbeg=0.d0
    xend=1.d0
  case(12)
! left- rho_pol boundary, right - inner boundary, 0<x<1
    xbeg=0.d0
    xend=1.d0-relmargin
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
    xbeg=relmargin
    xend=1.d0
  case(22)
! two inner boundaries, 0<x<1
    xbeg=relmargin
    xend=1.d0-relmargin
  case(23)
! left- inner boundary, right - X-point, 0<x<inf
    xbeg=relmargin
    xend=-0.5d0*log(relmargin/widthclass)
  case(24)
! left- inner boundary, right - X-point, 0<x<inf
    xbeg=relmargin
    xend=-0.25d0*log(relmargin/widthclass)
  case(31)
! left- X-point, right - rho_pol boundary, -inf<x<0
    xbeg=0.5d0*log(relmargin/widthclass)
    xend=0.d0
  case(32)
! left- X-point, right - inner boundary, -inf<x<0
    xbeg=0.5d0*log(relmargin/widthclass)
    xend=-relmargin
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
    xend=-relmargin
  case(43)
! two X-points, -inf<x<inf
    xend=-0.5d0*log(relmargin/widthclass)
    xbeg=-xend*0.5d0
  case(44)
! two X-points, -inf<x<inf
    xbeg=0.25d0*log(relmargin/widthclass)
    xend=-xbeg
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
  use sample_matrix_mod, only : n1,n2,x,amat
  use orbit_dim_mod,     only : neqm,next
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
!
  sigma=sigma_class(iclass)
  delta_R=R_class_end(iclass)-R_class_beg(iclass)
!
  call xi_func(ifuntype(iclass),x,xi,dxi_dx)
!
  Rst=R_class_beg(iclass)+delta_R*xi
!
  call starter_doublecount(toten,perpinv,sigma,Rst,   &
                           psiast,dpsiast_dRst,z,ierr)
!
  if(ierr.ne.0) then
    print *,'get_matrix: error in starter'
    return
  endif
!
  fullbounce=.true.
!  fullbounce=.false.
!
  if(next.eq.0) then
    if(fullbounce) then
!
      call find_bounce(next,velo,dtau,z,taub,delphi,extraset)
!
    else
!
      call first_return_map(sigma,Rst,sigma,Rst,taub,delphi,ierr)
      call first_return_map(sigma,Rst,sigma,Rst,tau_fr,dphi_fr,ierr)
!
      taub=taub+tau_fr
      delphi=delphi+dphi_fr
    endif
  else
    extraset=0.d0
!
    call find_bounce(next,velo_pphint,dtau,z,taub,delphi,extraset)
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
                                npoi,xarr,amat_arr
  use orbit_dim_mod,     only : next,numbasef
  use get_matrix_mod,    only : relerror,relmargin,iclass
  use global_invariants, only : dtau,toten,perpinv,sigma
  use form_classes_doublecount_mod, only : ifuntype,R_class_beg,R_class_end
  use cc_mod, only : dowrite
!
  implicit none
!
  integer :: iunit,ierr,i
  double precision :: psiastbeg,psiastend,widthclass
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
  ierr=0
!
  widthclass=abs(R_class_end(iclass)/R_class_beg(iclass)-1.d0)
!
  if(widthclass/relmargin.lt.1.d0) then
    print *,'ignore class'
    ierr=1
    return
  endif
!
  call classbounds(ifuntype(iclass),relmargin,widthclass,xbeg,xend)
!
  eps=relerror
!
  call sample_matrix(get_matrix_doublecount,ierr)
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
  subroutine interpolate_class_doublecount(x,vec,dvec)
!
! Interpolation on the adaptive grid over class parameter x of normalized toroidal
! moment, normalizaed bounce time, toroidal displacement per bounce time and bounce
! integrals (vec) and their derivatives (dvec)
!
  use sample_matrix_mod, only : nlagr,n1,npoi,xarr,amat_arr
  use get_matrix_mod,    only : iclass
  use form_classes_doublecount_mod, only : ifuntype
!
  implicit none
!
  integer, parameter :: nder=1
  integer            :: k,npoilag,nshift,ibeg,iend
  double precision   :: x,xin,xi,xib,dxi_dx,dxi_dxb,dpphi_dxib,xi_inf,a,b
!
  double precision, dimension(n1)             :: vec,dvec
  double precision, dimension(0:nder,nlagr+1) :: coef
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
  end subroutine interpolate_class_doublecount
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
  fmaxw=dens*exp(phi_elec-toten)/temp**1.5d0
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
! Finds maximum possible value of the normalized perpendicular invariant
! on the Poincare cut for given total energy. Needed for the upper limit
! of integration over perpendicular invariant
!
  use poicut_mod,        only : Rbou_hfs,Rbou_lfs
!
  implicit none
!
  double precision :: perpinv_max
!
  call find_minmax_bsc(.false.,Jperponcut,Rbou_hfs,Rbou_lfs,perpinv_max)
!
!------------
  contains
!------------
!
  subroutine Jperponcut(R,perpinv)
!
! Computes maximum possible value of the normalized perpendicular invariant
! for given total energy and value of cut parameter R_c
!
  use global_invariants, only : toten
!
  implicit none
!
  double precision :: R,phi_elec,Z,dZ_dR,bmod,perpinv
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
  perpinv=(toten-phi_elec)/bmod
!
  end subroutine Jperponcut
!
!------------
!
  end subroutine find_Jperpmax
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
