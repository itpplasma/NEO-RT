module bmod_pert_mod
  integer, parameter :: bmod_pert_success=0
  integer, parameter :: bmod_pert_outside_grid=1
  integer, parameter :: bmod_pert_load_failure=2
  integer, parameter :: bmod_pert_spline_failure=3
  integer, parameter :: bmod_pert_nonfinite=4
  logical :: prop=.true.
  integer :: bmod_pert_load_status=bmod_pert_success
  integer :: nrad,nzet,icp,iunit_bn=1265
  double precision :: hrad,hzet
  integer,          dimension(:),     allocatable :: imi,ima,jmi,jma
  integer,          dimension(:,:),   allocatable :: ipoint
  double precision, dimension(:),     allocatable :: rad,zet
  double precision, dimension(:,:,:), allocatable :: splbmod_re,splbmod_im
end module bmod_pert_mod
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine bmod_pert(R,Z,bmod_n)
!
! Historical evaluator.  It intentionally clamps to the rectangular spline
! grid because a few legacy diagnostics depend on that behavior.  Production
! resonance code must call bmod_pert_status instead.
!
  use bmod_pert_mod, only : nrad,nzet,hrad,hzet,rad,zet,bmod_pert_success
  implicit none
!
  integer :: ierr
  double precision :: R,Z,rrr,zzz
  complex(8) :: bmod_n
!
  call ensure_bmod_pert_loaded(ierr)
  if(ierr.ne.bmod_pert_success) then
    bmod_n=cmplx(0.d0,0.d0)
    return
  endif
  rrr=max(rad(1),min(rad(nrad),R))
  zzz=max(zet(1),min(zet(nzet),Z))
  call evaluate_bmod_pert(rrr,zzz,bmod_n,ierr)
  end subroutine bmod_pert
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine bmod_pert_status(R,Z,bmod_n,ierr)
!
! Strict evaluator for full-FOW production.  Grid boundaries are valid; an
! out-of-grid point is an error and never receives an extrapolated/clamped
! edge amplitude.
!
  use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
  use bmod_pert_mod, only : nrad,nzet,rad,zet,bmod_pert_success, &
                            bmod_pert_outside_grid,bmod_pert_nonfinite
  implicit none
!
  double precision, intent(in) :: R,Z
  complex(8), intent(out) :: bmod_n
  integer, intent(out) :: ierr
  double precision :: rlo,rhi,zlo,zhi
!
  bmod_n=cmplx(0.d0,0.d0)
  call ensure_bmod_pert_loaded(ierr)
  if(ierr.ne.bmod_pert_success) return
  if(.not.ieee_is_finite(R) .or. .not.ieee_is_finite(Z)) then
    ierr=bmod_pert_outside_grid
    return
  endif
  rlo=min(rad(1),rad(nrad))
  rhi=max(rad(1),rad(nrad))
  zlo=min(zet(1),zet(nzet))
  zhi=max(zet(1),zet(nzet))
  if(R.lt.rlo .or. R.gt.rhi .or. Z.lt.zlo .or. Z.gt.zhi) then
    ierr=bmod_pert_outside_grid
    return
  endif
  call evaluate_bmod_pert(R,Z,bmod_n,ierr)
  if(ierr.eq.bmod_pert_success) then
    if(.not.ieee_is_finite(real(bmod_n)) .or. &
       .not.ieee_is_finite(aimag(bmod_n))) ierr=bmod_pert_nonfinite
  endif
  end subroutine bmod_pert_status
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine ensure_bmod_pert_loaded(ierr)
  use bmod_pert_mod, only : prop,nrad,nzet,icp,iunit_bn,hrad,hzet,imi,ima,jmi,jma, &
                            ipoint,rad,zet,splbmod_re,splbmod_im, &
                            bmod_pert_success,bmod_pert_load_status, &
                            bmod_pert_load_failure
  implicit none
  integer, intent(out) :: ierr
  integer :: i,ios
  double precision, dimension(:,:), allocatable :: bmod_n_re,bmod_n_im
!
  ierr=bmod_pert_load_status
  if(.not.prop) return
  !$omp critical (bmod_pert_init)
  if(prop) then
    ios=0
    open(iunit_bn,form='unformatted',file='bmod_n.dat',status='old',iostat=ios)
    if(ios.eq.0) read(iunit_bn,iostat=ios) nrad,nzet
    if(ios.eq.0 .and. nrad.ge.2 .and. nzet.ge.2) then
      allocate(rad(nrad),zet(nzet),bmod_n_re(nrad,nzet),bmod_n_im(nrad,nzet))
      read(iunit_bn,iostat=ios) rad,zet
      if(ios.eq.0) read(iunit_bn,iostat=ios) bmod_n_re,bmod_n_im
    else
      ios=1
    endif
    close(iunit_bn)
    if(ios.eq.0) then
      print *,'bmod_n.dat read'
      hrad=rad(2)-rad(1)
      hzet=zet(2)-zet(1)
      allocate(imi(nzet),ima(nzet),jmi(nrad),jma(nrad))
      imi=1
      ima=nrad
      jmi=1
      jma=nzet
      icp=0
      do i=1,nzet
        if(imi(i).gt.0 .and. ima(i).gt.0) icp=icp+ima(i)-imi(i)+1
      enddo
      allocate(splbmod_re(6,6,icp),splbmod_im(6,6,icp),ipoint(nrad,nzet))
      call s2dcut(nrad,nzet,hrad,hzet,bmod_n_re,imi,ima,jmi,jma,icp, &
                  splbmod_re,ipoint)
      call s2dcut(nrad,nzet,hrad,hzet,bmod_n_im,imi,ima,jmi,jma,icp, &
                  splbmod_im,ipoint)
      deallocate(bmod_n_re,bmod_n_im)
      bmod_pert_load_status=bmod_pert_success
      prop=.false.
    else
      bmod_pert_load_status=bmod_pert_load_failure
      if(allocated(bmod_n_re)) deallocate(bmod_n_re,bmod_n_im)
    endif
  endif
  !$omp end critical (bmod_pert_init)
  ierr=bmod_pert_load_status
  end subroutine ensure_bmod_pert_loaded
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine evaluate_bmod_pert(R,Z,bmod_n,ierr)
  use bmod_pert_mod, only : nrad,nzet,icp,hrad,hzet,rad,zet,ipoint, &
                            splbmod_re,splbmod_im,bmod_pert_success, &
                            bmod_pert_spline_failure
  implicit none
  double precision, intent(in) :: R,Z
  complex(8), intent(out) :: bmod_n
  integer, intent(out) :: ierr
  integer :: spline_ierr
  double precision :: bmod_re,bmod_im,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
!
  call spline(nrad,nzet,rad,zet,hrad,hzet,icp,splbmod_re,ipoint,R,Z, &
              bmod_re,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,spline_ierr)
  ierr=spline_ierr
  call spline(nrad,nzet,rad,zet,hrad,hzet,icp,splbmod_im,ipoint,R,Z, &
              bmod_im,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,spline_ierr)
  if(ierr.eq.0 .and. spline_ierr.ne.0) ierr=spline_ierr
  bmod_n=cmplx(bmod_re,bmod_im)
  if(ierr.ne.0) ierr=bmod_pert_spline_failure
  if(ierr.eq.0) ierr=bmod_pert_success
  end subroutine evaluate_bmod_pert
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine test_bmodpert
!
  implicit none
!
  integer :: i,j
  double precision :: R,Z
  complex(8) :: bmod_n
  double precision, dimension(-150:150) :: bre,bim
!
  do i=0,300
    R=108.d0+116.d0*dble(i)/300.d0
    do j=-150,150
      Z=-7.3d0+194.d0*dble(j)/300.d0
      call bmod_pert(R,Z,bmod_n)
      bre(j)=real(bmod_n)
      bim(j)=aimag(bmod_n)
    enddo
    write(7001,*) bre
    write(7002,*) bim
  enddo
!
  stop
  end subroutine test_bmodpert
