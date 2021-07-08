module bmod_pert_mod
  logical :: prop = .true.  
  integer :: nrad,nzet,icp,iunit_bn=1265
  double precision :: hrad,hzet
  integer,          dimension(:),     allocatable  :: imi,ima,jmi,jma
  integer,          dimension(:,:),   allocatable  :: ipoint
  double precision, dimension(:),     allocatable  :: rad,zet
  double precision, dimension(:,:,:), allocatable  :: splbmod_re,splbmod_im
end module bmod_pert_mod
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine bmod_pert(R,Z,bmod_n)
!
! Toroidal Fourier amplitude of the perturbation of field module
!
  use bmod_pert_mod, only : prop,nrad,nzet,icp,iunit_bn,hrad,hzet,imi,ima,jmi,jma, &
                            ipoint,rad,zet,splbmod_re,splbmod_im
!
  implicit none
!
  integer :: i,ierr
  double precision :: R,Z,rrr,zzz
  double precision :: bmod_re,bmod_im,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2
  double complex   :: bmod_n
  double precision, dimension(:,:), allocatable :: bmod_n_re,bmod_n_im
!
  if(prop) then
    prop=.false.
!
    open(iunit_bn,form='unformatted',file='bmod_n.dat')
    read(iunit_bn) nrad,nzet
    allocate(rad(nrad),zet(nzet),bmod_n_re(nrad,nzet),bmod_n_im(nrad,nzet))
    read(iunit_bn) rad,zet
    read(iunit_bn) bmod_n_re,bmod_n_im
    close(iunit_bn)
    print *,'bmod_n.dat read'
!
    hrad = rad(2) - rad(1)
    hzet = zet(2) - zet(1)
!
! rectangular domain:
    allocate( imi(nzet),ima(nzet),jmi(nrad),jma(nrad) )
    imi = 1
    ima = nrad
    jmi = 1
    jma = nzet
!
!  Computation of the number of data in spline
    icp = 0
    do i=1,nzet
      if ( imi(i) .gt. 0 .and. ima(i) .gt. 0 ) then
         icp = icp + ima(i) - imi(i) + 1
      endif
    enddo
!
    allocate( splbmod_re(6,6,icp), splbmod_im(6,6,icp), ipoint(nrad,nzet) )
!
    call s2dcut(nrad,nzet,hrad,hzet,bmod_n_re,imi,ima,jmi,jma,icp,splbmod_re,ipoint)
    call s2dcut(nrad,nzet,hrad,hzet,bmod_n_im,imi,ima,jmi,jma,icp,splbmod_im,ipoint)
!
    deallocate(bmod_n_re,bmod_n_im)
  endif
!
  rrr=max(rad(1),min(rad(nrad),R))
  zzz=max(zet(1),min(zet(nzet),Z))
!
  call spline(nrad,nzet,rad,zet,hrad,hzet,icp,splbmod_re,ipoint,rrr,zzz, &
              bmod_re,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
  call spline(nrad,nzet,rad,zet,hrad,hzet,icp,splbmod_im,ipoint,rrr,zzz, &
              bmod_im,dpsidr,dpsidz,d2psidr2,d2psidrdz,d2psidz2,ierr)
!
  bmod_n=cmplx(bmod_re,bmod_im)
!
  end subroutine bmod_pert
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine test_bmodpert
!
  implicit none
!
  integer :: i,j
  double precision :: R,Z,rrr,zzz
  double complex   :: bmod_n
  double precision, dimension(-150:150) :: bre,bim
!
  do i=0,300
    R=108.d0+116.d0*dble(i)/300.d0
    do j=-150,150
      Z=-7.3d0+194.d0*dble(j)/300.d0
      call bmod_pert(R,Z,bmod_n)
      bre(j)=dble(bmod_n)
      bim(j)=dimag(bmod_n)
    enddo
    write(7001,*) bre
    write(7002,*) bim
  enddo
!
  stop
  end subroutine test_bmodpert
