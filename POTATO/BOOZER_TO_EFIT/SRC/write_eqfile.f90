  subroutine write_eqfile(nwEQD,nhEQD,psiAxis,psiSep,bt0,rzero,fpol,rad,zet,psiRZ)
!
  use mag_surf_module, only : iunit
!
  implicit none
!
  integer, intent(in) :: nwEQD, nhEQD
  real (kind=8), intent(in) :: bt0, rzero, psiAxis, psiSep
  real (kind=8), dimension(nwEQD), intent(in) :: fpol
  real (kind=8), intent(in) :: rad(nwEQD), zet(nhEQD)
  real (kind=8), dimension(nwEQD,nhEQD), intent(in) :: psiRZ

  integer :: gunit, idum
  character(len=10) :: dummy(6)
  integer :: i,j
  real (kind=8) :: xdim,zdim,r1,zmid,rmaxis,zmaxis,xdum
  real (kind=8) :: plas_cur
  real (kind=8), dimension(nwEQD) :: pres,ffprim,pprime,qpsi

  integer :: n_bndyxy,nlimEQD
  real (kind=8), dimension(:), allocatable :: LCFS, limEQD

  gunit=iunit
  dummy='          '
  idum=0
  xdim=rad(nwEQD)-rad(1)
  zdim=zet(nhEQD)-zet(1)
  r1=rad(1)
  zmid=0.5d0*(zet(nhEQD)+zet(1))
  rmaxis=rzero
  zmaxis=0.d0
  plas_cur=0.d0
  xdum=0.d0
  pres=0.d0
  ffprim=0.d0
  pprime=0.d0
  qpsi=0.d0

  open(unit=gunit,file='efit.eqdsk')

! Equilibrium Parameters
  write(gunit,2000) (dummy(i),i=1,6),idum,nwEQD,nhEQD
  write(gunit,2010) xdim,zdim,rzero,r1,zmid
  write(gunit,2010) rmaxis,zmaxis,psiAxis,psiSep,bt0
  write(gunit,2010) plas_cur,psiAxis,xdum,rmaxis,xdum
  write(gunit,2010) zmaxis,xdum,psiSep,xdum,xdum
  write(gunit,2010) (fpol(i),i=1,nwEQD)
  write(gunit,2010) (pres(i),i=1,nwEQD)
  write(gunit,2010) (ffprim(i),i=1,nwEQD)
  write(gunit,2010) (pprime(i),i=1,nwEQD)
  write(gunit,2010) ((psiRZ(i,j),i=1,nwEQD),j=1,nhEQD)
  write(gunit,2010) (qpsi(i),i=1,nwEQD)
!
  n_bndyxy=1
  nlimEQD=1
  write(gunit,'(2i5)') n_bndyxy,nlimEQD
  write(gunit,2010) xdum,xdum
  write(gunit,2010) xdum,xdum
!
  close(gunit)

  return

2000  format(6a8,3i4)
2010  format(5(e16.9))
!
  end subroutine write_eqfile
