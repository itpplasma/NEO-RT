  module sample_matrix_mod
    integer :: nlagr,n1,n2,npoi,itermax,nstiff,i_int
! >0 marks a matrix row whose value 0 flags an invalid (non-closing orbit) node.
! sample_matrix then refuses to refine across such nodes (see sample_matrix.f90).
    integer :: isentinel = 0
    double precision :: x,xbeg,xend,eps
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
  end module sample_matrix_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE sample_matrix(get_matrix,ierr)
!
  USE sample_matrix_mod
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: nder=0
  DOUBLE PRECISION, PARAMETER :: symm_break=0.01d0
! Floor each row's refinement scale at this fraction of its global maximum.  The
! adaptive split test compares the local interpolation error to eps times the
! row's scale.  Without the floor that scale is the local stencil maximum, so a
! row spanning many orders of magnitude (the bounce-integral weights, rows 4+)
! gets its negligible small-value tail refined to full relative accuracy -- the
! cause of the grid explosion (96k of 134k nodes in one #30835 class).  Flooring
! at refine_floor*global_max leaves moderate-range rows (psiast and the resonance
! frequencies, rows 1-3, whose global~local) untouched, so resonance positions
! are unchanged, while the wide-range weights stop refining where they are small.
  DOUBLE PRECISION, PARAMETER :: refine_floor=1.d-2
  INTEGER :: i,j,iter,npoi_old,iold,inew,ibeg,iend,nshift,npoilag,ierr
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: isplit
!
  DOUBLE PRECISION :: h,hh
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: xold
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: coef,amat_maxmod,amat_glob
  COMPLEX(8),   DIMENSION(:,:),   ALLOCATABLE :: amat1,amat2
  COMPLEX(8),   DIMENSION(:,:,:), ALLOCATABLE :: amat_old
!
  external :: get_matrix
!
  ierr=0
!
  npoilag=nlagr+1
  nshift=nlagr/2
  npoi=npoilag+1
!
  ALLOCATE(coef(0:nder,npoilag))
!
  h=(xend-xbeg)/npoilag/(1.d0+symm_break)
  hh=symm_break*h/npoilag
!
! amat is the shared sample_matrix_mod array; callers (e.g. the resonance class
! preclip) may have allocated it alone, so deallocate each array on its own
! status rather than gating all three on amat being allocated.
  if(allocated(amat))     DEALLOCATE(amat)
  if(allocated(xarr))     DEALLOCATE(xarr)
  if(allocated(amat_arr)) DEALLOCATE(amat_arr)
!
  ALLOCATE(amat(n1,n2),xarr(npoi),amat_arr(n1,n2,npoi))
!
  x=xbeg
  CALL get_matrix
!
  xarr(1)=x
  amat_arr(:,:,1)=amat
!
  x=xend
  CALL get_matrix
  xarr(npoi)=x
  amat_arr(:,:,npoi)=amat
!
  DO i=2,npoi-1
    x=xbeg+h*(i-1)+hh*(i-1)**2
    xarr(i)=x
  ENDDO
  DO i=2,npoi-1
    x=xarr(i)
    CALL get_matrix
    amat_arr(:,:,i)=amat
  ENDDO
!
  ALLOCATE(amat1(n1,n2),amat2(n1,n2),amat_maxmod(n1,n2),amat_glob(n1,n2))
!
! first check which intervals should be splitted
  ALLOCATE(isplit(npoi))
  isplit=0
  DO i=1,n1
    DO j=1,n2
      amat_glob(i,j)=MAXVAL(ABS(amat_arr(i,j,1:npoi)))
    ENDDO
  ENDDO
  DO inew=1,npoi-1
! Do not refine across a sentinel node (row isentinel = 0): the 0<->finite step
! is a find_bounce X-point-grazer skip, not a feature to resolve, and refining it
! explodes the grid.  sample_class_doublecount compacts these nodes out before
! the resonance search, so they never reach the interpolation/root finder.
    if(isentinel.gt.0) then
      if(amat_arr(isentinel,1,inew).eq.0.d0 .or. &
         amat_arr(isentinel,1,inew+1).eq.0.d0) cycle
    endif
    x=0.5d0*(xarr(inew)+xarr(inew+1))
    ibeg=MAX(1,MIN(npoi-nlagr-1,inew-nshift-1))
    iend=ibeg+nlagr
    CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
    DO i=1,n1
      amat1(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
      DO j=1,n2
        amat_maxmod(i,j)=MAXVAL(ABS(amat_arr(i,j,ibeg:iend)))
      ENDDO
    ENDDO
    ibeg=MAX(2,MIN(npoi-nlagr,inew-nshift+1))
    iend=ibeg+nlagr
    CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
    DO i=1,n1
      amat2(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
      DO j=1,n2
        amat_maxmod(i,j)=MAX(amat_maxmod(i,j),ABS(amat_arr(i,j,iend)))
      ENDDO
    ENDDO
    DO i=1,n1
      DO j=1,n2
        IF(ABS(amat1(i,j)-amat2(i,j)).GT.eps*max(amat_maxmod(i,j),refine_floor*amat_glob(i,j))) isplit(inew)=1
      ENDDO
    ENDDO
  ENDDO
  IF(MAXVAL(isplit).GT.0) THEN
    npoi_old=npoi
    ALLOCATE(xold(npoi),amat_old(n1,n2,npoi))
    xold=xarr
    amat_old=amat_arr
  ELSE
    RETURN
  ENDIF
!
  iter=0
  DO
    iter=iter+1
    IF(iter.GT.itermax) THEN
      ierr=2
      PRINT *,'sample_matrix : maximum number of iterations exceeded'
      RETURN
    ENDIF
!
! determine the dimension of new arrays and re-allocate them:
    DO iold=1,npoi_old-1
      IF(isplit(iold).EQ.1) npoi=npoi+1
    ENDDO
    IF(ALLOCATED(xarr)) THEN
      DEALLOCATE(xarr,amat_arr)
    ENDIF
    ALLOCATE(xarr(npoi),amat_arr(n1,n2,npoi))
!
! fill new arrays:
    inew=0
    DO iold=1,npoi_old-1
      inew=inew+1
      xarr(inew)=xold(iold)
      amat_arr(:,:,inew)=amat_old(:,:,iold)
      IF(isplit(iold).EQ.1) THEN
        inew=inew+1
        x=0.5d0*(xold(iold)+xold(iold+1))
        CALL get_matrix
        xarr(inew)=x
        amat_arr(:,:,inew)=amat
      ENDIF
    ENDDO
    inew=inew+1
    xarr(inew)=xold(npoi_old)
    amat_arr(:,:,inew)=amat_old(:,:,npoi_old)
    DEALLOCATE(isplit)
!
! check which intervals should be splitted
    ALLOCATE(isplit(npoi))
    isplit=0
    DO i=1,n1
      DO j=1,n2
        amat_glob(i,j)=MAXVAL(ABS(amat_arr(i,j,1:npoi)))
      ENDDO
    ENDDO
    DO inew=1,npoi-1
      if(isentinel.gt.0) then
        if(amat_arr(isentinel,1,inew).eq.0.d0 .or. &
           amat_arr(isentinel,1,inew+1).eq.0.d0) cycle
      endif
      x=0.5d0*(xarr(inew)+xarr(inew+1))
      ibeg=MAX(1,MIN(npoi-nlagr-1,inew-nshift-1))
      iend=ibeg+nlagr
      CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
      DO i=1,n1
        amat1(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
        DO j=1,n2
          amat_maxmod(i,j)=MAXVAL(ABS(amat_arr(i,j,ibeg:iend)))
        ENDDO
      ENDDO
      ibeg=MAX(2,MIN(npoi-nlagr,inew-nshift+1))
      iend=ibeg+nlagr
      CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
      DO i=1,n1
        amat2(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
        DO j=1,n2
          amat_maxmod(i,j)=MAX(amat_maxmod(i,j),ABS(amat_arr(i,j,iend)))
        ENDDO
      ENDDO
      DO i=1,n1
        DO j=1,n2
          IF(ABS(amat1(i,j)-amat2(i,j)).GT.eps*max(amat_maxmod(i,j),refine_floor*amat_glob(i,j))) isplit(inew)=1
        ENDDO
      ENDDO
    ENDDO
    IF(MAXVAL(isplit).GT.0) THEN
      npoi_old=npoi
      DEALLOCATE(xold,amat_old)
      ALLOCATE(xold(npoi),amat_old(n1,n2,npoi))
      xold=xarr
      amat_old=amat_arr
    ELSE
      EXIT
    ENDIF
  ENDDO
!
  END SUBROUTINE sample_matrix
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
