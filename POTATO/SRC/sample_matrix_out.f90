  module sample_matrix_out_mod
    integer :: nlagr,n1,n2,npoi=0,itermax=0,icount
    double precision :: x,xbeg,xend,eps
    integer,          dimension(:),     allocatable :: ind_hist
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
  end module sample_matrix_out_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE sample_matrix_out(get_matrix,ierr)
!
! Generates adaptive grid with values of matrix function on this grid.
! Matrix function a_ij(x) is computed by external function "get_matrix".
! For grid refinement, Lagrange polynomial interpolations at two shifted stencils
! are compared in the middle of the examined grid interval. If the difference is above
! given tolerance such an interval is split in two.
!
!  Formal input/output:
! get_matrix (in)  - external subroutine for computation of matrix function
! ierr       (out) - error code: 0 - normal work, 2 - maximum number of the refinements cylces is exceeded
!
!  Module input/output (via module "sample_matrix_out_mod"):
! nlagr                (in)    - order of Lagrange polynomial for sampling
! n1,n2                (in)    - a_ij(x) matrix size
! npoi                 (inout) - grid size. On input it is the initital equidistant grid size,
!                                on output - refined, non-equidistant grid size
! itermax              (in)    - maximum number of refinement cycles (at each cycle all intervals are checked)
! icount               (inout) - historic counter of grid points
! ind_hist             (out)   - mapper of points from the historic sequence to the increasing sequence in x
! x                            - argument of matrix function a_ij(x) (input for called routine "get_matrix")
! amat(n1,n2)                  - values of matrix function a_ij(x) (output of called routine "get_matrix")
! xarr(npoi)           (out)   - refined argument grid
! amat_arr(n1,n2,npoi) (out)   - matrix function on the refined grid
!
  USE sample_matrix_out_mod
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: nder=0
  DOUBLE PRECISION, PARAMETER :: symm_break=0.01d0
  INTEGER :: i,j,iter,npoi_old,iold,inew,ibeg,iend,nshift,npoilag,ierr
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: isplit,ind_hist_old
!
  DOUBLE PRECISION :: h,hh
  DOUBLE PRECISION, DIMENSION(:),     ALLOCATABLE :: xold
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: coef,amat_maxmod
  COMPLEX(8),   DIMENSION(:,:),   ALLOCATABLE :: amat1,amat2
  COMPLEX(8),   DIMENSION(:,:,:), ALLOCATABLE :: amat_old
!
  external :: get_matrix
!
  ierr=0
!
  npoilag=nlagr+1
  nshift=nlagr/2
  npoi=max(npoi,npoilag+1)
!
  ALLOCATE(coef(0:nder,npoilag))
!
  h=(xend-xbeg)/dble(npoi-1)/(1.d0+symm_break)
  hh=symm_break*h/dble(npoi-1)
!
  if(allocated(amat)) then
    DEALLOCATE(amat,xarr,amat_arr,ind_hist)
  endif
!
  ALLOCATE(amat(n1,n2),xarr(npoi),amat_arr(n1,n2,npoi),ind_hist(npoi))
  icount=0
!
  x=xbeg
  CALL get_matrix
  xarr(1)=x
  amat_arr(:,:,1)=amat
  ind_hist(1)=icount
!
  x=xend
  CALL get_matrix
  xarr(npoi)=x
  amat_arr(:,:,npoi)=amat
  ind_hist(npoi)=icount
!
  DO i=2,npoi-1
    x=xbeg+h*(i-1)+hh*(i-1)**2
    xarr(i)=x
  ENDDO
  DO i=2,npoi-1
    x=xarr(i)
    CALL get_matrix
    amat_arr(:,:,i)=amat
    ind_hist(i)=icount
  ENDDO
!
  ALLOCATE(amat1(n1,n2),amat2(n1,n2),amat_maxmod(n1,n2))
!
! first check which intervals should be splitted
  ALLOCATE(isplit(npoi))
  isplit=0
  DO inew=1,npoi-1
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
        IF(ABS(amat1(i,j)-amat2(i,j)).GT.eps*amat_maxmod(i,j)) isplit(inew)=1
      ENDDO
    ENDDO
  ENDDO
  IF(MAXVAL(isplit).GT.0) THEN
    npoi_old=npoi
    ALLOCATE(xold(npoi),amat_old(n1,n2,npoi),ind_hist_old(npoi))
    xold=xarr
    amat_old=amat_arr
    ind_hist_old=ind_hist
  ELSE
    RETURN
  ENDIF
!
  iter=0
  DO
    iter=iter+1
    IF(iter.GT.itermax) THEN
      ierr=2
      PRINT *,'sample_matrix_out : maximum number of iterations exceeded'
      RETURN
    ENDIF
!
! determine the dimension of new arrays and re-allocate them:
    DO iold=1,npoi_old-1
      IF(isplit(iold).EQ.1) npoi=npoi+1
    ENDDO
    IF(ALLOCATED(xarr)) THEN
      DEALLOCATE(xarr,amat_arr,ind_hist)
    ENDIF
    ALLOCATE(xarr(npoi),amat_arr(n1,n2,npoi),ind_hist(npoi))
!
! fill new arrays:
    inew=0
    DO iold=1,npoi_old-1
      inew=inew+1
      xarr(inew)=xold(iold)
      amat_arr(:,:,inew)=amat_old(:,:,iold)
      ind_hist(inew)=ind_hist_old(iold)
      IF(isplit(iold).EQ.1) THEN
        inew=inew+1
        x=0.5d0*(xold(iold)+xold(iold+1))
        CALL get_matrix
        xarr(inew)=x
        amat_arr(:,:,inew)=amat
        ind_hist(inew)=icount
      ENDIF
    ENDDO
    inew=inew+1
    xarr(inew)=xold(npoi_old)
    amat_arr(:,:,inew)=amat_old(:,:,npoi_old)
    ind_hist(inew)=ind_hist_old(npoi_old)
    DEALLOCATE(isplit)
!
! check which intervals should be splitted
    ALLOCATE(isplit(npoi))
    isplit=0
    DO inew=1,npoi-1
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
          IF(ABS(amat1(i,j)-amat2(i,j)).GT.eps*amat_maxmod(i,j)) isplit(inew)=1
        ENDDO
      ENDDO
    ENDDO
    IF(MAXVAL(isplit).GT.0) THEN
      npoi_old=npoi
      DEALLOCATE(xold,amat_old,ind_hist_old)
      ALLOCATE(xold(npoi),amat_old(n1,n2,npoi),ind_hist_old(npoi))
      xold=xarr
      amat_old=amat_arr
      ind_hist_old=ind_hist
    ELSE
      EXIT
    ENDIF
  ENDDO
!
  END SUBROUTINE sample_matrix_out
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
