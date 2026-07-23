  module sample_matrix_out_mod
    integer :: nlagr,n1,n2,npoi=0,itermax=0,icount
    double precision :: x,xbeg,xend,eps
    integer,          dimension(:),     allocatable :: ind_hist
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
! The adaptive J_perp grid is per-energy-slice scratch.  Energy slices may run
! concurrently, so each worker keeps its own grid and interpolation workspace.
    !$omp threadprivate(nlagr,n1,n2,npoi,itermax,icount,x,xbeg,xend,eps, &
    !$omp               ind_hist,xarr,amat,amat_arr)
  end module sample_matrix_out_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE sample_matrix_out(get_matrix,ierr)
!
! Generates adaptive grid with values of matrix function on this grid.
! Matrix function a_ij(x) is computed by external function "get_matrix".
! For grid refinement, Lagrange polynomial interpolations at two shifted stencils
! are compared in the middle of the examined grid interval. Their difference is
! integrated over the grid and compared with the absolute integral scale. This
! integral-error criterion can converge physical resonance-onset cusps that are
! not pointwise-polynomial; intervals carrying more than an equal share of the
! requested error budget are split.
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
  USE matrix_callback_status_mod, ONLY : matrix_callback_error, &
      matrix_callback_ok, reset_matrix_callback_error
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
  DOUBLE PRECISION, DIMENSION(:,:),   ALLOCATABLE :: coef,integral_scale
  DOUBLE PRECISION, DIMENSION(:,:,:), ALLOCATABLE :: interval_error
  COMPLEX(8),   DIMENSION(:,:),   ALLOCATABLE :: amat1,amat2
  COMPLEX(8),   DIMENSION(:,:,:), ALLOCATABLE :: amat_old
!
  external :: get_matrix
!
  ierr=0
  call reset_matrix_callback_error
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
  if(matrix_callback_error.ne.matrix_callback_ok) then
    ierr=matrix_callback_error
    return
  endif
  xarr(1)=x
  amat_arr(:,:,1)=amat
  ind_hist(1)=icount
!
  x=xend
  CALL get_matrix
  if(matrix_callback_error.ne.matrix_callback_ok) then
    ierr=matrix_callback_error
    return
  endif
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
    if(matrix_callback_error.ne.matrix_callback_ok) then
      ierr=matrix_callback_error
      return
    endif
    amat_arr(:,:,i)=amat
    ind_hist(i)=icount
  ENDDO
!
  ALLOCATE(amat1(n1,n2),amat2(n1,n2),integral_scale(n1,n2))
!
! first check which intervals should be splitted
  ALLOCATE(isplit(npoi))
  isplit=0
  ALLOCATE(interval_error(n1,n2,npoi-1))
  interval_error=0.d0
  integral_scale=0.d0
  DO inew=1,npoi-1
    DO i=1,n1
      DO j=1,n2
        integral_scale(i,j)=integral_scale(i,j) &
            +0.5d0*(ABS(amat_arr(i,j,inew))+ABS(amat_arr(i,j,inew+1))) &
            *(xarr(inew+1)-xarr(inew))
      ENDDO
    ENDDO
    x=0.5d0*(xarr(inew)+xarr(inew+1))
    ibeg=MAX(1,MIN(npoi-nlagr-1,inew-nshift-1))
    iend=ibeg+nlagr
    CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
    DO i=1,n1
      amat1(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
    ENDDO
    ibeg=MAX(2,MIN(npoi-nlagr,inew-nshift+1))
    iend=ibeg+nlagr
    CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
    DO i=1,n1
      amat2(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
      DO j=1,n2
        interval_error(i,j,inew)=ABS(amat1(i,j)-amat2(i,j)) &
            *(xarr(inew+1)-xarr(inew))
      ENDDO
    ENDDO
  ENDDO
  DO i=1,n1
    DO j=1,n2
      IF(SUM(interval_error(i,j,:)).GT.eps*integral_scale(i,j)) THEN
        DO inew=1,npoi-1
          IF(interval_error(i,j,inew).GE. &
              eps*integral_scale(i,j)/dble(npoi-1)) isplit(inew)=1
        ENDDO
      ENDIF
    ENDDO
  ENDDO
  DEALLOCATE(interval_error)
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
      PRINT *,'sample_matrix_out : npoi, unresolved intervals, eps = ', &
          npoi,COUNT(isplit.eq.1),eps
      PRINT *,'sample_matrix_out : unresolved x range = ', &
          MINVAL(xarr(1:npoi-1),MASK=isplit(1:npoi-1).eq.1), &
          MAXVAL(xarr(2:npoi),MASK=isplit(1:npoi-1).eq.1)
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
        if(matrix_callback_error.ne.matrix_callback_ok) then
          ierr=matrix_callback_error
          return
        endif
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
    ALLOCATE(interval_error(n1,n2,npoi-1))
    interval_error=0.d0
    integral_scale=0.d0
    DO inew=1,npoi-1
      DO i=1,n1
        DO j=1,n2
          integral_scale(i,j)=integral_scale(i,j) &
              +0.5d0*(ABS(amat_arr(i,j,inew))+ABS(amat_arr(i,j,inew+1))) &
              *(xarr(inew+1)-xarr(inew))
        ENDDO
      ENDDO
      x=0.5d0*(xarr(inew)+xarr(inew+1))
      ibeg=MAX(1,MIN(npoi-nlagr-1,inew-nshift-1))
      iend=ibeg+nlagr
      CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
      DO i=1,n1
        amat1(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
      ENDDO
      ibeg=MAX(2,MIN(npoi-nlagr,inew-nshift+1))
      iend=ibeg+nlagr
      CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
      DO i=1,n1
        amat2(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
        DO j=1,n2
          interval_error(i,j,inew)=ABS(amat1(i,j)-amat2(i,j)) &
              *(xarr(inew+1)-xarr(inew))
        ENDDO
      ENDDO
    ENDDO
    DO i=1,n1
      DO j=1,n2
        IF(SUM(interval_error(i,j,:)).GT.eps*integral_scale(i,j)) THEN
          DO inew=1,npoi-1
            IF(interval_error(i,j,inew).GE. &
                eps*integral_scale(i,j)/dble(npoi-1)) isplit(inew)=1
          ENDDO
        ENDIF
      ENDDO
    ENDDO
    DEALLOCATE(interval_error)
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
