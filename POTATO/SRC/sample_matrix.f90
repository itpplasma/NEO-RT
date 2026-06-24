!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  SUBROUTINE sample_matrix(get_matrix,ierr)
!
  USE sample_matrix_mod
  USE form_classes_doublecount_mod, only : ifuntype,R_class_beg,R_class_end,sigma_class
! next is threadprivate (the resonance mode loop writes it); copyin it into the
! grid-build regions, which never write it, so each worker starts from the master
! value.
  USE orbit_dim_mod, only : next
!
  IMPLICIT NONE
!
  INTEGER, PARAMETER :: nder=0
  DOUBLE PRECISION, PARAMETER :: symm_break=0.01d0
  INTEGER :: i,j,iter,npoi_old,iold,inew,ibeg,iend,nshift,npoilag,ierr,nnew,k
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: isplit,newslots
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
  npoi=npoilag+1
!
  ALLOCATE(coef(0:nder,npoilag))
!
  h=(xend-xbeg)/npoilag/(1.d0+symm_break)
  hh=symm_break*h/npoilag
!
  if(allocated(amat)) then
    DEALLOCATE(amat,xarr,amat_arr)
  endif
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
! Initial interior fill: one independent get_matrix (starter+find_bounce) per
! node, writing the disjoint slot amat_arr(:,:,i).  This is the find_bounce-heavy
! cost.  x and amat are threadprivate scratch; the form-class bounds are already
! threadprivate, so copyin seeds each worker with the master values.  toten,perpinv
! stay SHARED and read-only here (the grid build never calls pertham, so nothing
! writes them); next is threadprivate but unwritten here, so copyin seeds each
! worker with the master value to keep the prior shared behavior.
  !$omp parallel do default(shared) private(i) schedule(dynamic) &
  !$omp   copyin(ifuntype,R_class_beg,R_class_end,sigma_class,next)
  DO i=2,npoi-1
    if(.not.allocated(amat)) allocate(amat(n1,n2))
    x=xarr(i)
    CALL get_matrix
    amat_arr(:,:,i)=amat
  ENDDO
  !$omp end parallel do
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
! Serial layout pass: copy the retained old nodes into their new slots and, for
! each split interval, reserve the new node's slot inew with its midpoint xarr.
! Record the new-node slots in newslots so the find_bounce-heavy get_matrix calls
! run in the parallel fill below instead of here.  Layout (slot assignment, array
! re-allocation, the isplit decision) stays serial.
    ALLOCATE(newslots(npoi))
    nnew=0
    inew=0
    DO iold=1,npoi_old-1
      inew=inew+1
      xarr(inew)=xold(iold)
      amat_arr(:,:,inew)=amat_old(:,:,iold)
      IF(isplit(iold).EQ.1) THEN
        inew=inew+1
        xarr(inew)=0.5d0*(xold(iold)+xold(iold+1))
        nnew=nnew+1
        newslots(nnew)=inew
      ENDIF
    ENDDO
    inew=inew+1
    xarr(inew)=xold(npoi_old)
    amat_arr(:,:,inew)=amat_old(:,:,npoi_old)
    DEALLOCATE(isplit)
!
! Parallel fill pass: one independent get_matrix per new node, writing the
! disjoint slot amat_arr(:,:,newslots(k)).  Same per-thread scratch (x,amat) and
! copyin of the shared form-class bounds and next as the initial fill.
    !$omp parallel do default(shared) private(k) schedule(dynamic) &
    !$omp   copyin(ifuntype,R_class_beg,R_class_end,sigma_class,next)
    DO k=1,nnew
      if(.not.allocated(amat)) allocate(amat(n1,n2))
      x=xarr(newslots(k))
      CALL get_matrix
      amat_arr(:,:,newslots(k))=amat
    ENDDO
    !$omp end parallel do
    DEALLOCATE(newslots)
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
