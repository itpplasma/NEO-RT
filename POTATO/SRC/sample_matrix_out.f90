  module sample_matrix_out_mod
    integer, parameter :: sample_matrix_success=0
    integer, parameter :: sample_matrix_nonconverged=2
    integer, parameter :: sample_matrix_topology_transition=3
    integer, parameter :: sample_matrix_topology_unresolved=4
    integer :: nlagr,n1,n2,npoi=0,itermax=0,icount
    double precision :: x,xbeg,xend,eps
    integer :: topology_signature=0,topology_error=0
    double precision :: topology_context_h=0.d0
    logical :: topology_probe_only=.false.
    integer,          dimension(:),     allocatable :: ind_hist
    integer,          dimension(:),     allocatable :: topology_arr
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
! The adaptive J_perp grid is per-energy-slice scratch.  Energy slices may run
! concurrently, so each worker keeps its own grid and interpolation workspace.
    !$omp threadprivate(nlagr,n1,n2,npoi,itermax,icount,x,xbeg,xend,eps, &
    !$omp               topology_signature,topology_error,topology_context_h, &
    !$omp               topology_probe_only, &
    !$omp               ind_hist,topology_arr,xarr,amat,amat_arr)

  contains

    pure logical function topology_stencil_is_compatible(signatures,ibeg,iend)
      integer, intent(in) :: signatures(:),ibeg,iend

      topology_stencil_is_compatible=.false.
      if(ibeg.lt.1) return
      if(iend.gt.size(signatures)) return
      if(ibeg.gt.iend) return
      if(signatures(ibeg).eq.0) return
      topology_stencil_is_compatible=all(signatures(ibeg:iend).eq.signatures(ibeg))
    end function topology_stencil_is_compatible

    integer function topology_signature_of_classes(nclasses,ifuntype,sigma_class)
      integer, intent(in) :: nclasses
      integer, intent(in) :: ifuntype(:)
      double precision, intent(in) :: sigma_class(:)

      integer :: i,sign_code

      topology_signature_of_classes=17+nclasses
      do i=1,nclasses
        sign_code=1
        if(sigma_class(i).lt.0.d0) sign_code=-1
        ! Keep the rolling value small before multiplication so the default
        ! integer kind cannot overflow.  The signature intentionally ignores
        ! continuously moving class bounds; it identifies the discrete class
        ! topology (count, endpoint types, and passing sign).
        topology_signature_of_classes=mod( &
            31*mod(topology_signature_of_classes,1000000)+7*ifuntype(i)+ &
            sign_code+13*i,1000000007)
      enddo
      if(topology_signature_of_classes.eq.0) topology_signature_of_classes=1
    end function topology_signature_of_classes

    subroutine evaluate_matrix_callback(get_matrix,ierr_eval)
      integer, intent(out) :: ierr_eval

      external :: get_matrix

      topology_signature=0
      topology_error=0
      call get_matrix
      if(topology_error.ne.0) then
        ierr_eval=topology_error
      elseif(topology_signature.eq.0) then
        ierr_eval=sample_matrix_topology_transition
      else
        ierr_eval=sample_matrix_success
      endif
      if(ierr_eval.ne.0) then
        print *,'sample_matrix_out : callback failed at H, J, topology, ierr = ', &
            topology_context_h,x,topology_signature,ierr_eval
      endif
    end subroutine evaluate_matrix_callback

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
! ierr       (out) - error code: 0 - normal work, 2 - maximum number of the refinements cycles is exceeded, 3 - topology transition
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
  INTEGER,          DIMENSION(:),     ALLOCATABLE :: isplit,ind_hist_old,topology_old
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
    DEALLOCATE(amat,xarr,amat_arr,ind_hist,topology_arr)
  endif
  !
  ALLOCATE(amat(n1,n2),xarr(npoi),amat_arr(n1,n2,npoi),ind_hist(npoi), &
      topology_arr(npoi))
  icount=0
  !
  x=xbeg
  CALL evaluate_matrix_callback(get_matrix,ierr)
  IF(ierr.NE.0) RETURN
  xarr(1)=x
  amat_arr(:,:,1)=amat
  ind_hist(1)=icount
  topology_arr(1)=topology_signature
  !
  x=xend
  CALL evaluate_matrix_callback(get_matrix,ierr)
  IF(ierr.NE.0) RETURN
  xarr(npoi)=x
  amat_arr(:,:,npoi)=amat
  ind_hist(npoi)=icount
  topology_arr(npoi)=topology_signature
!
  DO i=2,npoi-1
    x=xbeg+h*(i-1)+hh*(i-1)**2
    xarr(i)=x
  ENDDO
  DO i=2,npoi-1
    x=xarr(i)
    CALL evaluate_matrix_callback(get_matrix,ierr)
    IF(ierr.NE.0) RETURN
    amat_arr(:,:,i)=amat
    ind_hist(i)=icount
    topology_arr(i)=topology_signature
  ENDDO
!
  ALLOCATE(amat1(n1,n2),amat2(n1,n2),amat_maxmod(n1,n2))
!
! first check which intervals should be splitted
  ALLOCATE(isplit(npoi))
  isplit=0
  DO inew=1,npoi-1
    IF(topology_arr(inew).NE.topology_arr(inew+1)) THEN
      CALL report_topology_transition(inew,inew+1)
      ierr=sample_matrix_topology_transition
      RETURN
    ENDIF
    x=0.5d0*(xarr(inew)+xarr(inew+1))
    ibeg=MAX(1,MIN(npoi-nlagr-1,inew-nshift-1))
    iend=ibeg+nlagr
    IF(.NOT.topology_stencil_is_compatible(topology_arr,ibeg,iend)) THEN
      CALL report_topology_transition(ibeg,iend)
      ierr=sample_matrix_topology_transition
      RETURN
    ENDIF
    CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
    DO i=1,n1
      amat1(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
      DO j=1,n2
        amat_maxmod(i,j)=MAXVAL(ABS(amat_arr(i,j,ibeg:iend)))
      ENDDO
    ENDDO
    ibeg=MAX(2,MIN(npoi-nlagr,inew-nshift+1))
    iend=ibeg+nlagr
    IF(.NOT.topology_stencil_is_compatible(topology_arr,ibeg,iend)) THEN
      CALL report_topology_transition(ibeg,iend)
      ierr=sample_matrix_topology_transition
      RETURN
    ENDIF
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
    ALLOCATE(xold(npoi),amat_old(n1,n2,npoi),ind_hist_old(npoi),topology_old(npoi))
    xold=xarr
    amat_old=amat_arr
    ind_hist_old=ind_hist
    topology_old=topology_arr
  ELSE
    RETURN
  ENDIF
!
  iter=0
  DO
    iter=iter+1
    IF(iter.GT.itermax) THEN
      ierr=sample_matrix_nonconverged
      PRINT *,'sample_matrix_out : maximum number of iterations exceeded'
      PRINT *,'sample_matrix_out : H, J probe, J domain, npoi, ierr = ', &
          topology_context_h,x,xarr(1),xarr(npoi),npoi,ierr
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
    ALLOCATE(xarr(npoi),amat_arr(n1,n2,npoi),ind_hist(npoi),topology_arr(npoi))
!
! fill new arrays:
    inew=0
    DO iold=1,npoi_old-1
      inew=inew+1
      xarr(inew)=xold(iold)
      amat_arr(:,:,inew)=amat_old(:,:,iold)
      ind_hist(inew)=ind_hist_old(iold)
      topology_arr(inew)=topology_old(iold)
      IF(isplit(iold).EQ.1) THEN
        inew=inew+1
        x=0.5d0*(xold(iold)+xold(iold+1))
        CALL evaluate_matrix_callback(get_matrix,ierr)
        IF(ierr.NE.0) RETURN
        xarr(inew)=x
        amat_arr(:,:,inew)=amat
        ind_hist(inew)=icount
        topology_arr(inew)=topology_signature
      ENDIF
    ENDDO
    inew=inew+1
    xarr(inew)=xold(npoi_old)
    amat_arr(:,:,inew)=amat_old(:,:,npoi_old)
    ind_hist(inew)=ind_hist_old(npoi_old)
    topology_arr(inew)=topology_old(npoi_old)
    DEALLOCATE(isplit)
!
! check which intervals should be splitted
    ALLOCATE(isplit(npoi))
    isplit=0
    DO inew=1,npoi-1
      IF(topology_arr(inew).NE.topology_arr(inew+1)) THEN
        CALL report_topology_transition(inew,inew+1)
        ierr=sample_matrix_topology_transition
        RETURN
      ENDIF
      x=0.5d0*(xarr(inew)+xarr(inew+1))
      ibeg=MAX(1,MIN(npoi-nlagr-1,inew-nshift-1))
      iend=ibeg+nlagr
      IF(.NOT.topology_stencil_is_compatible(topology_arr,ibeg,iend)) THEN
        CALL report_topology_transition(ibeg,iend)
        ierr=sample_matrix_topology_transition
        RETURN
      ENDIF
      CALL plag_coeff(npoilag,nder,x,xarr(ibeg:iend),coef)
      DO i=1,n1
        amat1(i,:)=MATMUL(amat_arr(i,:,ibeg:iend),coef(0,:))
        DO j=1,n2
          amat_maxmod(i,j)=MAXVAL(ABS(amat_arr(i,j,ibeg:iend)))
        ENDDO
      ENDDO
      ibeg=MAX(2,MIN(npoi-nlagr,inew-nshift+1))
      iend=ibeg+nlagr
      IF(.NOT.topology_stencil_is_compatible(topology_arr,ibeg,iend)) THEN
        CALL report_topology_transition(ibeg,iend)
        ierr=sample_matrix_topology_transition
        RETURN
      ENDIF
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
      DEALLOCATE(xold,amat_old,ind_hist_old,topology_old)
      ALLOCATE(xold(npoi),amat_old(n1,n2,npoi),ind_hist_old(npoi),topology_old(npoi))
      xold=xarr
      amat_old=amat_arr
      ind_hist_old=ind_hist
      topology_old=topology_arr
    ELSE
      EXIT
    ENDIF
  ENDDO
  !
contains

  subroutine report_topology_transition(ibeg_report,iend_report)
    integer, intent(in) :: ibeg_report,iend_report

    PRINT *,'sample_matrix_out : topology transition blocks interpolation'
    PRINT *,'sample_matrix_out : H, Jlo, Jhi, topology_lo, topology_hi = ', &
        topology_context_h,xarr(ibeg_report),xarr(iend_report), &
        topology_arr(ibeg_report),topology_arr(iend_report)
  end subroutine report_topology_transition

  END SUBROUTINE sample_matrix_out
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc

  SUBROUTINE sample_matrix_out_partitioned(get_matrix,ierr)
!
! Discover discrete topology changes in a bounded J_perp interval, bracket each
! change with valid samples, and run sample_matrix_out independently on each
! resulting open subinterval.  The bracket is deliberately retained as a gap:
! no interpolation or trapezoid crosses an uncertified topology boundary.
!
  USE sample_matrix_out_mod
!
  IMPLICIT NONE
!
  INTEGER, INTENT(OUT) :: ierr
  INTEGER, PARAMETER :: topology_seed_points=31
  INTEGER, PARAMETER :: topology_max_probes=4096
  INTEGER, PARAMETER :: topology_max_endpoint_attempts=16
  DOUBLE PRECISION, PARAMETER :: topology_rel_tol=1.d-8
  DOUBLE PRECISION, PARAMETER :: topology_abs_tol=1.d-12
  INTEGER :: i,attempt,npoi_request,nprobe,position,nsegments,n_total
  INTEGER :: ierr_local,nseg_points,segment_start,segment_end
  INTEGER :: n1_request,n2_request
  DOUBLE PRECISION :: xbeg_full,xend_full,topology_tol,scale,step,xmid
  DOUBLE PRECISION, ALLOCATABLE :: probe_x(:),probe_x_new(:)
  INTEGER, ALLOCATABLE :: probe_sig(:),probe_sig_new(:)
  DOUBLE PRECISION, ALLOCATABLE :: segment_lo(:),segment_hi(:)
  INTEGER, ALLOCATABLE :: segment_sig(:)
  DOUBLE PRECISION, ALLOCATABLE :: xout(:),amatout(:,:,:),xseg(:),amatseg(:,:,:)
  INTEGER, ALLOCATABLE :: indout(:),sigout(:),indseg(:),sigseg(:)
  LOGICAL :: endpoint_contracted
!
  EXTERNAL :: get_matrix
!
  ierr=sample_matrix_success
  xbeg_full=xbeg
  xend_full=xend
  npoi_request=npoi
  n1_request=n1
  n2_request=n2
  endpoint_contracted=.false.
!
  IF(xend_full.LE.xbeg_full) THEN
    PRINT *,'sample_matrix_out_partitioned: invalid J_perp domain = ', &
        xbeg_full,xend_full
    ierr=sample_matrix_topology_unresolved
    RETURN
  ENDIF
  scale=MAX(1.d0,ABS(xbeg_full),ABS(xend_full))
  topology_tol=MAX(topology_abs_tol*scale,topology_rel_tol*ABS(xend_full-xbeg_full))
  npoi_request=MAX(npoi_request,nlagr+2)
  nprobe=MAX(npoi_request,topology_seed_points)
  topology_probe_only=.true.
  IF(.NOT.ALLOCATED(amat)) ALLOCATE(amat(n1_request,n2_request))
!
! Topology-only endpoint and seed probes.  Invalid global endpoints may be
! contracted inward; any invalid interior probe fails closed.
  ALLOCATE(probe_x(nprobe),probe_sig(nprobe))
  probe_x(1)=xbeg_full
  x=probe_x(1)
  CALL evaluate_matrix_callback(get_matrix,ierr_local)
  IF(ierr_local.NE.sample_matrix_success) THEN
    CALL contract_endpoint(1,xbeg_full,1.d0,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      ierr=sample_matrix_topology_unresolved
      topology_probe_only=.false.
      RETURN
    ENDIF
    endpoint_contracted=.true.
  ENDIF
  probe_sig(1)=topology_signature
!
  probe_x(nprobe)=xend_full
  x=probe_x(nprobe)
  CALL evaluate_matrix_callback(get_matrix,ierr_local)
  IF(ierr_local.NE.sample_matrix_success) THEN
    CALL contract_endpoint(nprobe,xend_full,-1.d0,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      ierr=sample_matrix_topology_unresolved
      topology_probe_only=.false.
      RETURN
    ENDIF
    endpoint_contracted=.true.
  ENDIF
  probe_sig(nprobe)=topology_signature
  IF(probe_x(nprobe).LE.probe_x(1)) THEN
    ierr=sample_matrix_topology_unresolved
    topology_probe_only=.false.
    RETURN
  ENDIF
!
  step=(probe_x(nprobe)-probe_x(1))/DBLE(nprobe-1)
  DO i=2,nprobe-1
    probe_x(i)=probe_x(1)+step*DBLE(i-1)
    x=probe_x(i)
    CALL evaluate_matrix_callback(get_matrix,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      PRINT *,'sample_matrix_out_partitioned: invalid interior probe H,J = ', &
          topology_context_h,x
      ierr=ierr_local
      topology_probe_only=.false.
      RETURN
    ENDIF
    probe_sig(i)=topology_signature
  ENDDO
!
! Walk the ordered probe list.  Equal signatures receive a midpoint probe so
! an interior component is not silently dropped.  Unequal signatures are
! bisected until the valid left/right bracket is no wider than topology_tol.
  topology_probe_only=.true.
  position=1
  DO WHILE(position.LT.nprobe)
    IF(nprobe.GE.topology_max_probes) THEN
      PRINT *,'sample_matrix_out_partitioned: topology probe budget exceeded'
      ierr=sample_matrix_topology_unresolved
      topology_probe_only=.false.
      RETURN
    ENDIF
    IF(probe_sig(position).EQ.probe_sig(position+1)) THEN
      xmid=0.5d0*(probe_x(position)+probe_x(position+1))
      IF(xmid.LE.probe_x(position)) THEN
        ierr=sample_matrix_topology_unresolved
        topology_probe_only=.false.
        RETURN
      ENDIF
      IF(xmid.GE.probe_x(position+1)) THEN
        ierr=sample_matrix_topology_unresolved
        topology_probe_only=.false.
        RETURN
      ENDIF
      x=xmid
      CALL evaluate_matrix_callback(get_matrix,ierr_local)
      IF(ierr_local.NE.sample_matrix_success) THEN
        PRINT *,'sample_matrix_out_partitioned: invalid interior midpoint H,J = ', &
            topology_context_h,x
        ierr=ierr_local
        topology_probe_only=.false.
        RETURN
      ENDIF
      IF(topology_signature.EQ.probe_sig(position)) THEN
        position=position+1
      ELSE
        CALL insert_probe(position+1,xmid,topology_signature)
      ENDIF
    ELSE
      DO WHILE(probe_x(position+1)-probe_x(position).GT.topology_tol)
        xmid=0.5d0*(probe_x(position)+probe_x(position+1))
        IF(xmid.LE.probe_x(position)) THEN
          ierr=sample_matrix_topology_unresolved
          topology_probe_only=.false.
          RETURN
        ENDIF
        IF(xmid.GE.probe_x(position+1)) THEN
          ierr=sample_matrix_topology_unresolved
          topology_probe_only=.false.
          RETURN
        ENDIF
        x=xmid
        CALL evaluate_matrix_callback(get_matrix,ierr_local)
        IF(ierr_local.NE.sample_matrix_success) THEN
          PRINT *,'sample_matrix_out_partitioned: invalid transition probe H,J = ', &
              topology_context_h,x
          ierr=ierr_local
          topology_probe_only=.false.
          RETURN
        ENDIF
        IF(topology_signature.EQ.probe_sig(position)) THEN
          probe_x(position)=xmid
        ELSEIF(topology_signature.EQ.probe_sig(position+1)) THEN
          probe_x(position+1)=xmid
        ELSE
          CALL insert_probe(position+1,xmid,topology_signature)
        ENDIF
      ENDDO
      position=position+1
    ENDIF
  ENDDO
!
! Convert adjacent valid samples into open fixed-topology intervals.  A zero
! width component at a global endpoint is a measure-zero endpoint state; an
! interior component below the certification tolerance is an unresolved
! topology and is rejected rather than discarded.
  ALLOCATE(segment_lo(nprobe),segment_hi(nprobe),segment_sig(nprobe))
  nsegments=0
  segment_start=1
  DO i=1,nprobe-1
    IF(probe_sig(i).NE.probe_sig(i+1)) THEN
      segment_end=i
      IF(probe_x(segment_end)-probe_x(segment_start).GT.topology_tol) THEN
        nsegments=nsegments+1
        segment_lo(nsegments)=probe_x(segment_start)
        segment_hi(nsegments)=probe_x(segment_end)
        segment_sig(nsegments)=probe_sig(segment_start)
      ELSEIF(segment_start.NE.1) THEN
        PRINT *,'sample_matrix_out_partitioned: unresolved interior component H,J = ', &
            topology_context_h,probe_x(segment_start),probe_x(segment_end)
        ierr=sample_matrix_topology_unresolved
        topology_probe_only=.false.
        RETURN
      ENDIF
      segment_start=i+1
    ENDIF
  ENDDO
  segment_end=nprobe
  IF(probe_x(segment_end)-probe_x(segment_start).GT.topology_tol) THEN
    nsegments=nsegments+1
    segment_lo(nsegments)=probe_x(segment_start)
    segment_hi(nsegments)=probe_x(segment_end)
    segment_sig(nsegments)=probe_sig(segment_start)
  ELSEIF(segment_start.NE.nprobe) THEN
    PRINT *,'sample_matrix_out_partitioned: unresolved interior terminal component H,J = ', &
        topology_context_h,probe_x(segment_start),probe_x(segment_end)
    ierr=sample_matrix_topology_unresolved
    topology_probe_only=.false.
    RETURN
  ENDIF
  IF(nsegments.EQ.0) THEN
    ierr=sample_matrix_topology_unresolved
    topology_probe_only=.false.
    RETURN
  ENDIF
!
! The topology probes do not populate resonant data.  Each open interval now
! gets a fresh adaptive interpolation grid and its own callback history.
  topology_probe_only=.false.
  n_total=0
  DO i=1,nsegments
    xbeg=segment_lo(i)
    xend=segment_hi(i)
    npoi=npoi_request
    CALL sample_matrix_out(get_matrix,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      ierr=ierr_local
      xbeg=xbeg_full
      xend=xend_full
      RETURN
    ENDIF
    nseg_points=npoi
    ALLOCATE(xseg(nseg_points),amatseg(n1_request,n2_request,nseg_points), &
        indseg(nseg_points),sigseg(nseg_points))
    xseg=xarr
    amatseg=amat_arr
    indseg=ind_hist
    sigseg=topology_arr
    IF(.NOT.ALL(sigseg.EQ.segment_sig(i))) THEN
      PRINT *,'sample_matrix_out_partitioned: segment signature changed during sampling'
      ierr=sample_matrix_topology_transition
      xbeg=xbeg_full
      xend=xend_full
      RETURN
    ENDIF
    CALL append_segment(xseg,amatseg,indseg,sigseg,nseg_points)
    DEALLOCATE(xseg,amatseg,indseg,sigseg)
  ENDDO
!
  IF(ALLOCATED(xarr)) DEALLOCATE(xarr,amat_arr,ind_hist,topology_arr)
  ALLOCATE(xarr(n_total),amat_arr(n1_request,n2_request,n_total), &
      ind_hist(n_total),topology_arr(n_total))
  xarr=xout
  amat_arr=amatout
  ind_hist=indout
  topology_arr=sigout
  DEALLOCATE(xout,amatout,indout,sigout)
  npoi=n_total
  xbeg=xbeg_full
  xend=xend_full
  topology_probe_only=.false.
  IF(endpoint_contracted) THEN
    PRINT *,'sample_matrix_out_partitioned: contracted invalid global endpoint H,J = ', &
        topology_context_h,xbeg_full,xend_full
  ENDIF
!
contains

  SUBROUTINE contract_endpoint(index,boundary,direction,ierr_endpoint)
    INTEGER, INTENT(IN) :: index
    DOUBLE PRECISION, INTENT(IN) :: boundary,direction
    INTEGER, INTENT(OUT) :: ierr_endpoint
    INTEGER :: k
    DOUBLE PRECISION :: candidate,distance,invalid_distance,valid_distance
    DOUBLE PRECISION :: valid_candidate,mid_distance

    ierr_endpoint=sample_matrix_topology_unresolved
    distance=0.5d0*ABS(xend_full-xbeg_full)/DBLE(MAX(2,nprobe-1))
    DO k=1,topology_max_endpoint_attempts
      candidate=boundary+direction*distance
      IF(candidate.LE.xbeg_full) THEN
        EXIT
      ENDIF
      IF(candidate.GE.xend_full) THEN
        EXIT
      ENDIF
      x=candidate
      CALL evaluate_matrix_callback(get_matrix,ierr_local)
      IF(ierr_local.EQ.sample_matrix_success) THEN
        valid_candidate=candidate
        valid_distance=distance
        invalid_distance=0.d0
        DO WHILE(valid_distance-invalid_distance.GT.topology_tol)
          mid_distance=0.5d0*(valid_distance+invalid_distance)
          candidate=boundary+direction*mid_distance
          x=candidate
          CALL evaluate_matrix_callback(get_matrix,ierr_local)
          IF(ierr_local.EQ.sample_matrix_success) THEN
            valid_distance=mid_distance
            valid_candidate=candidate
          ELSE
            invalid_distance=mid_distance
          ENDIF
        ENDDO
        probe_x(index)=valid_candidate
        ierr_endpoint=sample_matrix_success
        RETURN
      ENDIF
      IF(distance.GE.0.5d0*ABS(xend_full-xbeg_full)) EXIT
      distance=MIN(2.d0*distance,0.5d0*ABS(xend_full-xbeg_full))
    ENDDO
  END SUBROUTINE contract_endpoint

  SUBROUTINE insert_probe(index,xinsert,siginsert)
    INTEGER, INTENT(IN) :: index,siginsert
    DOUBLE PRECISION, INTENT(IN) :: xinsert

    ALLOCATE(probe_x_new(nprobe+1),probe_sig_new(nprobe+1))
    IF(index.GT.1) THEN
      probe_x_new(1:index-1)=probe_x(1:index-1)
      probe_sig_new(1:index-1)=probe_sig(1:index-1)
    ENDIF
    probe_x_new(index)=xinsert
    probe_sig_new(index)=siginsert
    IF(index.LE.nprobe) THEN
      probe_x_new(index+1:nprobe+1)=probe_x(index:nprobe)
      probe_sig_new(index+1:nprobe+1)=probe_sig(index:nprobe)
    ENDIF
    DEALLOCATE(probe_x,probe_sig)
    CALL MOVE_ALLOC(probe_x_new,probe_x)
    CALL MOVE_ALLOC(probe_sig_new,probe_sig)
    nprobe=nprobe+1
  END SUBROUTINE insert_probe

  SUBROUTINE append_segment(xin,ain,iin,sin,nin)
    INTEGER, INTENT(IN) :: nin
    DOUBLE PRECISION, INTENT(IN) :: xin(:),ain(:,:,:)
    INTEGER, INTENT(IN) :: iin(:),sin(:)
    DOUBLE PRECISION, ALLOCATABLE :: xnew(:),anew(:,:,:)
    INTEGER, ALLOCATABLE :: inew(:),snew(:)

    ALLOCATE(xnew(n_total+nin),anew(n1_request,n2_request,n_total+nin), &
        inew(n_total+nin),snew(n_total+nin))
    IF(n_total.GT.0) THEN
      xnew(1:n_total)=xout
      anew(:,:,1:n_total)=amatout
      inew(1:n_total)=indout
      snew(1:n_total)=sigout
    ENDIF
    xnew(n_total+1:n_total+nin)=xin
    anew(:,:,n_total+1:n_total+nin)=ain
    inew(n_total+1:n_total+nin)=iin
    snew(n_total+1:n_total+nin)=sin
    IF(n_total.GT.0) THEN
      DEALLOCATE(xout,amatout,indout,sigout)
    ENDIF
    CALL MOVE_ALLOC(xnew,xout)
    CALL MOVE_ALLOC(anew,amatout)
    CALL MOVE_ALLOC(inew,indout)
    CALL MOVE_ALLOC(snew,sigout)
    n_total=n_total+nin
  END SUBROUTINE append_segment

  END SUBROUTINE sample_matrix_out_partitioned
