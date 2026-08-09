  module sample_matrix_out_mod
    integer, parameter :: sample_matrix_success=0
    integer, parameter :: sample_matrix_nonconverged=2
    integer, parameter :: sample_matrix_topology_transition=3
    integer, parameter :: sample_matrix_topology_unresolved=4
    integer, parameter :: sample_matrix_contribution_unresolved=5
    integer :: nlagr,n1,n2,npoi=0,itermax=0,icount
    double precision :: x,xbeg,xend,eps
    integer :: topology_signature=0,topology_error=0
    double precision :: topology_context_h=0.d0
    logical :: topology_probe_only=.false.
    logical :: sample_matrix_preserve_history=.false.
    integer :: topology_candidate_count=0,topology_transition_count=0
    double precision :: topology_gap_measure=0.d0
    double precision :: topology_gap_geometric_bound=0.d0
    ! Kept as a compatibility alias for existing diagnostics; it is geometric
    ! only and must never be read as a torque-error estimate.
    double precision :: topology_gap_bound=0.d0
    double precision :: topology_contribution_error_bound=huge(1.d0)
    logical :: topology_contribution_error_certified=.false.
    double precision :: topology_partition_tol=1.d-8
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
    !$omp               sample_matrix_preserve_history, &
    !$omp               topology_candidate_count,topology_transition_count, &
    !$omp               topology_gap_measure,topology_gap_geometric_bound, &
    !$omp               topology_gap_bound,topology_contribution_error_bound, &
    !$omp               topology_contribution_error_certified, &
    !$omp               topology_partition_tol, &
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
  USE potato_symbolic_kernel_mod, ONLY : potato_gap_contribution_kernel
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
  if(.not.sample_matrix_preserve_history) icount=0
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
! Deprecated observational entry point.  It is retained as a diagnostic
! baseline for the narrow-island regression, but it is never an acceptance
! path: without an explicit complete root certificate it fails closed.
!
  USE sample_matrix_out_mod
!
  IMPLICIT NONE
!
  INTEGER, INTENT(OUT) :: ierr
!
  EXTERNAL :: get_matrix
!
  PRINT *,'sample_matrix_out_partitioned: complete topology certificate required'
  ierr=sample_matrix_topology_unresolved
  topology_probe_only=.false.
  sample_matrix_preserve_history=.false.
  RETURN
  END SUBROUTINE sample_matrix_out_partitioned

  SUBROUTINE sample_matrix_out_partitioned_certified(get_matrix,get_boundaries,ierr, &
                                                     get_envelope)
!
! Certified replacement for the observational sampler above.  get_boundaries
! returns the complete root set of the physical class/fixpoint boundary
! equations for the current H.  The callback signature is
!
!   get_boundaries(nmax,nfound,boundaries,ierr)
!
! Equal signatures without that certificate are not accepted: finite midpoint
! sampling cannot exclude a narrow A-B-A island.  Certified roots are kept as
! open boundaries and the finite numerical brackets are reported as an explicit
! geometric integration gap.  topology_gap_geometric_bound bounds only its
! coordinate measure.  A get_envelope callback is required to turn that
! measure into a certified contribution-error bound; omission fails closed.
!
  USE sample_matrix_out_mod
  USE potato_symbolic_kernel_mod, ONLY : potato_gap_contribution_kernel
  IMPLICIT NONE
!
  INTEGER, INTENT(OUT) :: ierr
  EXTERNAL :: get_matrix,get_boundaries,get_envelope
  INTEGER, PARAMETER :: topology_max_candidates=4096
  INTEGER, PARAMETER :: topology_max_endpoint_attempts=16
  DOUBLE PRECISION, PARAMETER :: topology_abs_tol=256.d0*EPSILON(1.d0)
  INTEGER :: i,j,nfound,nunique,nactive,nsegments,n_total,nseg_points
  INTEGER :: npoi_request,n1_request,n2_request,ierr_local,ierr_certificate
  INTEGER :: left_endpoint_sig,right_endpoint_sig,sig_l,sig_r,sig_expected
  DOUBLE PRECISION :: xbeg_full,xend_full,scale,topology_tol
  DOUBLE PRECISION :: candidate_merge_tol
  DOUBLE PRECISION :: key,val,prev_bound,next_bound,delta,endpoint_tol, &
                      delta_resolution
  DOUBLE PRECISION :: left_open,right_open,left_gap,right_gap,covered
  DOUBLE PRECISION :: geometric_gap_bound
  DOUBLE PRECISION, ALLOCATABLE :: candidates(:),active_x(:),active_delta(:)
  INTEGER, ALLOCATABLE :: active_left_sig(:),active_right_sig(:)
  DOUBLE PRECISION, ALLOCATABLE :: segment_lo(:),segment_hi(:)
  INTEGER, ALLOCATABLE :: segment_sig(:)
  DOUBLE PRECISION, ALLOCATABLE :: xout(:),amatout(:,:,:),xseg(:),amatseg(:,:,:)
  INTEGER, ALLOCATABLE :: indout(:),sigout(:),indseg(:),sigseg(:)
!
  ierr=sample_matrix_success
  topology_candidate_count=0
  topology_transition_count=0
  topology_gap_measure=0.d0
  topology_gap_geometric_bound=0.d0
  topology_gap_bound=0.d0
  topology_contribution_error_bound=huge(1.d0)
  topology_contribution_error_certified=.false.
  xbeg_full=xbeg
  xend_full=xend
  npoi_request=MAX(npoi,nlagr+2)
  n1_request=n1
  n2_request=n2
  IF(xend_full.LE.xbeg_full) THEN
    PRINT *,'sample_matrix_out_partitioned_certified: invalid J_perp domain = ', &
        xbeg_full,xend_full
    ierr=sample_matrix_topology_unresolved
    RETURN
  ENDIF
  IF(topology_partition_tol.LE.0.d0) THEN
    PRINT *,'sample_matrix_out_partitioned_certified: non-positive partition tolerance'
    ierr=sample_matrix_topology_unresolved
    RETURN
  ENDIF
  scale=MAX(1.d0,ABS(xbeg_full),ABS(xend_full))
  topology_tol=MAX(topology_abs_tol*scale, &
                   topology_partition_tol*ABS(xend_full-xbeg_full))
! Distinct certified roots must not be merged at the partition tolerance:
! an arbitrarily narrow A-B-A component is still a real component.  Only
! roots indistinguishable at the floating-point coordinate scale are merged;
! transitions outside that coordinate-resolution window remain distinct.
  ALLOCATE(candidates(topology_max_candidates),active_x(topology_max_candidates), &
      active_delta(topology_max_candidates),active_left_sig(topology_max_candidates), &
      active_right_sig(topology_max_candidates))
!
! Ask the physics layer for a complete, bounded root set.  This is the only
! accepted source of transition locations.
  nfound=0
  ierr_certificate=sample_matrix_success
  CALL get_boundaries(topology_max_candidates,nfound,candidates,ierr_certificate)
  IF(ierr_certificate.NE.sample_matrix_success .OR. &
     nfound.LT.0 .OR. nfound.GT.topology_max_candidates) THEN
    PRINT *,'sample_matrix_out_partitioned_certified: invalid certificate H,n,ierr = ', &
        topology_context_h,nfound,ierr_certificate
    CALL fail_certified(sample_matrix_topology_unresolved)
    RETURN
  ENDIF
  DO i=1,nfound
    val=candidates(i)
    IF(val.NE.val .OR. ABS(val).GT.HUGE(val)) THEN
      PRINT *,'sample_matrix_out_partitioned_certified: non-finite root H,J = ', &
          topology_context_h,val
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
  ENDDO
  DO i=2,nfound
    key=candidates(i)
    j=i-1
    DO WHILE(j.GE.1)
      IF(candidates(j).LE.key) EXIT
      candidates(j+1)=candidates(j)
      j=j-1
    ENDDO
    candidates(j+1)=key
  ENDDO
  nunique=0
  DO i=1,nfound
    endpoint_tol=128.d0*EPSILON(1.d0)*MAX(1.d-300, &
        ABS(candidates(i)),ABS(xbeg_full))
    IF(candidates(i).LT.xbeg_full-endpoint_tol .OR. &
       candidates(i).GT.xend_full+endpoint_tol) CYCLE
    IF(candidates(i).LE.xbeg_full+endpoint_tol) THEN
      IF(ABS(candidates(i)-xbeg_full).LE.endpoint_tol) CYCLE
      PRINT *,'sample_matrix_out_partitioned_certified: unresolved endpoint root H,J,endpoint,tol = ', &
          topology_context_h,candidates(i),xbeg_full,endpoint_tol
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
    endpoint_tol=128.d0*EPSILON(1.d0)*MAX(1.d-300, &
        ABS(candidates(i)),ABS(xend_full))
    IF(candidates(i).GE.xend_full-endpoint_tol) THEN
      IF(ABS(candidates(i)-xend_full).LE.endpoint_tol) CYCLE
      PRINT *,'sample_matrix_out_partitioned_certified: unresolved endpoint root H,J,endpoint,tol = ', &
          topology_context_h,candidates(i),xend_full,endpoint_tol
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
    IF(nunique.GT.0) THEN
      IF(candidates(i).EQ.candidates(nunique)) CYCLE
      candidate_merge_tol=32.d0*EPSILON(1.d0)*MAX(1.d-300, &
          ABS(candidates(i)),ABS(candidates(nunique)))
      IF(ABS(candidates(i)-candidates(nunique)).LE.candidate_merge_tol) CYCLE
    ENDIF
    nunique=nunique+1
    candidates(nunique)=candidates(i)
  ENDDO
  topology_candidate_count=nunique
!
! Establish valid one-sided endpoint samples.  An exact endpoint is used only
! when its signature agrees with the open-side probe; otherwise it is contracted
! and the contraction is charged to the explicit gap.
  topology_probe_only=.true.
  x=xbeg_full
  CALL evaluate_matrix_callback(get_matrix,ierr_local)
  IF(ierr_local.EQ.sample_matrix_success) THEN
    left_endpoint_sig=topology_signature
    delta=MIN(topology_tol,0.25d0*(xend_full-xbeg_full))
    x=xbeg_full+delta
    CALL evaluate_matrix_callback(get_matrix,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
    IF(topology_signature.NE.left_endpoint_sig) THEN
      left_open=x
      left_gap=delta
      left_endpoint_sig=topology_signature
    ELSE
      left_open=xbeg_full
      left_gap=0.d0
    ENDIF
  ELSE
    CALL contract_certified_endpoint(xbeg_full,1.d0,left_open,left_gap, &
        left_endpoint_sig,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
  ENDIF
  x=xend_full
  CALL evaluate_matrix_callback(get_matrix,ierr_local)
  IF(ierr_local.EQ.sample_matrix_success) THEN
    right_endpoint_sig=topology_signature
    delta=MIN(topology_tol,0.25d0*(xend_full-xbeg_full))
    x=xend_full-delta
    CALL evaluate_matrix_callback(get_matrix,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
    IF(topology_signature.NE.right_endpoint_sig) THEN
      right_open=x
      right_gap=delta
      right_endpoint_sig=topology_signature
    ELSE
      right_open=xend_full
      right_gap=0.d0
    ENDIF
  ELSE
    CALL contract_certified_endpoint(xend_full,-1.d0,right_open,right_gap, &
        right_endpoint_sig,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
  ENDIF
!
! Probe both sides of each candidate.  The candidate itself is not evaluated:
! it is a certified topology boundary, and its two open sides are the physical
! interpolation domains.  A harmless candidate whose side signatures agree is
! retained in the certificate count but creates no omitted bracket.
  nactive=0
  DO i=1,nunique
    prev_bound=xbeg_full
    IF(i.GT.1) prev_bound=candidates(i-1)
    next_bound=xend_full
    IF(i.LT.nunique) next_bound=candidates(i+1)
    delta=MIN(topology_tol,0.25d0*(candidates(i)-prev_bound), &
              0.25d0*(next_bound-candidates(i)))
    delta_resolution=128.d0*EPSILON(1.d0)*MAX(1.d-300, &
        ABS(candidates(i)),ABS(prev_bound),ABS(next_bound))
    IF(delta.LE.delta_resolution) THEN
      PRINT *,'sample_matrix_out_partitioned_certified: unresolved root spacing H,J = ', &
          topology_context_h,candidates(i),delta,prev_bound,next_bound,delta_resolution
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
    x=candidates(i)-delta
    CALL evaluate_matrix_callback(get_matrix,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      PRINT *,'sample_matrix_out_partitioned_certified: invalid left side H,J = ', &
          topology_context_h,x
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
    sig_l=topology_signature
    x=candidates(i)+delta
    CALL evaluate_matrix_callback(get_matrix,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      PRINT *,'sample_matrix_out_partitioned_certified: invalid right side H,J = ', &
          topology_context_h,x
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
    sig_r=topology_signature
    IF(sig_l.NE.sig_r) THEN
      nactive=nactive+1
      active_x(nactive)=candidates(i)
      active_delta(nactive)=delta
      active_left_sig(nactive)=sig_l
      active_right_sig(nactive)=sig_r
    ENDIF
  ENDDO
  topology_transition_count=nactive
!
! Check that the certificate produces one consistent ordered sequence of open
! topology classes.  This catches a stale H/J root certificate before matrix
! sampling begins.
  IF(nactive.EQ.0) THEN
    IF(left_endpoint_sig.NE.right_endpoint_sig) THEN
      PRINT *,'sample_matrix_out_partitioned_certified: endpoint signatures disagree H,left,right,nunique,nactive = ', &
          topology_context_h,left_endpoint_sig,right_endpoint_sig,nunique,nactive
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
  ELSE
    IF(left_endpoint_sig.NE.active_left_sig(1)) THEN
      PRINT *,'sample_matrix_out_partitioned_certified: left certificate mismatch'
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
    DO i=1,nactive-1
      IF(active_right_sig(i).NE.active_left_sig(i+1)) THEN
        PRINT *,'sample_matrix_out_partitioned_certified: adjacent certificate mismatch'
        CALL fail_certified(sample_matrix_topology_unresolved)
        RETURN
      ENDIF
    ENDDO
    IF(active_right_sig(nactive).NE.right_endpoint_sig) THEN
      PRINT *,'sample_matrix_out_partitioned_certified: right certificate mismatch'
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
  ENDIF
!
  nsegments=nactive+1
  ALLOCATE(segment_lo(nsegments),segment_hi(nsegments),segment_sig(nsegments))
  DO i=1,nsegments
    IF(i.EQ.1) THEN
      segment_lo(i)=left_open
      sig_expected=left_endpoint_sig
    ELSE
      segment_lo(i)=active_x(i-1)+active_delta(i-1)
      sig_expected=active_right_sig(i-1)
    ENDIF
    IF(i.EQ.nsegments) THEN
      segment_hi(i)=right_open
    ELSE
      segment_hi(i)=active_x(i)-active_delta(i)
    ENDIF
    segment_sig(i)=sig_expected
    delta_resolution=128.d0*EPSILON(1.d0)*MAX(1.d-300, &
        ABS(segment_hi(i)),ABS(segment_lo(i)))
    IF(segment_hi(i)-segment_lo(i).LE.delta_resolution) THEN
      PRINT *,'sample_matrix_out_partitioned_certified: empty open component H,Jlo,Jhi = ', &
          topology_context_h,segment_lo(i),segment_hi(i)
      CALL fail_certified(sample_matrix_topology_unresolved)
      RETURN
    ENDIF
  ENDDO
!
! Account for, but never integrate across, the transition brackets.
  covered=0.d0
  DO i=1,nsegments
    covered=covered+segment_hi(i)-segment_lo(i)
  ENDDO
  topology_gap_measure=MAX(0.d0,(xend_full-xbeg_full)-covered)
  geometric_gap_bound=left_gap+right_gap
  DO i=1,nactive
    geometric_gap_bound=geometric_gap_bound+2.d0*active_delta(i)
  ENDDO
  topology_gap_geometric_bound=geometric_gap_bound
  topology_gap_bound=topology_gap_geometric_bound
  IF(topology_gap_measure.GT.topology_gap_geometric_bound+128.d0*EPSILON(1.d0)*scale) THEN
    PRINT *,'sample_matrix_out_partitioned_certified: gap accounting failed H,gap,bound = ', &
        topology_context_h,topology_gap_measure,topology_gap_geometric_bound
    CALL fail_certified(sample_matrix_topology_unresolved)
    RETURN
  ENDIF
  topology_contribution_error_bound=0.d0
  CALL accumulate_gap_envelope(xbeg_full,left_open,ierr_local)
  IF(ierr_local.EQ.sample_matrix_success) THEN
    CALL accumulate_gap_envelope(right_open,xend_full,ierr_local)
  ENDIF
  IF(ierr_local.EQ.sample_matrix_success) THEN
    DO i=1,nactive
      CALL accumulate_gap_envelope(candidates(i)-active_delta(i), &
          candidates(i)+active_delta(i),ierr_local)
      IF(ierr_local.NE.sample_matrix_success) EXIT
    ENDDO
  ENDIF
  IF(ierr_local.NE.sample_matrix_success) THEN
    PRINT *,'sample_matrix_out_partitioned_certified: contribution envelope failed H,ierr = ', &
        topology_context_h,ierr_local
    CALL fail_certified(sample_matrix_contribution_unresolved)
    RETURN
  ENDIF
  topology_contribution_error_certified=.true.
!
! Keep the callback history global across segments.  This is required for the
! respoints_all/ind_hist permutation consumed by resonant_torque.
  topology_probe_only=.false.
  sample_matrix_preserve_history=.true.
  icount=0
  n_total=0
  DO i=1,nsegments
    xbeg=segment_lo(i)
    xend=segment_hi(i)
    npoi=npoi_request
    CALL sample_matrix_out(get_matrix,ierr_local)
    IF(ierr_local.NE.sample_matrix_success) THEN
      CALL fail_certified(ierr_local)
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
      PRINT *,'sample_matrix_out_partitioned_certified: segment signature changed'
      DEALLOCATE(xseg,amatseg,indseg,sigseg)
      CALL fail_certified(sample_matrix_topology_transition)
      RETURN
    ENDIF
    CALL append_segment(xseg,amatseg,indseg,sigseg,nseg_points)
    DEALLOCATE(xseg,amatseg,indseg,sigseg)
  ENDDO
  IF(n_total.LE.0) THEN
    CALL fail_certified(sample_matrix_topology_unresolved)
    RETURN
  ENDIF
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
  sample_matrix_preserve_history=.false.
  DEALLOCATE(candidates,active_x,active_delta,active_left_sig,active_right_sig)
  DEALLOCATE(segment_lo,segment_hi,segment_sig)
  PRINT *,'sample_matrix_out_partitioned_certified: ', &
      'H,candidates,transitions,geometric_gap,geometric_bound,', &
      'contribution_bound,tol = ', &
      topology_context_h,topology_candidate_count,topology_transition_count, &
      topology_gap_measure,topology_gap_geometric_bound, &
      topology_contribution_error_bound,topology_tol
  RETURN
!
CONTAINS

  SUBROUTINE accumulate_gap_envelope(gap_lo,gap_hi,status)
    DOUBLE PRECISION, INTENT(IN) :: gap_lo,gap_hi
    INTEGER, INTENT(OUT) :: status
    DOUBLE PRECISION :: envelope,measured_gap,contribution

    status=sample_matrix_success
    IF(gap_hi.LE.gap_lo) RETURN
    CALL get_envelope(gap_lo,gap_hi,envelope,status)
    IF(status.NE.sample_matrix_success) RETURN
    IF(envelope.NE.envelope) THEN
      status=sample_matrix_contribution_unresolved
      RETURN
    ENDIF
    IF(envelope.LT.0.d0) THEN
      status=sample_matrix_contribution_unresolved
      RETURN
    ENDIF
    IF(ABS(envelope).GT.HUGE(envelope)) THEN
      status=sample_matrix_contribution_unresolved
      RETURN
    ENDIF
    CALL potato_gap_contribution_kernel(envelope,gap_lo,gap_hi,measured_gap, &
        contribution)
    IF(measured_gap.NE.measured_gap .OR. contribution.NE.contribution .OR. &
        measured_gap.LT.0.d0 .OR. contribution.LT.0.d0 .OR. &
        ABS(measured_gap).GT.HUGE(measured_gap) .OR. &
        ABS(contribution).GT.HUGE(contribution)) THEN
      status=sample_matrix_contribution_unresolved
      RETURN
    ENDIF
    topology_contribution_error_bound=topology_contribution_error_bound+contribution
    IF(topology_contribution_error_bound.NE.topology_contribution_error_bound) THEN
      status=sample_matrix_contribution_unresolved
    ENDIF
  END SUBROUTINE accumulate_gap_envelope

  SUBROUTINE fail_certified(status)
    INTEGER, INTENT(IN) :: status
  topology_probe_only=.false.
  sample_matrix_preserve_history=.false.
  if(allocated(candidates)) deallocate(candidates,active_x,active_delta, &
      active_left_sig,active_right_sig)
  if(allocated(segment_lo)) deallocate(segment_lo,segment_hi,segment_sig)
  if(allocated(xout)) deallocate(xout,amatout,indout,sigout)
  ierr=status
  END SUBROUTINE fail_certified

  SUBROUTINE contract_certified_endpoint(boundary,direction,valid_point,gap,sig_out,status)
    DOUBLE PRECISION, INTENT(IN) :: boundary,direction
    DOUBLE PRECISION, INTENT(OUT) :: valid_point,gap
    INTEGER, INTENT(OUT) :: sig_out,status
    INTEGER :: k
    DOUBLE PRECISION :: distance,invalid_distance,valid_distance
    DOUBLE PRECISION :: candidate,valid_candidate,mid_distance

    status=sample_matrix_topology_unresolved
    sig_out=0
    valid_point=boundary
    gap=0.d0
    distance=0.5d0*ABS(xend_full-xbeg_full)/31.d0
    DO k=1,topology_max_endpoint_attempts
      candidate=boundary+direction*distance
      IF(candidate.LE.xbeg_full .OR. candidate.GE.xend_full) EXIT
      x=candidate
      CALL evaluate_matrix_callback(get_matrix,ierr_local)
      IF(ierr_local.EQ.sample_matrix_success) THEN
        valid_candidate=candidate
        valid_distance=distance
        invalid_distance=0.d0
        sig_out=topology_signature
        DO WHILE(valid_distance-invalid_distance.GT.topology_tol)
          mid_distance=0.5d0*(valid_distance+invalid_distance)
          candidate=boundary+direction*mid_distance
          x=candidate
          CALL evaluate_matrix_callback(get_matrix,ierr_local)
          IF(ierr_local.EQ.sample_matrix_success) THEN
            valid_distance=mid_distance
            valid_candidate=candidate
            sig_out=topology_signature
          ELSE
            invalid_distance=mid_distance
          ENDIF
        ENDDO
        valid_point=valid_candidate
        gap=ABS(valid_candidate-boundary)
        status=sample_matrix_success
        RETURN
      ENDIF
      IF(distance.GE.0.5d0*ABS(xend_full-xbeg_full)) EXIT
      distance=MIN(2.d0*distance,0.5d0*ABS(xend_full-xbeg_full))
    ENDDO
  END SUBROUTINE contract_certified_endpoint

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
    IF(n_total.GT.0) DEALLOCATE(xout,amatout,indout,sigout)
    CALL MOVE_ALLOC(xnew,xout)
    CALL MOVE_ALLOC(anew,amatout)
    CALL MOVE_ALLOC(inew,indout)
    CALL MOVE_ALLOC(snew,sigout)
    n_total=n_total+nin
  END SUBROUTINE append_segment

  END SUBROUTINE sample_matrix_out_partitioned_certified
