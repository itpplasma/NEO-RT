  module sample_matrix_mod
    implicit none

    integer, parameter :: sample_matrix_success=0
    integer, parameter :: sample_matrix_nonconverged=2

    integer :: nlagr,n1,n2,npoi,itermax,nstiff,i_int
! Status of the last get_matrix call: nonzero when the orbit could not be
! evaluated (left the field domain, crossed the wall).  amat then carries no
! usable value and the node must be kept out of refinement and interpolation.
    integer :: ierr_get_matrix=0, ierr_interpolate=0
    logical, dimension(:), allocatable :: valid_arr
    double precision :: x,xbeg,xend,eps
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
! The class-grid state is per energy slice. Inner node-fill loops copy the
! current grid config into their workers before filling disjoint slots.
    !$omp threadprivate(nlagr,n1,n2,npoi,itermax,x,xbeg,xend,eps,xarr,amat,amat_arr, &
    !$omp               ierr_get_matrix,ierr_interpolate,valid_arr)

  contains

    pure logical function sample_matrix_grid_usable(ierr)
      integer, intent(in) :: ierr

      sample_matrix_grid_usable = ierr.eq.sample_matrix_success .or. &
          ierr.eq.sample_matrix_nonconverged
    end function sample_matrix_grid_usable

    logical function sample_window_valid(k)
      integer, intent(in) :: k
      integer :: ibeg,iend

      ibeg=max(1,k-nlagr/2)
      iend=min(npoi,ibeg+nlagr)
      ibeg=iend-nlagr
      sample_window_valid=all(valid_arr(ibeg:iend))
    end function sample_window_valid

    subroutine find_valid_class_roots(fun,xb,xe,ierr)
! Search contiguous ranges of valid interpolation stencils separately, so a
! sign change across an unevaluable orbit cannot manufacture a resonance.
      use find_all_roots_mod, only : nroots,roots
      external :: fun
      double precision, intent(in) :: xb,xe
      integer, intent(out) :: ierr
      integer :: k,kfirst,klast,istat
      double precision :: left,right
      double precision, allocatable :: saved_roots(:)

      if(all(valid_arr)) then
        call find_all_roots_bracketed(fun,xb,xe,ierr)
        return
      endif

      ierr=0
      allocate(saved_roots(0))
      k=2
      do while(k.le.npoi)
        if(.not.sample_window_valid(k)) then
          k=k+1
          cycle
        endif
        kfirst=k
        do while(k.le.npoi)
          if(.not.sample_window_valid(k)) exit
          k=k+1
        enddo
        klast=k-1
        left=xb
        if(kfirst.gt.2) left=max(xb,xarr(kfirst-1))
        right=xe
! binsrc chooses the right stencil at a grid node. Stay just inside the
! valid range when the next stencil is invalid.
        if(klast.lt.npoi) right=min(xe,nearest(xarr(klast),-1.d0))
        if(right.le.left) cycle
        call find_all_roots_bracketed(fun,left,right,istat)
        if(istat.ne.0) ierr=istat
        if(nroots.gt.0) saved_roots=[saved_roots,roots]
      enddo
      nroots=size(saved_roots)
      call move_alloc(saved_roots,roots)
    end subroutine find_valid_class_roots

  end module sample_matrix_mod
