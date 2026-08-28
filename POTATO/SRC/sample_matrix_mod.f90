  module sample_matrix_mod
    implicit none

    integer, parameter :: sample_matrix_success=0
    integer, parameter :: sample_matrix_nonconverged=2

    integer :: nlagr,n1,n2,npoi,itermax,nstiff,i_int
! Status of the last get_matrix call: nonzero when the orbit could not be
! evaluated (left the field domain, crossed the wall).  amat then carries no
! usable value and the node must be kept out of the refinement decision.
    integer :: ierr_get_matrix
    double precision :: x,xbeg,xend,eps
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
! The class-grid state is per energy slice. Inner node-fill loops copy the
! current grid config into their workers before filling disjoint slots.
    !$omp threadprivate(nlagr,n1,n2,npoi,itermax,x,xbeg,xend,eps,xarr,amat,amat_arr, &
    !$omp               ierr_get_matrix)

  contains

    pure logical function sample_matrix_grid_usable(ierr)
      integer, intent(in) :: ierr

      sample_matrix_grid_usable = ierr.eq.sample_matrix_success .or. &
          ierr.eq.sample_matrix_nonconverged
    end function sample_matrix_grid_usable
  end module sample_matrix_mod
