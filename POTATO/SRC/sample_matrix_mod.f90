  module sample_matrix_mod
    implicit none

    integer, parameter :: sample_matrix_success=0
    integer, parameter :: sample_matrix_nonconverged=2

    integer :: nlagr,n1,n2,npoi,itermax,nstiff,i_int
    double precision :: x,xbeg,xend,eps
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
! The class-grid state is per energy slice. Inner node-fill loops copy the
! current grid config into their workers before filling disjoint slots.
    !$omp threadprivate(nlagr,n1,n2,npoi,itermax,x,xbeg,xend,eps,xarr,amat,amat_arr)

  contains

    pure logical function sample_matrix_grid_usable(ierr)
      integer, intent(in) :: ierr

      sample_matrix_grid_usable = ierr.eq.sample_matrix_success .or. &
          ierr.eq.sample_matrix_nonconverged
    end function sample_matrix_grid_usable

  end module sample_matrix_mod
