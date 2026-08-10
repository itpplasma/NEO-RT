  module sample_matrix_mod
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
    !$omp threadprivate(nlagr,n1,n2,npoi,itermax,x,xbeg,xend,eps,xarr,amat,amat_arr)
    !$omp threadprivate(ierr_get_matrix)
  end module sample_matrix_mod
