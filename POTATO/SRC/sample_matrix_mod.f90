  module sample_matrix_mod
    integer :: nlagr,n1,n2,npoi,itermax,nstiff,i_int
! >0 marks a matrix row whose value 0 flags an invalid (non-closing orbit) node.
! sample_matrix then refuses to refine across such nodes (see sample_matrix.f90).
    integer :: isentinel = 0
    double precision :: x,xbeg,xend,eps
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
! The whole class-grid scratch is per energy slice.  The coarse energy loop may
! build several class grids at once; each thread needs its own grid metadata,
! arrays, and per-node scratch.
    !$omp threadprivate(nlagr,n1,n2,npoi,itermax,isentinel,x,xbeg,xend,eps)
    !$omp threadprivate(xarr,amat,amat_arr)
  end module sample_matrix_mod
