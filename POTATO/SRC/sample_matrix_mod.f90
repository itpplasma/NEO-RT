  module sample_matrix_mod
    integer :: nlagr,n1,n2,npoi,itermax,nstiff,i_int
    double precision :: x,xbeg,xend,eps
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
! x and amat are the per-node input/output of get_matrix_doublecount: each call
! reads x and writes amat for one grid node.  The grid-build node-fill loops in
! sample_matrix run those calls in parallel over disjoint nodes, so x and amat
! are genuine per-thread scratch and must be threadprivate.  The grid layout
! (npoi,xarr,amat_arr,n1,n2,...) stays shared: the adaptive split/iter logic is
! serial and amat_arr(:,:,i) writes are to disjoint slots.
    !$omp threadprivate(x,amat)
  end module sample_matrix_mod
