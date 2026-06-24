  module sample_matrix_mod
    integer :: nlagr,n1,n2,npoi,itermax,nstiff,i_int
    double precision :: x,xbeg,xend,eps
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
! The whole class-grid scratch is per-energy-slice state: the energy loop
! (resonant_torque) builds and interpolates one class grid per slice, and the
! slices run in parallel, so each thread needs its own grid (xarr,amat_arr,npoi,
! n1,n2,eps,xbeg,xend) and per-node scratch (x,amat).  nstiff,i_int are unused.
    !$omp threadprivate(nlagr,n1,n2,npoi,itermax,x,xbeg,xend,eps,xarr,amat,amat_arr)
  end module sample_matrix_mod
