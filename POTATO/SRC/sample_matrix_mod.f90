  module sample_matrix_mod
    integer, parameter :: matrix_eval_success=0
    integer, parameter :: matrix_eval_starter_failure=1
    integer, parameter :: matrix_eval_orbit_failure=2
    integer, parameter :: matrix_eval_wall_loss=3
    integer, parameter :: matrix_eval_nonfinite=4
    integer, parameter :: matrix_boundary_missing_limit=5
    integer :: nlagr,n1,n2,npoi,itermax,nstiff,i_int
    double precision :: x,xbeg,xend,eps
    logical :: matrix_eval_valid=.true.
    integer :: matrix_eval_error=matrix_eval_success
    integer :: matrix_boundary_error=matrix_eval_success
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
! The class-grid state is per energy slice. Inner node-fill loops copy the
! current grid config into their workers before filling disjoint slots.
    !$omp threadprivate(nlagr,n1,n2,npoi,itermax,x,xbeg,xend,eps, &
    !$omp&                matrix_eval_valid,matrix_eval_error,     &
    !$omp&                matrix_boundary_error,xarr,amat,amat_arr)

  contains

    pure logical function matrix_failure_is_open_boundary(error_code)
      integer, intent(in) :: error_code

      matrix_failure_is_open_boundary = error_code == matrix_eval_orbit_failure .or. &
                                        error_code == matrix_eval_wall_loss
    end function matrix_failure_is_open_boundary

  end module sample_matrix_mod
