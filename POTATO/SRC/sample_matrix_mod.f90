  module sample_matrix_mod
    integer :: nlagr,n1,n2,npoi,itermax,nstiff,i_int
    double precision :: x,xbeg,xend,eps
    double precision, dimension(:),     allocatable :: xarr
    double precision, dimension(:,:),   allocatable :: amat
    double precision, dimension(:,:,:), allocatable :: amat_arr
! The class-grid state is per energy slice. Inner node-fill loops copy the
! current grid config into their workers before filling disjoint slots.
    !$omp threadprivate(nlagr,n1,n2,npoi,itermax,x,xbeg,xend,eps,xarr,amat,amat_arr)
  end module sample_matrix_mod

  module matrix_callback_status_mod
    implicit none

    integer, parameter :: matrix_callback_ok = 0
    integer, parameter :: matrix_callback_starter_failed = 101
    integer, parameter :: matrix_callback_bounce_failed = 102
    integer, parameter :: matrix_callback_bounds_failed = 103
    integer, parameter :: matrix_callback_classes_failed = 104
    integer, parameter :: matrix_callback_class_sampling_failed = 105
    integer, parameter :: matrix_callback_resonance_search_failed = 106
    integer, parameter :: matrix_callback_resonance_orbit_failed = 107
    integer, parameter :: matrix_callback_nonfinite_weight = 108

    integer :: matrix_callback_error = matrix_callback_ok
    integer :: jperp_attempted = 0
    integer :: jperp_succeeded = 0
    integer :: jperp_failed = 0
    integer :: class_attempted = 0
    integer :: class_succeeded = 0
    integer :: class_failed = 0
    integer :: orbit_attempted = 0
    integer :: orbit_failed = 0
    integer :: resonance_search_attempted = 0
    integer :: resonance_search_failed = 0
    integer :: resonance_root_attempted = 0
    integer :: resonance_root_succeeded = 0
    integer :: resonance_root_failed = 0
    integer :: resonance_root_zeroed = 0
    !$omp threadprivate(matrix_callback_error,jperp_attempted,jperp_succeeded, &
    !$omp               jperp_failed,class_attempted,class_succeeded, &
    !$omp               class_failed,orbit_attempted,orbit_failed, &
    !$omp               resonance_search_attempted,resonance_search_failed, &
    !$omp               resonance_root_attempted,resonance_root_succeeded, &
    !$omp               resonance_root_failed,resonance_root_zeroed)

  contains

    subroutine reset_matrix_callback_error
      matrix_callback_error = matrix_callback_ok
    end subroutine reset_matrix_callback_error

    subroutine set_matrix_callback_error(code)
      integer, intent(in) :: code

      if (matrix_callback_error == matrix_callback_ok) matrix_callback_error = code
    end subroutine set_matrix_callback_error

    subroutine reset_torque_completeness
      matrix_callback_error = matrix_callback_ok
      jperp_attempted = 0
      jperp_succeeded = 0
      jperp_failed = 0
      class_attempted = 0
      class_succeeded = 0
      class_failed = 0
      orbit_attempted = 0
      orbit_failed = 0
      resonance_search_attempted = 0
      resonance_search_failed = 0
      resonance_root_attempted = 0
      resonance_root_succeeded = 0
      resonance_root_failed = 0
      resonance_root_zeroed = 0
    end subroutine reset_torque_completeness

  end module matrix_callback_status_mod
