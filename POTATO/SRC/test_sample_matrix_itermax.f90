program test_sample_matrix_itermax
! The sampler must retain a usable last grid when the refinement budget is
! exhausted, while reporting that the requested tolerance was not reached.
  use sample_matrix_mod, only : nlagr,n1,n2,npoi,itermax,x,xbeg,xend,eps, &
      xarr,amat,amat_arr,sample_matrix_success,sample_matrix_nonconverged
  implicit none

  integer :: ierr,i,npoi_initial
  complex(8) :: expected
  external :: sample_matrix,quartic_matrix,linear_matrix

  nlagr=2
  n1=1
  n2=1
  itermax=1
  eps=1.d-12
  xbeg=0.d0
  xend=1.d0

  call sample_matrix(quartic_matrix,ierr)
  call require(ierr.eq.sample_matrix_nonconverged, &
      'itermax cap did not report nonconvergence')
  npoi_initial=nlagr+2
  call require(npoi.gt.npoi_initial, &
      'itermax cap did not retain the completed refinement pass')
  call require(allocated(xarr).and.allocated(amat_arr), &
      'itermax cap did not retain the interpolation grid')
  do i=1,npoi
    expected=cmplx(xarr(i)**4,0.d0,kind=8)
    call require(abs(amat_arr(1,1,i)-expected).lt.1.d-13, &
        'retained grid contains a value not produced by the callback')
  enddo

  itermax=1
  eps=1.d-12
  call sample_matrix(linear_matrix,ierr)
  call require(ierr.eq.sample_matrix_success, &
      'converged sample was reported as nonconverged')
  call require(npoi.eq.npoi_initial, &
      'converged linear sample was refined unnecessarily')

  print *,'test_sample_matrix_itermax: OK'

contains

  subroutine require(ok,message)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: message

    if(.not.ok) error stop message
  end subroutine require

end program test_sample_matrix_itermax

subroutine quartic_matrix
  use sample_matrix_mod, only : x,amat
  implicit none

  amat(1,1)=cmplx(x**4,0.d0,kind=8)
end subroutine quartic_matrix

subroutine linear_matrix
  use sample_matrix_mod, only : x,amat
  implicit none

  amat(1,1)=cmplx(x,0.d0,kind=8)
end subroutine linear_matrix
