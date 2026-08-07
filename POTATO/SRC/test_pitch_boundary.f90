program test_pitch_boundary
  use pitch_boundary_mod, only : pitch_boundary_forbidden, &
                                 pitch_boundary_success,resolve_pitch_squared
  implicit none

  double precision :: pitch_squared
  integer :: ierr

  call resolve_pitch_squared(0.75d0,pitch_squared,ierr)
  call require(ierr.eq.pitch_boundary_success,'interior status')
  call require(abs(pitch_squared-0.25d0).lt.1.d-15,'interior value')

  call resolve_pitch_squared(1.d0+5.d-10,pitch_squared,ierr)
  call require(ierr.eq.pitch_boundary_success,'roundoff boundary status')
  call require(pitch_squared.eq.0.d0,'roundoff boundary clamp')

  call resolve_pitch_squared(1.d0+2.d-9,pitch_squared,ierr)
  call require(ierr.eq.pitch_boundary_forbidden,'physical forbidden status')

contains

  subroutine require(ok,msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg

    if(ok) return
    print *,trim(msg)
    error stop 1
  end subroutine require

end program test_pitch_boundary
