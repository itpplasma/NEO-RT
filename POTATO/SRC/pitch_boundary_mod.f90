module pitch_boundary_mod
  implicit none
  integer, parameter :: pitch_boundary_success=0
  integer, parameter :: pitch_boundary_forbidden=2
contains
  pure subroutine resolve_pitch_squared(ratio, pitch_squared, ierr)
    double precision, intent(in) :: ratio
    double precision, intent(out) :: pitch_squared
    integer, intent(out) :: ierr
    double precision, parameter :: relative_tolerance=1.d-9
    double precision :: tolerance

    pitch_squared=1.d0-ratio
    tolerance=relative_tolerance*max(1.d0,abs(ratio))
    if(pitch_squared.lt.-tolerance) then
      ierr=pitch_boundary_forbidden
      return
    endif
    pitch_squared=max(0.d0,pitch_squared)
    ierr=pitch_boundary_success
  end subroutine resolve_pitch_squared
end module pitch_boundary_mod
