program test_set_radial_boxes
  implicit none

  integer, parameter :: nbox = 100
  real(8) :: sbox(nbox)
  integer :: i

  call set_radial_boxes(nbox, .false., sbox)
  do i = 1, nbox
    call require(abs(sbox(i) - dble(i)/dble(nbox)) < 1.d-14, 'uniform s_pol edge')
  end do

  call set_radial_boxes(nbox, .true., sbox)
  do i = 1, nbox
    call require(abs(sqrt(sbox(i)) - dble(i)/dble(nbox)) < 1.d-14, &
                 'uniform rho_pol edge')
  end do
  call require(abs(sbox(nbox) - 1.d0) < 1.d-14, 'outermost edge is the separatrix')

contains

  subroutine require(ok, msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg

    if (ok) return
    print *, 'FAILED: ', msg
    error stop 1
  end subroutine require

end program test_set_radial_boxes
