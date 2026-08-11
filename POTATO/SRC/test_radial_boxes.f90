program test_radial_boxes
  implicit none

  integer, parameter :: nbox=4
  real(8) :: sbox(nbox)
  integer :: i

  call set_radial_boxes(nbox, .false., sbox)
  do i=1,nbox
    if (abs(sbox(i)-real(i,8)/real(nbox,8)) > 1d-14) then
      error stop 'non-uniform radial-box boundaries are wrong'
    endif
  enddo

  call set_radial_boxes(nbox, .true., sbox)
  do i=1,nbox
    if (abs(sbox(i)-(real(i,8)/real(nbox,8))**2) > 1d-14) then
      error stop 'uniform-rho_pol radial-box boundaries are wrong'
    endif
  enddo
end program test_radial_boxes
