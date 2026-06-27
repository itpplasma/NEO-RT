program test_wall_loss
    use wall_loss_mod, only: load_wall, outside_wall, wall_loaded
    implicit none

    integer :: iunit
    logical :: ok

    open (newunit=iunit, file='test_wall_square.dat', status='replace', &
          action='write')
    write (iunit, *) 0.0d0, 0.0d0
    write (iunit, *) 10.0d0, 0.0d0
    write (iunit, *) 10.0d0, 10.0d0
    write (iunit, *) 0.0d0, 10.0d0
    close (iunit)

    call load_wall('test_wall_square.dat')

    ok = .true.
    if (.not. wall_loaded) ok = .false.
    if (outside_wall(5.0d0, 5.0d0)) ok = .false.
    if (.not. outside_wall(15.0d0, 5.0d0)) ok = .false.
    if (.not. outside_wall(5.0d0, -1.0d0)) ok = .false.
    if (.not. outside_wall(-1.0d0, 5.0d0)) ok = .false.
    if (.not. outside_wall(5.0d0, 11.0d0)) ok = .false.

    if (ok) then
        print *, 'test_wall_loss: PASSED'
    else
        print *, 'test_wall_loss: FAILED'
        error stop 1
    end if
end program test_wall_loss
