program test_wall_loss
    use potato_input_mod, only: read_potato_input
    use wall_loss_mod, only: outside_wall, wall_loaded
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

    open (newunit=iunit, file='convexwall.dat', status='replace', &
          action='write')
    write (iunit, *) 100.0d0, 100.0d0
    write (iunit, *) 110.0d0, 100.0d0
    write (iunit, *) 110.0d0, 110.0d0
    write (iunit, *) 100.0d0, 110.0d0
    close (iunit)

    open (newunit=iunit, file='field_divB0.inp', status='replace', &
          action='write')
    write (iunit, '(A)') '0  ipert'
    write (iunit, '(A)') '1  iequil'
    write (iunit, '(A)') '1.0  ampl'
    write (iunit, '(A)') '1  ntor'
    write (iunit, '(A)') '0.99  cutoff'
    write (iunit, '(A)') '4  icftype'
    write (iunit, '(A)') "'unused.eqdsk'  gfile"
    write (iunit, '(A)') "'unused'  pfile"
    write (iunit, '(A)') "'test_wall_square.dat'  convexfile"
    write (iunit, '(A)') "'unused'  fluxdatapath"
    write (iunit, '(A)') '0  nwindow_r'
    write (iunit, '(A)') '0  nwindow_z'
    write (iunit, '(A)') '1  ieqfile'
    close (iunit)

    open (newunit=iunit, file='test_potato_wall.in', status='replace', &
          action='write')
    write (iunit, '(A)') '&potato_nml'
    write (iunit, '(A)') '/'
    close (iunit)

    call read_potato_input('test_potato_wall.in')

    ok = .true.
    if (.not. wall_loaded) ok = .false.
    if (outside_wall(5.0d0, 5.0d0)) ok = .false.
    if (outside_wall(9.0d0, 9.0d0)) ok = .false.
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
