program test_wall_config
    use neort_config, only: config_t, read_and_set_config, validate_frequency_model
    use neort_wall_io, only: WALL_IO_OK, load_wall_polygon, wall_polygon_t
    use driftorbit, only: FREQUENCY_MODEL_GC_FULL

    implicit none

    character(len=64) :: mode
    integer :: unit
    type(config_t) :: config
    type(wall_polygon_t) :: wall
    integer :: status

    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('missing')
        config = config_t()
        config%frequency_model = FREQUENCY_MODEL_GC_FULL
        config%inp_swi = 11
        config%wall_file = 'missing_wall_for_config.dat'
        call validate_frequency_model(FREQUENCY_MODEL_GC_FULL, 11, .true., &
            .false., .false., config%wall_file, 'm')
        error stop 'missing model-2 wall was not rejected'
    case ('malformed')
        call write_file('malformed_wall_for_config.dat', '1.0 0.0'//new_line('a')//'2.0')
        config = config_t()
        config%frequency_model = FREQUENCY_MODEL_GC_FULL
        config%inp_swi = 11
        config%wall_file = 'malformed_wall_for_config.dat'
        call validate_frequency_model(FREQUENCY_MODEL_GC_FULL, 11, .true., &
            .false., .false., config%wall_file, 'm')
        error stop 'malformed model-2 wall was not rejected'
    case ('namelist_missing')
        call write_file('missing_wall_for_config.in', &
            '&params'//new_line('a')//'    frequency_model = 2'//new_line('a')// &
            '    inp_swi = 11'//new_line('a')//'/')
        call read_and_set_config('missing_wall_for_config.in')
        error stop 'model-2 namelist without wall_file was not rejected'
    case default
        call write_file('valid_wall_for_config.dat', &
            '1.0 0.0'//new_line('a')//'2.0 0.0'//new_line('a')//'2.0 1.0')
        call load_wall_polygon('valid_wall_for_config.dat', wall, status)
        if (status /= WALL_IO_OK) error stop 'test wall could not be created'

        call validate_frequency_model(0, 0, .false., .false., .false., &
            'missing_wall_for_config.dat')
        call validate_frequency_model(1, 11, .true., .false., .false., &
            'missing_wall_for_config.dat')
        call validate_frequency_model(FREQUENCY_MODEL_GC_FULL, 11, .true., &
            .false., .false., 'valid_wall_for_config.dat', 'm')

        open (newunit=unit, file='valid_wall_for_config.dat', status='old')
        close (unit, status='delete')
    end select
end program test_wall_config

subroutine write_file(path, contents)
    character(len=*), intent(in) :: path, contents
    integer :: unit

    open (newunit=unit, file=path, status='replace', form='formatted')
    write (unit, '(A)', advance='no') contents
    close (unit)
end subroutine write_file
