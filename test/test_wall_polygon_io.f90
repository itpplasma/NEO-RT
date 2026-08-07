program test_wall_polygon_io
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_config, only: config_t, configured_wall_file, read_and_set_config, &
        set_config, validate_frequency_model
    use neort_wall_io, only: WALL_IO_INVALID_UNITS, WALL_IO_MALFORMED_ROW, &
        WALL_IO_MISSING_FILE, WALL_IO_NONFINITE, WALL_IO_NONPOSITIVE_RADIUS, &
        WALL_IO_OK, WALL_IO_SELF_INTERSECTION, WALL_IO_TOO_FEW_VERTICES, &
        WALL_IO_TOO_LONG_RECORD, WALL_IO_ZERO_AREA, &
        load_wall_polygon, wall_polygon_t

    implicit none

    type(wall_polygon_t) :: wall
    type(config_t) :: config
    integer :: status, unit, hash_index
    character(len=128) :: message

    call write_file('wall_valid.dat', &
        '  # a closed polygon is implicit'//new_line('a')// &
        '1.0 0.0'//new_line('a')// &
        '2.0 0.0 ! inline comment'//new_line('a')// &
        ''//new_line('a')// &
        '2.0 1.0'//new_line('a')// &
        '1.0 1.0')
    call load_wall_polygon('wall_valid.dat', wall, status, message, 'cm')
    if (status /= WALL_IO_OK) error stop 'valid wall polygon was rejected'
    if (.not. allocated(wall%vertices)) error stop 'wall vertices were not allocated'
    if (size(wall%vertices, 1) /= 2) error stop 'wall coordinate dimension changed'
    if (size(wall%vertices, 2) /= 4) error stop 'wall vertex count changed'
    if (abs(wall%vertices(1, 3) - 2.0_dp) > 1.0e-14_dp) then
        error stop 'wall R coordinate was not loaded'
    end if
    if (abs(wall%vertices(2, 4) - 1.0_dp) > 1.0e-14_dp) then
        error stop 'wall Z coordinate was not loaded'
    end if
    if (len_trim(wall%hash) /= 64 .or. any([(index('0123456789ABCDEF', &
            wall%hash(hash_index:hash_index)) == 0, hash_index = 1, 64)])) then
        error stop 'wall hash was not a 64-hex SHA-256 digest'
    end if
    block
        character(len=64) :: centimeter_hash
        centimeter_hash = wall%hash
        call load_wall_polygon('wall_valid.dat', wall, status, message, 'm')
        if (status /= WALL_IO_OK .or. trim(wall%hash) == trim(centimeter_hash)) then
            error stop 'wall hash did not distinguish normalized geometry'
        end if
    end block
    call load_wall_polygon('wall_valid.dat', wall, status, message, 'm')
    if (status /= WALL_IO_OK .or. abs(wall%vertices(1, 1) - 100.0_dp) > &
            1.0e-12_dp) error stop 'metre to centimetre wall conversion failed'

    call write_file('wall_malformed.dat', '1.0 0.0'//new_line('a')//'2.0')
    call load_wall_polygon('wall_malformed.dat', wall, status, input_units='cm')
    if (status /= WALL_IO_MALFORMED_ROW) error stop 'malformed wall row was accepted'

    call write_file('wall_extra_column.dat', '1.0 0.0 7.0'//new_line('a')// &
        '2.0 0.0'//new_line('a')//'2.0 1.0')
    call load_wall_polygon('wall_extra_column.dat', wall, status, input_units='cm')
    if (status /= WALL_IO_MALFORMED_ROW) error stop 'extra wall column was accepted'

    call write_file('wall_too_few.dat', '1.0 0.0'//new_line('a')//'2.0 0.0')
    call load_wall_polygon('wall_too_few.dat', wall, status, input_units='cm')
    if (status /= WALL_IO_TOO_FEW_VERTICES) error stop 'short polygon was accepted'

    call write_file('wall_nonfinite.dat', 'NaN 0.0'//new_line('a')//'2.0 0.0'// &
        new_line('a')//'2.0 1.0')
    call load_wall_polygon('wall_nonfinite.dat', wall, status, input_units='cm')
    if (status /= WALL_IO_NONFINITE) error stop 'nonfinite wall coordinate was accepted'

    call load_wall_polygon('wall_file_that_does_not_exist.dat', wall, status, &
        input_units='cm')
    if (status /= WALL_IO_MISSING_FILE) error stop 'missing wall file was accepted'

    call write_file('wall_invalid_units.dat', '1.0 0.0'//new_line('a')// &
        '2.0 0.0'//new_line('a')//'2.0 1.0')
    call load_wall_polygon('wall_invalid_units.dat', wall, status, &
        input_units='inch')
    if (status /= WALL_IO_INVALID_UNITS) error stop 'invalid wall units were accepted'

    call write_file('wall_duplicates.dat', '1.0 0.0'//new_line('a')// &
        '2.0 0.0'//new_line('a')//'2.0 1.0'//new_line('a')//'1.0 1.0'// &
        new_line('a')//'1.0 1.0'//new_line('a')//'1.0 0.0')
    call load_wall_polygon('wall_duplicates.dat', wall, status, &
        input_units='cm')
    if (status /= WALL_IO_OK .or. size(wall%vertices, 2) /= 4) then
        error stop 'wall duplicate/closing normalization failed'
    end if

    call write_file('wall_zero_area.dat', '1.0 0.0'//new_line('a')// &
        '2.0 0.0'//new_line('a')//'3.0 0.0')
    call load_wall_polygon('wall_zero_area.dat', wall, status, input_units='cm')
    if (status /= WALL_IO_ZERO_AREA) error stop 'zero-area wall was accepted'

    call write_file('wall_nonpositive_radius.dat', '0.0 0.0'//new_line('a')// &
        '2.0 0.0'//new_line('a')//'2.0 1.0')
    call load_wall_polygon('wall_nonpositive_radius.dat', wall, status, &
        input_units='cm')
    if (status /= WALL_IO_NONPOSITIVE_RADIUS) then
        error stop 'nonpositive-R wall was accepted'
    end if

    call write_file('wall_self_intersecting.dat', '1.0 0.0'//new_line('a')// &
        '3.0 2.0'//new_line('a')//'1.0 2.0'//new_line('a')//'3.0 0.0'// &
        new_line('a')//'2.0 3.0')
    call load_wall_polygon('wall_self_intersecting.dat', wall, status, &
        input_units='cm')
    if (status /= WALL_IO_SELF_INTERSECTION) then
        error stop 'self-intersecting wall was accepted'
    end if

    call write_file('wall_overlong.dat', repeat('x', 5000))
    call load_wall_polygon('wall_overlong.dat', wall, status, input_units='cm')
    if (status /= WALL_IO_TOO_LONG_RECORD) error stop 'overlong wall row was accepted'

    config = config_t()
    config%frequency_model = 0
    config%wall_file = 'wall_file_that_does_not_exist.dat'
    call set_config(config)

    call validate_frequency_model(1, 11, .true., .false., .false., &
        'wall_file_that_does_not_exist.dat')

    call write_file('wall_config.in', '&params'//new_line('a')// &
        '    frequency_model = 0'//new_line('a')// &
        '    wall_file = "wall_valid.dat"'//new_line('a')//'/')
    call read_and_set_config('wall_config.in')
    if (trim(configured_wall_file) /= 'wall_valid.dat') then
        error stop 'wall_file was not exposed through the params namelist'
    end if

    config = config_t()
    config%frequency_model = 2
    config%inp_swi = 11
    config%wall_file = 'wall_valid.dat'
    config%wall_units = 'm'
    call set_config(config)
    if (trim(configured_wall_file) /= 'wall_valid.dat') then
        error stop 'wall_file was not transferred from config_t'
    end if

    open (newunit=unit, file='wall_valid.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_malformed.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_extra_column.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_too_few.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_nonfinite.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_invalid_units.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_duplicates.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_zero_area.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_nonpositive_radius.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_self_intersecting.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_overlong.dat', status='old')
    close (unit, status='delete')
    open (newunit=unit, file='wall_config.in', status='old')
    close (unit, status='delete')
end program test_wall_polygon_io

subroutine write_file(path, contents)
    character(len=*), intent(in) :: path, contents
    integer :: unit

    open (newunit=unit, file=path, status='replace', form='formatted')
    write (unit, '(A)', advance='no') contents
    close (unit)
end subroutine write_file
