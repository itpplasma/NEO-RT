program neort_field_probe
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: Bphcov, Bthcov, bfac, do_magfie, &
        init_magfie_at_s, inp_swi, iota, magfie_thread_init, &
        psi_pr, read_boozer_file, s_outside_data, set_s
    use logger, only: set_log_level

    implicit none

    character(len=4096) :: bc_path, points_path, output_path, argument, line
    integer :: input_unit, output_unit, status, line_number, rows
    real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)

    if (command_argument_count() /= 4) then
        write (error_unit, '(a)') &
            'usage: neo_rt_field_probe.x BCF INPUT_SWITCH POINTS OUTPUT'
        error stop 2
    end if
    call get_command_argument(1, bc_path)
    call get_command_argument(2, argument)
    read (argument, *, iostat=status) inp_swi
    if (status /= 0) error stop 'invalid INPUT_SWITCH'
    if (inp_swi /= 8 .and. inp_swi /= 9) then
        error stop 'field probe requires Boozer input switch 8 or 9'
    end if
    call get_command_argument(3, points_path)
    call get_command_argument(4, output_path)
    call set_log_level(-1)
    bfac = 1.0_dp
    call magfie_thread_init()
    call read_boozer_file(trim(bc_path))

    open (newunit=input_unit, file=trim(points_path), status='old', &
        action='read', iostat=status)
    if (status /= 0) error stop 'cannot open point list'
    open (newunit=output_unit, file=trim(output_path), status='replace', &
        action='write', iostat=status)
    if (status /= 0) error stop 'cannot open output'
    write (output_unit, '(a)') '# bc_path='//trim(bc_path)
    write (output_unit, '(a,i0)') '# inp_swi=', inp_swi
    write (output_unit, '(a)') &
        '# Native coordinates (s_tor,phi_B,theta_B); angles rad; bfac=1.'
    write (output_unit, '(a)') &
        '# SI conversion: B/1e4, Bcov/1e6, sqrtg/1e6, psi_pr/1e8.'
    write (output_unit, '(a)') &
        's_tor,phi_B_rad,theta_B_rad,Bmod_G,Bphcov_Gcm,Bthcov_Gcm,'// &
        'iota,sqrtg_cm3,psi_pr_Gcm2'
    rows = 0
    line_number = 0
    do
        read (input_unit, '(a)', iostat=status) line
        if (status < 0) exit
        if (status /= 0) error stop 'cannot read point list'
        line_number = line_number + 1
        line = adjustl(line)
        if (len_trim(line) == 0) cycle
        if (line(1:1) == '#') cycle
        read (line, *, iostat=status) x
        if (status /= 0) then
            write (error_unit, '(a,i0)') 'invalid point at line ', line_number
            error stop 3
        end if
        if (.not. all(ieee_is_finite(x))) error stop 'nonfinite coordinate'
        if (s_outside_data(x(1))) error stop 'point outside Boozer radial range'
        call set_s(x(1))
        call init_magfie_at_s()
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        if (.not. all(ieee_is_finite([bmod, Bphcov, Bthcov, &
                iota, sqrtg, psi_pr]))) error stop 'nonfinite field result'
        if (bmod <= 0.0_dp) error stop 'nonpositive field magnitude'
        write (output_unit, '(8(es25.16e3,","),es25.16e3)') &
            x, bmod, Bphcov, Bthcov, iota, sqrtg, psi_pr
        rows = rows + 1
    end do
    close (input_unit)
    close (output_unit)
    if (rows == 0) error stop 'empty point list'
end program neort_field_probe
