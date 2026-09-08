program neort_diag
    use iso_fortran_env, only: real64
    use diag_bounce_nonlin, only: run_bounce_nonlin_diag
    use diag_atten_map,    only: run_atten_map_diag
    use diag_contrib_map,  only: run_contrib_diag
    use diag_bounce_debug, only: run_bounce_debug
    use diag_action_trace, only: run_action_trace_diag
    use diag_orbit_trace, only: run_orbit_trace_diag
    use diag_resonance_contour, only: run_resonance_contour_diag
    use diag_resonance_scan, only: run_resonance_scan_diag
    use diag_pitch_action, only: run_pitch_action_diag, run_pitch_coeff_diag
    implicit none
    character(len=256) :: diag, runname, ux_arg, eta_arg, nsteps_arg, &
        mth_arg, neta_arg, surface_file_arg
    real(real64) :: ux, eta
    integer :: ios, neta, nsteps, mth_value
    call get_command_argument(1, diag)
    call get_command_argument(2, runname)
    if (len_trim(diag) == 0 .or. len_trim(runname) == 0) then
        print *, "Usage: neo_rt_diag.x <diagnostic> <runname> [ux] [neta] [surface_file]"
        print *, "Diagnostics: bounce_nonlin, atten_map, contrib, bounce_debug, action_trace, orbit_trace, resonance_contour, resonance_scan"
        print *, "pitch_action <runname> <points_file>: off-root frequency and action"
        print *, "pitch_action_tight <runname> <points_file>: diagnostic tight-tolerance action"
        print *, "pitch_action_ultra <runname> <points_file>: diagnostic ultra-tight action"
        print *, "pitch_coeff <runname> <points_file>: native coefficient/event table"
        print *, "orbit_trace <runname> <ux> <eta> [nsteps] [mth]: ordered complex orbit packet"
        stop 1
    end if

    select case (trim(adjustl(diag)))
    case ("bounce_nonlin")
        call run_bounce_nonlin_diag(trim(adjustl(runname)))
    case ("atten_map")
        call run_atten_map_diag(trim(adjustl(runname)))
    case ("contrib")
        call run_contrib_diag(trim(adjustl(runname)))
    case ("bounce_debug")
        call run_bounce_debug(trim(adjustl(runname)))
    case ("action_trace")
        call run_action_trace_diag(trim(adjustl(runname)))
    case ("orbit_trace")
        ux = 1.0_real64
        eta = 1.0_real64
        nsteps = 257
        mth_value = 0
        ios = 0
        call get_command_argument(3, ux_arg)
        if (len_trim(ux_arg) > 0) read(ux_arg, *, iostat=ios) ux
        if (len_trim(ux_arg) > 0 .and. ios /= 0) error stop "UX must be a real number"
        ios = 0
        call get_command_argument(4, eta_arg)
        if (len_trim(eta_arg) > 0) read(eta_arg, *, iostat=ios) eta
        if (len_trim(eta_arg) > 0 .and. ios /= 0) error stop "ETA must be a real number"
        ios = 0
        call get_command_argument(5, nsteps_arg)
        if (len_trim(nsteps_arg) > 0) read(nsteps_arg, *, iostat=ios) nsteps
        if (len_trim(nsteps_arg) > 0 .and. ios /= 0) error stop "NSTEPS must be an integer"
        ios = 0
        call get_command_argument(6, mth_arg)
        if (len_trim(mth_arg) > 0) read(mth_arg, *, iostat=ios) mth_value
        if (len_trim(mth_arg) > 0 .and. ios /= 0) error stop "MTH must be an integer"
        call run_orbit_trace_diag(trim(adjustl(runname)), ux, eta, nsteps, mth_value)
    case ("resonance_contour")
        call run_resonance_contour_diag(trim(adjustl(runname)))
    case ("resonance_scan")
        ux = 1.5_real64
        ios = 0
        call get_command_argument(3, ux_arg)
        if (len_trim(ux_arg) > 0) read(ux_arg, *, iostat=ios) ux
        if (len_trim(ux_arg) > 0 .and. ios /= 0) error stop "UX must be a real number"
        neta = 180
        ios = 0
        call get_command_argument(4, neta_arg)
        if (len_trim(neta_arg) > 0) read(neta_arg, *, iostat=ios) neta
        if (len_trim(neta_arg) > 0 .and. ios /= 0) error stop "NETA must be an integer"
        call get_command_argument(5, surface_file_arg)
        if (len_trim(surface_file_arg) > 0) then
            call run_resonance_scan_diag(trim(adjustl(runname)), ux, neta, trim(adjustl(surface_file_arg)))
        else
            call run_resonance_scan_diag(trim(adjustl(runname)), ux, neta)
        end if
    case ("pitch_action")
        call get_command_argument(3, surface_file_arg)
        if (len_trim(surface_file_arg) == 0) &
            error stop "pitch_action requires a point-file argument"
        call run_pitch_action_diag(trim(adjustl(runname)), trim(surface_file_arg))
    case ("pitch_action_tight")
        call get_command_argument(3, surface_file_arg)
        if (len_trim(surface_file_arg) == 0) &
            error stop "pitch_action_tight requires a point-file argument"
        call run_pitch_action_diag(trim(adjustl(runname)), trim(surface_file_arg), .true.)
    case ("pitch_action_ultra")
        call get_command_argument(3, surface_file_arg)
        if (len_trim(surface_file_arg) == 0) &
            error stop "pitch_action_ultra requires a point-file argument"
        call run_pitch_action_diag(trim(adjustl(runname)), trim(surface_file_arg), .false., .true.)
    case ("pitch_coeff")
        call get_command_argument(3, surface_file_arg)
        if (len_trim(surface_file_arg) == 0) &
            error stop "pitch_coeff requires a point-file argument"
        call run_pitch_coeff_diag(trim(adjustl(runname)), trim(surface_file_arg))
    case default
        print *, "Unknown diagnostic:", trim(diag)
        stop 2
    end select

end program neort_diag
