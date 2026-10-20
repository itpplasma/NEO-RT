program neort_diag
    use iso_fortran_env, only: real64
    use diag_bounce_nonlin, only: run_bounce_nonlin_diag
    use diag_atten_map,    only: run_atten_map_diag
    use diag_contrib_map,  only: run_contrib_diag
    use diag_bounce_debug, only: run_bounce_debug
    use diag_resonance_contour, only: run_resonance_contour_diag
    use diag_resonance_scan, only: run_resonance_scan_diag
    implicit none
    character(len=256) :: diag, runname, ux_arg
    real(real64) :: ux
    integer :: ios
    call get_command_argument(1, diag)
    call get_command_argument(2, runname)
    if (len_trim(diag) == 0 .or. len_trim(runname) == 0) then
        print *, "Usage: neo_rt_diag.x <diagnostic> <runname>"
        print *, "Diagnostics: bounce_nonlin, atten_map, contrib, bounce_debug, resonance_contour, resonance_scan"
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
    case ("resonance_contour")
        call run_resonance_contour_diag(trim(adjustl(runname)))
    case ("resonance_scan")
        ux = 1.5_real64
        ios = 0
        call get_command_argument(3, ux_arg)
        if (len_trim(ux_arg) > 0) read(ux_arg, *, iostat=ios) ux
        if (len_trim(ux_arg) > 0 .and. ios /= 0) error stop "UX must be a real number"
        call run_resonance_scan_diag(trim(adjustl(runname)), ux)
    case default
        print *, "Unknown diagnostic:", trim(diag)
        stop 2
    end select

end program neort_diag
