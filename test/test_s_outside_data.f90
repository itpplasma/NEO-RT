program test_s_outside_data
    ! examples/circ.bc tabulates 63 surfaces, s = 2.3438e-2 ... 9.9219e-1 as printed.
    ! do_magfie clamps outside that range, so a request below the first
    ! surface used to return the first surface's field silently (identical
    ! eps, iota, psi_pr for every s < 0.0234375). It must be rejected.
    ! Without arguments the range predicate is checked; with argument "init"
    ! init_magfie_at_s is called at s below the data, which must error stop.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: init_magfie_at_s, inp_swi, magfie_thread_init, &
        read_boozer_file, s_outside_data, set_s
    use logger, only: set_log_level

    implicit none

    real(dp), parameter :: s_first = 2.3438e-2_dp, s_last = 9.9219e-1_dp
    character(len=16) :: arg

    call set_log_level(-1)
    inp_swi = 8
    call magfie_thread_init()
    call read_boozer_file("circ.bc")

    if (command_argument_count() > 0) then
        call get_command_argument(1, arg)
        if (trim(arg) == "init") then
            call set_s(0.0117_dp)
            call init_magfie_at_s()
            print *, "FAIL: init_magfie_at_s accepted s below the data"
            stop 0
        end if
    end if

    call expect(s_outside_data(0.0117_dp), "s=0.0117 below data")
    call expect(s_outside_data(0.0195_dp), "s=0.0195 below data")
    call expect(s_outside_data(0.9999_dp), "s=0.9999 above data")
    call expect(.not. s_outside_data(s_first), "first surface is inside")
    call expect(.not. s_outside_data(s_last), "last surface is inside")
    call expect(.not. s_outside_data(0.039063_dp), "s=0.039063 inside")
    call set_s(0.039063_dp)
    call init_magfie_at_s()
    print *, "PASS"

contains

    subroutine expect(ok, label)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: label
        if (.not. ok) then
            print *, "FAIL: ", label
            error stop 1
        end if
    end subroutine expect

end program test_s_outside_data
