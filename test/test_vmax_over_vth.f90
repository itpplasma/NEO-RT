program test_vmax_over_vth
    use iso_fortran_env, only: dp => real64
    use neort, only: configured => vmax_over_vth
    use neort_config, only: config_t, set_config, read_and_set_config

    implicit none

    type(config_t) :: config
    integer :: unit

    ! config_t default preserves the historical hard-coded cutoff
    if (config%vmax_over_vth /= 3.0_dp) error stop "config_t default changed"

    ! set_config transfers an explicit value to the module variable
    config%vmax_over_vth = 4.5_dp
    call set_config(config)
    if (configured /= 4.5_dp) error stop "config_t did not set vmax_over_vth"

    config = config_t()
    call set_config(config)
    if (configured /= 3.0_dp) error stop "config_t default not restored"

    ! namelist parsing picks up an explicit value
    call write_namelist("vmax.in", "    vmax_over_vth = 5.0")
    call read_and_set_config("vmax.in")
    if (configured /= 5.0_dp) error stop "namelist value not read"

    ! namelist parsing falls back to the default when the field is absent
    call write_namelist("vmax.in", "")
    call read_and_set_config("vmax.in")
    if (configured /= 3.0_dp) error stop "namelist default changed"

    open (newunit=unit, file="vmax.in", status="old")
    close (unit, status="delete")

contains

    subroutine write_namelist(path, extra_line)
        character(len=*), intent(in) :: path, extra_line
        integer :: u

        open (newunit=u, file=path, status="replace", form="formatted")
        write (u, '(A)') "&params"
        write (u, '(A)') "    qs = 1.0"
        write (u, '(A)') "    ms = 2.0"
        if (len_trim(extra_line) > 0) write (u, '(A)') trim(extra_line)
        write (u, '(A)') "/"
        close (u)
    end subroutine write_namelist

end program test_vmax_over_vth
