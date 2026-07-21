program test_input_switches
    use neort_config, only: config_t, set_config, read_and_set_config
    use do_magfie_mod, only: axisymmetric_switch => inp_swi
    use do_magfie_pert_mod, only: perturbation_switch => inp_swi_pert

    implicit none

    type(config_t) :: config
    integer :: unit

    if (config%inp_swi_pert /= -1) error stop "perturbation switch default changed"

    config%inp_swi = 10
    config%inp_swi_pert = 9
    config%pertfile = .true.
    call set_config(config)
    if (axisymmetric_switch /= 10) error stop "axisymmetric switch not set"
    if (perturbation_switch /= 9) error stop "perturbation switch not set"

    config = config_t()
    config%inp_swi = 9
    config%pertfile = .true.
    call set_config(config)
    if (perturbation_switch /= 9) error stop "inheritance from inp_swi failed"

    call write_namelist("input_switches.in", .true.)
    call read_and_set_config("input_switches.in")
    if (axisymmetric_switch /= 10) error stop "namelist axisymmetric switch not read"
    if (perturbation_switch /= 9) error stop "namelist perturbation switch not read"

    call write_namelist("input_switches.in", .false.)
    call read_and_set_config("input_switches.in")
    if (axisymmetric_switch /= 9) error stop "namelist inherited axisymmetric switch changed"
    if (perturbation_switch /= 9) error stop "namelist inheritance failed"

    open (newunit=unit, file="input_switches.in", status="old")
    close (unit, status="delete")

contains

    subroutine write_namelist(path, mixed)
        character(len=*), intent(in) :: path
        logical, intent(in) :: mixed
        integer :: u

        open (newunit=u, file=path, status="replace", form="formatted")
        write (u, '(A)') "&params"
        write (u, '(A)') "    qs = 1.0"
        write (u, '(A)') "    ms = 2.0"
        write (u, '(A)') "    pertfile = .true."
        if (mixed) then
            write (u, '(A)') "    inp_swi = 10"
            write (u, '(A)') "    inp_swi_pert = 9"
        else
            write (u, '(A)') "    inp_swi = 9"
        end if
        write (u, '(A)') "/"
        close (u)
    end subroutine write_namelist

end program test_input_switches
