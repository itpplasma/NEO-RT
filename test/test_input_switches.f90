program test_input_switches
    use neort_config, only: config_t, set_config, read_and_set_config, &
        perturbation_chart_is_compatible, frequency_model_requires_direct_eqdsk
    use do_magfie_mod, only: axisymmetric_switch => inp_swi
    use do_magfie_pert_mod, only: perturbation_switch => inp_swi_pert
    use driftorbit, only: FREQUENCY_MODEL_BOOZER_THIN, FREQUENCY_MODEL_GC_FULL, &
        configured_frequency_model => frequency_model

    implicit none

    type(config_t) :: config
    character(len=32) :: mode
    integer :: unit

    call get_command_argument(1, mode)
    if (trim(mode) == "reject_boozer_rz") then
        config%inp_swi = 9
        config%inp_swi_pert = 11
        config%pertfile = .true.
        call set_config(config)
        stop
    end if
    if (trim(mode) == "reject_full_boozer") then
        config%frequency_model = FREQUENCY_MODEL_GC_FULL
        config%inp_swi = 10
        call set_config(config)
        stop
    end if
    if (trim(mode) == "reject_full_supban") then
        config%frequency_model = FREQUENCY_MODEL_GC_FULL
        config%inp_swi = 11
        config%supban = .true.
        call set_config(config)
        stop
    end if
    if (trim(mode) == "reject_full_nonlin") then
        config%frequency_model = FREQUENCY_MODEL_GC_FULL
        config%inp_swi = 11
        config%nonlin = .true.
        call set_config(config)
        stop
    end if
    if (trim(mode) == "reject_model1") then
        config%frequency_model = 1
        call set_config(config)
        stop
    end if

    if (config%inp_swi_pert /= -1) error stop "perturbation switch default changed"
    if (perturbation_chart_is_compatible(11, 9)) then
        error stop "direct GEQDSK accepted a Boozer-angle perturbation"
    end if
    if (.not. perturbation_chart_is_compatible(11, 11)) then
        error stop "direct GEQDSK rejected an R-Z perturbation"
    end if
    if (.not. perturbation_chart_is_compatible(10, 9)) then
        error stop "Boozer chartmap rejected a Boozer perturbation"
    end if
    if (frequency_model_requires_direct_eqdsk(FREQUENCY_MODEL_BOOZER_THIN)) then
        error stop "Boozer thin frequency model unexpectedly requires direct GEQDSK"
    end if
    if (frequency_model_requires_direct_eqdsk(1)) then
        error stop "removed frequency model unexpectedly requires direct GEQDSK"
    end if
    if (.not. frequency_model_requires_direct_eqdsk(FREQUENCY_MODEL_GC_FULL)) then
        error stop "GC full frequency model lacks its direct-GEQDSK gate"
    end if

    call write_wall("input_switches_wall.dat")

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

    config = config_t()
    config%frequency_model = FREQUENCY_MODEL_GC_FULL
    config%inp_swi = 11
    config%wall_file = "input_switches_wall.dat"
    config%wall_units = "m"
    call set_config(config)
    if (configured_frequency_model /= FREQUENCY_MODEL_GC_FULL) then
        error stop "GC full frequency model was not accepted for direct GEQDSK"
    end if

    call write_namelist("input_switches.in", .true.)
    call read_and_set_config("input_switches.in")
    if (axisymmetric_switch /= 10) error stop "namelist axisymmetric switch not read"
    if (perturbation_switch /= 9) error stop "namelist perturbation switch not read"

    call write_namelist("input_switches.in", .false., FREQUENCY_MODEL_GC_FULL)
    call read_and_set_config("input_switches.in")
    if (axisymmetric_switch /= 11) error stop "namelist direct axisymmetric switch not read"
    if (perturbation_switch /= 11) error stop "namelist direct perturbation switch not inherited"
    if (configured_frequency_model /= FREQUENCY_MODEL_GC_FULL) then
        error stop "namelist GC full frequency model was not read"
    end if

    open (newunit=unit, file="input_switches.in", status="old")
    close (unit, status="delete")
    open (newunit=unit, file="input_switches_wall.dat", status="old")
    close (unit, status="delete")

contains

    subroutine write_namelist(path, mixed, selected_model)
        character(len=*), intent(in) :: path
        logical, intent(in) :: mixed
        integer, intent(in), optional :: selected_model
        integer :: u

        open (newunit=u, file=path, status="replace", form="formatted")
        write (u, '(A)') "&params"
        write (u, '(A)') "    qs = 1.0"
        write (u, '(A)') "    ms = 2.0"
        if (present(selected_model)) then
            if (selected_model == FREQUENCY_MODEL_GC_FULL) then
                write (u, '(A)') "    pertfile = .false."
            else
                write (u, '(A)') "    pertfile = .true."
            end if
        else
            write (u, '(A)') "    pertfile = .true."
        end if
        if (mixed) then
            write (u, '(A)') "    inp_swi = 10"
            write (u, '(A)') "    inp_swi_pert = 9"
        else if (present(selected_model)) then
            write (u, '(A)') "    inp_swi = 11"
        else
            write (u, '(A)') "    inp_swi = 9"
        end if
        if (present(selected_model)) then
            write (u, '(A,I0)') "    frequency_model = ", selected_model
            if (selected_model == FREQUENCY_MODEL_GC_FULL) then
                write (u, '(A)') '    wall_file = "input_switches_wall.dat"'
                write (u, '(A)') '    wall_units = "m"'
            end if
        end if
        write (u, '(A)') "/"
        close (u)
    end subroutine write_namelist

    subroutine write_wall(path)
        character(len=*), intent(in) :: path
        integer :: u

        open (newunit=u, file=path, status="replace", form="formatted")
        write (u, '(A)') "1.0 0.0"
        write (u, '(A)') "2.0 0.0"
        write (u, '(A)') "2.0 1.0"
        write (u, '(A)') "1.0 1.0"
        close (u)
    end subroutine write_wall

end program test_input_switches
