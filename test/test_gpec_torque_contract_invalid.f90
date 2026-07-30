program test_gpec_torque_contract_invalid
    use gpec_torque_contract, only: gpec_control_torque, &
        read_gpec_control_torque, require_gpec_volume_current

    implicit none

    character(16) :: case_name
    character(512) :: fixture
    type(gpec_control_torque) :: result

    call get_command_argument(1, case_name)
    select case (trim(case_name))
    case ("volume")
        call require_gpec_volume_current
    case ("zero_mode")
        call read_gpec_control_torque("unused", 0, result)
    case ("missing_torque")
        call get_command_argument(2, fixture)
        call read_gpec_control_torque(trim(fixture), 3, result)
    case default
        error stop "unknown invalid GPEC torque case"
    end select
end program test_gpec_torque_contract_invalid
