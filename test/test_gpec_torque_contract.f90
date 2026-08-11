program test_gpec_torque_contract
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gpec_torque_contract, only: gpec_control_torque, read_gpec_control_torque

    implicit none

    type(gpec_control_torque) :: result
    character(512) :: fixture

    call get_command_argument(1, fixture)
    if (len_trim(fixture) == 0) error stop "missing GPEC_CONTROL fixture"

    call read_gpec_control_torque(trim(fixture), 3, result)
    if (result%toroidal_mode /= 3) error stop "wrong GPEC toroidal mode"
    if (abs(result%torque_nm - (-1.23093339_dp)) > 1.0e-12_dp) then
        error stop "wrong native GPEC toroidal torque"
    end if
end program test_gpec_torque_contract
