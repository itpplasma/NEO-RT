program test_response_jxb_invalid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use response_jxb, only: cylindrical_toroidal_torque, integrate_mars_profile

    implicit none

    character(16) :: case_name
    complex(dp) :: one(1), two(2)
    real(dp) :: edges(2), radius(1), weight(1), value

    call get_command_argument(1, case_name)
    one = cmplx(1.0_dp, 0.0_dp, dp)
    two = cmplx(1.0_dp, 0.0_dp, dp)
    radius = 1.0_dp
    weight = 1.0_dp
    edges = [0.0_dp, 1.0_dp]

    select case (trim(case_name))
    case ("shape")
        value = cylindrical_toroidal_torque(radius, weight, two, one, one, one)
    case ("profile")
        value = integrate_mars_profile(edges, [1.0_dp, 2.0_dp], 1.0_dp)
    case default
        error stop "Unknown invalid-input test"
    end select
    print *, value
    error stop "Invalid input was accepted"

end program test_response_jxb_invalid
