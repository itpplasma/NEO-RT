program test_electromagnetic_torque_invalid
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use electromagnetic_torque, only: staggered_cross_contraction

    implicit none

    character(8) :: case_name
    complex(dp) :: half_one(1, 1), half_two(2, 1), full_two(2, 1), full_three(3, 1)
    real(dp), allocatable :: value(:)

    call get_command_argument(1, case_name)
    half_one = cmplx(1.0_dp, 0.0_dp, dp)
    half_two = cmplx(1.0_dp, 0.0_dp, dp)
    full_two = cmplx(1.0_dp, 0.0_dp, dp)
    full_three = cmplx(1.0_dp, 0.0_dp, dp)

    select case (trim(case_name))
    case ("half")
        value = staggered_cross_contraction(half_one, half_two, full_two, full_two)
    case ("full")
        value = staggered_cross_contraction(half_one, half_one, full_three, full_two)
    case default
        error stop "unknown invalid electromagnetic-torque case"
    end select
    print *, value
    error stop "invalid electromagnetic-torque input was accepted"

end program test_electromagnetic_torque_invalid
