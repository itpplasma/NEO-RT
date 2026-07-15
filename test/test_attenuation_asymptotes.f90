program test_attenuation_asymptotes
    use iso_fortran_env, only: dp => real64
    use thetadata_mod, only: init_attenuation_data
    implicit none

    real(dp) :: theta

    call init_attenuation_data()

    call attenuation_factor(0.0_dp, theta)
    call assert_close("zero-D limit", theta, 0.0_dp, 0.0_dp)

    call attenuation_factor(1.0e-8_dp, theta)
    call assert_close("low-D asymptote", theta, 1.7555267e-8_dp, 1.0e-14_dp)

    call attenuation_factor(1.0e8_dp, theta)
    call assert_close("high-D asymptote", theta, 1.0_dp, 0.0_dp)

contains

    subroutine assert_close(name, actual, expected, tolerance)
        character(*), intent(in) :: name
        real(dp), intent(in) :: actual, expected, tolerance

        if (abs(actual - expected) > tolerance) then
            print *, trim(name), " failed: ", actual, expected
            error stop 1
        end if
    end subroutine assert_close

end program test_attenuation_asymptotes
