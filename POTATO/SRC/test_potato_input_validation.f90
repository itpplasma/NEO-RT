program test_potato_input_validation
    use potato_input_mod, only : itest_type, rho_pol, rho_pol_max, &
        potato_input_is_valid
    implicit none

    itest_type = 3
    rho_pol = 0.9d0
    rho_pol_max = 0.9d0
    call require(.not. potato_input_is_valid(), &
        'torque topology domain equal to deposition boundary was accepted')

    rho_pol_max = 0.99d0
    call require(potato_input_is_valid(), &
        'extended torque topology domain was rejected')

    rho_pol_max = 1.d0
    call require(.not. potato_input_is_valid(), &
        'rho_pol_max at the separatrix was accepted')

    itest_type = 1
    rho_pol_max = rho_pol
    call require(potato_input_is_valid(), &
        'non-torque orbit cut at rho_pol was rejected')

contains

    subroutine require(ok, message)
        logical, intent(in) :: ok
        character(len=*), intent(in) :: message

        if (ok) return
        print *, trim(message)
        error stop 1
    end subroutine require

end program test_potato_input_validation
