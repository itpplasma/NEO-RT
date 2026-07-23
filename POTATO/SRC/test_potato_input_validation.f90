program test_potato_input_validation
    use potato_input_mod, only : itest_type, rho_pol, rho_pol_max, &
        eps_sampling, itermax_sampling, class_eps_sampling, &
        class_itermax_sampling, potato_input_is_valid
    implicit none

    itest_type = 3
    rho_pol = 0.9d0
    rho_pol_max = 0.9d0
    call require(.not. potato_input_is_valid(), &
        'torque topology domain equal to deposition boundary was accepted')

    rho_pol_max = 0.99d0
    call require(potato_input_is_valid(), &
        'extended torque topology domain was rejected')

    class_itermax_sampling = 0
    call require(.not. potato_input_is_valid(), &
        'nonpositive class sampler iteration limit was accepted')
    class_itermax_sampling = 20

    class_eps_sampling = 0.d0
    call require(.not. potato_input_is_valid(), &
        'nonpositive class sampler tolerance was accepted')
    class_eps_sampling = 1.d-3

    itermax_sampling = 0
    call require(.not. potato_input_is_valid(), &
        'nonpositive outer sampler iteration limit was accepted')
    itermax_sampling = 5

    eps_sampling = 0.d0
    call require(.not. potato_input_is_valid(), &
        'nonpositive outer sampler tolerance was accepted')
    eps_sampling = 1.d-2

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
