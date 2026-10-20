program test_mixture_collisions
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use neort_species, only: mixture_t, species_state_t, configure_mixture, &
        evaluate_mixture_at_s, selected_species
    use neort_mixture_collisions
    use neort_krook, only: krook_collisions, krook_nu, krook_model
    implicit none
    type(mixture_t) :: mixture
    character(len=32) :: mode
    real(dp) :: thermal
    call get_command_argument(1, mode)
    if (mode == 'missing_background') then
        krook_rate_source = 1
        call validate_collision_background(.true.)
        stop 0
    end if
    if (mode == 'nonfinite_scale') then
        krook_rate_scale = ieee_value(0.0_dp, ieee_quiet_nan)
        call validate_collision_background()
        stop 0
    end if
    if (mode == 'duplicate_shape') then
        krook_rate_source = 2
        krook_model = 2
        call validate_collision_background()
        stop 0
    end if
    mixture%coulomb_log = 20.0_dp
    allocate(mixture%species(3))
    mixture%species(1) = species_state_t('D', 3.3435837724e-27_dp, 1.0_dp, &
        5.0e15_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
    mixture%species(2) = species_state_t('T', 5.0073567446e-27_dp, 1.0_dp, &
        5.0e15_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
    mixture%species(3) = species_state_t('e', 9.1093837139e-31_dp, -1.0_dp, &
        1.0e16_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
    call configure_mixture(mixture)
    call evaluate_mixture_at_s(0.0_dp)
    selected_species = 3
    call initialize_collision_background()
    thermal = selected_thermal_deflection_rate()
    krook_collisions = .true.
    krook_rate_source = 1
    krook_rate_scale = 2.0_dp
    krook_model = 2
    call close(effective_krook_rate(2.0_dp, 0.0_dp, 3.0_dp), 0.75_dp*thermal)
    krook_model = 1
    call close(effective_krook_rate(2.0_dp, 2.0_dp, 3.0_dp), 1.5_dp*thermal)
    krook_rate_source = 2
    krook_model = 0
    call close(effective_krook_rate(1.0_dp, 9.0_dp, 3.0_dp), 6.0_dp*thermal)
    krook_rate_scale = 0.0_dp
    if (finite_krook_active()) error stop 'zero rate scale did not restore collisionless dispatch'
    krook_rate_source = 0
    krook_nu = 7.0_dp
    call close(effective_krook_rate(1.0_dp, 0.0_dp, 3.0_dp), 21.0_dp)
    krook_collisions = .false.
    if (finite_krook_active()) error stop 'disabled collisions dispatch is not exact'
    print *, 'Mixture rate selection and legacy/zero controls pass'
contains
    subroutine close(actual, expected)
        real(dp), intent(in) :: actual, expected
        if (abs(actual-expected) > 1.0e-12_dp*abs(expected)) &
            error stop 'effective mixture width has wrong scaling or shape'
    end subroutine close
end program test_mixture_collisions
