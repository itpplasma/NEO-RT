program test_species_mixture_run
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_is_finite
    use neort_lib, only: config_t, transport_data_t, mixture_t, species_state_t, &
        mixture_result_t, neort_init, neort_prepare_mixture, neort_compute_mixture_at_s, &
        neort_compute_at_s
    use neort_species, only: selected_species
    use neort_mixture_collisions, only: krook_rate_source, krook_rate_scale
    use neort_krook, only: krook_nu
    use neort_output, only: write_mixture_output
    use omp_lib, only: omp_set_num_threads, omp_set_dynamic
    implicit none
    type(config_t) :: config
    type(mixture_t) :: mixture
    type(mixture_result_t) :: result_
    type(transport_data_t) :: serial_result
    character(len=1024) :: output_path
    real(dp) :: summed(3), expected(3), actual(3)
    integer :: k, rate_source
    character(len=32) :: rate_arg
    call get_command_argument(1, output_path)
    if (len_trim(output_path) == 0) error stop 'test requires an output prefix'
    call get_command_argument(2, rate_arg)
    rate_source = 0
    if (len_trim(rate_arg) > 0) read(rate_arg, *) rate_source
    call setup_config(config)
    config%krook_rate_source = rate_source
    config%krook_rate_scale = 1.0e4_dp
    call neort_init(config, 'in_file')
    call setup_mixture(mixture)
    call neort_prepare_mixture(mixture)
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(3)
    call neort_compute_mixture_at_s(config%s, result_)
    summed = 0.0_dp
    do k = 1, 3
        if (result_%species(k)%species%name /= mixture%species(k)%name) &
            error stop 'parallel species order changed'
        call assert_close('native temperature and common physical rotation', &
            [result_%species(k)%magfie%params%T_in_eV, &
             result_%species(k)%magfie%params%Om_tE], [2000.0_dp, -700.0_dp])
        selected_species = k
        ! Independent rate module tests supply the physical oracle. This
        ! check verifies that native parallel dispatch uses that species rate.
        if (rate_source == 1) then
            krook_rate_source = 0
            krook_nu = result_%species(k)%collision_thermal_deflection_rate*krook_rate_scale
        end if
        call neort_compute_at_s(config%s, serial_result)
        actual = [result_%species(k)%transport%torque%Tco, &
            result_%species(k)%transport%torque%Tctr, result_%species(k)%transport%torque%Tt]
        expected = [serial_result%torque%Tco, serial_result%torque%Tctr, serial_result%torque%Tt]
        call assert_close('parallel versus isolated native species solve', actual, expected)
        if (.not. all(ieee_is_finite(actual))) error stop 'native species torque is nonfinite'
        if (abs(actual(3)) <= tiny(1.0_dp)) error stop 'finite-width trapped torque is absent'
        summed = summed + actual
    end do
    call assert_close('total physical torque', result_%torque_total, summed)
    call write_mixture_output(result_, trim(output_path), 'both')
    print *, 'Native parallel D/T/e mixture matches isolated species solves'
contains
    subroutine setup_config(value)
        type(config_t), intent(out) :: value
        value%s = 0.09_dp
        value%epsmn = 1.0e-3_dp
        value%m0 = 0
        value%mph = 3
        value%comptorque = .true.
        value%supban = .true.
        value%nopassing = .true.
        value%inp_swi = 8
        value%vsteps = 2
        value%mth_max_abs = 0
        value%vmax_over_vth = 1.0_dp
        value%log_level = -1
        value%krook_collisions = .true.
        value%krook_nu = 1.0e3_dp
        value%krook_eta_steps = 4
        value%krook_rtol = 1.0e-4_dp
    end subroutine setup_config

    subroutine setup_mixture(value)
        type(mixture_t), intent(out) :: value
        value%reference_s = 0.09_dp
        value%omega_e = -700.0_dp
        value%coulomb_log = 20.0_dp
        allocate(value%species(3))
        value%species(1) = species_state_t('D', 3.3435837724e-27_dp, 1.0_dp, &
            5.0e15_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
        value%species(2) = species_state_t('T', 5.0073567446e-27_dp, 1.0_dp, &
            5.0e15_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
        value%species(3) = species_state_t('e', 9.1093837139e-31_dp, -1.0_dp, &
            1.0e16_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
    end subroutine setup_mixture

    subroutine assert_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:), expected(:)
        if (any(abs(actual - expected) > &
            1.0e-11_dp*max(abs(expected), 1.0e-100_dp))) then
            print *, trim(label), actual, expected
            error stop 'native mixture regression mismatch'
        end if
    end subroutine assert_close
end program test_species_mixture_run
