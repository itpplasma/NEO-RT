program test_resonance_ledger
    use iso_fortran_env, only: dp => real64
    use neort_lib, only: config_t, transport_data_t, mixture_t, species_state_t, &
        mixture_result_t, neort_init, neort_prepare_splines, neort_compute_at_s, &
        neort_prepare_mixture, neort_compute_mixture_at_s
    use omp_lib, only: omp_set_dynamic, omp_set_num_threads, omp_get_num_threads
    implicit none
    type(config_t) :: config
    type(transport_data_t) :: results(2)
    type(mixture_t) :: mixture
    type(mixture_result_t) :: mixture_result
    character(32) :: mode
    character(1024) :: output
    integer :: k, unit, team
    call get_command_argument(1, mode)
    call get_command_argument(2, output)
    config%s = 0.09_dp
    config%M_t = 0.036_dp
    config%magdrift = .true.
    config%epsmn = 1.0e-3_dp
    config%m0 = 0
    config%mph = 3
    config%comptorque = mode /= 'no-torque'
    config%supban = .false.
    config%nopassing = .false.
    config%inp_swi = 8
    config%vsteps = 16
    config%mth_max_abs = 8
    config%vmax_over_vth = 4.0_dp
    config%log_level = -1
    if (mode == 'mono') config%monoenergetic_x = 2.25_dp
    if (mode == 'krook') then
        config%krook_collisions = .true.
        config%krook_nu = 1.0e3_dp
        config%krook_eta_steps = 4
    end if
    call neort_init(config, 'in_file')
    call neort_prepare_splines('plasma.in', 'profile.in')
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(2)
    team = 0
    !$omp parallel
    !$omp single
    team = omp_get_num_threads()
    !$omp end single
    !$omp end parallel
    if (team /= 2) error stop 'ledger test requires two threads'
    open(newunit=unit, file=trim(output), status='replace')
    if (mode == 'species') then
        mixture%reference_s = config%s
        mixture%omega_e = -700.0_dp
        mixture%coulomb_log = 20.0_dp
        allocate(mixture%species(3))
        mixture%species(1) = species_state_t('D', 3.3435837724e-27_dp, 1.0_dp, &
            5.0e15_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
        mixture%species(2) = species_state_t('T', 5.0073567446e-27_dp, 1.0_dp, &
            5.0e15_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
        mixture%species(3) = species_state_t('e', 9.1093837139e-31_dp, -1.0_dp, &
            1.0e16_dp, 2000.0_dp, 0.0_dp, 0.0_dp)
        call neort_prepare_mixture(mixture)
        call neort_compute_mixture_at_s(config%s, mixture_result)
        do k = 1, 3
            call write_expected(unit, k, mixture_result%species(k)%transport)
        end do
    else if (mode == 'parallel') then
        !$omp parallel do private(k)
        do k = 1, 2
            call neort_compute_at_s(config%s + (k - 1)*0.03_dp, results(k))
        end do
        !$omp end parallel do
        do k = 1, 2
            call write_expected(unit, 1, results(k))
        end do
    else
        call neort_compute_at_s(config%s, results(1))
        call write_expected(unit, 1, results(1))
    end if
    close(unit)
contains
    subroutine write_expected(unit, species, result_)
        integer, intent(in) :: unit, species
        type(transport_data_t), intent(in) :: result_
        integer :: j
        do j = 1, size(result_%harmonics)
            write(unit, *) result_%torque%s, species, result_%harmonics(j)%mth, &
                result_%harmonics(j)%Tresco, result_%harmonics(j)%Tresctr, &
                result_%harmonics(j)%Trest
        end do
    end subroutine write_expected
end program test_resonance_ledger
