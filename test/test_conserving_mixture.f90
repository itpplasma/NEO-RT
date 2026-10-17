program test_conserving_mixture
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_lib, only: config_t, mixture_t, species_state_t, mixture_result_t, &
        neort_init, neort_prepare_mixture, neort_compute_mixture_at_s
    use neort_output, only: write_mixture_output
    use omp_lib, only: omp_set_dynamic, omp_set_num_threads
    implicit none
    type(config_t) :: config
    type(mixture_t) :: mixture, permuted
    type(mixture_result_t) :: result_, serial_result, permuted_result, zero_result
    type(mixture_result_t) :: combined_result
    real(dp) :: summed(3), torque(3), reference(3), scale, work_scale
    complex(dp) :: work(3)
    character(len=1024) :: prefix
    integer :: species, driver, index
    integer, parameter :: permutation(3) = [3, 1, 2]

    config%s = 0.5_dp
    config%epsmn = 1.0e-3_dp
    config%m0 = -6
    config%mph = 3
    config%comptorque = .true.
    config%magdrift = .false.
    config%noshear = .true.
    config%inp_swi = 8
    config%log_level = -1
    config%conserving_fp_collisions = .true.
    config%conserving_speed_steps = 2
    config%conserving_pitch_steps = 2
    config%conserving_theta_order = 8
    config%conserving_velocity_order = 4
    config%conserving_cm_order = 2
    config%conserving_rel_order = 2
    config%conserving_polar_order = 3
    config%conserving_azimuth_order = 6
    mixture%reference_s = config%s
    mixture%omega_e = -700.0_dp
    mixture%coulomb_log = 20.0_dp
    allocate(mixture%species(3))
    mixture%species(1) = species_state_t('D', 3.3435837724e-27_dp, 1.0_dp, &
        5.0e18_dp, 2000.0_dp, -2.0e18_dp, 400.0_dp)
    mixture%species(2) = species_state_t('T', 5.0073567446e-27_dp, 1.0_dp, &
        5.0e18_dp, 2000.0_dp, -2.0e18_dp, 400.0_dp)
    mixture%species(3) = species_state_t('e', 9.1093837139e-31_dp, -1.0_dp, &
        1.0e19_dp, 2000.0_dp, -4.0e18_dp, 400.0_dp)
    call neort_init(config, 'lorentz_circle.bc')
    call neort_prepare_mixture(mixture)
    call omp_set_dynamic(.false.)
    call omp_set_num_threads(4)
    call neort_compute_mixture_at_s(config%s, result_)
    summed = 0.0_dp
    do species = 1, 3
        associate(data => result_%species(species)%transport)
            call require(data%conserving_active, 'native conserving feature')
            call require(.not. data%susceptibility_valid, 'no auxiliary D claim')
            call require(all(data%conserving_cross_work_qualified), &
                'nonzero-detuning driver qualification')
            call require(data%conserving_rank == 51, 'joint 51-DOF operator')
            call require(data%conserving_residual < 1.0e-10_dp, &
                'physical coupled solve residual')
            work = sum(data%conserving_cross_work, dim=2)
            work_scale = maxval(abs(data%conserving_work))
            call require(work_scale > 0.0_dp, 'native EM work is nonzero')
            call require(maxval(abs(work - data%conserving_work)) < &
                2.0e-10_dp*work_scale, 'all driver work sums to physical work')
            torque = [data%torque%Tt, data%torque%Tco, data%torque%Tctr]
            reference = -real(config%mph, dp)*acos(-1.0_dp)* &
                aimag(data%conserving_work)
            call close(torque, reference, 2.0e-13_dp, 'direct physical work torque')
            call require(data%conserving_entropy >= 0.0_dp, &
                'physical collision entropy')
            ! A diagonal background-rate solve could not create a receiving
            ! response to another species' source. Observe the actual output.
            driver = 1 + mod(species, 3)
            call require(maxval(abs(data%conserving_cross_work(:, driver))) > &
                1.0e-15_dp*work_scale, 'cross-species field-particle response')
            summed = summed + [torque(2), torque(3), torque(1)]
        end associate
    end do
    call close(result_%torque_total, summed, 1.0e-13_dp, 'physical branch sum')
    call omp_set_num_threads(1)
    call neort_compute_mixture_at_s(config%s, serial_result)
    call close(serial_result%torque_total, result_%torque_total, 2.0e-10_dp, &
        'one versus four threads')
    do species = 1, 3
        scale = maxval(abs(result_%species(species)%transport%conserving_work))
        call require(maxval(abs(serial_result%species(species)%transport% &
            conserving_work - result_%species(species)%transport% &
            conserving_work)) < 2.0e-10_dp*scale, 'per-species parallel setup')
    end do
    permuted = mixture
    do species = 1, 3
        permuted%species(species) = mixture%species(permutation(species))
    end do
    call neort_prepare_mixture(permuted)
    call neort_compute_mixture_at_s(config%s, permuted_result)
    do species = 1, 3
        index = permutation(species)
        work_scale = maxval(abs(result_%species(index)%transport%conserving_work))
        call require(maxval(abs(permuted_result%species(species)%transport% &
            conserving_work - result_%species(index)%transport% &
            conserving_work)) < 2.0e-9_dp*work_scale, 'D/T/e permutation covariance')
    end do
    call get_command_argument(1, prefix)
    if (len_trim(prefix) > 0) &
        call write_mixture_output(result_, trim(prefix), 'both')
    mixture%omega_e = 0.0_dp
    mixture%species%density_gradient_m3 = 0.0_dp
    mixture%species%temperature_gradient_ev = 0.0_dp
    call neort_prepare_mixture(mixture)
    call neort_compute_mixture_at_s(config%s, zero_result)
    call require(all(zero_result%torque_total == 0.0_dp), &
        'zero physical force at undamped conserved nulls')
    call require(zero_result%species(1)%transport%conserving_nullity >= 5, &
        'number/common-flow/common-temperature nulls survive')
    ! At zero detuning choose dn/n=-dT/T. Independent Maxwellian moments
    ! <u^2>=3/2 and <u^4>=15/4 then cancel each species density source:
    ! A1=-5*A2/2. Quasineutrality cancels the SUMMED energy source, although
    ! each driver individually forces that conserved common-temperature mode.
    ! Resolve source moments separately from the collision/basis quadrature.
    config%conserving_velocity_order = 12
    call neort_init(config, 'lorentz_circle.bc')
    do species = 1, 3
        mixture%species(species)%density_gradient_m3 = &
            -0.2_dp*mixture%species(species)%density_m3
        mixture%species(species)%temperature_gradient_ev = 400.0_dp
    end do
    call neort_prepare_mixture(mixture)
    call neort_compute_mixture_at_s(config%s, combined_result)
    do species = 1, 3
        associate(data => combined_result%species(species)%transport)
            call require(data%conserving_compatibility < 1.0e-9_dp, &
                'nonzero combined physical source compatibility')
            call require(.not. any(data%conserving_cross_work_qualified), &
                'individual energy-forcing driver responses unqualified')
            call require(maxval(abs(data%conserving_work)) > 0.0_dp, &
                'compatible combined physical response is nonzero')
        end associate
    end do
    print '(a)', 'Native coupled D/T/e Landau response/work/permutation checks passed'

contains

    subroutine close(actual, expected, tolerance, message)
        real(dp), intent(in) :: actual(:), expected(:), tolerance
        character(len=*), intent(in) :: message
        call require(maxval(abs(actual - expected)) <= &
            tolerance*max(maxval(abs(expected)), tiny(1.0_dp)), message)
    end subroutine close

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message
        if (.not. condition) then
            print '(a)', message
            error stop 'conserving mixture oracle failed'
        end if
    end subroutine require
end program test_conserving_mixture
