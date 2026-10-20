program test_species_mixture
    use iso_fortran_env, only: dp => real64
    use neort_species, only: species_state_t, mixture_t, configure_mixture, &
        evaluate_mixture_at_s, mixture_at_s, mixture_enabled, reset_mixture, &
        selected_species, common_omega_e, common_domega_e_ds, read_species_input
    use neort_profiles, only: init_plasma_at_s, init_profile_at_s, &
        init_thermodynamic_forces, ns, ts, dnsds, dtsds, vth, dvthds, &
        M_t, dM_tds, Om_tE, dOm_tEds, A1, A2
    use do_magfie_mod, only: s, sign_theta
    use util, only: mi, qi, qe, ev, c
    implicit none
    type(mixture_t) :: mixture
    character(len=32) :: mode
    real(dp) :: expected(3), thermal_speed, expected_a1, expected_a2
    real(dp), parameter :: radius = 620.0_dp, psi_pr = 7.0e7_dp, safety = 1.7_dp
    integer :: k

    call get_command_argument(1, mode)
    call make_mixture(mixture)
    if (mode == 'charge-imbalance') mixture%species(3)%density_m3 = 3.0e19_dp
    if (mode == 'gradient-imbalance') mixture%species(3)%density_gradient_m3 = 0.0_dp
    if (mode == 'negative-mass') mixture%species(1)%mass_kg = -1.0_dp
    call configure_mixture(mixture)
    if (len_trim(mode) > 0) stop
    s = 0.43_dp
    call evaluate_mixture_at_s(s)
    call assert_close('common rotation', [common_omega_e, common_domega_e_ds], &
        [-12000.0_dp + 800.0_dp*(s - mixture%reference_s), 800.0_dp])
    do k = 1, size(mixture%species)
        selected_species = k
        call init_plasma_at_s()
        expected = [mixture%species(k)%mass_kg*1.0e3_dp, &
            mixture%species(k)%charge_e*qe, mixture_at_s(k)%density_m3*1.0e-6_dp]
        call assert_close('explicit species', [mi, qi, ns], expected)
        call assert_close('explicit derivatives', [dnsds, dtsds], &
            [mixture%species(k)%density_gradient_m3*1.0e-6_dp, &
             mixture%species(k)%temperature_gradient_ev])
        thermal_speed = sqrt(2.0_dp*ev*mixture_at_s(k)%temperature_ev/expected(1))
        call assert_close('species thermal speed', [vth], [thermal_speed])
        call init_profile_at_s(radius, 1.0_dp, 1.0_dp)
        call assert_close('physical rotation held across masses', [Om_tE, dOm_tEds], &
            [common_omega_e, common_domega_e_ds])
        call assert_close('species Mach', [M_t, dM_tds], &
            [radius*common_omega_e/thermal_speed, &
             radius*(common_domega_e_ds/thermal_speed - &
                     common_omega_e*dvthds/thermal_speed**2)])
        call init_thermodynamic_forces(psi_pr, safety)
        expected_a2 = mixture%species(k)%temperature_gradient_ev/ts
        expected_a1 = dnsds/ns - expected(2)*sign_theta*psi_pr*common_omega_e/ &
            (safety*c*ts*ev) - 1.5_dp*expected_a2
        call assert_close('signed affine Maxwellian force', [A1, A2], &
            [expected_a1, expected_a2])
    end do
    call check_file_input()
    call reset_mixture()
    if (mixture_enabled) error stop 'reset must restore isolated-species mode'
    print *, 'General D/T/e mixture, physical rotation and forces pass'
contains
    subroutine make_mixture(value)
        type(mixture_t), intent(out) :: value
        value%reference_s = 0.2_dp
        value%omega_e = -12000.0_dp
        value%domega_e_ds = 800.0_dp
        value%coulomb_log = 20.0_dp
        allocate(value%species(3))
        value%species(1) = species_state_t('D', 3.3435837724e-27_dp, 1.0_dp, &
            8.0e19_dp, 8000.0_dp, -1.6e19_dp, -2400.0_dp)
        value%species(2) = species_state_t('T', 5.0073567446e-27_dp, 1.0_dp, &
            2.0e19_dp, 7000.0_dp, -2.0e18_dp, -1400.0_dp)
        value%species(3) = species_state_t('e', 9.1093837139e-31_dp, -1.0_dp, &
            1.0e20_dp, 9000.0_dp, -1.8e19_dp, -2250.0_dp)
    end subroutine make_mixture

    subroutine check_file_input()
        integer :: unit
        open(newunit=unit, file='species-test.in', status='replace')
        write(unit, '(a)') '&species_count n_species=3 /'
        write(unit, '(a)') '&species reference_s=.2, omega_e=-12000, domega_e_ds=800, coulomb_log=20,'
        write(unit, '(a)') "species_name='D','T','e', charge_e=1,1,-1,"
        write(unit, '(a)') 'mass_kg=3.3435837724e-27,5.0073567446e-27,9.1093837139e-31,'
        write(unit, '(a)') 'density_m3=8e19,2e19,1e20, temperature_ev=8000,7000,9000,'
        write(unit, '(a)') 'density_gradient_m3=-1.6e19,-2e18,-1.8e19,'
        write(unit, '(a)') 'temperature_gradient_ev=-2400,-1400,-2250 /'
        close(unit)
        call read_species_input('species-test.in')
        call evaluate_mixture_at_s(s)
        call assert_close('explicit file profiles', &
            [mixture_at_s(1)%density_m3, mixture_at_s(2)%temperature_ev], &
            [8.0e19_dp - 1.6e19_dp*(s - 0.2_dp), &
             7000.0_dp - 1400.0_dp*(s - 0.2_dp)])
        open(newunit=unit, file='species-test.in', status='old')
        close(unit, status='delete')
    end subroutine check_file_input

    subroutine assert_close(label, actual, expected)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: actual(:), expected(:)
        if (any(abs(actual - expected) > 1.0e-12_dp*max(abs(expected), tiny(1.0_dp)))) then
            print *, trim(label), actual, expected
            error stop 'mixture oracle mismatch'
        end if
    end subroutine assert_close
end program test_species_mixture
