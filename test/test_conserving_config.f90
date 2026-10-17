program test_conserving_config
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use neort_config, only: config_t, set_config, read_and_set_config
    use neort_fp_mode, only: fp_physical_mode
    use neort_fp_controls, only: validate_fp_background
    use neort_species, only: mixture_t, configure_mixture, selected_species
    use neort_profiles, only: init_plasma_at_s, init_thermodynamic_forces, vth, &
        A1, A2, Om_tE
    use do_magfie_mod, only: sign_theta
    use util, only: qi, c
    use neort_fp_controls, only: conserving_fp_collisions, conserving_speed_steps, &
        conserving_pitch_steps, conserving_umax, conserving_theta_order, &
        conserving_rate_scale, &
        conserving_velocity_order, conserving_cm_order, conserving_rel_order, &
        conserving_polar_order, conserving_azimuth_order, conserving_rank_rtol, &
        conserving_compatibility_rtol
    implicit none
    type(config_t) :: config
    type(mixture_t) :: mixture
    real(dp), parameter :: elementary_charge = 1.602176634e-19_dp
    real(dp), parameter :: mass(3) = &
        [3.345228e-27_dp, 5.017842e-27_dp, 9.1091e-31_dp]
    real(dp), parameter :: charge(3) = [1.0_dp, 1.0_dp, -1.0_dp]
    real(dp), parameter :: density(3) = [5.0e18_dp, 5.0e18_dp, 1.0e19_dp]
    real(dp), parameter :: temperature = 2000.0_dp, psi_cgs = 1.5e6_dp
    real(dp), parameter :: omega = -700.0_dp, safety = 2.0_dp
    character(1), parameter :: species_name(3) = ['D', 'T', 'e']
    real(dp) :: expected_speed, force_scale
    character(64) :: mode
    integer :: unit, species

    call get_command_argument(1, mode)
    config%conserving_fp_collisions = .true.
    select case (trim(mode))
    case ('rate')
        config%conserving_rate_scale = 0.0_dp
    case ('speed')
        config%conserving_speed_steps = 1
    case ('pitch')
        config%conserving_pitch_steps = 1
    case ('umax')
        config%conserving_umax = 0.0_dp
    case ('quadrature')
        config%conserving_rel_order = 1
    case ('azimuth')
        config%conserving_azimuth_order = 3
    case ('rank')
        config%conserving_rank_rtol = 0.0_dp
    case ('nan-rank')
        config%conserving_rank_rtol = ieee_value(0.0_dp, ieee_quiet_nan)
    case ('compatibility')
        config%conserving_compatibility_rtol = 1.0_dp
    case ('tp')
        config%tp_fp_collisions = .true.
    case ('lorentz')
        config%lorentz_pitch_collisions = .true.
    case ('krook')
        config%krook_collisions = .true.
    case ('drift-phase')
        config%drift_phase = .true.
    case ('fow')
        config%fow_radial = .true.
    case ('nopassing')
        config%nopassing = .true.
    end select
    call set_config(config)
    if (len_trim(mode) > 0) error stop 'invalid conserving controls accepted'
    if (.not. conserving_fp_collisions .or. .not. fp_physical_mode()) &
        error stop 'active conserving mode is not physical FP'
    ! Independent physical units and signed electrostatic-equilibrium oracle.
    ! Wb = 1e-8 G cm^2; SI dn/n = (q/T_J)*(Psi_Wb/q_safety)*Omega.
    allocate(mixture%species(3))
    mixture%coulomb_log = 20.0_dp
    do species = 1, 3
        mixture%species(species)%name = species_name(species)
        mixture%species(species)%mass_kg = mass(species)
        mixture%species(species)%charge_e = charge(species)
        mixture%species(species)%density_m3 = density(species)
        mixture%species(species)%temperature_ev = temperature
    end do
    do species = 1, 3
        ! Balance the chosen species independently; a companion carries the
        ! neutralizing gradient. A common equal-T mixture cannot balance all
        ! signed species simultaneously at finite electric field.
        mixture%species%density_gradient_m3 = 0.0_dp
        mixture%species(species)%density_gradient_m3 = density(species)* &
            charge(species)/temperature*(sign_theta*psi_cgs*1.0e-8_dp/safety)*omega
        if (species == 3) then
            mixture%species(1)%density_gradient_m3 = &
                mixture%species(species)%density_gradient_m3
        else
            mixture%species(3)%density_gradient_m3 = &
                mixture%species(species)%density_gradient_m3
        end if
        call configure_mixture(mixture)
        selected_species = species
        call init_plasma_at_s()
        call validate_fp_background()
        expected_speed = 100.0_dp* &
            sqrt(2.0_dp*temperature*elementary_charge/mass(species))
        if (abs(vth/expected_speed - 1.0_dp) > 1.0e-14_dp) &
            error stop 'conserving thermal speed does not represent the SI Maxwellian'
        if (abs(10.0_dp*qi/c/(charge(species)*elementary_charge) - &
                1.0_dp) > 1.0e-14_dp) &
            error stop 'conserving electromagnetic charge is not physical SI'
        Om_tE = omega
        call init_thermodynamic_forces(psi_cgs, safety)
        force_scale = abs(charge(species)/temperature* &
            (psi_cgs*1.0e-8_dp/safety)*omega)
        if (abs(A1) > 1.0e-13_dp*force_scale .or. A2 /= 0.0_dp) &
            error stop 'conserving signed electrostatic equilibrium has nonzero drive'
    end do
    config%conserving_fp_collisions = .false.
    call set_config(config)
    if (fp_physical_mode()) error stop 'conserving mode leaked into disabled path'

    ! A disabled namelist still reads and retains every distinct control.
    open(newunit=unit, file='conserving-config.in', status='replace')
    write(unit, '(a)') '&params qs=1,ms=2,conserving_fp_collisions=F,'
    write(unit, '(a)') 'conserving_speed_steps=3,conserving_pitch_steps=4,'
    write(unit, '(a)') 'conserving_umax=5,conserving_rate_scale=.7,'
    write(unit, '(a)') 'conserving_theta_order=6,'
    write(unit, '(a)') 'conserving_velocity_order=7,conserving_cm_order=8,'
    write(unit, '(a)') 'conserving_rel_order=9,conserving_polar_order=10,'
    write(unit, '(a)') 'conserving_azimuth_order=12,conserving_rank_rtol=1e-10,'
    write(unit, '(a)') 'conserving_compatibility_rtol=1e-8 /'
    close(unit)
    call read_and_set_config('conserving-config.in')
    if (conserving_speed_steps /= 3 .or. conserving_pitch_steps /= 4 .or. &
            conserving_umax /= 5.0_dp .or. conserving_rate_scale /= 0.7_dp .or. &
            conserving_theta_order /= 6 .or. &
            conserving_velocity_order /= 7 .or. conserving_cm_order /= 8 .or. &
            conserving_rel_order /= 9 .or. conserving_polar_order /= 10 .or. &
            conserving_azimuth_order /= 12 .or. &
            conserving_rank_rtol /= 1.0e-10_dp .or. &
            conserving_compatibility_rtol /= 1.0e-8_dp) &
        error stop 'distinct conserving namelist controls were not retained'
    open(newunit=unit, file='conserving-config.in', status='replace')
    write(unit, '(a)') '&params qs=1,ms=2 /'
    close(unit)
    call read_and_set_config('conserving-config.in')
    if (fp_physical_mode() .or. conserving_speed_steps /= 2 .or. &
            conserving_theta_order /= 8) error stop 'conserving defaults leaked'
    open(newunit=unit, file='conserving-config.in', status='old')
    close(unit, status='delete')
    print *, 'conserving FP config and disabled reset pass'
end program test_conserving_config
