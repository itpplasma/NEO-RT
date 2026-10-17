program test_lorentz_dispatch
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_lib, only: config_t, mixture_t, mixture_result_t, neort_init, &
        neort_prepare_mixture, neort_compute_mixture_at_s, neort_setup_at_s
    use neort_lorentz_native, only: compute_lorentz_transport
    use neort_species, only: selected_species, species_mixture, mixture_enabled
    use neort_config, only: set_config
    use neort_lorentz_controls, only: lorentz_pitch_collisions, lorentz_nu, &
        lorentz_rate_source, lorentz_rate_scale, lorentz_pitch_steps
    use neort_profiles, only: vth
    use neort, only: mth_max_abs
    use omp_lib, only: omp_set_num_threads
    implicit none
    real(dp), parameter :: mass(3) = [3.345228e-27_dp, 5.017842e-27_dp, 9.1091e-31_dp]
    real(dp), parameter :: temperature(3) = [1800.0_dp, 2600.0_dp, 2100.0_dp]
    real(dp), parameter :: density(3) = [5.0e15_dp, 5.0e15_dp, 1.0e16_dp]
    real(dp), parameter :: charge(3) = [1.0_dp, 1.0_dp, -1.0_dp]
    type(config_t) :: config
    type(mixture_t) :: in_memory_mixture
    type(mixture_result_t) :: result_, parallel_result
    real(dp) :: expected_d(2, 3), expected_t(3), actual_t(3), scale
    integer :: species, unit, info, harmonic

    open(newunit=unit, file='lorentz-species.in', status='replace')
    write(unit, '(a)') '&species_count n_species=3 /'
    write(unit, '(a)') &
        '&species reference_s=.5,omega_e=-700,domega_e_ds=0,coulomb_log=20,'
    write(unit, '(a)') "species_name='D','T','e',charge_e=1,1,-1,"
    write(unit, '(a,3(es24.16e3,","))') 'mass_kg=', mass
    write(unit, '(a,3(es24.16e3,","))') 'density_m3=', density
    write(unit, '(a,3(es24.16e3,","))') 'temperature_ev=', temperature
    write(unit, '(a)') 'density_gradient_m3=0,0,0,temperature_gradient_ev=0,0,0 /'
    close(unit)
    config%s = 0.5_dp
    config%species_file = 'lorentz-species.in'
    config%epsmn = 1.0e-3_dp
    config%m0 = -6
    config%mph = 3
    config%comptorque = .true.
    config%magdrift = .false.
    config%noshear = .true.
    config%inp_swi = 8
    config%vsteps = 1
    config%vmax_over_vth = 2.0_dp
    config%lorentz_pitch_collisions = .true.
    config%lorentz_rate_source = 1
    config%lorentz_rate_scale = 0.7_dp
    config%lorentz_pitch_steps = 3
    config%mth_max_abs = 0
    call neort_init(config, 'lorentz_circle.bc')
    call omp_set_num_threads(1)
    call neort_compute_mixture_at_s(0.5_dp, result_)
    do species = 1, 3
        selected_species = species
        call neort_setup_at_s(0.5_dp)
        call compute_lorentz_transport(1.0e-6_dp*vth, 2.0_dp*vth, 1, 3, &
            independent_rate, expected_d, expected_t, info)
        if (info /= 0) error stop 'independent rate callback transport failed'
        actual_t = [result_%species(species)%transport%torque%Tt, &
            result_%species(species)%transport%torque%Tco, &
            result_%species(species)%transport%torque%Tctr]
        scale = maxval(abs(expected_t))
        if (scale <= 0.0_dp) error stop 'all-background Lorentz oracle vanished'
        if (maxval(abs(actual_t - expected_t)) > 5.0e-11_dp*scale) &
            error stop 'dispatch lost physical speed, '// &
                'all-background shape, or scale'
        if (size(result_%species(species)%transport%harmonics) /= 2) &
            error stop 'Lorentz solve repeated in independent harmonic loop'
        do harmonic = 1, 2
            associate(row => result_%species(species)%transport%harmonics(harmonic))
                if (harmonic == 1) then
                    if (row%mth /= 0 .or. row%Tresco /= 0.0_dp .or. &
                            row%Tresctr /= 0.0_dp) &
                        error stop 'trapped harmonic label/content'
                else
                    if (row%mth /= -6 .or. row%Trest /= 0.0_dp) &
                        error stop 'passing common-phase harmonic label/content'
                end if
            end associate
        end do
    end do
    ! A wider legacy harmonic bound cannot duplicate the graph solve. At the
    ! same time, changing the species schedule must preserve physical rates.
    mth_max_abs = 20
    call omp_set_num_threads(3)
    call neort_compute_mixture_at_s(0.5_dp, parallel_result)
    do species = 1, 3
        actual_t = [parallel_result%species(species)%transport%torque%Tco, &
            parallel_result%species(species)%transport%torque%Tctr, &
            parallel_result%species(species)%transport%torque%Tt]
        expected_t = [result_%species(species)%transport%torque%Tco, &
            result_%species(species)%transport%torque%Tctr, &
            result_%species(species)%transport%torque%Tt]
        if (maxval(abs(actual_t - expected_t)) > &
                5.0e-11_dp*maxval(abs(expected_t))) &
            error stop 'Lorentz graph depends on thread count or legacy harmonic bound'
    end do

    ! The documented library sequence prepares an in-memory mixture after init.
    in_memory_mixture = species_mixture
    config%species_file = ''
    call neort_init(config, 'lorentz_circle.bc')
    call neort_prepare_mixture(in_memory_mixture)
    call neort_compute_mixture_at_s(0.5_dp, parallel_result)
    if (maxval(abs(parallel_result%torque_total - result_%torque_total)) > &
            5.0e-11_dp*maxval(abs(result_%torque_total))) &
        error stop 'in-memory Lorentz mixture changed physical transport'
    call set_config(config_t())
    if (lorentz_pitch_collisions .or. lorentz_nu /= 0.0_dp .or. &
            lorentz_rate_source /= 0 .or. lorentz_rate_scale /= 1.0_dp .or. &
            lorentz_pitch_steps /= 32 .or. mixture_enabled) &
        error stop 'struct initialization retained Lorentz or mixture state'
    print *, 'Lorentz D/T/e dispatch matches independent SI all-background rate oracle'

contains

    function independent_rate(speed_cm_s) result(rate)
        real(dp), intent(in) :: speed_cm_s
        real(dp) :: rate, speed, background_speed, x, chandrasekhar, shape, gamma
        real(dp), parameter :: elementary_charge = 1.602176634e-19_dp
        real(dp), parameter :: epsilon0 = 8.8541878188e-12_dp
        real(dp), parameter :: pi = acos(-1.0_dp)
        integer :: background
        speed = speed_cm_s/100.0_dp
        rate = 0.0_dp
        do background = 1, 3
            background_speed = sqrt(2*elementary_charge*temperature(background)/ &
                mass(background))
            x = speed/background_speed
            chandrasekhar = (erf(x) - 2*x*exp(-x*x)/sqrt(pi))/(2*x*x)
            shape = erf(x) - chandrasekhar
            gamma = density(background)*(charge(selected_species)*charge(background)* &
                elementary_charge**2)**2*20.0_dp/ &
                (4*pi*epsilon0**2*mass(selected_species)**2)
            rate = rate + gamma*shape/speed**3
        end do
        rate = 0.7_dp*rate
    end function independent_rate

end program test_lorentz_dispatch
