program test_krook_config
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_value, ieee_quiet_nan, ieee_positive_inf
    use neort_config, only: config_t, set_config, read_and_set_config
    use neort_krook, only: krook_collisions, krook_nu, krook_model, &
                       krook_harmonic, krook_geometry_trapped, krook_geometry_passing, &
                           krook_eta_steps, krook_rtol, krook_atol, krook_max_depth
    use neort_profiles, only: profile_species
    implicit none
    type(config_t) :: config
    character(len=64) :: mode
    integer :: unit

    call get_command_argument(1, mode)
    select case (trim(mode))
    case ('negative-rate')
        config%krook_nu = -1.0_dp
    case ('invalid-model')
        config%krook_model = 3
    case ('invalid-harmonic')
        config%krook_harmonic = 2
    case ('invalid-geometry')
        config%krook_geometry_passing = 0.0_dp
    case ('invalid-steps')
        config%krook_eta_steps = 1
    case ('invalid-tolerance')
        config%krook_rtol = 0.0_dp
    case ('nan-rate')
        config%krook_nu = ieee_value(0.0_dp, ieee_quiet_nan)
    case ('infinite-rate')
        config%krook_nu = ieee_value(0.0_dp, ieee_positive_inf)
    case ('nan-geometry')
        config%krook_geometry_trapped = ieee_value(0.0_dp, ieee_quiet_nan)
    case ('infinite-tolerance')
        config%krook_rtol = ieee_value(0.0_dp, ieee_positive_inf)
    case ('electron-nonlinear')
        config%profile_species = 0
        config%nonlin = .true.
    case ('krook-nonlinear')
        config%nonlin = .true.
        config%krook_collisions = .true.
        config%krook_nu = 1.0_dp
    end select
    call set_config(config)
    if (len_trim(mode) > 0) stop
    if (krook_collisions .or. krook_nu /= 0.0_dp) error stop 'default not collisionless'
    config%krook_collisions = .true.
    config%krook_nu = 7.0_dp
    config%krook_model = 1
    config%krook_harmonic = 1
    config%krook_geometry_trapped = 2.0_dp
    config%krook_geometry_passing = 3.0_dp
    config%krook_eta_steps = 64
    config%krook_rtol = 2.0e-6_dp
    config%krook_atol = 1.0e-12_dp
    config%krook_max_depth = 12
    config%profile_species = 2
    call set_config(config)
    call assert_settings()
    open (newunit=unit, file='krook-config.in', status='replace')
    write (unit, '(a)') '&params qs=1,ms=3,krook_collisions=T,krook_nu=7,'
    write (unit, '(a)') 'krook_model=1,krook_harmonic=1,krook_geometry_trapped=2,'
    write (unit, '(a)') 'krook_geometry_passing=3,krook_eta_steps=64,krook_rtol=2e-6,'
    write (unit, '(a)') 'krook_atol=1e-12,krook_max_depth=12,profile_species=2 /'
    close (unit)
    call read_and_set_config('krook-config.in')
    call assert_settings()
    open (newunit=unit, file='krook-config.in', status='replace')
    write (unit, '(a)') '&params qs=1,ms=2 /'
    close (unit)
    call read_and_set_config('krook-config.in')
    if (krook_collisions .or. krook_nu /= 0.0_dp) error stop 'namelist defaults leaked'
    if (profile_species /= 1) error stop 'species default leaked'
    open (newunit=unit, file='krook-config.in', status='old')
    close (unit, status='delete')
    print *, 'test_krook_config PASSED'

contains

    subroutine assert_settings()
      if (.not. krook_collisions .or. krook_nu /= 7.0_dp) error stop 'Krook enable/rate'
        if (krook_model /= 1 .or. krook_harmonic /= 1) error stop 'Krook model'
        if (krook_geometry_trapped /= 2.0_dp) error stop 'Krook trapped geometry'
        if (krook_geometry_passing /= 3.0_dp) error stop 'Krook passing geometry'
        if (krook_eta_steps /= 64 .or. krook_max_depth /= 12) error stop 'Krook grid'
 if (krook_rtol /= 2.0e-6_dp .or. krook_atol /= 1.0e-12_dp) error stop 'Krook tolerance'
        if (profile_species /= 2) error stop 'species selector'
    end subroutine assert_settings
end program test_krook_config
