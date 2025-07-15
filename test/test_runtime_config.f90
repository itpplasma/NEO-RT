program test_runtime_config
    ! Test program for runtime configuration switching between thin/thick orbits
    use iso_fortran_env, only: real64
    use runtime_config, only: get_use_thick_orbits, set_use_thick_orbits
    use orbit_interface, only: orbit_calculator_t, orbit_calculator_factory
    use field_interface, only: field_evaluator_t, field_evaluator_factory
    implicit none
    
    integer, parameter :: dp = real64
    
    ! Test variables
    logical :: use_thick_orbits
    real(dp) :: v, eta, taub, delphi, extraset(7)
    real(dp) :: s_flux, theta_boozer, phi_boozer
    logical :: success
    class(orbit_calculator_t), allocatable :: calculator
    class(field_evaluator_t), allocatable :: evaluator
    real(dp) :: psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2
    real(dp) :: R, Z
    
    print *, '=== Runtime Configuration Test ==='
    
    ! Test 1: Default configuration (should be thin orbits)
    print *, 'Test 1: Default configuration'
    use_thick_orbits = get_use_thick_orbits()
    print *, 'Default use_thick_orbits:', use_thick_orbits
    
    ! Test 2: Set to thin orbits explicitly
    print *, 'Test 2: Set to thin orbits'
    call set_use_thick_orbits(.false.)
    use_thick_orbits = get_use_thick_orbits()
    print *, 'After setting to thin orbits:', use_thick_orbits
    
    ! Test thin orbit calculation
    calculator = orbit_calculator_factory(use_thick_orbits)
    v = 1.0d6
    eta = 0.5d0
    s_flux = 0.5d0
    theta_boozer = 0.0d0
    phi_boozer = 0.0d0
    
    call calculator%find_bounce(v, eta, s_flux, theta_boozer, phi_boozer, &
                               taub, delphi, extraset, success)
    print *, 'Thin orbit bounce time:', taub
    call calculator%cleanup()
    
    ! Test thin orbit field evaluation
    R = 1.5d0
    Z = 0.0d0
    evaluator = field_evaluator_factory(use_thick_orbits)
    call evaluator%evaluate_field(R, Z, psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2)
    print *, 'Thin orbit field psi:', psif
    call evaluator%cleanup()
    
    ! Test 3: Set to thick orbits
    print *, 'Test 3: Set to thick orbits'
    call set_use_thick_orbits(.true.)
    use_thick_orbits = get_use_thick_orbits()
    print *, 'After setting to thick orbits:', use_thick_orbits
    
    ! Test thick orbit calculation
    calculator = orbit_calculator_factory(use_thick_orbits)
    call calculator%find_bounce(v, eta, s_flux, theta_boozer, phi_boozer, &
                               taub, delphi, extraset, success)
    print *, 'Thick orbit bounce time:', taub
    call calculator%cleanup()
    
    ! Test thick orbit field evaluation
    evaluator = field_evaluator_factory(use_thick_orbits)
    call evaluator%evaluate_field(R, Z, psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2)
    print *, 'Thick orbit field psi:', psif
    call evaluator%cleanup()
    
    print *, '=== Runtime Configuration Test Complete ==='
    
end program test_runtime_config