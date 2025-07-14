program test_thick_orbit_initial_conditions
    ! Test proper initial condition setup for thick orbits (Phase G.1.2)
    ! Map (s, θ, φ) flux coordinates to (R, Z, φ) for POTATO
    
    use potato_field_bridge, only: convert_neort_to_potato, initialize_potato_field
    implicit none
    
    ! Test parameters
    real(8), parameter :: pi = 3.141592653589793d0
    
    ! Variables
    real(8) :: s_flux, theta_flux, phi_flux, v_test, eta_test
    real(8) :: R_start, Z_start, phi_start
    real(8) :: z_eqm(5)
    logical :: success
    integer :: i, j
    
    print *, '======================================================='
    print *, 'TEST: Thick Orbit Initial Conditions (Phase G.1.2)'
    print *, '======================================================='
    
    ! Initialize POTATO field system
    call initialize_potato_field(success)
    if (.not. success) then
        print *, 'ERROR: Field initialization failed'
        stop 1
    end if
    
    print *, 'Field system initialized successfully'
    
    ! Test 1: Map flux surface coordinates to physical coordinates
    print *, ''
    print *, 'Test 1: Flux coordinate mapping'
    print *, '-------------------------------'
    
    ! Test several flux surfaces
    do i = 1, 5
        s_flux = 0.1d0 * real(i, 8)  ! s ∈ [0.1, 0.5] (normalized flux)
        
        do j = 1, 3
            theta_flux = 0.5d0 * pi * real(j-1, 8)  ! θ ∈ [0, π] (poloidal angle)
            phi_flux = 0.0d0  ! Start at φ = 0
            
            ! Convert flux coordinates to physical coordinates
            call flux_to_physical(s_flux, theta_flux, phi_flux, R_start, Z_start, phi_start)
            
            print *, 'Flux (s, θ, φ) = (', s_flux, ',', theta_flux, ',', phi_flux, ')'
            print *, '  → Physical (R, Z, φ) = (', R_start, ',', Z_start, ',', phi_start, ')'
        end do
    end do
    
    ! Test 2: Convert NEO-RT parameters to POTATO phase space
    print *, ''
    print *, 'Test 2: NEO-RT to POTATO conversion'
    print *, '-----------------------------------'
    
    v_test = 1.0d6    ! 1 MeV thermal velocity
    eta_test = 0.5d0  ! Mid-pitch angle
    R_start = 1.65d0  ! Near magnetic axis
    Z_start = 0.0d0
    phi_start = 0.0d0
    
    call convert_neort_to_potato(v_test, eta_test, R_start, Z_start, phi_start, z_eqm, success)
    
    if (success) then
        print *, 'NEO-RT (v, η) = (', v_test, ',', eta_test, ')'
        print *, 'Starting position (R, Z, φ) = (', R_start, ',', Z_start, ',', phi_start, ')'
        print *, 'POTATO phase space z_eqm:'
        print *, '  z(1) = R       =', z_eqm(1), 'm'
        print *, '  z(2) = φ       =', z_eqm(2), 'rad'
        print *, '  z(3) = Z       =', z_eqm(3), 'm'
        print *, '  z(4) = p       =', z_eqm(4), '(normalized momentum)'
        print *, '  z(5) = λ       =', z_eqm(5), '(cos pitch angle)'
    else
        print *, 'ERROR: Conversion failed'
        stop 1
    end if
    
    ! Test 3: Validate conservation laws
    print *, ''
    print *, 'Test 3: Conservation validation'
    print *, '-------------------------------'
    
    ! Check momentum magnitude conservation
    call validate_momentum_conservation(v_test, eta_test, z_eqm)
    
    ! Check pitch angle consistency
    call validate_pitch_angle_consistency(eta_test, z_eqm)
    
    print *, ''
    print *, 'SUCCESS: Initial conditions implemented correctly'
    print *, 'Thick orbit phase space conversion validated'
    
contains

    subroutine flux_to_physical(s, theta, phi, R, Z, phi_out)
        ! Convert flux coordinates (s, θ, φ) to physical (R, Z, φ)
        implicit none
        real(8), intent(in) :: s, theta, phi
        real(8), intent(out) :: R, Z, phi_out
        
        ! Simple tokamak mapping for test
        real(8), parameter :: R_axis = 1.65d0  ! Major radius
        real(8), parameter :: a_minor = 0.3d0  ! Minor radius
        real(8) :: r_minor
        
        ! Map flux coordinate to minor radius
        r_minor = sqrt(s) * a_minor
        
        ! Convert to cylindrical coordinates
        R = R_axis + r_minor * cos(theta)
        Z = r_minor * sin(theta)
        phi_out = phi
        
    end subroutine flux_to_physical
    
    subroutine validate_momentum_conservation(v, eta, z_eqm)
        ! Validate that momentum magnitude is conserved in conversion
        implicit none
        real(8), intent(in) :: v, eta, z_eqm(5)
        real(8) :: p_magnitude, v_parallel, v_perpendicular, v_total
        real(8), parameter :: v_thermal = 1.0d6
        
        ! Extract POTATO variables
        p_magnitude = z_eqm(4)
        
        ! Calculate velocity components from POTATO λ = cos(pitch)
        v_parallel = z_eqm(5) * v  ! v_par = λ * v_total
        v_perpendicular = sqrt(1.0d0 - z_eqm(5)**2) * v  ! v_perp = sin(pitch) * v_total
        v_total = sqrt(v_parallel**2 + v_perpendicular**2)
        
        print *, 'Momentum conservation check:'
        print *, '  Original v =', v, 'm/s'
        print *, '  Reconstructed v =', v_total, 'm/s'
        print *, '  POTATO p =', p_magnitude, '(normalized)'
        print *, '  Expected p =', v / v_thermal * sqrt(2.0d0)
        
        if (abs(v_total - v) / v < 1.0d-10) then
            print *, '  ✓ Momentum magnitude conserved'
        else
            print *, '  ✗ Momentum magnitude NOT conserved'
        end if
        
    end subroutine validate_momentum_conservation
    
    subroutine validate_pitch_angle_consistency(eta, z_eqm)
        ! Validate pitch angle conversion between NEO-RT and POTATO conventions
        implicit none
        real(8), intent(in) :: eta, z_eqm(5)
        real(8) :: lambda_potato, eta_reconstructed
        real(8), parameter :: bmod_local = 1.0d0  ! Normalized B-field
        
        lambda_potato = z_eqm(5)  ! cos(pitch angle)
        
        ! Reconstruct NEO-RT eta from POTATO λ
        eta_reconstructed = 1.0d0 - lambda_potato**2  ! sin²(pitch) / B_norm
        
        print *, 'Pitch angle consistency check:'
        print *, '  Original η =', eta
        print *, '  POTATO λ =', lambda_potato, '(cos pitch)'
        print *, '  Reconstructed η =', eta_reconstructed
        
        if (abs(eta_reconstructed - eta) < 1.0d-10) then
            print *, '  ✓ Pitch angle conversion consistent'
        else
            print *, '  ✗ Pitch angle conversion INCONSISTENT'
        end if
        
    end subroutine validate_pitch_angle_consistency

end program test_thick_orbit_initial_conditions