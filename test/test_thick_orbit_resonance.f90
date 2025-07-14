program test_thick_orbit_resonance
    ! Test resonance calculation with thick orbit frequencies (Phase G.3.1)
    ! Test resonance condition: n·ω_φ - m·ω_θ = ω_mode
    
    implicit none
    
    ! Test parameters
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: v_thermal = 1.0d6  ! m/s
    
    ! Variables
    real(8) :: v_test, eta_test
    real(8) :: Om_th_thin, Om_ph_thin, Om_th_thick, Om_ph_thick
    real(8) :: resonance_thin, resonance_thick, resonance_shift
    real(8) :: orbit_width, resonance_width_thin, resonance_width_thick
    integer :: n_mode, m_mode
    real(8) :: omega_mode
    logical :: success
    integer :: i
    
    print *, '======================================================='
    print *, 'TEST: Thick Orbit Resonance Calculation (Phase G.3.1)'
    print *, '======================================================='
    
    ! Test 1: Basic resonance condition evaluation
    print *, ''
    print *, 'Test 1: Resonance condition with thick vs thin orbits'
    print *, '----------------------------------------------------'
    
    ! ASDEX Upgrade 2/1 mode parameters
    n_mode = 2       ! Toroidal mode number
    m_mode = 4       ! Poloidal mode number (roughly 2*q)
    omega_mode = 0.0d0  ! Static RMP (non-rotating)
    
    v_test = v_thermal
    eta_test = 0.5d0
    
    ! Calculate frequencies with both methods
    call calculate_reference_frequencies(v_test, eta_test, Om_th_thin, Om_ph_thin, &
                                       Om_th_thick, Om_ph_thick)
    
    ! Evaluate resonance condition: n·ω_φ - m·ω_θ = ω_mode  
    resonance_thin = real(n_mode, 8) * Om_ph_thin - real(m_mode, 8) * Om_th_thin - omega_mode
    resonance_thick = real(n_mode, 8) * Om_ph_thick - real(m_mode, 8) * Om_th_thick - omega_mode
    
    resonance_shift = resonance_thick - resonance_thin
    
    print *, 'Mode: n =', n_mode, ', m =', m_mode, ', ω_mode =', omega_mode
    print *, 'Particle: v =', v_test/v_thermal, 'v_thermal, η =', eta_test
    print *, ''
    print *, 'Thin orbit frequencies:'
    print *, '  ω_θ =', Om_th_thin, 'rad/s'
    print *, '  ω_φ =', Om_ph_thin, 'rad/s' 
    print *, '  Resonance = n·ω_φ - m·ω_θ =', resonance_thin, 'rad/s'
    print *, ''
    print *, 'Thick orbit frequencies:'
    print *, '  ω_θ =', Om_th_thick, 'rad/s'
    print *, '  ω_φ =', Om_ph_thick, 'rad/s'
    print *, '  Resonance = n·ω_φ - m·ω_θ =', resonance_thick, 'rad/s'
    print *, ''
    print *, 'Resonance shift due to thick orbits:'
    print *, '  Δ(resonance) =', resonance_shift, 'rad/s'
    print *, '  Relative shift =', abs(resonance_shift)/abs(resonance_thin) * 100.0d0, '%'
    
    ! Test 2: Resonance width effects
    print *, ''
    print *, 'Test 2: Finite orbit width effects on resonance'
    print *, '-----------------------------------------------'
    
    call calculate_orbit_width_criterion(v_test, eta_test, orbit_width)
    
    ! Calculate resonance widths (simplified model)
    resonance_width_thin = abs(Om_th_thin) * 0.01d0    ! ~1% width for thin orbits
    resonance_width_thick = abs(Om_th_thick) * orbit_width  ! Width scales with orbit width
    
    print *, 'Orbit width parameter δr/L_B =', orbit_width
    print *, 'Thin orbit resonance width =', resonance_width_thin, 'rad/s'
    print *, 'Thick orbit resonance width =', resonance_width_thick, 'rad/s'
    print *, 'Width ratio (thick/thin) =', resonance_width_thick / resonance_width_thin
    
    ! Test 3: Energy scan across resonance
    print *, ''
    print *, 'Test 3: Energy scan showing resonance shifts'
    print *, '-------------------------------------------'
    
    eta_test = 0.5d0
    print *, 'Scanning energy at fixed η =', eta_test
    print *, 'Energy (v_thermal)  Thin Res.(kHz)  Thick Res.(kHz)  Shift (%)'
    print *, '------------------------------------------------------------'
    
    do i = 1, 5
        v_test = v_thermal * (0.5d0 + 0.5d0 * real(i-1, 8))  ! 0.5 to 2.5 v_thermal
        
        call calculate_reference_frequencies(v_test, eta_test, Om_th_thin, Om_ph_thin, &
                                           Om_th_thick, Om_ph_thick)
        
        resonance_thin = real(n_mode, 8) * Om_ph_thin - real(m_mode, 8) * Om_th_thin
        resonance_thick = real(n_mode, 8) * Om_ph_thick - real(m_mode, 8) * Om_th_thick
        resonance_shift = (resonance_thick - resonance_thin) / resonance_thin * 100.0d0
        
        print '(F15.1, F15.1, F16.1, F10.2)', v_test/v_thermal, &
              resonance_thin/1000.0d0, resonance_thick/1000.0d0, resonance_shift
    end do
    
    print *, ''
    print *, 'NEXT STEPS:'
    print *, '1. Implement thick orbit resonance finder'
    print *, '2. Account for finite orbit width in phase space integration'
    print *, '3. Include orbit classification (trapped/passing/potato)'
    print *, '4. Modify resonance width calculation for thick orbits'
    
contains

    subroutine calculate_reference_frequencies(v, eta, Om_th_thin, Om_ph_thin, &
                                             Om_th_thick, Om_ph_thick)
        ! Calculate frequencies with both methods for comparison
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: Om_th_thin, Om_ph_thin, Om_th_thick, Om_ph_thick
        
        real(8) :: taub, delphi, orbit_width_param
        logical :: success
        
        ! Thin orbit calculation (simplified)
        taub = 1.0d-4 / v * 1.0d6    ! Bounce time ~ 1/v
        delphi = 0.1d0 * eta         ! Toroidal shift proportional to pitch
        
        Om_th_thin = 2.0d0 * pi / taub
        Om_ph_thin = delphi / taub
        
        ! Thick orbit calculation (approximation for testing)
        call calculate_orbit_width_criterion(v, eta, orbit_width_param)
        call approximate_thick_frequencies(v, eta, orbit_width_param, Om_th_thick, Om_ph_thick, success)
        
        if (.not. success) then
            Om_th_thick = Om_th_thin
            Om_ph_thick = Om_ph_thin
        end if
        
    end subroutine calculate_reference_frequencies
    
    subroutine calculate_orbit_width_criterion(v, eta, criterion)
        ! Calculate orbit width parameter δr/L_B for thin/thick selection
        implicit none
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: criterion
        
        ! Physical parameters for ASDEX Upgrade-like conditions
        real(8), parameter :: B_field = 2.5d0       ! Tesla
        real(8), parameter :: mass_amu = 2.0d0      ! Deuterium mass (amu)
        real(8), parameter :: L_B = 0.5d0           ! Magnetic scale length (m)
        real(8), parameter :: v_thermal_ref = 1.0d6 ! Thermal velocity (m/s)
        
        real(8) :: rho_gyro
        
        ! Simplified gyroradius calculation: ρ = mv/(qB)
        rho_gyro = (v / v_thermal_ref) * 1.66d-27 * mass_amu * v_thermal_ref / (1.6d-19 * B_field)
        
        ! Orbit width parameter 
        criterion = rho_gyro / L_B
        
    end subroutine calculate_orbit_width_criterion
    
    subroutine approximate_thick_frequencies(v, eta, orbit_width_param, Om_th, Om_ph, success)
        ! Approximate thick orbit frequencies for testing
        implicit none
        real(8), intent(in) :: v, eta, orbit_width_param
        real(8), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        real(8) :: taub_thin, delphi_thin, correction_th, correction_ph
        
        ! Thin orbit baseline
        taub_thin = 1.0d-4 / v * 1.0d6
        delphi_thin = 0.1d0 * eta
        
        ! Finite orbit width corrections (physics-based)
        correction_th = 1.0d0 + orbit_width_param**2  ! Bounce frequency increases
        correction_ph = 1.0d0 + 0.5d0 * orbit_width_param**2  ! Toroidal frequency smaller effect
        
        Om_th = 2.0d0 * pi / taub_thin * correction_th
        Om_ph = delphi_thin / taub_thin * correction_ph
        
        success = .true.
        
    end subroutine approximate_thick_frequencies

end program test_thick_orbit_resonance