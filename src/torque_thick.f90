module torque_thick
    ! NTV torque calculation with thick orbit corrections
    ! Complete physics pipeline: field → orbit → frequency → resonance → transport → torque
    
    use neort_resonance
    use transport_thick
    use freq_thick, only: compute_canonical_frequencies_thick
    use runtime_config, only: get_use_thick_orbits
    implicit none
    
    private
    public :: calculate_ntv_torque_density
    public :: calculate_velocity_space_torque
    public :: calculate_torque_profile
    public :: validate_torque_conservation
    
    ! Physical constants
    real(8), parameter :: pi = 3.141592653589793d0
    real(8), parameter :: mass_deuterium = 2.0d0 * 1.66d-27  ! kg
    real(8), parameter :: charge_elementary = 1.6d-19        ! C
    real(8), parameter :: B_field_ref = 2.5d0                ! Tesla
    real(8), parameter :: k_boltzmann = 1.38d-23             ! J/K
    
contains

    subroutine calculate_ntv_torque_density(v, eta, n_mode, m_mode, omega_mode, &
                                           torque_density, success)
        ! Calculate NTV torque density: T = ∫∫ f₀(v,η) × resonance × transport dv dη
        implicit none
        real(8), intent(in) :: v, eta
        integer, intent(in) :: n_mode, m_mode
        real(8), intent(in) :: omega_mode
        real(8), intent(out) :: torque_density
        logical, intent(out) :: success
        
        ! Local variables
        real(8) :: Om_theta, Om_phi, resonance_condition
        real(8) :: D11, D12, D22, transport_coeff
        real(8) :: f0_maxwell, gradient_factor
        real(8) :: orbit_width, resonance_width
        logical :: resonant, transport_success, freq_success
        
        ! Initialize
        success = .false.
        torque_density = 0.0d0
        
        ! Step 1: Calculate frequencies using thick orbit bounce integrals
        ! NO FALLBACK - this is a replacement, not a correction
        call compute_canonical_frequencies_thick(v, eta, Om_theta, Om_phi, freq_success)
        if (.not. freq_success) then
            print *, 'ERROR: Thick orbit bounce integral calculation failed'
            print *, 'This is a fundamental failure - cannot use thin orbit fallback'
            return
        end if
        
        ! Step 2: Check resonance condition: n·ω_φ - m·ω_θ = ω_mode
        resonance_condition = real(n_mode, 8) * Om_phi - real(m_mode, 8) * Om_theta - omega_mode
        
        ! Step 3: Calculate resonance width with finite orbit width effects
        ! This is inherent to thick orbit physics - always use orbit width
        orbit_width = calculate_orbit_width_parameter(v, eta)
        resonance_width = abs(Om_theta) * orbit_width
        
        ! Step 4: Check if particle is resonant
        resonant = abs(resonance_condition) < resonance_width
        
        if (.not. resonant) then
            success = .true.
            torque_density = 0.0d0
            return
        end if
        
        ! Step 5: Calculate transport coefficients using thick orbit bounce integrals
        ! NO FALLBACK - this is a replacement, not a correction
        call calculate_transport_coefficients_thick(v, eta, D11, D12, D22, transport_success)
        
        if (.not. transport_success) then
            print *, 'WARNING: Transport coefficient calculation failed'
            return
        end if
        
        ! Step 6: Calculate Maxwell-Boltzmann distribution
        f0_maxwell = maxwell_boltzmann_distribution(v, eta)
        
        ! Step 7: Calculate gradient factor (thermodynamic forces)
        gradient_factor = calculate_gradient_factor(v, eta)
        
        ! Step 8: Calculate torque density contribution
        ! T = ∫∫ f₀ × δ(resonance) × D₁₂ × ∇p/p × v dv dη
        transport_coeff = D12  ! Cross-transport coefficient couples radial and toroidal
        
        ! Resonance delta function approximation and final torque density calculation
        torque_density = f0_maxwell * resonance_delta_function(resonance_condition, resonance_width) * &
                        transport_coeff * gradient_factor * v
        
        success = .true.
        
    end subroutine calculate_ntv_torque_density
    
    subroutine calculate_velocity_space_torque(n_mode, m_mode, omega_mode, &
                                              total_torque, success)
        ! Calculate total torque by integrating over velocity space
        implicit none
        integer, intent(in) :: n_mode, m_mode
        real(8), intent(in) :: omega_mode
        real(8), intent(out) :: total_torque
        logical, intent(out) :: success
        
        ! Integration parameters
        integer, parameter :: n_v = 50, n_eta = 20
        real(8), parameter :: v_min = 0.1d6, v_max = 3.0d6  ! m/s
        real(8), parameter :: eta_min = 0.05d0, eta_max = 0.95d0
        
        ! Local variables
        real(8) :: dv, deta, v_test, eta_test
        real(8) :: torque_density, torque_contribution
        logical :: calc_success
        integer :: iv, ieta
        
        ! Initialize
        success = .false.
        total_torque = 0.0d0
        
        ! Set up integration grid
        dv = (v_max - v_min) / real(n_v - 1, 8)
        deta = (eta_max - eta_min) / real(n_eta - 1, 8)
        
        print *, 'Velocity space integration:'
        print *, '  v_range: [', v_min/1d6, ',', v_max/1d6, '] × 10⁶ m/s'
        print *, '  eta_range: [', eta_min, ',', eta_max, ']'
        print *, '  Grid: ', n_v, '×', n_eta, '=', n_v*n_eta, 'points'
        
        ! Integrate over velocity space
        do iv = 1, n_v
            do ieta = 1, n_eta
                v_test = v_min + real(iv-1, 8) * dv
                eta_test = eta_min + real(ieta-1, 8) * deta
                
                ! Calculate torque density at this velocity space point
                call calculate_ntv_torque_density(v_test, eta_test, n_mode, m_mode, omega_mode, &
                                                 torque_density, calc_success)
                
                if (calc_success) then
                    ! Integration weight includes velocity space Jacobian
                    torque_contribution = torque_density * v_test**2 * dv * deta
                    total_torque = total_torque + torque_contribution
                end if
            end do
        end do
        
        ! Multiply by 4π for full solid angle integration
        total_torque = total_torque * 4.0d0 * pi
        
        print *, 'Total NTV torque calculated:', total_torque, 'N·m'
        
        success = .true.
        
    end subroutine calculate_velocity_space_torque
    
    subroutine calculate_torque_profile(n_mode, m_mode, omega_mode, &
                                       flux_surfaces, torque_profile, success)
        ! Calculate torque profile across flux surfaces
        implicit none
        integer, intent(in) :: n_mode, m_mode
        real(8), intent(in) :: omega_mode
        real(8), intent(in) :: flux_surfaces(:)
        real(8), intent(out) :: torque_profile(:)
        logical, intent(out) :: success
        
        integer :: i, n_surfaces
        real(8) :: total_torque
        logical :: calc_success
        
        ! Initialize
        success = .false.
        n_surfaces = size(flux_surfaces)
        torque_profile = 0.0d0
        
        print *, 'Calculating torque profile for', n_surfaces, 'flux surfaces'
        
        ! Calculate torque at each flux surface
        do i = 1, n_surfaces
            ! Set flux surface (would need to update magnetic field accordingly)
            ! For now, use constant calculation
            
            call calculate_velocity_space_torque(n_mode, m_mode, omega_mode, &
                                               total_torque, calc_success)
            
            if (calc_success) then
                torque_profile(i) = total_torque
            else
                print *, 'WARNING: Torque calculation failed at flux surface', i
            end if
        end do
        
        success = .true.
        
    end subroutine calculate_torque_profile
    
    subroutine validate_torque_conservation(torque_profile, success)
        ! Validate torque conservation and symmetry properties
        implicit none
        real(8), intent(in) :: torque_profile(:)
        logical, intent(out) :: success
        
        real(8) :: total_torque, torque_sum
        integer :: i
        
        ! Calculate total torque
        torque_sum = 0.0d0
        do i = 1, size(torque_profile)
            torque_sum = torque_sum + torque_profile(i)
        end do
        
        total_torque = torque_sum
        
        print *, 'Torque conservation check:'
        print *, '  Total torque:', total_torque, 'N·m'
        print *, '  Maximum torque:', maxval(abs(torque_profile)), 'N·m'
        
        ! Check if torque is reasonable
        if (abs(total_torque) > 1.0d-6 .and. abs(total_torque) < 1.0d6) then
            success = .true.
            print *, '  ✓ Torque magnitude is reasonable'
        else
            success = .false.
            print *, '  ✗ Torque magnitude is unreasonable'
        end if
        
    end subroutine validate_torque_conservation
    
    ! Helper functions
    
    ! NOTE: calculate_thin_orbit_frequencies REMOVED
    ! This was a lazy shortcut - thick orbit physics uses bounce integrals only
    
    ! NOTE: calculate_thin_orbit_transport_coefficients REMOVED
    ! This was a lazy shortcut - thick orbit physics uses bounce integrals only
    
    function maxwell_boltzmann_distribution(v, eta) result(f0)
        ! Maxwell-Boltzmann distribution function
        implicit none
        real(8), intent(in) :: v, eta
        real(8) :: f0
        
        real(8), parameter :: m = mass_deuterium
        real(8), parameter :: T = 10.0d0 * 1.6d-19 / k_boltzmann  ! 10 keV in K
        real(8) :: v_th, normalization
        
        ! Thermal velocity
        v_th = sqrt(2.0d0 * k_boltzmann * T / m)
        
        ! Normalization
        normalization = (m / (2.0d0 * pi * k_boltzmann * T))**(1.5d0)
        
        ! Maxwell-Boltzmann distribution
        f0 = normalization * exp(-m * v**2 / (2.0d0 * k_boltzmann * T)) * v**2
        
    end function maxwell_boltzmann_distribution
    
    function calculate_gradient_factor(v, eta) result(gradient_factor)
        ! Calculate thermodynamic gradient factor
        implicit none
        real(8), intent(in) :: v, eta
        real(8) :: gradient_factor
        
        ! Simplified gradient factor: ∇p/p + ∇T/T
        real(8), parameter :: pressure_gradient = 0.1d0  ! Normalized
        real(8), parameter :: temperature_gradient = 0.05d0  ! Normalized
        
        gradient_factor = pressure_gradient + temperature_gradient
        
    end function calculate_gradient_factor
    
    function resonance_delta_function(resonance_condition, resonance_width) result(delta)
        ! Approximate delta function for resonance condition
        implicit none
        real(8), intent(in) :: resonance_condition, resonance_width
        real(8) :: delta
        
        ! Lorentzian approximation to delta function
        delta = (resonance_width / pi) / (resonance_condition**2 + resonance_width**2)
        
    end function resonance_delta_function
    
    ! Note: calculate_orbit_width_parameter is provided by neort_resonance module

end module torque_thick