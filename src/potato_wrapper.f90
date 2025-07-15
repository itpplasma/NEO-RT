module potato_wrapper
    ! Wrapper module for POTATO integration
    ! This module provides the interface between NEO-RT and POTATO
    ! handling coordinate transformations and initialization
    
    ! NOTE: Actual POTATO integration currently disabled due to missing
    ! field evaluation functions (psif, dpsidr, dpsidz) in POTATO field_eq_mod
    ! 
    ! To complete integration:
    ! 1. Implement field evaluation functions connecting to NEO-RT magnetic field
    ! 2. Enable POTATO compilation in CMakeLists.txt
    ! 3. Replace stub calls with actual POTATO find_bounce calls
    
    implicit none
    
    private
    public :: potato_wrapper_init, potato_wrapper_find_bounce, &
              potato_wrapper_calculate_frequencies
    
    ! POTATO configuration parameters
    logical :: potato_initialized = .false.
    real(8) :: potato_dtau_default = 1.0d-10
    
contains
    
    subroutine potato_wrapper_init()
        ! Initialize POTATO modules
        ! Currently using stub implementation until actual POTATO integration completed
        ! 
        ! For actual POTATO integration, this would:
        ! 1. Initialize POTATO field_eq_mod with NEO-RT magnetic field interface
        ! 2. Set up orbit_dim_mod parameters
        ! 3. Configure POTATO orbit integration tolerances
        
        potato_initialized = .true.
        
        print *, 'POTATO wrapper initialized (stub mode - see CMakeLists.txt for integration status)'
        
    end subroutine potato_wrapper_init
    
    subroutine potato_wrapper_find_bounce(v, eta, taub, delphi, extraset)
        ! Wrapper for POTATO find_bounce with NEO-RT interface
        ! Handles coordinate transformation from NEO-RT to POTATO
        use potato_field_bridge, only: convert_neort_to_potato
        use runtime_config, only: get_use_thick_orbits
        
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: taub, delphi
        real(8), intent(inout) :: extraset(:)
        
        ! Local variables for POTATO interface
        integer :: next
        real(8) :: z_eqm(5)
        real(8) :: dtau_in
        logical :: success
        
        if (.not. potato_initialized) then
            call potato_wrapper_init()
        endif
        
        ! Transform NEO-RT (v, eta) to POTATO phase space coordinates z_eqm
        ! z_eqm(1:3) = spatial coordinates (R, phi, Z)  
        ! z_eqm(4) = momentum module p (normalized to thermal velocity)
        ! z_eqm(5) = lambda = cos(pitch angle)
        
        ! Use proper coordinate conversion through potato_field_bridge
        call convert_neort_to_potato(v, eta, 0.5d0, 0.0d0, 0.0d0, z_eqm, success)
        
        if (.not. success) then
            print *, 'ERROR: convert_neort_to_potato failed in potato_wrapper_find_bounce'
            print *, 'v=', v, 'eta=', eta
            stop
        end if
        
        next = size(extraset)
        dtau_in = potato_dtau_default
        
        ! Call actual POTATO find_bounce when thick orbits enabled
        if (get_use_thick_orbits()) then
            ! Use real POTATO find_bounce with velo subroutine
            external :: velo
            call find_bounce(next, velo, dtau_in, z_eqm, taub, delphi, extraset)
        else
            ! Use stub implementation for thin orbits
            call potato_stub_find_bounce(v, eta, taub, delphi, extraset)
        end if
        
    end subroutine potato_wrapper_find_bounce
    
    subroutine potato_wrapper_calculate_frequencies(taub, delphi, omega_bounce, omega_toroidal)
        ! Calculate canonical frequencies from POTATO bounce results with PROPER time normalization
        ! 
        ! From POTATO find_bounce output:
        ! - taub: bounce time in dimensionless units (tau = sqrt(2*T/m)*t)
        ! - delphi: toroidal shift per bounce time (Deltaphi_bounce)
        !
        ! Canonical frequencies:
        ! - omega_bounce: bounce frequency = 2*pi*sqrt(2*T/m)/taub
        ! - omega_toroidal: toroidal precession frequency = delphi*sqrt(2*T/m)/taub
        !   (Converted to physical frequencies using proper time normalization)
        
        use time_normalization, only: calculate_thermal_velocity, convert_frequency_to_physical
        
        real(8), intent(in) :: taub, delphi
        real(8), intent(out) :: omega_bounce, omega_toroidal
        
        ! Physical parameters for time normalization
        real(8), parameter :: temperature_keV = 1.0d0  ! 1 keV typical temperature
        real(8), parameter :: mass_number = 2.0d0      ! Deuterium mass number
        
        real(8) :: v_thermal
        
        ! Calculate thermal velocity for proper time normalization
        v_thermal = calculate_thermal_velocity(temperature_keV, mass_number)
        
        ! Convert dimensionless bounce time to physical bounce frequency
        omega_bounce = convert_frequency_to_physical(taub, v_thermal)
        
        ! Convert dimensionless toroidal precession to physical toroidal frequency
        ! omega_toroidal = delphi * v_thermal / taub
        omega_toroidal = delphi * v_thermal / taub
        
    end subroutine potato_wrapper_calculate_frequencies
    
    ! Temporary stub implementation until POTATO build fully integrated
    subroutine potato_stub_find_bounce(v, eta, taub, delphi, extraset)
        use potato_stub, only: potato_find_bounce
        
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: taub, delphi
        real(8), intent(inout) :: extraset(:)
        
        call potato_find_bounce(v, eta, taub, delphi, extraset)
        
    end subroutine potato_stub_find_bounce
    
end module potato_wrapper