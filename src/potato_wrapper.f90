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
        
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: taub, delphi
        real(8), intent(inout) :: extraset(:)
        
        ! Local variables for POTATO interface
        integer :: next
        real(8) :: z_eqm(5)
        real(8) :: dtau_in
        
        if (.not. potato_initialized) then
            call potato_wrapper_init()
        endif
        
        ! Transform NEO-RT (v, eta) to POTATO phase space coordinates z_eqm
        ! z_eqm(1:3) = spatial coordinates (R, phi, Z)
        ! z_eqm(4) = momentum module p (normalized)
        ! z_eqm(5) = lambda = cos(pitch angle)
        
        ! For stub implementation, use simplified mapping
        z_eqm(1) = 1.0d0  ! R - would come from orbit starting point
        z_eqm(2) = 0.0d0  ! phi
        z_eqm(3) = 0.0d0  ! Z
        z_eqm(4) = v / 1.0d6  ! Normalized momentum (rough scaling)
        z_eqm(5) = sqrt(1.0d0 - eta)  ! lambda from eta
        
        next = size(extraset)
        dtau_in = potato_dtau_default
        
        ! TODO: Call actual POTATO find_bounce when integration complete
        ! call find_bounce(next, velo, dtau_in, z_eqm, taub, delphi, extraset)
        
        ! Currently using stub implementation
        call potato_stub_find_bounce(v, eta, taub, delphi, extraset)
        
    end subroutine potato_wrapper_find_bounce
    
    subroutine potato_wrapper_calculate_frequencies(taub, delphi, omega_bounce, omega_toroidal)
        ! Calculate canonical frequencies from POTATO bounce results
        
        real(8), intent(in) :: taub, delphi
        real(8), intent(out) :: omega_bounce, omega_toroidal
        
        real(8), parameter :: pi = 3.14159265358979323846d0
        
        ! Bounce frequency: omega_b = 2*pi / taub
        omega_bounce = 2.0d0 * pi / taub
        
        ! Toroidal precession frequency: omega_phi = delphi / taub
        omega_toroidal = delphi / taub
        
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