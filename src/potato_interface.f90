module potato_interface
    use orbit_types, only: orbit_type_t
    implicit none
    
    private
    public :: thick_orbit_type_t
    
    type, extends(orbit_type_t) :: thick_orbit_type_t
    contains
        procedure :: calculate_bounce_time => thick_calculate_bounce_time
        procedure :: calculate_frequencies => thick_calculate_frequencies
    end type thick_orbit_type_t
    
contains

    subroutine thick_calculate_bounce_time(this, v, eta, taub, bounceavg)
#ifndef USE_THICK_ORBITS
        use neort_orbit, only: bounce
#endif
        
        class(thick_orbit_type_t), intent(in) :: this
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: taub
        real(8), intent(out) :: bounceavg(:)
        
#ifdef USE_THICK_ORBITS
        ! TODO: Call POTATO's find_bounce() here
        ! For now, use placeholder values to test interface
        taub = 1.0d-6
        bounceavg = 1.0d0
#else
        ! Fallback to thin orbit if thick orbits not enabled
        call bounce(v, eta, taub, bounceavg)
#endif
        
    end subroutine thick_calculate_bounce_time
    
    subroutine thick_calculate_frequencies(this, eta, omega_theta, omega_phi)
        class(thick_orbit_type_t), intent(in) :: this
        real(8), intent(in) :: eta
        real(8), intent(out) :: omega_theta, omega_phi
        
#ifdef USE_THICK_ORBITS
        ! TODO: Call POTATO's frequency calculation here
        ! For now, use placeholder values to test interface
        omega_theta = 2.0d0
        omega_phi = 3.0d0
#else
        ! Fallback to thin orbit if thick orbits not enabled
        omega_theta = 1.0d0
        omega_phi = 1.0d0
#endif
        
    end subroutine thick_calculate_frequencies
    
end module potato_interface