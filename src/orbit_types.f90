module orbit_types
    implicit none
    
    private
    public :: orbit_type_t, thin_orbit_type_t
    
    type, abstract :: orbit_type_t
    contains
        procedure(calculate_bounce_time_interface), deferred :: calculate_bounce_time
        procedure(calculate_frequencies_interface), deferred :: calculate_frequencies
    end type orbit_type_t
    
    type, extends(orbit_type_t) :: thin_orbit_type_t
    contains
        procedure :: calculate_bounce_time => thin_calculate_bounce_time
        procedure :: calculate_frequencies => thin_calculate_frequencies
    end type thin_orbit_type_t
    
    abstract interface
        subroutine calculate_bounce_time_interface(this, v, eta, taub, bounceavg)
            import :: orbit_type_t
            class(orbit_type_t), intent(in) :: this
            real(8), intent(in) :: v, eta
            real(8), intent(out) :: taub
            real(8), intent(out) :: bounceavg(:)
        end subroutine calculate_bounce_time_interface
        
        subroutine calculate_frequencies_interface(this, eta, omega_theta, omega_phi)
            import :: orbit_type_t
            class(orbit_type_t), intent(in) :: this
            real(8), intent(in) :: eta
            real(8), intent(out) :: omega_theta, omega_phi
        end subroutine calculate_frequencies_interface
    end interface
    
contains

    subroutine thin_calculate_bounce_time(this, v, eta, taub, bounceavg)
        use neort_orbit, only: bounce, nvar
        
        class(thin_orbit_type_t), intent(in) :: this
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: taub
        real(8), intent(out) :: bounceavg(:)
        
        call bounce(v, eta, taub, bounceavg)
    end subroutine thin_calculate_bounce_time
    
    subroutine thin_calculate_frequencies(this, eta, omega_theta, omega_phi)
        class(thin_orbit_type_t), intent(in) :: this
        real(8), intent(in) :: eta
        real(8), intent(out) :: omega_theta, omega_phi
        
        omega_theta = 1.0d0
        omega_phi = 1.0d0
    end subroutine thin_calculate_frequencies
    
end module orbit_types