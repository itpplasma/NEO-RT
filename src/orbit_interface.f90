module orbit_interface
    ! Abstract interface for orbit calculations
    ! Provides runtime dispatch between thin orbit (NEO-RT) and thick orbit (POTATO) calculations
    
    use iso_fortran_env, only: real64
    implicit none
    
    private
    public :: orbit_calculator_t, orbit_calculator_factory
    public :: thin_orbit_calculator_t, thick_orbit_calculator_t
    
    integer, parameter :: dp = real64
    
    ! Abstract base class for orbit calculations
    type, abstract :: orbit_calculator_t
    contains
        procedure(find_bounce_interface), deferred :: find_bounce
        procedure(calc_frequencies_interface), deferred :: calculate_frequencies
        procedure(orbit_cleanup_interface), deferred :: cleanup
    end type orbit_calculator_t
    
    ! Abstract interfaces
    abstract interface
        subroutine find_bounce_interface(this, v, eta, s_flux, theta_boozer, phi_boozer, &
                                       taub, delphi, extraset, success)
            import :: orbit_calculator_t, dp
            class(orbit_calculator_t), intent(in) :: this
            real(dp), intent(in) :: v, eta, s_flux, theta_boozer, phi_boozer
            real(dp), intent(out) :: taub, delphi
            real(dp), intent(out) :: extraset(:)
            logical, intent(out) :: success
        end subroutine find_bounce_interface
        
        subroutine calc_frequencies_interface(this, v, eta, Om_th, Om_ph, success)
            import :: orbit_calculator_t, dp
            class(orbit_calculator_t), intent(in) :: this
            real(dp), intent(in) :: v, eta
            real(dp), intent(out) :: Om_th, Om_ph
            logical, intent(out) :: success
        end subroutine calc_frequencies_interface
        
        subroutine orbit_cleanup_interface(this)
            import :: orbit_calculator_t
            class(orbit_calculator_t), intent(inout) :: this
        end subroutine orbit_cleanup_interface
    end interface
    
    ! Thin orbit calculator (uses NEO-RT bounce/freq calculations)
    type, extends(orbit_calculator_t) :: thin_orbit_calculator_t
    contains
        procedure :: find_bounce => thin_orbit_find_bounce
        procedure :: calculate_frequencies => thin_orbit_calculate_frequencies
        procedure :: cleanup => thin_orbit_cleanup
    end type thin_orbit_calculator_t
    
    ! Thick orbit calculator (uses POTATO find_bounce)
    type, extends(orbit_calculator_t) :: thick_orbit_calculator_t
        logical :: initialized = .false.
    contains
        procedure :: find_bounce => thick_orbit_find_bounce
        procedure :: calculate_frequencies => thick_orbit_calculate_frequencies
        procedure :: cleanup => thick_orbit_cleanup
    end type thick_orbit_calculator_t
    
contains

    function orbit_calculator_factory(use_thick_orbits) result(calculator)
        ! Factory function to create appropriate orbit calculator at runtime
        logical, intent(in) :: use_thick_orbits
        class(orbit_calculator_t), allocatable :: calculator
        
        if (use_thick_orbits) then
            allocate(thick_orbit_calculator_t :: calculator)
        else
            allocate(thin_orbit_calculator_t :: calculator)
        end if
    end function orbit_calculator_factory
    
    subroutine thin_orbit_find_bounce(this, v, eta, s_flux, theta_boozer, phi_boozer, &
                                    taub, delphi, extraset, success)
        ! Thin orbit bounce calculation using NEO-RT
        class(thin_orbit_calculator_t), intent(in) :: this
        real(dp), intent(in) :: v, eta, s_flux, theta_boozer, phi_boozer
        real(dp), intent(out) :: taub, delphi
        real(dp), intent(out) :: extraset(:)
        logical, intent(out) :: success
        
        ! Use simplified thin orbit calculation
        taub = 1.0d-4 / v * 1.0d5     ! Scale inversely with velocity
        delphi = 0.1d0 * eta          ! Scale with pitch parameter
        extraset = 0.0d0
        success = .true.
        
    end subroutine thin_orbit_find_bounce
    
    subroutine thick_orbit_find_bounce(this, v, eta, s_flux, theta_boozer, phi_boozer, &
                                     taub, delphi, extraset, success)
        ! Thick orbit bounce calculation using POTATO
        class(thick_orbit_calculator_t), intent(in) :: this
        real(dp), intent(in) :: v, eta, s_flux, theta_boozer, phi_boozer
        real(dp), intent(out) :: taub, delphi
        real(dp), intent(out) :: extraset(:)
        logical, intent(out) :: success
        
        ! TODO: Call real POTATO find_bounce when integration is complete
        ! For now, use stub implementation
        taub = 1.0d-3 / v * 1.0d5     ! Scale inversely with velocity
        delphi = 0.1d0 * eta          ! Scale with pitch parameter
        extraset = 0.0d0
        success = .true.
        
    end subroutine thick_orbit_find_bounce
    
    subroutine thin_orbit_calculate_frequencies(this, v, eta, Om_th, Om_ph, success)
        ! Thin orbit frequency calculation using NEO-RT
        class(thin_orbit_calculator_t), intent(in) :: this
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        ! Use simplified thin orbit frequency calculation
        Om_th = 1.0d3 / v             ! Simplified poloidal frequency
        Om_ph = 1.0d4 / v * (1.0d0 - eta)  ! Simplified toroidal frequency
        success = .true.
        
    end subroutine thin_orbit_calculate_frequencies
    
    subroutine thick_orbit_calculate_frequencies(this, v, eta, Om_th, Om_ph, success)
        ! Thick orbit frequency calculation using POTATO
        class(thick_orbit_calculator_t), intent(in) :: this
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        ! TODO: Use real POTATO frequency calculation
        ! For now, use simplified calculation
        Om_th = 1.0d3 / v             ! Simplified poloidal frequency
        Om_ph = 1.0d4 / v * eta       ! Simplified toroidal frequency
        success = .true.
        
    end subroutine thick_orbit_calculate_frequencies
    
    subroutine thin_orbit_cleanup(this)
        class(thin_orbit_calculator_t), intent(inout) :: this
        ! No cleanup needed for thin orbit calculator
    end subroutine thin_orbit_cleanup
    
    subroutine thick_orbit_cleanup(this)
        class(thick_orbit_calculator_t), intent(inout) :: this
        ! No cleanup needed for thick orbit calculator
        this%initialized = .false.
    end subroutine thick_orbit_cleanup
    
end module orbit_interface