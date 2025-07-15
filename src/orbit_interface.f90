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
        ! Thin orbit bounce calculation using REAL NEO-RT bounce calculation
        use neort_orbit, only: bounce_time, bounce
        use driftorbit, only: init_done
        class(thin_orbit_calculator_t), intent(in) :: this
        real(dp), intent(in) :: v, eta, s_flux, theta_boozer, phi_boozer
        real(dp), intent(out) :: taub, delphi
        real(dp), intent(out) :: extraset(:)
        logical, intent(out) :: success
        
        ! Local variables for real NEO-RT calculation
        real(dp) :: bounceavg(7)
        
        ! Check if driftorbit is initialized
        if (.not. init_done) then
            success = .false.
            taub = 0.0d0
            delphi = 0.0d0
            extraset = 0.0d0
            return
        end if
        
        ! Call REAL NEO-RT bounce calculation
        call bounce(v, eta, taub, bounceavg)
        
        ! Extract toroidal shift from bounce averages
        ! bounceavg(3) contains the toroidal velocity v_ph
        ! delphi = <v_ph> * taub gives toroidal shift per bounce
        delphi = bounceavg(3) * taub
        
        ! Set extraset to bounce averages for compatibility
        if (size(extraset) >= 7) then
            extraset(1:7) = bounceavg(1:7)
        else
            extraset = bounceavg(1:size(extraset))
        end if
        
        success = .true.
        
    end subroutine thin_orbit_find_bounce
    
    subroutine thick_orbit_find_bounce(this, v, eta, s_flux, theta_boozer, phi_boozer, &
                                     taub, delphi, extraset, success)
        ! Thick orbit bounce calculation using POTATO
        use potato_stub, only: potato_find_bounce
        class(thick_orbit_calculator_t), intent(in) :: this
        real(dp), intent(in) :: v, eta, s_flux, theta_boozer, phi_boozer
        real(dp), intent(out) :: taub, delphi
        real(dp), intent(out) :: extraset(:)
        logical, intent(out) :: success
        
        ! Call POTATO find_bounce directly (stub for now)
        call potato_find_bounce(v, eta, taub, delphi, extraset)
        
        ! POTATO stub always succeeds for now
        success = .true.
        
    end subroutine thick_orbit_find_bounce
    
    subroutine thin_orbit_calculate_frequencies(this, v, eta, Om_th, Om_ph, success)
        ! Thin orbit frequency calculation using REAL NEO-RT frequency calculation
        use neort_freq, only: om_th_func => om_th, om_ph_func => om_ph
        use driftorbit, only: init_done, sign_vpar_htheta
        class(thin_orbit_calculator_t), intent(in) :: this
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        ! Local variables for real NEO-RT calculation
        real(dp) :: dOmthdv, dOmthdeta, dOmphdv, dOmphdeta
        
        ! Check if driftorbit is initialized
        if (.not. init_done) then
            success = .false.
            Om_th = 0.0d0
            Om_ph = 0.0d0
            return
        end if
        
        ! Call REAL NEO-RT frequency calculations
        call om_th_func(v, eta, Om_th, dOmthdv, dOmthdeta)
        call om_ph_func(v, eta, Om_ph, dOmphdv, dOmphdeta)
        
        ! Apply sign convention for consistency
        Om_th = Om_th * sign_vpar_htheta
        
        success = .true.
        
    end subroutine thin_orbit_calculate_frequencies
    
    subroutine thick_orbit_calculate_frequencies(this, v, eta, Om_th, Om_ph, success)
        ! Thick orbit frequency calculation using POTATO
        use potato_stub, only: potato_find_bounce, potato_calculate_frequencies
        class(thick_orbit_calculator_t), intent(in) :: this
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: Om_th, Om_ph
        logical, intent(out) :: success
        
        ! Local variables
        real(dp) :: taub, delphi, extraset(7)
        real(dp) :: omega_bounce, omega_toroidal
        
        ! Get bounce time from POTATO
        call potato_find_bounce(v, eta, taub, delphi, extraset)
        
        ! Calculate frequencies from bounce results
        call potato_calculate_frequencies(taub, delphi, omega_bounce, omega_toroidal)
        
        ! Map to NEO-RT frequency convention
        Om_th = omega_bounce      ! Poloidal bounce frequency
        Om_ph = omega_toroidal    ! Toroidal precession frequency
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