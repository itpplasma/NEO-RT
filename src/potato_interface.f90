module potato_interface
    use orbit_types, only: orbit_type_t
    implicit none
    
    private
    public :: thick_orbit_type_t, thin_orbit_type_t
    
    type, extends(orbit_type_t) :: thick_orbit_type_t
    contains
        procedure :: calculate_bounce_time => thick_calculate_bounce_time
        procedure :: calculate_frequencies => thick_calculate_frequencies
    end type thick_orbit_type_t
    
    type, extends(orbit_type_t) :: thin_orbit_type_t
    contains
        procedure :: calculate_bounce_time => thin_calculate_bounce_time
        procedure :: calculate_frequencies => thin_calculate_frequencies
    end type thin_orbit_type_t
    
contains

    subroutine thick_calculate_bounce_time(this, v, eta, taub, bounceavg)
        use potato_wrapper, only: potato_wrapper_find_bounce
        
        class(thick_orbit_type_t), intent(in) :: this
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: taub
        real(8), intent(out) :: bounceavg(:)
        
        ! Call POTATO through wrapper interface
        real(8) :: delphi, extraset(size(bounceavg))
        extraset = bounceavg
        call potato_wrapper_find_bounce(v, eta, taub, delphi, extraset)
        bounceavg = extraset
        
    end subroutine thick_calculate_bounce_time
    
    subroutine thick_calculate_frequencies(this, eta, omega_theta, omega_phi)
        use potato_wrapper, only: potato_wrapper_find_bounce, potato_wrapper_calculate_frequencies
        
        class(thick_orbit_type_t), intent(in) :: this
        real(8), intent(in) :: eta
        real(8), intent(out) :: omega_theta, omega_phi
        
        ! Calculate frequencies using POTATO wrapper
        real(8) :: v_dummy, taub, delphi, extraset(7)
        real(8) :: omega_bounce, omega_toroidal
        
        v_dummy = 1.0d6  ! Dummy velocity for frequency calculation
        extraset = 0.0d0
        call potato_wrapper_find_bounce(v_dummy, eta, taub, delphi, extraset)
        call potato_wrapper_calculate_frequencies(taub, delphi, omega_bounce, omega_toroidal)
        
        omega_theta = omega_bounce
        omega_phi = omega_toroidal
        
    end subroutine thick_calculate_frequencies
    
    subroutine thin_calculate_bounce_time(this, v, eta, taub, bounceavg)
        use neort_orbit, only: bounce
        
        class(thin_orbit_type_t), intent(in) :: this
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: taub
        real(8), intent(out) :: bounceavg(:)
        
        call bounce(v, eta, taub, bounceavg)
        
    end subroutine thin_calculate_bounce_time
    
    subroutine thin_calculate_frequencies(this, eta, omega_theta, omega_phi)
        use neort_freq, only: Om_th, Om_ph
        use neort_profiles, only: vth
        
        class(thin_orbit_type_t), intent(in) :: this
        real(8), intent(in) :: eta
        real(8), intent(out) :: omega_theta, omega_phi
        
        real(8) :: dOmthdv, dOmthdeta, dOmphdv, dOmphdeta
        
        call Om_th(vth, eta, omega_theta, dOmthdv, dOmthdeta)
        call Om_ph(vth, eta, omega_phi, dOmphdv, dOmphdeta)
        
    end subroutine thin_calculate_frequencies
    
end module potato_interface