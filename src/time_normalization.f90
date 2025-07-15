module time_normalization
    ! Time normalization module for NEO-RT/POTATO interface
    ! Handles conversion between NEO-RT physical time and POTATO dimensionless time
    ! 
    ! POTATO uses dimensionless time: tau = sqrt(2*T/m) * t
    ! where T is temperature and m is particle mass
    
    use util, only: mi, qi, pi
    implicit none
    
    private
    public :: calculate_thermal_velocity, convert_time_to_physical, &
              convert_time_to_dimensionless, convert_frequency_to_physical
    
    ! Physical constants (SI units)
    real(8), parameter :: electron_mass = 9.1094d-31    ! kg
    real(8), parameter :: proton_mass = 1.6726d-27      ! kg
    real(8), parameter :: elementary_charge = 1.6022d-19 ! C
    real(8), parameter :: keV_to_joule = 1.6022d-16     ! J/keV
    
contains

    function calculate_thermal_velocity(temperature_keV, mass_number) result(v_thermal)
        ! Calculate thermal velocity sqrt(2*T/m) for given temperature and mass
        real(8), intent(in) :: temperature_keV  ! Temperature in keV
        real(8), intent(in) :: mass_number      ! Mass number (2.0 for deuterium)
        real(8) :: v_thermal                    ! Thermal velocity in m/s
        
        real(8) :: temperature_J, mass_kg
        
        ! Convert temperature from keV to Joules
        temperature_J = temperature_keV * keV_to_joule
        
        ! Calculate particle mass in kg
        mass_kg = mass_number * proton_mass
        
        ! Calculate thermal velocity: v_thermal = sqrt(2*T/m)
        v_thermal = sqrt(2.0d0 * temperature_J / mass_kg)
        
    end function calculate_thermal_velocity
    
    function convert_time_to_physical(tau_dimensionless, v_thermal) result(t_physical)
        ! Convert POTATO dimensionless time to physical time
        ! tau = sqrt(2*T/m) * t  =>  t = tau / sqrt(2*T/m)
        real(8), intent(in) :: tau_dimensionless  ! Dimensionless time from POTATO
        real(8), intent(in) :: v_thermal          ! Thermal velocity sqrt(2*T/m)
        real(8) :: t_physical                     ! Physical time in seconds
        
        t_physical = tau_dimensionless / v_thermal
        
    end function convert_time_to_physical
    
    function convert_time_to_dimensionless(t_physical, v_thermal) result(tau_dimensionless)
        ! Convert physical time to POTATO dimensionless time
        ! tau = sqrt(2*T/m) * t
        real(8), intent(in) :: t_physical         ! Physical time in seconds
        real(8), intent(in) :: v_thermal          ! Thermal velocity sqrt(2*T/m)
        real(8) :: tau_dimensionless              ! Dimensionless time for POTATO
        
        tau_dimensionless = t_physical * v_thermal
        
    end function convert_time_to_dimensionless
    
    function convert_frequency_to_physical(tau_bounce_dimensionless, v_thermal) result(omega_physical)
        ! Convert POTATO dimensionless bounce time to physical frequency
        ! omega = 2*pi*sqrt(2*T/m) / tau_bounce
        real(8), intent(in) :: tau_bounce_dimensionless  ! Dimensionless bounce time
        real(8), intent(in) :: v_thermal                 ! Thermal velocity sqrt(2*T/m)
        real(8) :: omega_physical                         ! Physical frequency in rad/s
        
        omega_physical = 2.0d0 * pi * v_thermal / tau_bounce_dimensionless
        
    end function convert_frequency_to_physical
    
end module time_normalization