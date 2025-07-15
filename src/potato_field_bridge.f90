module potato_field_bridge
    use iso_fortran_env, only: real64
    implicit none

    private
    public :: psif_neo_rt, dpsidr_neo_rt, dpsidz_neo_rt
    public :: field_eq, initialize_potato_field
    public :: convert_neort_to_potato
    public :: real_find_bounce_calculation

    integer, parameter :: dp = real64
    
    ! External POTATO spline routines
    interface
        subroutine s2dcut(nx, ny, hx, hy, f, imi, ima, jmi, jma, icount, spl, ipoint)
            import :: real64
            integer, intent(in) :: nx, ny
            real(real64), intent(in) :: hx, hy
            real(real64), intent(in) :: f(nx, ny)
            integer, intent(in) :: imi(ny), ima(ny), jmi(nx), jma(nx)
            integer, intent(out) :: icount
            real(real64), intent(out) :: spl(6, 6, *)
            integer, intent(out) :: ipoint(nx, ny)
        end subroutine s2dcut
        
        subroutine spline(nx, ny, x, y, hx, hy, icount, spl, ipoint, xb, yb, &
                         u, ux, uy, uxx, uxy, uyy, ierr)
            import :: real64
            integer, intent(in) :: nx, ny, icount
            real(real64), intent(in) :: x(nx), y(ny), hx, hy
            real(real64), intent(in) :: spl(6, 6, *)
            integer, intent(in) :: ipoint(nx, ny)
            real(real64), intent(in) :: xb, yb
            real(real64), intent(out) :: u, ux, uy, uxx, uxy, uyy
            integer, intent(out) :: ierr
        end subroutine spline
    end interface

contains

    subroutine psif_neo_rt(R, Z, psi_result)
        use field_interface, only: field_evaluator_t, field_evaluator_factory
        use runtime_config, only: get_use_thick_orbits
        implicit none
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: psi_result
        
        ! Local variables
        class(field_evaluator_t), allocatable :: evaluator
        real(dp) :: dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2
        logical :: use_thick_orbits
        
        ! Get runtime configuration
        use_thick_orbits = get_use_thick_orbits()
        
        ! Create appropriate field evaluator
        evaluator = field_evaluator_factory(use_thick_orbits)
        
        ! Evaluate field
        call evaluator%evaluate_field(R, Z, psi_result, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2)
        
        ! Cleanup
        call evaluator%cleanup()
    end subroutine psif_neo_rt

    subroutine dpsidr_neo_rt(R, Z, dpsidr_result)
        use field_interface, only: field_evaluator_t, field_evaluator_factory
        use runtime_config, only: get_use_thick_orbits
        implicit none
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: dpsidr_result
        
        ! Local variables
        class(field_evaluator_t), allocatable :: evaluator
        real(dp) :: psif, dpsidz, d2psidr2, d2psidrdz, d2psidz2
        logical :: use_thick_orbits
        
        ! Get runtime configuration
        use_thick_orbits = get_use_thick_orbits()
        
        ! Create appropriate field evaluator
        evaluator = field_evaluator_factory(use_thick_orbits)
        
        ! Evaluate field
        call evaluator%evaluate_field(R, Z, psif, dpsidr_result, dpsidz, d2psidr2, d2psidrdz, d2psidz2)
        
        ! Cleanup
        call evaluator%cleanup()
    end subroutine dpsidr_neo_rt

    subroutine dpsidz_neo_rt(R, Z, dpsidz_result)
        use field_interface, only: field_evaluator_t, field_evaluator_factory
        use runtime_config, only: get_use_thick_orbits
        implicit none
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: dpsidz_result
        
        ! Local variables
        class(field_evaluator_t), allocatable :: evaluator
        real(dp) :: psif, dpsidr, d2psidr2, d2psidrdz, d2psidz2
        logical :: use_thick_orbits
        
        ! Get runtime configuration
        use_thick_orbits = get_use_thick_orbits()
        
        ! Create appropriate field evaluator
        evaluator = field_evaluator_factory(use_thick_orbits)
        
        ! Evaluate field
        call evaluator%evaluate_field(R, Z, psif, dpsidr, dpsidz_result, d2psidr2, d2psidrdz, d2psidz2)
        
        ! Cleanup
        call evaluator%cleanup()
    end subroutine dpsidz_neo_rt

    subroutine field_eq(R, Z)
        use field_interface, only: field_evaluator_t, field_evaluator_factory
        use runtime_config, only: get_use_thick_orbits
        implicit none
        real(dp), intent(in) :: R, Z
        
        ! This is a legacy interface - should be replaced with direct field_interface usage
        ! For now, just use the field interface with runtime configuration
        class(field_evaluator_t), allocatable :: evaluator
        real(dp) :: psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2
        logical :: use_thick_orbits
        
        ! Get runtime configuration
        use_thick_orbits = get_use_thick_orbits()
        
        ! Create appropriate field evaluator
        evaluator = field_evaluator_factory(use_thick_orbits)
        
        ! Evaluate field
        call evaluator%evaluate_field(R, Z, psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2)
        
        ! Cleanup
        call evaluator%cleanup()
        
    end subroutine field_eq

    subroutine initialize_potato_field(success)
        implicit none
        logical, intent(out) :: success
        
        ! Initialize success flag
        success = .true.
        
        ! With the abstract interface approach, field initialization
        ! is handled by the field evaluator implementations
        ! No need for explicit grid setup here
        
    end subroutine initialize_potato_field
    
    subroutine convert_neort_to_potato(v, eta, s_flux, theta_boozer, phi_boozer, z_eqm, success)
        ! Convert NEO-RT parameters to POTATO phase space using REAL coordinate transformation
        use do_magfie_mod, only: booz_to_cyl, s, inp_swi
        use util, only: mi, qi, c
        use time_normalization, only: calculate_thermal_velocity
        implicit none
        real(dp), intent(in) :: v, eta, s_flux, theta_boozer, phi_boozer
        real(dp), intent(out) :: z_eqm(5)
        logical, intent(out) :: success
        
        ! Local variables for coordinate conversion
        real(dp) :: x_boozer(3), r_cyl(3)
        real(dp) :: s_saved, bmod_local, v_thermal, lambda_potato
        
        ! Save current flux surface coordinate
        s_saved = s
        
        ! Set flux surface for coordinate conversion
        s = s_flux
        
        ! Use REAL Boozer coordinate transformation
        x_boozer(1) = s_flux        ! Flux surface coordinate s ∈ [0,1]
        x_boozer(2) = phi_boozer    ! Boozer toroidal angle
        x_boozer(3) = theta_boozer  ! Boozer poloidal angle
        
        ! Convert Boozer coordinates to cylindrical coordinates
        call booz_to_cyl(x_boozer, r_cyl)
        
        ! Set POTATO spatial coordinates (NEO-RT uses cm, POTATO uses m)
        z_eqm(1) = r_cyl(1) / 100.0d0  ! R coordinate (cm -> m)
        z_eqm(2) = r_cyl(2)            ! phi coordinate (radians)
        z_eqm(3) = r_cyl(3) / 100.0d0  ! Z coordinate (cm -> m)
        
        success = .true.
        
        ! Convert NEO-RT (v, eta) to POTATO velocity space
        ! NEO-RT: eta = v_perp²/(v_perp² + v_parallel²), v = |v|
        ! POTATO: lambda = v_parallel/v, normalized momentum
        
        ! Convert to POTATO's momentum and pitch angle convention
        lambda_potato = sqrt(1.0d0 - eta)  ! cos(pitch_angle) for eta = μB/E
        
        ! Use proper thermal velocity calculation from time_normalization module
        v_thermal = calculate_thermal_velocity(1.0d0, 2.0d0)  ! 1 keV deuterium
        
        ! POTATO momentum normalization: p = v/v_thermal * sqrt(2)
        z_eqm(4) = v / v_thermal * sqrt(2.0d0)  ! Normalized momentum
        z_eqm(5) = lambda_potato                ! Cosine of pitch angle
        
        ! Restore original flux surface coordinate
        s = s_saved
        
        ! Validate physical constraints
        success = (v > 0.0d0) .and. (eta >= 0.0d0) .and. (eta <= 1.0d0) .and. &
                  (s_flux >= 0.0d0) .and. (s_flux <= 1.0d0) .and. &
                  (z_eqm(1) > 0.0d0) .and. (abs(lambda_potato) <= 1.0d0)
        
    end subroutine convert_neort_to_potato
    
    subroutine real_find_bounce_calculation(v, eta, taub, delphi, success)
        ! Real POTATO integration for thick orbit bounce calculation
        use orbit_interface, only: orbit_calculator_t, orbit_calculator_factory
        use runtime_config, only: get_use_thick_orbits
        implicit none
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: taub, delphi
        logical, intent(out) :: success
        
        ! Local variables
        class(orbit_calculator_t), allocatable :: calculator
        real(dp) :: extraset(7)
        logical :: use_thick_orbits
        real(dp) :: s_flux, theta_boozer, phi_boozer
        
        ! Get runtime configuration
        use_thick_orbits = get_use_thick_orbits()
        
        ! Use default test Boozer coordinates
        s_flux = 0.5d0
        theta_boozer = 0.0d0
        phi_boozer = 0.0d0
        
        ! Create appropriate orbit calculator
        calculator = orbit_calculator_factory(use_thick_orbits)
        
        ! Calculate bounce time
        call calculator%find_bounce(v, eta, s_flux, theta_boozer, phi_boozer, &
                                   taub, delphi, extraset, success)
        
        ! Cleanup
        call calculator%cleanup()
        
    end subroutine real_find_bounce_calculation
    
end module potato_field_bridge
