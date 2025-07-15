module field_interface
    ! Abstract interface for magnetic field evaluation
    ! Provides runtime dispatch between thin orbit (NEO-RT) and thick orbit (POTATO) field evaluation
    
    use iso_fortran_env, only: real64
    implicit none
    
    private
    public :: field_evaluator_t, field_evaluator_factory
    public :: thin_orbit_field_evaluator_t, thick_orbit_field_evaluator_t
    
    integer, parameter :: dp = real64
    
    ! Abstract base class for field evaluation
    type, abstract :: field_evaluator_t
    contains
        procedure(field_eval_interface), deferred :: evaluate_field
        procedure(field_cleanup_interface), deferred :: cleanup
    end type field_evaluator_t
    
    ! Abstract interfaces
    abstract interface
        subroutine field_eval_interface(this, R, Z, psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2)
            import :: field_evaluator_t, dp
            class(field_evaluator_t), intent(in) :: this
            real(dp), intent(in) :: R, Z
            real(dp), intent(out) :: psif, dpsidr, dpsidz
            real(dp), intent(out) :: d2psidr2, d2psidrdz, d2psidz2
        end subroutine field_eval_interface
        
        subroutine field_cleanup_interface(this)
            import :: field_evaluator_t
            class(field_evaluator_t), intent(inout) :: this
        end subroutine field_cleanup_interface
    end interface
    
    ! Thin orbit field evaluator (uses NEO-RT magfie)
    type, extends(field_evaluator_t) :: thin_orbit_field_evaluator_t
    contains
        procedure :: evaluate_field => thin_orbit_evaluate_field
        procedure :: cleanup => thin_orbit_cleanup
    end type thin_orbit_field_evaluator_t
    
    ! Thick orbit field evaluator (uses POTATO field_eq)
    type, extends(field_evaluator_t) :: thick_orbit_field_evaluator_t
        logical :: initialized = .false.
    contains
        procedure :: evaluate_field => thick_orbit_evaluate_field
        procedure :: cleanup => thick_orbit_cleanup
    end type thick_orbit_field_evaluator_t
    
contains

    function field_evaluator_factory(use_thick_orbits) result(evaluator)
        ! Factory function to create appropriate field evaluator at runtime
        logical, intent(in) :: use_thick_orbits
        class(field_evaluator_t), allocatable :: evaluator
        
        if (use_thick_orbits) then
            allocate(thick_orbit_field_evaluator_t :: evaluator)
        else
            allocate(thin_orbit_field_evaluator_t :: evaluator)
        end if
    end function field_evaluator_factory
    
    subroutine thin_orbit_evaluate_field(this, R, Z, psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2)
        ! Thin orbit field evaluation using NEO-RT magfie
        use do_magfie_mod, only: do_magfie, psi_pr
        class(thin_orbit_field_evaluator_t), intent(in) :: this
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: psif, dpsidr, dpsidz
        real(dp), intent(out) :: d2psidr2, d2psidrdz, d2psidz2
        
        ! Use REAL NEO-RT field evaluation through do_magfie call
        real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp), parameter :: delta = 1.0d-6  ! Small step for finite differences
        real(dp) :: psi_plus, psi_minus
        
        ! Set up coordinate input for magfie (Boozer coordinates)
        x(1) = R
        x(2) = 0.0d0  ! phi
        x(3) = Z
        
        ! Call magfie to get field information
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        
        ! Extract psi from magfie result
        psif = psi_pr  ! This is the poloidal flux from magfie
        
        ! Calculate dpsidr
        x(1) = R + delta
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        psi_plus = psi_pr
        
        x(1) = R - delta
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        psi_minus = psi_pr
        
        dpsidr = (psi_plus - psi_minus) / (2.0d0 * delta)
        
        ! Calculate dpsidz  
        x(1) = R
        x(3) = Z + delta
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        psi_plus = psi_pr
        
        x(3) = Z - delta
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        psi_minus = psi_pr
        
        dpsidz = (psi_plus - psi_minus) / (2.0d0 * delta)
        
        ! Second derivatives (simplified - set to zero for now)
        d2psidr2 = 0.0d0
        d2psidrdz = 0.0d0
        d2psidz2 = 0.0d0
        
    end subroutine thin_orbit_evaluate_field
    
    subroutine thick_orbit_evaluate_field(this, R, Z, psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2)
        ! Thick orbit field evaluation using NEO-RT field for now
        use do_magfie_mod, only: do_magfie, psi_pr
        class(thick_orbit_field_evaluator_t), intent(in) :: this
        real(dp), intent(in) :: R, Z
        real(dp), intent(out) :: psif, dpsidr, dpsidz
        real(dp), intent(out) :: d2psidr2, d2psidrdz, d2psidz2
        
        ! Use same NEO-RT field evaluation as thin orbits for now
        ! TODO: Replace with real POTATO field evaluation once circular dependency is resolved
        real(dp) :: x(3), bmod, sqrtg, bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp), parameter :: delta = 1.0d-6  ! Small step for finite differences
        real(dp) :: psi_plus, psi_minus
        
        ! Set up coordinate input for magfie (Boozer coordinates)
        x(1) = R
        x(2) = 0.0d0  ! phi
        x(3) = Z
        
        ! Call magfie to get field information
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        
        ! Extract psi from magfie result
        psif = psi_pr  ! This is the poloidal flux from magfie
        
        ! Calculate dpsidr
        x(1) = R + delta
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        psi_plus = psi_pr
        
        x(1) = R - delta
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        psi_minus = psi_pr
        
        dpsidr = (psi_plus - psi_minus) / (2.0d0 * delta)
        
        ! Calculate dpsidz  
        x(1) = R
        x(3) = Z + delta
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        psi_plus = psi_pr
        
        x(3) = Z - delta
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        psi_minus = psi_pr
        
        dpsidz = (psi_plus - psi_minus) / (2.0d0 * delta)
        
        ! Second derivatives (simplified - set to zero for now)
        d2psidr2 = 0.0d0
        d2psidrdz = 0.0d0
        d2psidz2 = 0.0d0
        
    end subroutine thick_orbit_evaluate_field
    
    subroutine thin_orbit_cleanup(this)
        class(thin_orbit_field_evaluator_t), intent(inout) :: this
        ! No cleanup needed for thin orbit evaluator
    end subroutine thin_orbit_cleanup
    
    subroutine thick_orbit_cleanup(this)
        class(thick_orbit_field_evaluator_t), intent(inout) :: this
        ! No cleanup needed for thick orbit evaluator
        this%initialized = .false.
    end subroutine thick_orbit_cleanup
    
end module field_interface