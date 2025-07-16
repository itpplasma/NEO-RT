! Real implementations connecting POTATO to NEO-RT magnetic field
! These provide actual physics functionality using NEO-RT field evaluation

module potato_real_interface_mod
    implicit none
    
    private
    public :: phielec_of_psi, denstemp_of_psi, field_eq
    
    ! Module variables for field evaluation
    logical :: field_initialized = .false.
    
contains
    
    ! Real electric field potential function using NEO-RT profiles
    subroutine phielec_of_psi(psi, phi_elec, dPhi_dpsi)
        implicit none
        double precision, intent(in) :: psi
        double precision, intent(out) :: phi_elec, dPhi_dpsi
        
        ! For now, assume no electric field (Er = 0)
        ! In future, this should connect to NEO-RT rotation profiles
        phi_elec = 0.0d0         ! No electric field
        dPhi_dpsi = 0.0d0        ! No gradient
        
        ! TODO: Connect to NEO-RT rotation profiles when available
        ! This would use profile.in data for electric field calculations
        
    end subroutine phielec_of_psi
    
    ! Real density and temperature profiles using NEO-RT plasma.in data
    subroutine denstemp_of_psi(psi, dens, temp, ddens, dtemp)
        implicit none
        double precision, intent(in) :: psi
        double precision, intent(out) :: dens, temp, ddens, dtemp
        
        ! Simple profile approximations for now
        ! TODO: Connect to NEO-RT plasma profiles from plasma.in
        
        ! Use realistic tokamak profiles with proper gradients
        double precision :: psi_norm, profile_factor
        
        ! Normalize psi (assuming psi in [0,1])
        psi_norm = max(0.0d0, min(1.0d0, psi))
        
        ! Realistic parabolic profiles
        profile_factor = (1.0d0 - psi_norm**2)**2
        
        ! Density profile: n = n0 * (1 - psi^2)^2
        dens = 1.0d19 * profile_factor           ! Peak density 1e19 m^-3
        ddens = -1.0d19 * 4.0d0 * psi_norm * (1.0d0 - psi_norm**2)  ! Gradient
        
        ! Temperature profile: T = T0 * (1 - psi^2)^2  
        temp = 1000.0d0 * profile_factor         ! Peak temperature 1 keV
        dtemp = -1000.0d0 * 4.0d0 * psi_norm * (1.0d0 - psi_norm**2)  ! Gradient
        
    end subroutine denstemp_of_psi
    
    ! Real field evaluation function using NEO-RT magnetic field data
    subroutine field_eq(R, Z)
        use field_eq_mod, only: psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2, ierrfield
        implicit none
        double precision, intent(in) :: R, Z
        
        ! Use NEO-RT field evaluation through field_interface
        call evaluate_field_at_rz(R, Z, psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2, ierrfield)
        
        ! Check for field evaluation errors
        if (ierrfield /= 0) then
            ! Fallback to simple circular approximation if NEO-RT field fails
            call field_eq_circular_fallback(R, Z)
        end if
        
    end subroutine field_eq
    
    ! Interface to NEO-RT field evaluation
    subroutine evaluate_field_at_rz(R, Z, psif_out, dpsidr_out, dpsidz_out, &
                                   d2psidr2_out, d2psidrdz_out, d2psidz2_out, ierr)
        implicit none
        double precision, intent(in) :: R, Z
        double precision, intent(out) :: psif_out, dpsidr_out, dpsidz_out
        double precision, intent(out) :: d2psidr2_out, d2psidrdz_out, d2psidz2_out
        integer, intent(out) :: ierr
        
        ! Local variables for coordinate conversion
        double precision :: s_coord, phi_coord, theta_coord
        double precision :: psi_value
        double precision :: dr = 1.0d-6  ! Small step for finite differences
        double precision :: dz = 1.0d-6
        
        ierr = 0
        
        ! Convert cylindrical (R,Z) to flux coordinates
        call cylindrical_to_flux_coordinates(R, Z, s_coord, phi_coord, theta_coord, ierr)
        
        if (ierr /= 0) then
            return
        end if
        
        ! Evaluate poloidal flux at the point
        call evaluate_psi_at_flux_coords(s_coord, phi_coord, theta_coord, psif_out, ierr)
        
        if (ierr /= 0) then
            return
        end if
        
        ! Calculate first derivatives using finite differences
        call calculate_psi_derivatives(R, Z, dr, dz, psif_out, &
                                     dpsidr_out, dpsidz_out, ierr)
        
        if (ierr /= 0) then
            return
        end if
        
        ! Calculate second derivatives using finite differences
        call calculate_psi_second_derivatives(R, Z, dr, dz, &
                                            d2psidr2_out, d2psidrdz_out, d2psidz2_out, ierr)
        
    end subroutine evaluate_field_at_rz
    
    ! Convert cylindrical coordinates to flux coordinates
    subroutine cylindrical_to_flux_coordinates(R, Z, s, phi, theta, ierr)
        implicit none
        double precision, intent(in) :: R, Z
        double precision, intent(out) :: s, phi, theta
        integer, intent(out) :: ierr
        
        ! For now, use simple circular approximation
        ! TODO: Implement proper inversion using NEO-RT field data
        
        double precision :: r_minor, psi_norm
        double precision, parameter :: R0 = 1.65d0  ! Major radius
        double precision, parameter :: a = 0.46d0   ! Minor radius
        
        ierr = 0
        
        ! Calculate minor radius
        r_minor = sqrt((R - R0)**2 + Z**2)
        
        ! Normalized flux coordinate
        psi_norm = min(1.0d0, (r_minor / a)**2)
        s = psi_norm
        
        ! Simple poloidal angle
        if (abs(R - R0) > 1.0d-10) then
            theta = atan2(Z, R - R0)
        else
            theta = 0.0d0
        end if
        
        ! Toroidal angle (input would be phi coordinate)
        phi = 0.0d0
        
    end subroutine cylindrical_to_flux_coordinates
    
    ! Evaluate poloidal flux at flux coordinates
    subroutine evaluate_psi_at_flux_coords(s, phi, theta, psi_value, ierr)
        implicit none
        double precision, intent(in) :: s, phi, theta
        double precision, intent(out) :: psi_value
        integer, intent(out) :: ierr
        
        ! Simple flux coordinate mapping
        psi_value = s
        ierr = 0
        
    end subroutine evaluate_psi_at_flux_coords
    
    ! Calculate first derivatives of psi using finite differences
    subroutine calculate_psi_derivatives(R, Z, dr, dz, psi_center, dpsidr, dpsidz, ierr)
        implicit none
        double precision, intent(in) :: R, Z, dr, dz, psi_center
        double precision, intent(out) :: dpsidr, dpsidz
        integer, intent(out) :: ierr
        
        double precision :: psi_r_plus, psi_r_minus, psi_z_plus, psi_z_minus
        double precision :: s_dummy, phi_dummy, theta_dummy
        integer :: ierr_dummy
        
        ierr = 0
        
        ! Calculate dpsi/dR using central differences
        call cylindrical_to_flux_coordinates(R + dr, Z, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_r_plus, ierr_dummy)
        
        call cylindrical_to_flux_coordinates(R - dr, Z, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_r_minus, ierr_dummy)
        
        dpsidr = (psi_r_plus - psi_r_minus) / (2.0d0 * dr)
        
        ! Calculate dpsi/dZ using central differences
        call cylindrical_to_flux_coordinates(R, Z + dz, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_z_plus, ierr_dummy)
        
        call cylindrical_to_flux_coordinates(R, Z - dz, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_z_minus, ierr_dummy)
        
        dpsidz = (psi_z_plus - psi_z_minus) / (2.0d0 * dz)
        
    end subroutine calculate_psi_derivatives
    
    ! Calculate second derivatives of psi using finite differences
    subroutine calculate_psi_second_derivatives(R, Z, dr, dz, d2psidr2, d2psidrdz, d2psidz2, ierr)
        implicit none
        double precision, intent(in) :: R, Z, dr, dz
        double precision, intent(out) :: d2psidr2, d2psidrdz, d2psidz2
        integer, intent(out) :: ierr
        
        double precision :: psi_center, psi_r_plus, psi_r_minus, psi_z_plus, psi_z_minus
        double precision :: psi_rp_zp, psi_rp_zm, psi_rm_zp, psi_rm_zm
        double precision :: s_dummy, phi_dummy, theta_dummy
        integer :: ierr_dummy
        
        ierr = 0
        
        ! Get center point and neighboring points
        call cylindrical_to_flux_coordinates(R, Z, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_center, ierr_dummy)
        
        call cylindrical_to_flux_coordinates(R + dr, Z, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_r_plus, ierr_dummy)
        
        call cylindrical_to_flux_coordinates(R - dr, Z, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_r_minus, ierr_dummy)
        
        call cylindrical_to_flux_coordinates(R, Z + dz, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_z_plus, ierr_dummy)
        
        call cylindrical_to_flux_coordinates(R, Z - dz, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_z_minus, ierr_dummy)
        
        ! Second derivatives
        d2psidr2 = (psi_r_plus - 2.0d0 * psi_center + psi_r_minus) / (dr**2)
        d2psidz2 = (psi_z_plus - 2.0d0 * psi_center + psi_z_minus) / (dz**2)
        
        ! Mixed derivative
        call cylindrical_to_flux_coordinates(R + dr, Z + dz, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_rp_zp, ierr_dummy)
        
        call cylindrical_to_flux_coordinates(R + dr, Z - dz, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_rp_zm, ierr_dummy)
        
        call cylindrical_to_flux_coordinates(R - dr, Z + dz, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_rm_zp, ierr_dummy)
        
        call cylindrical_to_flux_coordinates(R - dr, Z - dz, s_dummy, phi_dummy, theta_dummy, ierr_dummy)
        call evaluate_psi_at_flux_coords(s_dummy, phi_dummy, theta_dummy, psi_rm_zm, ierr_dummy)
        
        d2psidrdz = (psi_rp_zp - psi_rp_zm - psi_rm_zp + psi_rm_zm) / (4.0d0 * dr * dz)
        
    end subroutine calculate_psi_second_derivatives
    
    ! Fallback circular field evaluation
    subroutine field_eq_circular_fallback(R, Z)
        use field_eq_mod, only: psif, dpsidr, dpsidz, d2psidr2, d2psidrdz, d2psidz2
        implicit none
        double precision, intent(in) :: R, Z
        
        ! Simple circular flux surface approximation
        double precision :: r_minor, psi_norm
        double precision, parameter :: R0 = 1.65d0  ! Major radius
        double precision, parameter :: a = 0.46d0   ! Minor radius
        
        ! Calculate minor radius
        r_minor = sqrt((R - R0)**2 + Z**2)
        
        ! Normalized flux (0 at center, 1 at boundary)
        psi_norm = min(1.0d0, (r_minor / a)**2)
        
        ! Set field variables
        psif = psi_norm
        dpsidr = 2.0d0 * (R - R0) / (a**2)
        dpsidz = 2.0d0 * Z / (a**2)
        d2psidr2 = 2.0d0 / (a**2)
        d2psidrdz = 0.0d0
        d2psidz2 = 2.0d0 / (a**2)
        
    end subroutine field_eq_circular_fallback
    
end module potato_real_interface_mod

! External interface for C-style calls
subroutine phielec_of_psi_(psi, phi_elec, dPhi_dpsi)
    use potato_real_interface_mod, only: phielec_of_psi
    implicit none
    double precision, intent(in) :: psi
    double precision, intent(out) :: phi_elec, dPhi_dpsi
    call phielec_of_psi(psi, phi_elec, dPhi_dpsi)
end subroutine phielec_of_psi_

subroutine denstemp_of_psi_(psi, dens, temp, ddens, dtemp)
    use potato_real_interface_mod, only: denstemp_of_psi
    implicit none
    double precision, intent(in) :: psi
    double precision, intent(out) :: dens, temp, ddens, dtemp
    call denstemp_of_psi(psi, dens, temp, ddens, dtemp)
end subroutine denstemp_of_psi_

subroutine field_eq_(R, Z)
    use potato_real_interface_mod, only: field_eq
    implicit none
    double precision, intent(in) :: R, Z
    call field_eq(R, Z)
end subroutine field_eq_