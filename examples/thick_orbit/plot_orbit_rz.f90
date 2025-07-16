program plot_orbit_rz
    ! Visualize and compare thin vs thick orbit trajectories in R-Z plane
    ! This program demonstrates the orbit width differences between the two approaches
    ! Currently uses synthetic data to illustrate the expected physics
    
    use fortplot
    
    implicit none
    
    ! Parameters
    integer, parameter :: npts = 200       ! number of points along orbit
    real(8), parameter :: R0 = 1.65d0      ! major radius [m]
    real(8), parameter :: a = 0.5d0        ! minor radius [m]
    real(8), parameter :: s = 0.6d0        ! normalized flux surface
    real(8), parameter :: eta = 0.7d0      ! pitch parameter (trapped)
    real(8), parameter :: q = 3.5d0        ! safety factor
    real(8), parameter :: eps = a * sqrt(s) / R0  ! local inverse aspect ratio
    
    ! Orbit arrays
    real(8), dimension(npts) :: r_thin, z_thin, r_thick, z_thick
    real(8), dimension(npts) :: theta
    
    ! Local variables
    integer :: i
    real(8) :: dtheta, r_flux, banana_width, radial_shift
    real(8) :: kappa  ! elongation parameter for trapped orbit
    
    print *, '========================================'
    print *, 'Orbit R-Z Comparison: Thin vs Thick'
    print *, '========================================'
    print *, 'Parameters:'
    print '(A,F6.3)', '  s (flux surface) = ', s
    print '(A,F6.3)', '  eta (pitch)      = ', eta
    print '(A,F6.3)', '  q (safety factor)= ', q
    print '(A,F6.3)', '  epsilon          = ', eps
    print *, ''
    
    ! Calculate orbit parameters
    r_flux = R0 + a * sqrt(s)  ! flux surface radius
    kappa = sqrt(1.0d0 - eta)   ! measure of trappedness
    
    ! Estimate banana width for thick orbit
    ! Banana width ~ q * rho_pol where rho_pol ~ rho_tor/sqrt(eps)
    ! For demonstration, use ~ 2-3% of minor radius
    banana_width = 0.025d0 * a * sqrt(s) / sqrt(eps)
    
    print *, 'Generating orbit trajectories...'
    print '(A,F6.3,A)', '  Flux surface radius: ', r_flux, ' m'
    print '(A,F6.3,A)', '  Estimated banana width: ', banana_width, ' m'
    
    ! Generate poloidal angle array (full bounce)
    dtheta = 4.0d0 * asin(kappa) / real(npts-1, 8)  ! trapped orbit range
    
    do i = 1, npts
        ! Poloidal angle for trapped orbit (symmetric around theta=0)
        theta(i) = -2.0d0 * asin(kappa) + real(i-1, 8) * dtheta
        
        ! Thin orbit - standard banana on flux surface
        r_thin(i) = r_flux + a * sqrt(s) * cos(theta(i))
        z_thin(i) = a * sqrt(s) * sin(theta(i))
        
        ! Thick orbit - includes finite orbit width effects
        ! 1. Radial excursion varies along orbit (maximum at banana tips)
        radial_shift = banana_width * sin(theta(i)/2.0d0)**2
        
        ! 2. Orbit center is shifted outward (grad-B drift)
        r_thick(i) = r_flux + radial_shift + a * sqrt(s) * cos(theta(i) * 0.98d0)
        z_thick(i) = a * sqrt(s) * sin(theta(i) * 0.98d0)
    end do
    
    ! Create comparison plot
    print *, ''
    print *, 'Creating orbit comparison plot...'
    
    call figure(1000, 800)
    
    ! Plot thin orbit
    call plot(r_thin, z_thin, label="Thin orbit", linestyle="b-")
    
    ! Plot thick orbit  
    call plot(r_thick, z_thick, label="Thick orbit", linestyle="r--")
    
    call xlabel("R [m]")
    call ylabel("Z [m]")
    call title("Orbit Trajectory Comparison: Thin vs Thick (s=0.6, η=0.7)")
    call legend()
    
    call savefig("orbit_rz_comparison.png")
    
    print *, 'Plot saved as orbit_rz_comparison.png'
    print *, ''
    print *, 'Orbit characteristics (synthetic data):'
    print '(A,F6.3,A)', '  Banana width:        ', &
            maxval(sqrt((r_thick - R0)**2 + z_thick**2)) - &
            minval(sqrt((r_thick - R0)**2 + z_thick**2)), ' m'
    print '(A,F6.3,A)', '  Tip radial shift:   ', &
            r_thick(1) - r_thin(1), ' m'
    print '(A,F6.3,A)', '  Orbit center shift: ', &
            sum(r_thick)/npts - sum(r_thin)/npts, ' m'
    print *, ''
    print *, 'NOTE: This demonstration uses synthetic data to illustrate'
    print *, '      the expected differences between thin and thick orbits.'
    print *, '      Real POTATO integration will provide accurate physics.'
    
end program plot_orbit_rz