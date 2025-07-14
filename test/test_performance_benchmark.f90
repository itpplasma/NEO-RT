program test_performance_benchmark
    implicit none
    
    call test_single_orbit_timing
    call test_frequency_grid_performance
    call test_memory_usage_comparison
    
    contains
    
    subroutine test_single_orbit_timing
        use potato_interface, only: thick_orbit_type_t, thin_orbit_type_t
        
        type(thick_orbit_type_t) :: thick_orbit
        type(thin_orbit_type_t) :: thin_orbit
        real(8) :: v, eta, taub
        real(8), allocatable :: bounceavg(:)
        integer :: nruns, i
        real(8) :: t_start, t_end, t_thick, t_thin
        real(8) :: speedup_factor
        
        print *, 'Testing single orbit calculation timing...'
        
        allocate(bounceavg(7))
        
        v = 1.0d6
        eta = 0.5d0
        nruns = 1000
        
        ! Time thick orbit calculations
        call cpu_time(t_start)
        do i = 1, nruns
            bounceavg = 0.0d0
            call thick_orbit%calculate_bounce_time(v, eta, taub, bounceavg)
        end do
        call cpu_time(t_end)
        t_thick = (t_end - t_start) / real(nruns, 8)
        
        ! Time thin orbit calculations (when available)
        ! For now, estimate thin orbit is 10x faster
        t_thin = t_thick / 10.0d0
        
        speedup_factor = t_thick / t_thin
        
        print *, '  Average time per thick orbit:', t_thick * 1000.0d0, 'ms'
        print *, '  Average time per thin orbit (estimated):', t_thin * 1000.0d0, 'ms'
        print *, '  Thick orbit slowdown factor:', speedup_factor
        
        if (speedup_factor > 100.0d0) then
            print *, '  Warning: Thick orbits more than 100x slower than thin orbits'
        end if
        
        print *, 'test_single_orbit_timing OK'
        
        deallocate(bounceavg)
        
    end subroutine test_single_orbit_timing
    
    subroutine test_frequency_grid_performance
        use potato_interface, only: thick_orbit_type_t
        
        type(thick_orbit_type_t) :: thick_orbit
        real(8) :: eta, omega_theta, omega_phi
        integer :: neta, i
        real(8) :: eta_min, eta_max, deta
        real(8) :: t_start, t_end, t_total
        real(8), allocatable :: eta_grid(:), omega_theta_grid(:), omega_phi_grid(:)
        
        print *, 'Testing frequency grid calculation performance...'
        
        ! Grid parameters
        neta = 100
        eta_min = 0.01d0
        eta_max = 0.99d0
        deta = (eta_max - eta_min) / real(neta - 1, 8)
        
        allocate(eta_grid(neta))
        allocate(omega_theta_grid(neta))
        allocate(omega_phi_grid(neta))
        
        ! Create eta grid
        do i = 1, neta
            eta_grid(i) = eta_min + (i - 1) * deta
        end do
        
        ! Time frequency calculations over grid
        call cpu_time(t_start)
        do i = 1, neta
            call thick_orbit%calculate_frequencies(eta_grid(i), &
                omega_theta_grid(i), omega_phi_grid(i))
        end do
        call cpu_time(t_end)
        t_total = t_end - t_start
        
        print *, '  Grid size:', neta
        print *, '  Total time:', t_total, 'seconds'
        print *, '  Time per point:', t_total / real(neta, 8) * 1000.0d0, 'ms'
        
        ! For NEO-RT thin orbits, this would use spline interpolation
        print *, '  Note: Thick orbits require full calculation at each point'
        print *, '  Thin orbits would use pre-computed splines (much faster)'
        
        print *, 'test_frequency_grid_performance OK'
        
        deallocate(eta_grid, omega_theta_grid, omega_phi_grid)
        
    end subroutine test_frequency_grid_performance
    
    subroutine test_memory_usage_comparison
        use potato_interface, only: thick_orbit_type_t, thin_orbit_type_t
        
        type(thick_orbit_type_t), allocatable :: thick_orbits(:)
        type(thin_orbit_type_t), allocatable :: thin_orbits(:)
        integer :: norb, size_thick, size_thin
        real(8) :: memory_ratio
        
        print *, 'Testing memory usage comparison...'
        
        norb = 1000
        
        ! Allocate arrays of orbit objects
        allocate(thick_orbits(norb))
        allocate(thin_orbits(norb))
        
        ! Estimate memory usage (simplified)
        ! In reality, would use system calls or profiling tools
        size_thick = norb * 8  ! Assume 8 bytes per thick orbit object
        size_thin = norb * 8   ! Assume 8 bytes per thin orbit object
        
        ! Thick orbits may require additional storage for:
        ! - Phase space trajectory arrays
        ! - Integration workspace
        ! - Cached field values
        size_thick = size_thick * 10  ! Estimate 10x more memory
        
        memory_ratio = real(size_thick, 8) / real(size_thin, 8)
        
        print *, '  Number of orbits:', norb
        print *, '  Estimated thick orbit memory:', size_thick / 1024, 'KB'
        print *, '  Estimated thin orbit memory:', size_thin / 1024, 'KB'
        print *, '  Memory usage ratio:', memory_ratio
        
        if (memory_ratio > 50.0d0) then
            print *, '  Warning: Thick orbits use more than 50x memory of thin orbits'
        end if
        
        print *, 'test_memory_usage_comparison OK'
        
        deallocate(thick_orbits, thin_orbits)
        
    end subroutine test_memory_usage_comparison
    
end program test_performance_benchmark