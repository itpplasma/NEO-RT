program test_cmake_thick_enabled
    implicit none
    
    call test_thick_orbit_enabled_flag
    
    contains
    
    subroutine test_thick_orbit_enabled_flag
        print *, 'Testing CMake thick orbit enabled flag...'
        
#ifdef USE_THICK_ORBITS
        print *, 'test_cmake_thick_enabled OK - USE_THICK_ORBITS correctly defined when enabled'
#else
        print *, 'test_cmake_thick_enabled FAILED - USE_THICK_ORBITS should be defined when enabled'
        error stop
#endif
        
    end subroutine test_thick_orbit_enabled_flag
    
end program test_cmake_thick_enabled