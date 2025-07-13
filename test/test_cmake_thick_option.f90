program test_cmake_thick_option
    implicit none
    
    call test_thick_orbit_compilation_flag
    
    contains
    
    subroutine test_thick_orbit_compilation_flag
        print *, 'Testing CMake thick orbit compilation flag...'
        
#ifdef USE_THICK_ORBITS
        print *, 'test_cmake_thick_option FAILED - USE_THICK_ORBITS should not be defined by default'
        error stop
#else
        print *, 'test_cmake_thick_option OK - USE_THICK_ORBITS correctly undefined by default'
#endif
        
    end subroutine test_thick_orbit_compilation_flag
    
end program test_cmake_thick_option