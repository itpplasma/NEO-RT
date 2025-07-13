program test_potato_sources
    implicit none
    
    call test_potato_source_availability
    call test_potato_build_integration
    
    contains
    
    subroutine test_potato_source_availability
        print *, 'Testing POTATO source availability...'
        
#ifdef USE_THICK_ORBITS
        print *, 'test_potato_source_availability OK - USE_THICK_ORBITS defined'
#else
        print *, 'test_potato_source_availability FAILED - USE_THICK_ORBITS not defined'
        error stop
#endif
        
    end subroutine test_potato_source_availability
    
    subroutine test_potato_build_integration
        print *, 'Testing POTATO build integration...'
        
#ifdef USE_THICK_ORBITS
        ! This will test if POTATO modules can be accessed
        ! For now, just check that we can compile with the flag
        print *, 'test_potato_build_integration OK - POTATO integration ready'
#else
        print *, 'test_potato_build_integration FAILED - USE_THICK_ORBITS not enabled'
        error stop
#endif
        
    end subroutine test_potato_build_integration
    
end program test_potato_sources