program test_init
    ! Minimal test to verify NEO-RT physics initialization
    ! This tests the basic initialization without complex dependencies
    
    use do_magfie_mod, only: do_magfie_init
    use driftorbit, only: s
    implicit none
    
    ! Test parameters
    logical :: success = .true.
    
    print *, '========================================'
    print *, 'Test: Basic NEO-RT Physics Initialization'
    print *, '========================================'
    print *, ''
    
    ! Set flux surface
    s = 0.6d0
    
    ! Test magnetic field initialization
    print *, 'Testing magnetic field initialization...'
    
    ! This should read the in_file in the current directory
    call do_magfie_init
    
    print *, 'Magnetic field initialization completed successfully'
    print *, 'Flux surface s = ', s
    
    print *, ''
    print *, 'Test completed successfully!'
    
end program test_init