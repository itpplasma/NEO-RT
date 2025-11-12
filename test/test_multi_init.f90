program test_multi_init
    ! Test that we can call initialization functions multiple times without memory errors
    ! This test verifies that the deallocation-before-reallocation pattern works correctly
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie_init, nmode, nflux
    use neort, only: read_and_set_control

    implicit none

    integer, parameter :: N_ITERATIONS = 5
    integer :: iter
    integer :: nmode_prev, nflux_prev

    print *, "============================================"
    print *, "Testing multiple sequential init calls..."
    print *, "============================================"
    print *, ""

    ! Initialize once to load data
    call read_and_set_control("driftorbit")
    call do_magfie_init("in_file")

    nmode_prev = nmode
    nflux_prev = nflux

    print *, "Initial call successful"
    print *, "  nmode = ", nmode
    print *, "  nflux = ", nflux
    print *, ""

    ! Now call multiple times to test reallocation
    do iter = 1, N_ITERATIONS
        print *, "Iteration ", iter, " of ", N_ITERATIONS

        ! Re-initialize without cleanup - should handle deallocation internally
        call do_magfie_init("in_file")

        ! Verify data is still consistent
        if (nmode /= nmode_prev .or. nflux /= nflux_prev) then
            print *, "ERROR: Data changed between calls!"
            print *, "  Previous: nmode=", nmode_prev, ", nflux=", nflux_prev
            print *, "  Current:  nmode=", nmode, ", nflux=", nflux
            error stop "Test FAILED: Inconsistent data"
        end if

        print *, "  Re-initialization successful (nmode=", nmode, ", nflux=", nflux, ")"
        print *, ""
    end do

    print *, "============================================"
    print *, "All ", N_ITERATIONS, " iterations successful!"
    print *, "TEST PASSED: Multiple init calls work"
    print *, "============================================"

end program test_multi_init
