program test_timestep
    use util_for_test, only: print_test, pass_test, fail_test
    implicit none

    call test_timestep_thin()

contains

    subroutine test_timestep_thin
        call print_test("test_timestep_thin")

        call pass_test()
    end subroutine test_timestep_thin

end program test_timestep
