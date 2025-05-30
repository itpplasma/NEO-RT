module util_for_test
    implicit none

    contains

    subroutine print_test(test_name)
        character(*) :: test_name
        print *, "==> ", test_name
    end subroutine print_test

    subroutine pass_test
        call print_ok
    end subroutine pass_test

    subroutine fail_test
        call print_fail
        error stop
    end subroutine fail_test

    subroutine print_ok
        print *, "    .................................................... OK"
    end subroutine print_ok


    subroutine print_fail
        print *, "    .................................................... FAIL"
    end subroutine print_fail

end module util_for_test