program test_torque_prog
    use neort, only: runname, read_control, init_profile, init_plasma, init_test, &
                     test_magfie, compute_torque, compute_transport_coeffs
    use do_magfie_mod, only: do_magfie_init
    implicit none

    call test_torque

    contains

    subroutine test_torque
        runname = "driftorbit64.new"
        call read_control
        call do_magfie_init
        call init_profile
        call init_plasma
        call init_test(use_thermodynamic_profiles=.True.)
        call test_magfie
        call compute_transport_coeffs
        call compute_torque
    end subroutine test_torque

    subroutine compare_output
    end subroutine compare_output
end program test_torque_prog
