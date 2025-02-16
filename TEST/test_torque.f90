program test_torque_prog
    use neort, only: runname, read_control, init_profile, init_plasma, init_test, &
                     test_magfie, compute_torque, compute_transport_coeffs
    implicit none

    call test_torque

    contains

    subroutine test_torque
        runname = "driftorbit64.new"
        call read_control
        call init_profile
        call init_plasma
        call init_test(use_thermodynamic_profiles=.True.)
        call test_magfie
        call compute_torque
        call compute_transport_coeffs
    end subroutine test_torque
end program test_torque_prog
