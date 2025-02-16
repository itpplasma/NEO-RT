program test_torque_prog
    use neort, only: runname, read_control, init_profile, init_plasma, init_test, &
                     test_magfie, compute_torque_harmonic, compute_transport_coeff_harmonic
    use do_magfie_mod, only: do_magfie_init
    implicit none

    real(8) :: Dco(2), Dctr(2), Dt(2)
    real(8) :: Tco, Tctr, Tt
    integer, parameter :: MTH = -3

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
        call compute_transport_coeff_harmonic(MTH, Dco, Dctr, Dt)
        call compute_torque_harmonic(MTH, Tco, Tctr, Tt)
    end subroutine test_torque

    subroutine compare_output
    end subroutine compare_output
end program test_torque_prog
