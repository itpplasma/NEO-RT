program test_torque_prog
    use neort, only: runname, read_control, init_profile, init_plasma, init_test, M_t, &
                     test_magfie, compute_torque_harmonic, &
                     compute_transport_coeff_harmonic
    use driftorbit, only: A1, A2, ni1, vth, B0
    use do_magfie_mod, only: do_magfie_init, R0, iota, psi_pr
    use util
    implicit none

    real(8) :: Dco(2), Dctr(2), Dt(2)
    real(8) :: Tco, Tctr, Tt
    integer, parameter :: MTH = -3

    call test_torque
    call print_torque_from_transport

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

    subroutine print_torque_from_transport
        real(8) :: TTco, TTctr, TTt

        TTco = torque_from_transport(Dco)
        TTctr = torque_from_transport(Dctr)
        TTt = torque_from_transport(Dt)

        print *, ""
        print *, "torque_from_transport: Mt = ", M_t, ", mth = ", MTH
        write (*, "(4ES12.2,2F12.2)") TTco, TTctr, TTt, TTco + TTctr + TTt
    end subroutine print_torque_from_transport

    function torque_from_transport(D) result(Tphi)
        real(8), intent(in) :: D(2)
        real(8) :: Tphi
        real(8) :: Dp

        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
        Tphi = -(iota*psi_pr/c)*qe*(-ni1*Dp*(D(1)*A1 + D(2)*A2))
    end function torque_from_transport
end program test_torque_prog
