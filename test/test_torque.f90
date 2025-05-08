program test_torque_prog
    use neort, only: read_control, check_magfie, init_profiles, init_profile_input, &
                      init_plasma_input, init, compute_transport_harmonic, &
                      runname, s, M_t
    use driftorbit, only: A1, A2, ni1, vth, B0, a, efac
    use do_magfie_mod, only: do_magfie_init, do_magfie, R0, iota, bfac
    use util
    implicit none

    real(8) :: Dco(2), Dctr(2), Dt(2)
    real(8) :: Tco, Tctr, Tt

    integer, parameter :: MTH = -3
    real(8), parameter :: AVG_NABLA_S = 0.02162413513580247d0
    real(8), parameter :: DV_OVER_DS = 12621249.179676009d0
    real(8), parameter :: FLUX_SURFACE_AREA = DV_OVER_DS*AVG_NABLA_S
    real(8), parameter :: TOL = 1d-7
    character(1024) :: TEST_RUN = "driftorbit64.new"

    call setup
    call test_torque

contains

    subroutine setup
        runname = TEST_RUN
        call read_control
        call do_magfie_init
        call init_profiles(R0)
        call init_plasma_input(s)
        call init_profile_input(s, R0, efac, bfac)
        call init
        call check_magfie

        Dco = 0d0
        Dctr = 0d0
        Dt = 0d0
        Tco = 0d0
        Tctr = 0d0
        Tt = 0d0
        call compute_transport_harmonic(MTH, Dco, Dctr, Dt, Tco, Tctr, Tt)
        call correct_transport_coeff(Dco)
        call correct_transport_coeff(Dctr)
        call correct_transport_coeff(Dt)
        call correct_torque(Tco)
        call correct_torque(Tctr)
        call correct_torque(Tt)
    end subroutine setup

    subroutine test_torque
        real(8) :: TTco, TTctr, TTt

        TTco = torque_from_transport(Dco)
        TTctr = torque_from_transport(Dctr)
        TTt = torque_from_transport(Dt)

        print *, ""
        print *, "torque_from_transport: Mt = ", M_t, ", mth = ", MTH
        write (*, "(4ES12.2,2F12.2)") TTco, TTctr, TTt, TTco + TTctr + TTt
        print *, ""
        print *, "Ratios to direct calculation:"
        print *, "TTctr/Tctr = ", TTctr/Tctr
        print *, "TTt/Tt = ", TTt/Tt

        if (abs(TTctr/Tctr - 1d0) > TOL) then
            error stop
        end if

        if (abs(TTt/Tt - 1d0) > TOL) then
            error stop
        end if

        print *, "test_torque ... OK"

    end subroutine test_torque

    function torque_from_transport(D) result(Tphi)
        real(8), intent(in) :: D(2)
        real(8) :: Tphi
        real(8) :: Dp
        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(8) :: sqrtgBth

        x(1) = s
        x(2) = 0d0
        x(3) = 0d0
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        sqrtgBth = sqrtg*hctrvr(3)*bmod*AVG_NABLA_S
        print *, "sqrtgBth = ", sqrtgBth

        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
        Tphi = (sqrtgBth/c)*qe*(-ni1*Dp*(D(1)*A1 + D(2)*A2))*AVG_NABLA_S
        Tphi = Tphi
    end function torque_from_transport

    ! Correct ad-hoc approximation for AVG_NABLA_S
    ! This then finally cancels out in torque calculation
    subroutine correct_transport_coeff(D)
        real(8), intent(inout) :: D(2)
        D = D*(2d0/a*sqrt(s))**2/AVG_NABLA_S**2
    end subroutine correct_transport_coeff


    ! Correct torque calculation
    subroutine correct_torque(T)
        real(8), intent(inout) :: T
        T = T/DV_OVER_DS
    end subroutine correct_torque
end program test_torque_prog
