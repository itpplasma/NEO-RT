program test_torque_prog
    use iso_fortran_env, only: dp => real64
    use neort, only: check_magfie, init_profiles, &
                     init, compute_transport_harmonic, write_magfie_data_to_files
    use neort_config, only: read_and_set_config
    use neort_main, only: runname
    use neort_datatypes
    use neort_profiles, only: read_and_init_plasma_input, read_and_init_profile_input, M_t
    use driftorbit, only: A1, A2, ni1, vth, B0, a, efac
    use do_magfie_mod, only: do_magfie_init, do_magfie, R0, iota, bfac, s
    use util
    implicit none

    real(dp) :: Dco(2), Dctr(2), Dt(2)
    real(dp) :: Tco, Tctr, Tt

    integer, parameter :: MTH = -3
    real(dp), parameter :: AVG_NABLA_S = 0.02162413513580247_dp
    real(dp), parameter :: DV_OVER_DS = 12621249.179676009_dp
    real(dp), parameter :: FLUX_SURFACE_AREA = DV_OVER_DS * AVG_NABLA_S
    real(dp), parameter :: TOL = 1.0e-7_dp
    character(1024) :: TEST_RUN = "driftorbit64.new"
    type(magfie_data_t) :: magfie_data
    type(transport_harmonic_t) :: harmonic_data

    call setup
    call test_torque

contains

    subroutine setup
        runname = TEST_RUN
        call read_and_set_config(trim(runname)//".in")
        call do_magfie_init("in_file")
        call init_profiles(R0)
        call read_and_init_plasma_input("plasma.in", s)
        call read_and_init_profile_input("profile.in", s, R0, efac, bfac)
        call init
        call check_magfie(magfie_data)
        call write_magfie_data_to_files(magfie_data, runname)

        Dco = 0.0_dp
        Dctr = 0.0_dp
        Dt = 0.0_dp
        Tco = 0.0_dp
        Tctr = 0.0_dp
        Tt = 0.0_dp
        call compute_transport_harmonic(MTH, Dco, Dctr, Dt, Tco, Tctr, Tt, harmonic_data)
        call correct_transport_coeff(Dco)
        call correct_transport_coeff(Dctr)
        call correct_transport_coeff(Dt)
        call correct_torque(Tco)
        call correct_torque(Tctr)
        call correct_torque(Tt)
    end subroutine setup

    subroutine test_torque
        real(dp) :: TTco, TTctr, TTt

        TTco = torque_from_transport(Dco)
        TTctr = torque_from_transport(Dctr)
        TTt = torque_from_transport(Dt)

        print *, ""
        print *, "torque_from_transport: Mt = ", M_t, ", mth = ", MTH
        write (*, "(4ES12.2,2F12.2)") TTco, TTctr, TTt, TTco + TTctr + TTt
        print *, ""
        print *, "Ratios to direct calculation:"
        print *, "TTctr/Tctr = ", TTctr / Tctr
        print *, "TTt/Tt = ", TTt / Tt

        if (abs(TTctr / Tctr - 1.0_dp) > TOL) then
            error stop
        end if

        if (abs(TTt / Tt - 1.0_dp) > TOL) then
            error stop
        end if

        print *, "test_torque ... OK"

    end subroutine test_torque

    function torque_from_transport(D) result(Tphi)
        real(dp), intent(in) :: D(2)
        real(dp) :: Tphi
        real(dp) :: D_plateau
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: sqrtgBth

        x(1) = s
        x(2) = 0.0_dp
        x(3) = 0.0_dp
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        sqrtgBth = sqrtg * hctrvr(3) * bmod * AVG_NABLA_S
        print *, "sqrtgBth = ", sqrtgBth

        D_plateau = pi * vth**3 / (16.0_dp * R0 * iota * (qi * B0 / (mi * c))**2)
        Tphi = (sqrtgBth / c) * qe * (-ni1 * D_plateau * (D(1) * A1 + D(2) * A2)) * AVG_NABLA_S
        Tphi = Tphi
    end function torque_from_transport

    ! Correct ad-hoc approximation for AVG_NABLA_S
    ! This then finally cancels out in torque calculation
    subroutine correct_transport_coeff(D)
        real(dp), intent(inout) :: D(2)
        D = D * (2.0_dp / a * sqrt(s))**2 / AVG_NABLA_S**2
    end subroutine correct_transport_coeff

    ! Correct torque calculation
    subroutine correct_torque(T)
        real(dp), intent(inout) :: T
        T = T / DV_OVER_DS
    end subroutine correct_torque
end program test_torque_prog
