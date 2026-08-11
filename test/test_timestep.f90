program test_timestep
    use iso_fortran_env, only: dp => real64
    use util_for_test, only: print_test, pass_test, fail_test
    use neort_config, only: read_and_set_config
    use do_magfie_mod, only: read_boozer_file, init_magfie_at_s, do_magfie, s
    use neort_orbit, only: timestep_poloidal_motion, th0, vpar
    use driftorbit, only: vth

    implicit none

    call setup
    call test_timestep_poloidal()

contains

    ! Reads the same configuration and magnetic-field input as the
    ! ripple-plateau integration test (test/ripple_plateau): the
    ! driftorbit_test.in namelist and the in_file Boozer file.
    subroutine setup
        character(len=*), parameter :: config_file = "driftorbit_test.in"
        character(len=*), parameter :: boozer_file = "in_file"

        call read_and_set_config(config_file)
        call read_boozer_file(boozer_file)
        call init_magfie_at_s()
    end subroutine setup

    ! Verifies that timestep_poloidal_motion computes the poloidal
    ! equations of motion consistently with the local magnetic-field
    ! quantities returned by do_magfie.
    subroutine test_timestep_poloidal
        integer, parameter :: neq = 2
        real(dp), parameter :: tolerance = 1.0e-12_dp
        real(dp), parameter :: eta_value = 1.0e-6_dp

        real(dp) :: v, t, y(neq), ydot(neq), ydot_expected(neq)
        real(dp) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: hthctr, hderth, v_par

        call print_test("test_timestep_poloidal")

        v = vth
        t = 0.0_dp

        ! Replicate the field evaluation done inside timestep_poloidal_motion.
        x = [s, 0.0_dp, th0]
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        hthctr = hctrvr(3)
        hderth = bder(3)
        v_par = vpar(v, eta_value, bmod)

        y(1) = th0
        y(2) = v_par
        call timestep_poloidal_motion(v, eta_value, neq, t, y, ydot)

        ydot_expected(1) = v_par * hthctr
        ydot_expected(2) = -v**2 * eta_value / 2.0_dp * hthctr * hderth * bmod

        if (maxval(abs(ydot - ydot_expected)) > tolerance) then
            print *, "ydot =", ydot
            print *, "expected =", ydot_expected
            call fail_test()
        end if
        call pass_test()
    end subroutine test_timestep_poloidal

end program test_timestep
