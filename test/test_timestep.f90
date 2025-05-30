program test_timestep
    use util_for_test, only: print_test, pass_test, fail_test
    use util, only: qe, mu
    implicit none

    real(8), parameter :: tol = 1e-13
    real(8) :: v, eta

    call setup
    call test_timestep_poloidal

contains

    subroutine setup
        use neort, only: init
        use driftorbit, only: do_magfie_init, s, inp_swi, etamin, etamax, &
            Om_tE, dOm_tEds, etatp, etadt, epst, sign_vpar, vth, M_t, dM_tds
        use do_magfie_mod, only: R0

        inp_swi = 8
        s = 3.9063d-2
        vth = 1.0d8
        M_t = 0.1d0

        call do_magfie_init
        call init

        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0

        etamin = (1+epst)*etatp
        etamax = (1-epst)*etadt
        sign_vpar = 1

        v = vth
        eta = 0.5d0*(etamin + etamax)
    end subroutine setup

    subroutine test_timestep_poloidal()
        use neort_orbit, only: timestep_poloidal_motion
        use do_magfie_mod, only: do_magfie, s
        use util, only: mu

        integer, parameter :: neq = 2
        real(8) :: y(neq), ydot(neq)
        real(8) :: t, dt
        real(8) :: bmod, sqrtg
        real(8), dimension(3) :: x, hder, hcovar, hctrvr, hcurl
        real(8) :: vpar_val, vperp_squared
        real(8) :: mag_moment, mirror_force, expected_accel
        real(8) :: energy_initial, energy_final
        logical :: test_passed

        call print_test("test_timestep_poloidal")

        x(1) = s
        x(2) = 0d0
        x(3) = 1d-1
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        vpar_val = sqrt(2d0*v**2*eta*bmod)
        vperp_squared = 2d0*v**2*(1d0-eta*bmod)

        y(1) = 0d0
        y(2) = vpar_val
        t = 0d0
        dt = 1d-6

        call timestep_poloidal_motion(v, eta, neq, t, y, ydot)

        test_passed = .true.

        if (.not. (abs(ydot(1)) > 0d0 .and. abs(ydot(1)) < 1d10)) then
            call fail_test("Poloidal velocity out of reasonable range: " // trim(adjustl(str(ydot(1)))))
        end if

        if (.not. (abs(ydot(2)) < 1d20)) then
            call fail_test("Parallel acceleration too large: " // trim(adjustl(str(ydot(2)))))
        end if

        if (sign(1.0d0, ydot(1)) /= sign(1.0d0, y(2)*hctrvr(3))) then
            call fail_test("Poloidal velocity direction inconsistent with parallel velocity: " // &
                          "ydot(1)=" // trim(adjustl(str(ydot(1)))) // &
                          ", Expected direction=" // trim(adjustl(str(sign(1.0d0, y(2)*hctrvr(3))))))
        end if

        if (hder(3) /= 0d0) then
            if (sign(1.0d0, ydot(2)) == sign(1.0d0, hder(3))) then
                call fail_test("Acceleration direction inconsistent with mirror force: " // &
                              "ydot(2)=" // trim(adjustl(str(ydot(2)))) // &
                              ", B gradient=" // trim(adjustl(str(hder(3)))))
            end if
        end if

        call pass_test()
    end subroutine test_timestep_poloidal

    ! Helper function to convert real to string
    function str(val)
        real(8), intent(in) :: val
        character(len=30) :: str
        write(str, '(ES14.6)') val
    end function str

end program test_timestep
