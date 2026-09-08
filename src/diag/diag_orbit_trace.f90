module diag_orbit_trace
    use ieee_arithmetic, only: ieee_is_finite
    use iso_fortran_env, only: dp => real64
    use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use neort, only: init, check_magfie
    use neort_config, only: read_and_set_config
    use neort_main, only: runname
    use neort_profiles, only: init_profiles, read_and_init_plasma_input, &
        read_and_init_profile_input, vth
    use neort_freq, only: Om_th, Om_ph
    use neort_resonance, only: valid_resonance_jacobian
    use neort_transport, only: evaluate_hamiltonian, timestep_transport, &
        transport_Omth => Omth, transport_dOmthdv => dOmthdv, &
        transport_dOmthdeta => dOmthdeta
    use neort_datatypes, only: magfie_data_t
    use neort_orbit, only: nvar, th0
    use driftorbit, only: mph, mth, pertfile, sign_vpar
    use do_magfie_mod, only: do_magfie, do_magfie_init, R0, s, q
    use do_magfie_pert_mod, only: do_magfie_pert_init
    implicit none

contains

    subroutine run_orbit_trace_diag(arg_runname, ux_target, eta_target, nsteps, mth_target)
        character(*), intent(in) :: arg_runname
        real(dp), intent(in) :: ux_target, eta_target
        integer, intent(in) :: nsteps, mth_target

        logical :: file_exists
        integer :: i, unit, istate, orientation
        real(dp) :: v, taub, dt, target_time, theta, phi
        real(dp) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: y0(nvar), atol(nvar)
        real(dp), allocatable :: yout(:)
        real(dp) :: omph, domphdv, domphdeta, residual, jacobian
        real(dp) :: x(3), H_action_re, H_action_im
        complex(dp) :: Hn
        type(magfie_data_t) :: magfie_data
        type(vode_state_t) :: vstate
        type(fortnum_status_t) :: status
        character(len=160) :: orbit_id

        if (.not. ieee_is_finite(ux_target) .or. ux_target <= 0.0_dp) then
            error stop "orbit_trace requires finite positive ux"
        end if
        if (.not. ieee_is_finite(eta_target) .or. eta_target <= 0.0_dp) then
            error stop "orbit_trace requires finite positive eta"
        end if
        if (nsteps < 2) error stop "orbit_trace requires at least two samples"

        runname = trim(arg_runname)
        call read_and_set_config(trim(runname)//".in")
        call do_magfie_init("in_file")
        if (pertfile) call do_magfie_pert_init("in_file_pert")
        call init_profiles(R0)

        inquire (file="plasma.in", exist=file_exists)
        if (file_exists) call read_and_init_plasma_input("plasma.in", s)
        inquire (file="profile.in", exist=file_exists)
        if (file_exists) call read_and_init_profile_input("profile.in", s, R0, &
            1.0_dp, 1.0_dp)

        call init
        call check_magfie(magfie_data)
        mth = mth_target

        if (s < 0.0_dp) error stop "orbit_trace requires nonnegative s_tor"

        v = ux_target * vth
        x(1) = s
        x(2) = 0.0_dp
        x(3) = th0
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        if (.not. ieee_is_finite(bmod) .or. bmod <= 0.0_dp) then
            error stop "orbit_trace encountered invalid local field"
        end if
        if (eta_target * bmod >= 1.0_dp) then
            error stop "orbit_trace eta is outside the local trapped/passing speed domain"
        end if

        sign_vpar = 1.0_dp
        call Om_th(v, eta_target, transport_Omth, transport_dOmthdv, &
            transport_dOmthdeta)
        call Om_ph(v, eta_target, omph, domphdv, domphdeta)
        residual = real(mth, dp) * transport_Omth + real(mph, dp) * omph
        jacobian = real(mth, dp) * transport_dOmthdeta + &
            real(mph, dp) * domphdeta
        if (.not. ieee_is_finite(transport_Omth) .or. &
            transport_Omth == 0.0_dp) then
            error stop "orbit_trace encountered invalid poloidal frequency"
        end if
        if (.not. ieee_is_finite(jacobian) .or. jacobian == 0.0_dp) then
            error stop "orbit_trace encountered invalid resonance Jacobian"
        end if
        if (.not. valid_resonance_jacobian(jacobian)) then
            error stop "orbit_trace encountered unrepresentable resonance Jacobian"
        end if

        taub = 2.0_dp * acos(-1.0_dp) / abs(transport_Omth)
        dt = taub / real(nsteps - 1, dp)
        y0 = 1.0e-15_dp
        y0(1) = th0
        y0(2) = sign(1.0_dp, hctrvr(3)) * &
            v * sqrt(1.0_dp - eta_target * bmod)
        y0(3:) = 0.0_dp
        atol = 1.0e-10_dp
        call vode_init(vstate, nvar, 0.0_dp, y0)

        orbit_id = trim(runname)//"_orbit_0001"
        open (newunit=unit, file=trim(runname)//"_orbit_trace.dat", &
            status="replace", action="write")
        write (unit, '(A)') "# schema: iter-tc24-neort-common-orbit-trace-v1"
        write (unit, '(A,F18.10)') "# s_tor = ", s
        write (unit, '(A,F18.10)') "# rho_tor = ", sqrt(s)
        write (unit, '(A)') "# position_coordinates = Boozer(s_tor,phi,theta)"
        write (unit, '(A)') "# phase_gauge = t=0 at theta=th0 and phi=0"
        write (unit, '(A)') "# orientation_convention = sign(v_parallel)"
        write (unit, '(A,F18.10)') "# ux = ", ux_target
        write (unit, '(A,F18.10)') "# eta = ", eta_target
        write (unit, '(A,I0)') "# mth = ", mth
        write (unit, '(A,I0)') "# mph = ", mph
        write (unit, '(A,ES24.16)') "# residual = ", residual
        write (unit, '(A,ES24.16)') "# jacobian_dres_deta = ", jacobian
        write (unit, '(A,ES24.16)') "# taub = ", taub
        write (unit, '(A,A)') "# orbit_id = ", trim(orbit_id)
        write (unit, '(A)') "# columns: orbit_id sample_index time time_fraction "// &
            "bounce_angle theta phi s_tor rho_tor ux eta vpar bmod "// &
            "H_inst_re H_inst_im H_action_re H_action_im residual "// &
            "jacobian_dres_deta orientation istate"

        do i = 1, nsteps
            target_time = dt * real(i - 1, dp)
            if (i == 1) then
                yout = y0
            else
                call vode_integrate_to(vode_rhs, vstate, target_time, 1.0e-9_dp, &
                    atol, yout, status)
                if (status%code /= FORTNUM_OK) then
                    error stop "orbit_trace ODE integration failed"
                end if
            end if

            theta = yout(1)
            phi = q * (theta - th0)
            x(1) = s
            x(2) = phi
            x(3) = theta
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            call evaluate_hamiltonian(v, eta_target, target_time, theta, bmod, &
                transport_Omth, Hn)
            ! timestep_transport stores the integrated complex Hamiltonian in
            ! y(3:4); y(5:6) are reserved for nonlinear attenuation moments.
            H_action_re = yout(3)
            H_action_im = yout(4)
            if (i == 1) then
                orientation = merge(1, -1, yout(2) >= 0.0_dp)
            else if (abs(yout(2)) > 1.0e-12_dp * max(v, 1.0_dp)) then
                orientation = merge(1, -1, yout(2) >= 0.0_dp)
            end if
            write (unit, '(A,1X,I0,1X,17(ES24.16,1X),I0,1X,I0)') &
                trim(orbit_id), i - 1, target_time, target_time / taub, &
                target_time * abs(transport_Omth), theta, phi, s, sqrt(s), &
                ux_target, eta_target, yout(2), bmod, real(Hn), aimag(Hn), &
                H_action_re, H_action_im, residual, jacobian, orientation, 2
        end do
        close (unit)

    contains

        subroutine vode_rhs(t_, y_, dydt_, ctx_)
            real(dp), intent(in) :: t_
            real(dp), intent(in) :: y_(:)
            real(dp), intent(out) :: dydt_(:)
            class(*), intent(in), optional :: ctx_
            real(dp) :: y_fixed(nvar), dydt_fixed(nvar)
            associate (dummy => ctx_)
            end associate
            y_fixed = y_
            call timestep_transport(v, eta_target, nvar, t_, y_fixed, dydt_fixed)
            dydt_ = dydt_fixed
        end subroutine vode_rhs

    end subroutine run_orbit_trace_diag

end module diag_orbit_trace
