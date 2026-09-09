module diag_orbit_trace
    use ieee_arithmetic, only: ieee_is_finite
    use iso_fortran_env, only: dp => real64
    use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use neort, only: init, check_magfie
    use neort_config, only: read_and_set_config
    use neort_main, only: runname
    use neort_profiles, only: init_profiles, read_and_init_plasma_input, &
        read_and_init_profile_input, vth, Om_tE
    use neort_freq, only: Om_th, Om_ph
    use neort_resonance, only: valid_resonance_jacobian
    use neort_transport, only: evaluate_hamiltonian, timestep_transport, &
        transport_Omth => Omth, transport_dOmthdv => dOmthdv, &
        transport_dOmthdeta => dOmthdeta
    use neort_datatypes, only: magfie_data_t
    use neort_orbit, only: nvar, th0, magnetic_toroidal_drift_per_v2
    use driftorbit, only: mph, mth, pertfile, sign_vpar, etatp, magdrift, &
        magdrift_passing, supban
    use do_magfie_mod, only: do_magfie, do_magfie_init, R0, s, q
    use do_magfie_pert_mod, only: do_magfie_pert_init
    implicit none

contains

    pure function physical_orientation(vpar_state, hctrvr_theta) result(orientation)
        ! The second orbit state is the signed coordinate velocity state.  The
        ! physical parallel direction has one additional chart factor because
        ! ydot(theta) = vpar_state * hctrvr(3).  Return zero when either factor
        ! is exactly zero (the direction is then undefined).  NaNs also fall
        ! through to zero instead of being admitted as a sign.
        real(dp), intent(in) :: vpar_state, hctrvr_theta
        integer :: orientation

        if (vpar_state == 0.0_dp .or. hctrvr_theta == 0.0_dp) then
            orientation = 0
        else if ((vpar_state > 0.0_dp .and. hctrvr_theta > 0.0_dp) .or. &
                (vpar_state < 0.0_dp .and. hctrvr_theta < 0.0_dp)) then
            orientation = 1
        else if ((vpar_state > 0.0_dp .and. hctrvr_theta < 0.0_dp) .or. &
                (vpar_state < 0.0_dp .and. hctrvr_theta > 0.0_dp)) then
            orientation = -1
        else
            orientation = 0
        end if
    end function physical_orientation

    pure function toroidal_velocity_from_components(vpar, hctrvr_phi, &
            v2_omtb, omte) result(phi_dot)
        ! Thin-orbit guiding-centre toroidal rate in the native chart.
        ! ``vpar*hctrvr_phi`` is the field-line contribution, ``v2_omtb`` is
        ! the magnetic drift rate already multiplied by v**2, and ``omte`` is
        ! the electric drift rate.  Keeping the sum as a pure helper gives the
        ! diagnostic an independent algebraic oracle without changing the
        ! production transport callback.
        real(dp), intent(in) :: vpar, hctrvr_phi, v2_omtb, omte
        real(dp) :: phi_dot

        phi_dot = vpar * hctrvr_phi + v2_omtb + omte
    end function toroidal_velocity_from_components

    subroutine run_orbit_trace_diag(arg_runname, ux_target, eta_target, nsteps, &
            mth_target, emit_physical_phi)
        character(*), intent(in) :: arg_runname
        real(dp), intent(in) :: ux_target, eta_target
        integer, intent(in) :: nsteps, mth_target
        logical, intent(in), optional :: emit_physical_phi

        logical :: file_exists, trapped_orbit
        logical :: physical_phi
        integer :: i, unit, istate, orientation_state, orientation_vpar
        real(dp) :: v, taub, dt, target_time, theta, phi
        real(dp) :: bmod, sqrtg, hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: hctrvr_theta, hctrvr_phi, omtb_v
        real(dp) :: vpar_physical, phi_gc, phi_gc_dot
        real(dp) :: phi_field_dot, phi_magnetic_dot, phi_electric_dot
        real(dp) :: phi_period_avg, phi_period_error
        real(dp) :: y0(nvar), atol(nvar)
        real(dp), allocatable :: yout(:)
        real(dp) :: omph, domphdv, domphdeta, residual, jacobian
        real(dp) :: x(3), H_action_re, H_action_im
        complex(dp) :: Hn
        type(magfie_data_t) :: magfie_data
        type(vode_state_t) :: vstate
        type(fortnum_status_t) :: status
        character(len=160) :: orbit_id

        physical_phi = .false.
        if (present(emit_physical_phi)) physical_phi = emit_physical_phi

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
        trapped_orbit = eta_target > etatp

        if (physical_phi) then
            if (supban) error stop &
                "physical phi diagnostic does not support supban frequency"
        end if

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
        if (physical_phi) then
            write (unit, '(A)') "# schema: iter-tc24-neort-common-orbit-trace-v3"
        else
            write (unit, '(A)') "# schema: iter-tc24-neort-common-orbit-trace-v2"
        end if
        write (unit, '(A,F18.10)') "# s_tor = ", s
        write (unit, '(A,F18.10)') "# rho_tor = ", sqrt(s)
        if (physical_phi) then
            write (unit, '(A)') &
                "# position_coordinates = Boozer(s_tor,phi_gc,theta)"
            write (unit, '(A)') &
                "# toroidal_coordinate = thin_orbit_guiding_center_candidate"
            write (unit, '(A)') &
                "# toroidal_phase_definition = phi_gc(0)=0; dphi_gc/dt = "// &
                "vpar*hctrvr_phi + v**2*Om_tB_over_v2 + Om_tE"
            write (unit, '(A)') &
                "# toroidal_phase_excludes = canonical_phi_H_periodic_Delta_phi "// &
                "and_finite_orbit_chart_terms"
            write (unit, '(A)') &
                "# toroidal_velocity_model = native_thin_orbit; supban=false"
            write (unit, '(A,ES24.16)') "# omph_native = ", omph
            write (unit, '(A)') &
                "# phi_period_identity = mean(phi_gc_dot) = omph_native "// &
                "only on this native thin_orbit model"
        else
            ! The displayed toroidal coordinate follows the unperturbed field
            ! line.  It is not the full guiding-centre toroidal position:
            ! canonical toroidal angle and secular precession are absent.
            write (unit, '(A)') "# position_coordinates = Boozer(s_tor,phi_trace,theta)"
            write (unit, '(A)') "# toroidal_coordinate = field_line_phase_only"
            write (unit, '(A)') "# toroidal_phase_definition = phi_trace=q*(theta-th0)"
            write (unit, '(A)') "# toroidal_phase_excludes = "// &
                "canonical_toroidal_angle_and_Omega_t_secular_drift"
        end if
        ! NEO-RT starts at the local minimum-field point and closes one native
        ! period.  For trapped input this is a full bounce; passing and
        ! separatrix inputs are labelled separately.  MARS' KJPCOEFF trace is
        ! a half-bounce between turning points; recording the class and span
        ! prevents a consumer from pairing equal-looking fractions blindly.
        if (trapped_orbit) then
            write (unit, '(A)') "# orbit_class = trapped"
            write (unit, '(A)') "# orbit_span = full_bounce"
        else if (eta_target < etatp) then
            write (unit, '(A)') "# orbit_class = passing"
            write (unit, '(A)') "# orbit_span = full_transit"
        else
            write (unit, '(A)') "# orbit_class = separatrix"
            write (unit, '(A)') "# orbit_span = separatrix"
        end if
        write (unit, '(A)') "# start_point = local_Bmin"
        write (unit, '(A)') "# end_point = local_Bmin"
        write (unit, '(A)') "# endpoint_bounce_angle = 2*pi"
        write (unit, '(A)') "# time_orientation = increasing_native_time"
        write (unit, '(A)') "# phase_gauge = t=0 at theta=th0 and phi_trace=0"
        write (unit, '(A)') "# state_velocity_convention = vpar_state = "// &
            "sign(hctrvr_theta)*v_parallel"
        write (unit, '(A)') "# orientation_convention = sign(v_parallel) = "// &
            "sign(vpar_state*hctrvr_theta)"
        write (unit, '(A)') "# orientation_zero = 0 when vpar_state or "// &
            "hctrvr_theta is zero"
        write (unit, '(A,F18.10)') "# ux = ", ux_target
        write (unit, '(A,F18.10)') "# eta = ", eta_target
        write (unit, '(A,I0)') "# mth = ", mth
        write (unit, '(A,I0)') "# mph = ", mph
        write (unit, '(A,ES24.16)') "# residual = ", residual
        write (unit, '(A,ES24.16)') "# jacobian_dres_deta = ", jacobian
        write (unit, '(A,ES24.16)') "# taub = ", taub
        write (unit, '(A,A)') "# orbit_id = ", trim(orbit_id)
        if (physical_phi) then
            write (unit, '(A)') "# columns: orbit_id sample_index time time_fraction "// &
                "bounce_angle theta phi_trace s_tor rho_tor ux eta vpar_state "// &
                "bmod hctrvr_theta vpar_physical hctrvr_phi phi_gc phi_gc_dot "// &
                "phi_field_dot phi_magnetic_dot phi_electric_dot phi_period_avg "// &
                "phi_period_error H_inst_re H_inst_im H_action_re H_action_im "// &
                "residual jacobian_dres_deta orientation_state orientation_vpar istate"
        else
            write (unit, '(A)') "# columns: orbit_id sample_index time time_fraction "// &
                "bounce_angle theta phi_trace s_tor rho_tor ux eta vpar_state bmod "// &
                "hctrvr_theta H_inst_re H_inst_im H_action_re H_action_im residual "// &
                "jacobian_dres_deta orientation_state orientation_vpar istate"
        end if

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
            ! This is the drift-free field-line phase used to evaluate the
            ! axisymmetric equilibrium and its toroidal Fourier amplitude.  It
            ! is intentionally not a physical guiding-centre phi trajectory.
            phi = q * (theta - th0)
            x(1) = s
            x(2) = phi
            x(3) = theta
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            hctrvr_theta = hctrvr(3)
            call evaluate_hamiltonian(v, eta_target, target_time, theta, bmod, &
                transport_Omth, Hn)
            ! timestep_transport stores the integrated complex Hamiltonian in
            ! y(3:4); y(5:6) are reserved for nonlinear attenuation moments.
            H_action_re = yout(3)
            H_action_im = yout(4)

            if (physical_phi) then
                if (.not. ieee_is_finite(hctrvr_theta) .or. hctrvr_theta == 0.0_dp) then
                    error stop "physical phi diagnostic encountered zero theta chart factor"
                end if
                hctrvr_phi = hctrvr(2)
                vpar_physical = yout(2) * sign(1.0_dp, hctrvr_theta)
                omtb_v = 0.0_dp
                if (magdrift) then
                    if (trapped_orbit) then
                        omtb_v = magnetic_toroidal_drift_per_v2(eta_target, bmod, &
                            hder(1), hctrvr_theta)
                    else if (magdrift_passing > 0) then
                        omtb_v = magnetic_toroidal_drift_per_v2(eta_target, bmod, &
                            hder(1), hctrvr_theta)
                    end if
                end if
                phi_field_dot = vpar_physical * hctrvr_phi
                phi_magnetic_dot = v**2 * omtb_v
                phi_electric_dot = Om_tE
                phi_gc_dot = toroidal_velocity_from_components(vpar_physical, &
                    hctrvr_phi, phi_magnetic_dot, phi_electric_dot)
                phi_gc = yout(7)
                if (target_time > 0.0_dp) then
                    phi_period_avg = phi_gc / target_time
                    phi_period_error = phi_period_avg - omph
                else
                    phi_period_avg = 0.0_dp
                    phi_period_error = 0.0_dp
                end if
            end if
            if (yout(2) == 0.0_dp) then
                ! An exactly sampled turning point has no signed
                ! orientation.  Do not retain the previous leg's value;
                ! near-zero nonzero values remain one-sided signs.
                orientation_state = 0
            else
                orientation_state = merge(1, -1, yout(2) >= 0.0_dp)
            end if
            orientation_vpar = physical_orientation(yout(2), hctrvr_theta)
            if (physical_phi) then
                write (unit, '(A,1X,I0,1X,27(ES24.16,1X),I0,1X,I0,1X,I0)') &
                    trim(orbit_id), i - 1, target_time, target_time / taub, &
                    target_time * abs(transport_Omth), theta, phi, s, sqrt(s), &
                    ux_target, eta_target, yout(2), bmod, hctrvr_theta, &
                    vpar_physical, hctrvr_phi, phi_gc, phi_gc_dot, phi_field_dot, &
                    phi_magnetic_dot, phi_electric_dot, phi_period_avg, &
                    phi_period_error, real(Hn), aimag(Hn), H_action_re, H_action_im, &
                    residual, jacobian, orientation_state, orientation_vpar, 2
            else
                write (unit, '(A,1X,I0,1X,18(ES24.16,1X),I0,1X,I0,1X,I0)') &
                    trim(orbit_id), i - 1, target_time, target_time / taub, &
                    target_time * abs(transport_Omth), theta, phi, s, sqrt(s), &
                    ux_target, eta_target, yout(2), bmod, hctrvr_theta, &
                    real(Hn), aimag(Hn), H_action_re, H_action_im, residual, &
                    jacobian, orientation_state, orientation_vpar, 2
            end if
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
            if (physical_phi) then
                x(1) = s
                x(2) = 0.0_dp
                x(3) = y_fixed(1)
                call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
                if (.not. ieee_is_finite(hctrvr(3)) .or. hctrvr(3) == 0.0_dp) then
                    error stop "physical phi RHS encountered zero theta chart factor"
                end if
                vpar_physical = y_fixed(2) * sign(1.0_dp, hctrvr(3))
                omtb_v = 0.0_dp
                if (magdrift) then
                    if (trapped_orbit) then
                        omtb_v = magnetic_toroidal_drift_per_v2(eta_target, bmod, &
                            hder(1), hctrvr(3))
                    else if (magdrift_passing > 0) then
                        omtb_v = magnetic_toroidal_drift_per_v2(eta_target, bmod, &
                            hder(1), hctrvr(3))
                    end if
                end if
                dydt_fixed(7) = toroidal_velocity_from_components(vpar_physical, &
                    hctrvr(2), v**2 * omtb_v, Om_tE)
            end if
            dydt_ = dydt_fixed
        end subroutine vode_rhs

    end subroutine run_orbit_trace_diag

end module diag_orbit_trace
