module diag_action_trace
    use ieee_arithmetic, only: ieee_is_finite
    use iso_fortran_env, only: dp => real64
    use neort, only: init, check_magfie, write_magfie_data_to_files, &
        set_to_passing_region, set_to_trapped_region, harmonic_bounds, &
        vsteps, vmax_over_vth, mth_max_abs
    use neort_config, only: read_and_set_config
    use neort_main, only: runname
    use neort_datatypes, only: magfie_data_t
    use neort_profiles, only: read_and_init_profile_input, read_and_init_plasma_input, &
        init_profiles, vth, Om_tE
    use neort_nonlin, only: nonlinear_attenuation
    use neort_freq, only: Om_th, Om_ph
    use neort_transport, only: timestep_transport, Tphi_int, transport_Omth => Omth, &
        transport_dOmthdv => dOmthdv, transport_dOmthdeta => dOmthdeta
    use neort_orbit, only: bounce_fast, nvar
    use neort_resonance, only: driftorbit_coarse, driftorbit_root
    use neort_action_trace_contract, only: action_trace_weight
    use driftorbit, only: mth, mph, mi, nlev, pertfile, nonlin, etamin, etamax, &
        sign_vpar, nopassing
    use do_magfie_mod, only: R0, s, q, do_magfie_init
    use do_magfie_pert_mod, only: do_magfie_pert_init
    implicit none

contains

    subroutine run_action_trace_diag(arg_runname)
        character(*), intent(in) :: arg_runname
        logical :: file_exists
        integer :: mth_min, mth_max, j, unit
        real(dp) :: vminp, vmaxp, vmint, vmaxt
        type(magfie_data_t) :: magfie_data

        runname = trim(arg_runname)
        call read_and_set_config(trim(runname)//".in")
        if (vsteps <= 0) error stop "action trace requires positive vsteps"
        call do_magfie_init("in_file")
        if (pertfile) call do_magfie_pert_init("in_file_pert")
        call init_profiles(R0)

        inquire (file="plasma.in", exist=file_exists)
        if (file_exists) call read_and_init_plasma_input("plasma.in", s)
        inquire (file="profile.in", exist=file_exists)
        if (file_exists) call read_and_init_profile_input("profile.in", s, R0, 1.0_dp, 1.0_dp)

        call init
        call check_magfie(magfie_data)
        call write_magfie_data_to_files(magfie_data, runname)
        call harmonic_bounds(mph, q, mth_max_abs, mth_min, mth_max)

        vminp = 1.0e-6_dp * vth
        vmaxp = vmax_over_vth * vth
        vmint = vminp
        vmaxt = vmaxp

        open (newunit=unit, file=trim(runname)//"_action_trace.dat", status="replace", action="write")
        write (unit, '(A,F18.10)') "# s_tor = ", s
        write (unit, '(A,I0)') "# mph = ", mph
        write (unit, '(A,I0,A,I0)') "# mth range = ", mth_min, "..", mth_max
        write (unit, '(A,I0,A,ES18.10)') "# vsteps = ", vsteps, " vmax_over_vth = ", vmax_over_vth
        write (unit, '(A)') "# columns: branch mth ux v eta dR_deta root_lo root_hi Omth Omph residual dOmth_deta dOmph_deta taub bounce_re bounce_im Hmn2 attenuation Tphi_int du weighted_T istate"

        do j = mth_min, mth_max
            mth = j
            if (.not. nopassing) then
                sign_vpar = 1.0_dp
                call set_to_passing_region(etamin, etamax)
                call trace_branch("passing_co", vminp, vmaxp, unit)
                sign_vpar = -1.0_dp
                call set_to_passing_region(etamin, etamax)
                call trace_branch("passing_ctr", vminp, vmaxp, unit)
            end if
            sign_vpar = 1.0_dp
            call set_to_trapped_region(etamin, etamax)
            call trace_branch("trapped", vmint, vmaxt, unit)
        end do
        close (unit)
    end subroutine run_action_trace_diag

    subroutine trace_branch(label, vmin, vmax, unit)
        character(*), intent(in) :: label
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: unit
        integer :: ku, kr, nroots, istate_dv
        real(dp) :: du, ux, v, eta
        real(dp) :: roots(nlev, 3), eta_res(2)
        real(dp) :: omph, domphdv, domphdeta, residual
        real(dp) :: taub, bounceavg(nvar), hmn2, attenuation
        real(dp) :: tphi, weighted

        du = (vmax - vmin) / (real(vsteps, dp) * vth)
        ux = vmin / vth + du / 2.0_dp
        do ku = 1, vsteps
            v = ux * vth
            call driftorbit_coarse(v, etamin, etamax, roots, nroots)
            do kr = 1, nroots
                eta_res = driftorbit_root(v, 1.0e-8_dp * abs(Om_tE), roots(kr, 1), roots(kr, 2))
                if (.not. ieee_is_finite(eta_res(1))) error stop "nonfinite resonance root"
                if (eta_res(1) < 0.0_dp) then
                    if (eta_res(1) == -1.0_dp) cycle
                    error stop "unconverged resonance root"
                end if
                if (.not. ieee_is_finite(eta_res(2)) .or. eta_res(2) == 0.0_dp) then
                    error stop "nonfinite or zero resonance Jacobian"
                end if
                eta = eta_res(1)
                ! The production callback reads neort_transport's threadprivate
                ! Omth while bounce_fast advances the orbit.  Bind the same
                ! module state here; a local frequency gives a numerically
                ! plausible but non-equivalent action trace.
                call Om_th(v, eta, transport_Omth, transport_dOmthdv, transport_dOmthdeta)
                call Om_ph(v, eta, omph, domphdv, domphdeta)
                residual = real(mth, dp) * transport_Omth + real(mph, dp) * omph
                taub = 2.0_dp * acos(-1.0_dp) / abs(transport_Omth)
                call bounce_fast(v, eta, taub, bounceavg, timestep_transport, istate_dv)
                if (istate_dv /= 2) error stop "non-success bounce status"
                if (.not. all(ieee_is_finite(bounceavg))) error stop "nonfinite bounce average"
                hmn2 = (bounceavg(3)**2 + bounceavg(4)**2) * &
                    (mi * (ux * vth)**2 / 2.0_dp)**2
                attenuation = nonlinear_attenuation(ux, eta, bounceavg, transport_Omth, &
                    transport_dOmthdv, transport_dOmthdeta, hmn2)
                tphi = Tphi_int(ux, taub, hmn2)
                weighted = action_trace_weight(du, tphi, eta_res(2), attenuation)
                write (unit, '(A,1X,I0,1X,19(ES24.16,1X),I0)') trim(label), mth, &
                    ux, v, eta, eta_res(2), roots(kr, 1), roots(kr, 2), transport_Omth, omph, &
                    residual, transport_dOmthdeta, domphdeta, taub, bounceavg(3), bounceavg(4), &
                    hmn2, attenuation, tphi, du, weighted, istate_dv
            end do
            ux = ux + du
        end do
    end subroutine trace_branch

end module diag_action_trace
