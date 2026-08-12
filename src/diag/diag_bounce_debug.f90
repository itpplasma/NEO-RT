module diag_bounce_debug
    use iso_fortran_env, only: dp => real64
    use neort, only: init, check_magfie, write_magfie_data_to_files, &
        set_to_passing_region, set_to_trapped_region, vsteps, &
        vmax_over_vth
    use neort_config, only: read_and_set_config
    use neort_main, only: runname
    use neort_datatypes, only: magfie_data_t
    use neort_profiles, only: read_and_init_profile_input, read_and_init_plasma_input, init_profiles, vth, Om_tE
    use neort_freq, only: Om_th, Om_ph
    use neort_transport, only: timestep_transport
    use neort_orbit, only: nvar
    use neort_resonance, only: driftorbit_coarse, driftorbit_root
    use driftorbit, only: mth, mph, pertfile, etamin, etamax, sign_vpar, &
        etatp, etadt, efac
    use do_magfie_mod, only: R0, s, q, bfac, do_magfie_init
    use do_magfie_pert_mod, only: do_magfie_pert_init
    use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t
    use fortnum_ode_dop853, only: ode_integrate_dop
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    implicit none

contains

    subroutine run_bounce_debug(arg_runname)
        character(*), intent(in) :: arg_runname
        character(32) :: arg_mth, arg_branch
        logical :: file_exists
        integer :: j, mmin, mmax, nm, u, target_mth
        logical :: filter_mth
        real(dp) :: vminp, vmaxp, vmint, vmaxt
        real(dp) :: du, ux, v
        real(dp) :: roots(100,3), eta_res(2), eta
        integer :: nroots, kr
        real(dp) :: Omth, dOmthdv, dOmthdeta, taub
        integer :: istate
        integer :: istats(50)
        real(dp) :: rstats(50)
        type(magfie_data_t) :: magfie_data

        call get_command_argument(3, arg_mth)
        call get_command_argument(4, arg_branch)
        filter_mth = len_trim(arg_mth) > 0
        if (filter_mth) then
            read(arg_mth, *, err=900) target_mth
        else
            target_mth = 0
        end if
        if (len_trim(arg_branch) == 0) arg_branch = 'all'
        if (all(trim(arg_branch) /= ['all ', 'cop ', 'ctr ', 'trap'])) then
            error stop 'bounce_debug branch must be all, cop, ctr, or trap'
        end if

        ! Initialize environment just like main
        runname = trim(arg_runname)
        call read_and_set_config(trim(runname)//".in")
        call do_magfie_init("in_file")
        if (pertfile) call do_magfie_pert_init("in_file_pert")
        call init_profiles(R0)

        inquire(file="plasma.in", exist=file_exists)
        if (file_exists) call read_and_init_plasma_input("plasma.in", s)

        inquire(file="profile.in", exist=file_exists)
        if (file_exists) call read_and_init_profile_input("profile.in", s, R0, efac, bfac)

        call init
        call check_magfie(magfie_data)
        call write_magfie_data_to_files(magfie_data, runname)

        vminp = 1.0e-6_dp*vth
        vmaxp = vmax_over_vth*vth
        vmint = vminp
        vmaxt = vmaxp

        mmin = -ceiling(2.0_dp*abs(mph*q))
        mmax =  ceiling(2.0_dp*abs(mph*q))
        nm = mmax - mmin + 1

        open(newunit=u, file=trim(arg_runname)//'_bounce_debug.txt', status='replace', action='write')
        write(u,'(A)') '# bounce debug: logs where integration stalls (MXSTEP)'
        write(u,'(A,F10.6)') '# s = ', s
        write(u,'(A)') &
            '# cols: branch mth ux eta Omth Omph residual dres_deta taub istate_probe nst_probe last_h_probe'

        do j = 1, nm
            mth = mmin + (j - 1)
            if (filter_mth .and. mth /= target_mth) cycle

            ! Passing co
            if (trim(arg_branch) == 'all' .or. trim(arg_branch) == 'cop') then
                sign_vpar = 1.0_dp
                call set_to_passing_region(etamin, etamax)
                call scan_branch('cop', vminp, vmaxp, vsteps, u)
            end if

            ! Passing ctr
            if (trim(arg_branch) == 'all' .or. trim(arg_branch) == 'ctr') then
                sign_vpar = -1.0_dp
                call set_to_passing_region(etamin, etamax)
                call scan_branch('ctr', vminp, vmaxp, vsteps, u)
            end if

            ! Trapped
            if (trim(arg_branch) == 'all' .or. trim(arg_branch) == 'trap') then
                sign_vpar = 1.0_dp
                call set_to_trapped_region(etamin, etamax)
                call scan_branch('trap', vmint, vmaxt, vsteps, u)
            end if
        end do
        close(u)
        return

        900 error stop 'bounce_debug mth must be an integer'
    end subroutine run_bounce_debug

    subroutine scan_branch(label, vmin, vmax, vsteps_loc, u)
        character(*), intent(in) :: label
        real(dp), intent(in) :: vmin, vmax
        integer, intent(in) :: vsteps_loc
        integer, intent(in) :: u
        real(dp) :: du, ux, v
        real(dp) :: roots(100,3), eta_res(2), eta
        integer :: nroots, kr
        real(dp) :: Omth, dOmthdv, dOmthdeta, Omph, dOmphdv, dOmphdeta
        real(dp) :: residual, taub
        integer :: istate
        integer :: istats(50)
        real(dp) :: rstats(50)

        du = (vmax - vmin)/(real(vsteps_loc, dp)*vth)
        ux = vmin/vth + du/2.0_dp
        do while (ux*vth <= vmax + 1.0e-12_dp)
            v = ux*vth
            call driftorbit_coarse(v, etamin, etamax, roots, nroots)
            if (nroots > 0) then
                do kr = 1, nroots
                    eta_res = driftorbit_root(v, 1.0e-8_dp*abs(Om_tE), roots(kr,1), roots(kr,2))
                    if (eta_res(1) < 0.0_dp) cycle ! bracket-failure sentinel
                    eta = eta_res(1)
                    call Om_th(v, eta, Omth, dOmthdv, dOmthdeta)
                    call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
                    residual = real(mth, dp)*Omth + real(mph, dp)*Omph
                    taub = 2.0_dp*acos(-1.0_dp)/abs(Omth)
                    call probe_bounce(v, eta, taub, istate, istats, rstats)
                    write(u,'(A,1X,I4,1X,F12.8,1X,6(ES22.14E3,1X),I4,1X,I9,1X,ES14.6E3)') &
                        trim(label), mth, ux, eta, Omth, Omph, residual, eta_res(2), &
                        taub, istate, istats(1), rstats(6)
                end do
            end if
            ux = ux + du
        end do
    end subroutine scan_branch

    subroutine probe_bounce(v, eta, taub, istate, istats, rstats)
        real(dp), intent(in) :: v, eta, taub
        integer, intent(out) :: istate
        integer, intent(out) :: istats(:)
        real(dp), intent(out) :: rstats(:)

        real(dp) :: y0(nvar)
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status

        istats = 0
        rstats = 0.0_dp

        y0 = 1.0e-15_dp
        y0(1) = 0.0_dp
        y0(2) = 0.0_dp
        y0(3:6) = 0.0_dp

        problem%rhs => probe_rhs
        problem%t0 = 0.0_dp
        problem%t1 = taub
        problem%y0 = y0
        problem%rtol = 1.0e-8_dp
        problem%atol = 1.0e-9_dp

        call ode_integrate_dop(problem, workspace, solution, status)

        ! Report the integrator diagnostics in the slots scan_branch prints:
        ! istats(1) = accepted steps, rstats(6) = last step size.
        if (status%code == FORTNUM_OK) then
            istate = 2
        else
            istate = -1
        end if
        istats(1) = solution%nsteps
        istats(2) = solution%nrejected
        if (allocated(solution%h) .and. solution%nsteps > 0) then
            rstats(6) = solution%h(solution%nsteps)
        end if

    contains
        subroutine probe_rhs(t_, y_, dydt_, ctx_)
            real(dp), intent(in) :: t_
            real(dp), intent(in) :: y_(:)
            real(dp), intent(out) :: dydt_(:)
            class(*), intent(in), optional :: ctx_
            associate (dummy => ctx_)
            end associate
            call timestep_transport(v, eta, size(y_), t_, y_, dydt_)
        end subroutine probe_rhs
    end subroutine probe_bounce

end module diag_bounce_debug
