module neort_orbit
    use iso_fortran_env, only: dp => real64
    use logger, only: debug, trace, get_log_level, LOG_TRACE, error
    use util, only: imun, pi, mi, qi, c
    use spline, only: spline_coeff, spline_val_0
    use do_magfie_mod, only: do_magfie, s, iota, R0, eps, psi_pr, &
        bphcov, dbthcovds, dbphcovds, q, dqds, sign_theta
    use do_magfie_pert_mod, only: do_magfie_pert_amp
    use neort_profiles, only: vth, Om_tE, dOm_tEds
    use driftorbit, only: etatp, etadt, etamin, etamax, epsmn, mth, mph, m0, mth, &
        pertfile, magdrift, nonlin, epsst_spl, epssp_spl, epst_spl, epsp_spl, &
        sign_vpar, sign_vpar_htheta

    implicit none

    integer, parameter :: nvar = 7
    real(dp) :: th0 = 0.0_dp

    logical :: noshear = .false.      ! neglect magnetic shear

    !$omp threadprivate (th0)

    interface
        subroutine timestep_i(v, eta, neq, t, y, ydot)
            import :: dp
            real(dp), intent(in) :: v, eta
            integer, intent(in) :: neq
            real(dp), intent(in) :: t
            real(dp), intent(in) :: y(neq)
            real(dp), intent(out) :: ydot(neq)
        end subroutine timestep_i
    end interface

contains

    subroutine dvode_error_context(where, v_in, eta_in, tcur, tout, ist)
        use do_magfie_mod, only: s, iota, R0, q, psi_pr, eps
        use driftorbit, only: etatp, etadt, etamin, etamax, mth, mph, sign_vpar
        use neort_profiles, only: vth, Om_tE
        character(*), intent(in) :: where
        real(dp), intent(in) :: v_in, eta_in, tcur, tout
        integer, intent(in) :: ist
        character(len=512) :: msg
        character(len=64) :: reg
        if (eta_in < etatp) then
            reg = 'passing'
        else
            reg = 'trapped'
        end if
        write(0,'(A,1X,A,1X,A)') '[ERROR] DVODE MXSTEP in', trim(where), trim(reg)
        write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  v=', v_in, 'eta=', eta_in
        write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,I0)') '  tcur=', tcur, 'tout=', tout, 'istate=', ist
        write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  vth=', vth, 'Om_tE=', Om_tE
        write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') &
             '  s=', s, 'R0=', R0, 'q=', q, 'iota=', iota
        write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') &
             '  etatp=', etatp, 'etadt=', etadt, 'etamin=', etamin, 'etamax=', etamax
        write(0,'(A,1X,I0,2X,A,1X,I0,2X,A,1X,ES12.5)') '  mth=', mth, 'mph=', mph, 'sign_vpar=', dble(sign_vpar)
        write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  eps=', eps, 'psi_pr=', psi_pr
        call error('DVODE MXSTEP')
    end subroutine dvode_error_context

    pure function to_es(x) result(sout)
        real(dp), intent(in) :: x
        character(len=24) :: sout
        write(sout,'(ES12.5)') x
    end function to_es

    pure function to_i(i) result(sout)
        integer, intent(in) :: i
        character(len=12) :: sout
        write(sout,'(I0)') i
        sout = adjustl(sout)
    end function to_i

    subroutine bounce(v, eta, taub, bounceavg, taub_estimate)
        ! calculate all bounce averages
        real(dp), intent(in) :: v, eta
        real(dp), intent(out) :: taub, bounceavg(nvar)
        real(dp), optional :: taub_estimate  ! estimated bounce time (user input)
        real(dp) :: findroot_res(nvar + 1)
        real(dp) :: bmod, htheta
        real(dp) :: y0(nvar)

        ! Initialize bounce-averated quantities y0. Their meaning
        ! is defined inside subroutine timestep (thin orbit integration)
        call evaluate_bfield_local(bmod, htheta)
        sign_vpar_htheta = sign(1.0_dp, htheta)*sign_vpar
        y0 = 1.0e-15_dp
        y0(1) = th0  ! poloidal angle theta
        y0(2) = sign_vpar_htheta * vpar(v, eta, bmod)  ! parallel velocity vpar
        y0(3) = 0.0_dp  ! toroidal velocity v_ph for drift frequency Om_ph
        y0(4) = 0.0_dp  ! perturbed Hamiltonian real part
        y0(5) = 0.0_dp  ! perturbed Hamiltonian imaginary part
        y0(6) = 0.0_dp  ! 1/abs(B)
        ! y0(7) = 0.0_dp       ! abs(B)

        ! If bounce time estimate exists (elliptic integrals),
        ! initialize taub with it, owtherwise estimate here.
        if (present(taub_estimate)) then
            taub = taub_estimate
        else
            taub = 2.0 * pi / abs(vperp(v, eta, bmod) * iota / R0 * sqrt(eps / 2.0_dp))
        end if

        ! Look for exactly one orbit turn via root-finding.
        ! Start by looking for 5 points per turn.
        findroot_res = bounce_integral(v, eta, nvar, y0, taub / 5.0_dp, timestep)

        taub = findroot_res(1)
        bounceavg = findroot_res(2:) / taub
    end subroutine bounce

    subroutine evaluate_bfield_local(bmod, htheta)
        real(dp), intent(out) :: bmod, htheta
        real(dp) :: sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0.0_dp
        x(3) = th0
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        htheta = hctrvr(3)
    end subroutine evaluate_bfield_local

    pure function vpar(v, eta, bmod)
        !   parallel velocity
        real(dp) :: vpar
        real(dp), intent(in) :: v, eta, bmod
        vpar = v * sqrt(1.0_dp - eta * bmod)
        if (isnan(vpar)) then
            vpar = 0.0_dp
        end if
    end function vpar

    pure function vperp(v, eta, bmod)
        !   perpendicular velocity
        real(dp) :: vperp
        real(dp), intent(in) :: v, eta, bmod
        vperp = v*sqrt(eta*bmod)
        if (isnan(vperp)) then
            vperp = 0.0_dp
        end if
    end function vperp

    subroutine bounce_fast(v, eta, taub, bounceavg, ts, istate_out)
        use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t
        use fortnum_ode_dop853, only: ode_integrate_dop
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
            FORTNUM_CONVERGENCE_ERROR

        real(dp), intent(in) :: v, eta, taub
        real(dp), intent(out) :: bounceavg(nvar)
        procedure(timestep_i) :: ts
        integer, intent(out), optional :: istate_out

        ! Resolve the bounce period into at least this many DOP853 steps by
        ! capping hmax at taub/max_steps_per_turn. The bounceavg(3:4) integrands
        ! carry the perturbed Hamiltonian, which oscillates at the resonant
        ! harmonic (up to ~mth + q*mph cycles per bounce). Those components are
        ! driven, not fed back into the orbit dynamics, so the rtol error
        ! estimate barely constrains them and the unconstrained adaptive step
        ! grows too coarse, leaving the oscillatory integral under-resolved. The
        ! cap recovers the DVODE result to four significant figures.
        integer, parameter :: max_steps_per_turn = 200

        real(dp) :: t1, t2, bmod, htheta
        real(dp) :: y0(nvar)
        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(fortnum_status_t) :: status
        integer :: istate

        call trace('bounce_fast')

        t1 = 0.0_dp
        t2 = taub

        call evaluate_bfield_local(bmod, htheta)
        sign_vpar_htheta = sign(1.0_dp, htheta) * sign_vpar
        y0 = 1.0e-15_dp
        y0(1) = th0
        y0(2) = sign_vpar_htheta * vpar(v, eta, bmod)
        y0(3:6) = 0.0_dp

        problem%rhs => bounce_rhs
        problem%t0 = t1
        problem%t1 = t2
        problem%y0 = y0
        problem%rtol = 1.0e-9_dp
        problem%atol = 1.0e-10_dp
        problem%hmax = abs(taub) / real(max_steps_per_turn, dp)

        call ode_integrate_dop(problem, workspace, solution, status)

        ! Map the fortnum status onto the legacy istate convention the callers
        ! already branch on: 2 = success, -1 = step budget exhausted (the old
        ! DVODE MXSTEP case), anything else = unexpected failure.
        if (status%code == FORTNUM_OK) then
            istate = 2
        else if (status%code == FORTNUM_CONVERGENCE_ERROR) then
            istate = -1
        else
            istate = 0
        end if
        if (istate == -1) then
            call dvode_error_context('bounce_fast', v, eta, t1, t2, istate)
        end if

        bounceavg = solution%y(:, size(solution%t)) / taub
        if (present(istate_out)) istate_out = istate

        call trace('bounce_fast complete')

    contains

        subroutine bounce_rhs(t_, y_, dydt_, rhs_ctx)
            real(dp), intent(in) :: t_
            real(dp), intent(in) :: y_(:)
            real(dp), intent(out) :: dydt_(:)
            class(*), intent(in), optional :: rhs_ctx
            associate (dummy => rhs_ctx)
            end associate
            call ts(v, eta, size(y_), t_, y_, dydt_)
        end subroutine bounce_rhs
    end subroutine bounce_fast

    function bounce_time(v, eta, taub_estimate) result(taub)

        real(dp), intent(in) :: v, eta
        real(dp), intent(in), optional :: taub_estimate
        real(dp) :: taub

        integer, parameter :: neq = 2
        real(dp) :: y0(neq), roots(neq + 1)
        real(dp) :: bmod, htheta
        call trace('bounce_time')

        call evaluate_bfield_local(bmod, htheta)
        sign_vpar_htheta = sign(1.0_dp, htheta) * sign_vpar

        y0(1) = th0  ! poloidal angle theta
        y0(2) = sign_vpar_htheta * vpar(v, eta, bmod)  ! parallel velocity vpar

        if (present(taub_estimate)) then
            taub = taub_estimate
        else
            taub = 2.0 * pi / abs(vperp(v, eta, bmod) * iota / R0 * sqrt(eps / 2.0_dp))
        end if

        roots = bounce_integral(v, eta, neq, y0, taub, timestep_poloidal_motion)
        taub = roots(1)
        call trace('bounce_time complete')
    end function bounce_time

    subroutine timestep_poloidal_motion(v, eta, neq, t, y, ydot)
        real(dp), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(dp), intent(in) :: t, y(neq)
        real(dp), intent(out) :: ydot(neq)

        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)

        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        call poloidal_velocity(v, eta, bmod, hctrvr(3), hder(3), y(2), ydot)
    end subroutine timestep_poloidal_motion

    pure subroutine poloidal_velocity(v, eta, bmod, hthctr, hderth, v_par, ydot)
        real(dp), intent(in) :: v, eta, bmod, hthctr, hderth
        real(dp), intent(in) :: v_par
        real(dp), intent(out) :: ydot(2)

        ydot(1) = v_par * hthctr  ! theta
        ydot(2) = -v**2 * eta / 2.0_dp * hthctr * hderth * bmod  ! v_par
    end subroutine poloidal_velocity

    function bounce_integral(v, eta, neq, y0, dt, ts)
        !
        !  Finds the root of an orbit after the first turn
        !
        use fortnum_ode, only: ode_problem_t, ode_workspace_t, ode_solution_t, &
            ODE_EVENT_RISING, ODE_EVENT_FALLING
        use fortnum_ode_dop853, only: ode_integrate_dop
        use fortnum_ode_events, only: ode_event_scan, ode_event_result_t
        use fortnum_status, only: fortnum_status_t, FORTNUM_CONVERGENCE_ERROR

        real(dp) :: bounce_integral(neq + 1)
        real(dp), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(dp), intent(in) :: y0(neq), dt
        procedure(timestep_i) :: ts

        ! Number of dt-sized turns the search window spans. Matches the old
        ! chunked DVODE search, which advanced in dt steps up to this many turns.
        integer, parameter :: n_turns = 500

        type(ode_problem_t) :: problem
        type(ode_workspace_t) :: workspace
        type(ode_solution_t) :: solution
        type(ode_event_result_t) :: event_res
        type(fortnum_status_t) :: status
        logical :: passing
        integer :: direction

        passing = (eta < etatp)

        ! Passing orbits return after a full 2*pi turn in the sign_vpar_htheta
        ! direction; trapped orbits return to th0 with theta rising from below
        ! (the original (yold(1)-th0) < 0 accept filter). Pick the event branch
        ! and crossing direction accordingly.
        if (passing) then
            if (sign_vpar_htheta > 0.0_dp) then
                direction = ODE_EVENT_RISING
            else
                direction = ODE_EVENT_FALLING
            end if
        else
            direction = ODE_EVENT_RISING
        end if

        problem%rhs => bounce_int_rhs
        problem%t0 = 0.0_dp
        problem%t1 = real(n_turns - 1, dp) * dt
        problem%y0 = y0
        problem%rtol = 1.0e-9_dp
        problem%atol = 1.0e-10_dp
        problem%event => bounce_event
        problem%event_direction = direction
        problem%terminal_event = .true.

        if (get_log_level() >= LOG_TRACE) then
            write(*,'(A,2ES12.5,2A)') '[TRACE] bounce_integral start v,eta=', v, eta, &
                ' pass=', merge('T','F', passing)
        end if

        call ode_integrate_dop(problem, workspace, solution, status)
        if (status%code == FORTNUM_CONVERGENCE_ERROR) then
            call dvode_error_context('bounce_integral', v, eta, problem%t0, problem%t1, -1)
        end if

        call ode_event_scan(problem%rhs, problem%event, solution, &
            problem%event_direction, problem%event_tol, event_res, status)

        if (.not. event_res%found) then
            write(0,'(A)') '[ERROR] bounce_integral: no bounce event located'
            write(0,'(A,1X,A)') '  region =', merge('passing','trapped', passing)
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  v =', v, 'eta =', eta
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  t1 =', problem%t1, 'dt =', dt
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') '  etamin =', etamin, 'etamax =', etamax, 'etatp =', etatp
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  theta(y1) =', solution%y(1, size(solution%t)), 'th0 =', th0
            write(0,'(A,1X,I0,2X,A,1X,I0,2X,A,1X,ES12.5)') '  mth =', mth, 'mph =', mph, 'sign_vpar =', dble(sign_vpar)
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') '  s =', s, 'R0 =', R0, 'q =', q, 'iota =', iota
            bounce_integral(1) = solution%t(size(solution%t))
            bounce_integral(2:) = solution%y(:, size(solution%t))
            return
        end if

        bounce_integral(1) = event_res%t_event
        bounce_integral(2:) = event_res%y_event

    contains

        subroutine bounce_int_rhs(t_, y_, dydt_, ctx_)
            real(dp), intent(in) :: t_
            real(dp), intent(in) :: y_(:)
            real(dp), intent(out) :: dydt_(:)
            class(*), intent(in), optional :: ctx_
            associate (dummy => ctx_)
            end associate
            call ts(v, eta, size(y_), t_, y_, dydt_)
        end subroutine bounce_int_rhs

        function bounce_event(t_, y_, ctx_) result(g)
            real(dp), intent(in) :: t_
            real(dp), intent(in) :: y_(:)
            class(*), intent(in), optional :: ctx_
            real(dp) :: g
            associate (dummy_t => t_, dummy_c => ctx_)
            end associate
            if (passing) then
                ! Full 2*pi turn in the sign_vpar_htheta direction.
                g = (y_(1) - th0) - sign_vpar_htheta * 2.0_dp * pi
            else
                ! Return to the starting poloidal angle.
                g = y_(1) - th0
            end if
        end function bounce_event
    end function bounce_integral

    subroutine timestep(v, eta, neq, t, y, ydot)
        !
        !  Timestep function for orbit integration.
        !  Includes poloidal angle theta and parallel velocity.
        !  More integrands may be added starting from y(3)
        !

        real(dp), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(dp), intent(in) :: t
        real(dp), intent(in) :: y(neq)
        real(dp), intent(out) :: ydot(neq)

        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: Om_tB_v
        real(dp) :: shearterm

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        shearterm = Bphcov * dqds
        if (noshear) then
            shearterm = 0
        end if

        Om_tB_v = mi * c * q / (2.0_dp * qi * sign_theta * psi_pr * bmod) * ( &  ! Om_tB/v**2
                  -(2.0_dp - eta * bmod) * bmod * hder(1) &
                  + 2.0_dp * (1.0_dp - eta * bmod) * hctrvr(3) * &
                  (dBthcovds + q * dBphcovds + shearterm))

        ydot(1) = y(2) * hctrvr(3)  ! theta
        ydot(2) = -0.5_dp * v**2 * eta * hctrvr(3) * hder(3) * bmod  ! v_par
        ydot(3) = Om_tB_v  ! for bounce average of Om_tB/v**2
        ydot(4:) = 0.0_dp  ! remaining integrands not computed here
    end subroutine timestep

end module neort_orbit
