module neort_orbit
    use iso_fortran_env, only: dp => real64
    use logger, only: debug, trace, get_log_level, LOG_TRACE, error
    use util, only: imun, pi, mi, qi, c
    use spline, only: spline_coeff, spline_val_0
    use do_magfie_mod, only: do_magfie, s, iota, R0, eps, psi_pr, inp_swi, &
        bphcov, dbthcovds, dbphcovds, q, dqds, sign_theta
    use do_magfie_pert_mod, only: do_magfie_pert_amp
    use neort_profiles, only: vth, Om_tE, dOm_tEds
    use driftorbit, only: etatp, etadt, etamin, etamax, epsmn, mth, mph, m0, mth, &
        pertfile, magdrift, nonlin, epsst_spl, epssp_spl, epst_spl, epsp_spl, &
        sign_vpar, sign_vpar_htheta

    implicit none

    integer, parameter :: nvar = 7
    real(dp) :: th0 = 0.0_dp

    logical :: noshear = .false. ! neglect magnetic shear

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
        real(dp), optional :: taub_estimate ! estimated bounce time (user input)
        real(dp) :: findroot_res(nvar + 1)
        real(dp) :: bmod, htheta
        real(dp) :: y0(nvar)

        ! Initialize bounce-averated quantities y0. Their meaning
        ! is defined inside subroutine timestep (thin orbit integration)
        call evaluate_bfield_local(bmod, htheta)
        sign_vpar_htheta = sign(1.0_dp, htheta)*sign_vpar
        y0 = 1.0e-15_dp
        y0(1) = th0 ! poloidal angle theta
        y0(2) = sign_vpar_htheta * vpar(v, eta, bmod) ! parallel velocity vpar
        y0(3) = 0.0_dp ! toroidal velocity v_ph for drift frequency Om_ph
        y0(4) = 0.0_dp ! perturbed Hamiltonian real part
        y0(5) = 0.0_dp ! perturbed Hamiltonian imaginary part
        y0(6) = 0.0_dp ! 1/abs(B)
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
        use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
        use fortnum_status, only: fortnum_status_t, FORTNUM_OK, &
            FORTNUM_CONVERGENCE_ERROR

        real(dp), intent(in) :: v, eta, taub
        real(dp), intent(out) :: bounceavg(nvar)
        procedure(timestep_i) :: ts
        integer, intent(out), optional :: istate_out

        ! Nonstiff variable-order Adams (fortnum vode, DVODE MF=10) integrated
        ! over the single bounce span [0, taub], relative tolerance 1e-9 and a
        ! per-component absolute tolerance 1e-10 (DVODE ITOL=2). This is the
        ! same method NEO-RT drove through DVODE before the migration, so the
        ! variable-order controller resolves the oscillatory bounceavg(3:4)
        ! Hamiltonian integrands without an artificial step cap.
        real(dp), parameter :: rtol = 1.0e-9_dp
        real(dp), parameter :: atol_val = 1.0e-10_dp

        real(dp) :: t1, t2, bmod, htheta
        real(dp) :: y0(nvar), atol(nvar)
        real(dp), allocatable :: yend(:)
        type(vode_state_t) :: vstate
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

        atol = atol_val
        call vode_init(vstate, nvar, t1, y0)
        call vode_integrate_to(bounce_rhs, vstate, t2, rtol, atol, yend, status)

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

        bounceavg = yend / taub
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

        y0(1) = th0 ! poloidal angle theta
        y0(2) = sign_vpar_htheta * vpar(v, eta, bmod) ! parallel velocity vpar

        if (present(taub_estimate)) then
            taub = taub_estimate
        else
            taub = 2.0 * pi / abs(vperp(v, eta, bmod) * iota / R0 * sqrt(eps / 2.0_dp))
        end if

        roots = bounce_integral(v, eta, neq, y0, taub / 5.0_dp, timestep_poloidal_motion)
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

        ydot(1) = v_par * hthctr ! theta
        ydot(2) = -v**2 * eta / 2.0_dp * hthctr * hderth * bmod ! v_par
    end subroutine poloidal_velocity

    function bounce_integral(v, eta, neq, y0, dt, ts)
        !
        !  Finds the root of an orbit after the first turn
        !
        use fortnum_ode, only: ODE_EVENT_ANY
        use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
        use fortnum_status, only: fortnum_status_t, FORTNUM_CONVERGENCE_ERROR

        real(dp) :: bounce_integral(neq + 1)
        real(dp), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(dp), intent(in) :: y0(neq), dt
        procedure(timestep_i) :: ts

        ! Number of dt-sized chunks the search window spans, matching the old
        ! chunked DVODE search that advanced in dt steps up to this many turns.
        integer, parameter :: n_turns = 500
        real(dp), parameter :: rtol = 1.0e-9_dp
        real(dp), parameter :: atol_val = 1.0e-10_dp

        type(vode_state_t) :: vstate
        type(fortnum_status_t) :: status
        real(dp), allocatable :: y_out(:)
        real(dp) :: atol(neq), t_now, t_root, theta_before
        logical :: passing, found
        integer :: chunk

        passing = (eta < etatp)
        atol = atol_val

        if (get_log_level() >= LOG_TRACE) then
            write(*,'(A,2ES12.5,2A)') '[TRACE] bounce_integral start v,eta=', v, eta, &
                ' pass=', merge('T','F', passing)
        end if

        ! Reproduce the DVODE bounceroots search: two event functions monitored
        ! with NEVENTS=2, advancing in dt-sized chunks, stop at the first root
        ! of either that satisfies the turn acceptance test. fortnum vode
        ! locates the root on its own Nordsieck interpolant (relerr 1e-9,
        ! per-component abserr 1e-10, ITOL=2), so taub is the located root.
        !   g1 = theta - th0          (trapped: return to th0)
        !   g2 = 2*pi - (theta - th0) (passing: full +2*pi turn)
        ! DVODE accepted a root when passing, or when theta entered from below
        ! th0 (the old (yold(1)-th0) < 0 filter); otherwise it kept integrating.
        call vode_init(vstate, neq, 0.0_dp, y0)
        t_now = 0.0_dp
        theta_before = y0(1)
        found = .false.

        do chunk = 1, n_turns - 1
            ! Step one dt-sized window. The integration is never re-initialised:
            ! after a turn root vode stops at the root and continues from there
            ! on the next call, so the next window ends dt past wherever this
            ! one stopped. theta_before is the poloidal angle at the end of the
            ! previous window (DVODE's yold), used by the acceptance filter.
            call vode_integrate_to(bounce_int_rhs, vstate, t_now + dt, &
                rtol, atol, y_out, status, &
                event=root_theta, event_dir=ODE_EVENT_ANY, &
                event2=root_turn, event_dir2=ODE_EVENT_ANY, &
                t_root=t_root, root_found=found)
            if (status%code == FORTNUM_CONVERGENCE_ERROR) then
                call dvode_error_context('bounce_integral', v, eta, &
                    t_now, t_now + dt, -1)
            end if

            if (found) then
                if (passing .or. (theta_before - th0) < 0.0_dp) exit
                ! Reject this turning point: keep integrating, the next window
                ! ends dt past the located root (DVODE istate=2 continuation).
                found = .false.
                theta_before = y_out(1)
                t_now = t_root
                cycle
            end if

            theta_before = y_out(1)
            t_now = t_now + dt
        end do

        if (.not. found) then
            write(0,'(A)') '[ERROR] bounce_integral: no bounce event located'
            write(0,'(A,1X,A)') '  region =', merge('passing','trapped', passing)
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  v =', v, 'eta =', eta
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  t =', t_now, 'dt =', dt
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') '  etamin =', etamin, 'etamax =', etamax, 'etatp =', etatp
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  theta(y1) =', vstate%yh(1, 1), 'th0 =', th0
            write(0,'(A,1X,I0,2X,A,1X,I0,2X,A,1X,ES12.5)') '  mth =', mth, 'mph =', mph, 'sign_vpar =', dble(sign_vpar)
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') '  s =', s, 'R0 =', R0, 'q =', q, 'iota =', iota
            bounce_integral(1) = t_now
            bounce_integral(2:) = vstate%yh(:, 1)
            return
        end if

        bounce_integral(1) = t_root
        bounce_integral(2:) = y_out

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

        ! DVODE bounceroots GOUT(1): theta returns to th0 (trapped turn).
        function root_theta(t_, y_, ctx_) result(g)
            real(dp), intent(in) :: t_
            real(dp), intent(in) :: y_(:)
            class(*), intent(in), optional :: ctx_
            real(dp) :: g
            associate (dummy_t => t_, dummy_c => ctx_)
            end associate
            g = y_(1) - th0
        end function root_theta

        ! DVODE bounceroots GOUT(2): theta advances by 2*pi (passing turn).
        function root_turn(t_, y_, ctx_) result(g)
            real(dp), intent(in) :: t_
            real(dp), intent(in) :: y_(:)
            class(*), intent(in), optional :: ctx_
            real(dp) :: g
            associate (dummy_t => t_, dummy_c => ctx_)
            end associate
            g = 2.0_dp * pi - (y_(1) - th0)
        end function root_turn
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
        real(dp) :: Om_tB_v, cross_gradB_phi, cross_gradB_theta
        real(dp) :: curl_parallel, curvature_phi, curvature_theta
        real(dp) :: drift_phi_v, drift_theta_v
        real(dp) :: shearterm

        x(1) = s
        x(2) = 0.0_dp
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        if (inp_swi == 11) then
            ! Coordinate-independent first-order magnetic drift, evaluated from
            ! the direct cylindrical EQDSK field and transformed to
            ! (s_tor,phi,theta_geo).  hcurl is curl(h), and the parallel piece
            ! is removed to obtain h cross (h dot grad h).
            cross_gradB_phi = (hcovar(3)*hder(1) - hcovar(1)*hder(3))/sqrtg
            cross_gradB_theta = (hcovar(1)*hder(2) - hcovar(2)*hder(1))/sqrtg
            curl_parallel = sum(hcovar*hcurl)
            curvature_phi = hcurl(2) - hctrvr(2)*curl_parallel
            curvature_theta = hcurl(3) - hctrvr(3)*curl_parallel
            drift_phi_v = mi*c/(qi*bmod)*((1.0_dp - eta*bmod)*curvature_phi &
                + 0.5_dp*eta*bmod*cross_gradB_phi)
            drift_theta_v = mi*c/(qi*bmod)*((1.0_dp - eta*bmod)*curvature_theta &
                + 0.5_dp*eta*bmod*cross_gradB_theta)
            ! Canonical toroidal precession follows the field-line label
            ! alpha=phi-q*theta, not the small cylindrical phi component alone.
            Om_tB_v = drift_phi_v - q*drift_theta_v
        else
            shearterm = Bphcov * dqds
            if (noshear) then
                shearterm = 0
            end if

            Om_tB_v = mi * c * q / (2.0_dp * qi * sign_theta * psi_pr * bmod) * ( & ! Om_tB/v**2
                -(2.0_dp - eta * bmod) * bmod * hder(1) &
                + 2.0_dp * (1.0_dp - eta * bmod) * hctrvr(3) * &
                (dBthcovds + q * dBphcovds + shearterm))
        end if

        ydot(1) = y(2) * hctrvr(3) ! theta
        ydot(2) = -0.5_dp * v**2 * eta * hctrvr(3) * hder(3) * bmod ! v_par
        ydot(3) = Om_tB_v ! for bounce average of Om_tB/v**2
        ydot(4:) = 0.0_dp ! remaining integrands not computed here
        if (inp_swi == 11 .and. neq >= 7) ydot(7) = y(2)*hctrvr(2)
    end subroutine timestep

end module neort_orbit
