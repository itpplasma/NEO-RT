module neort_orbit
    use logger, only: trace, get_log_level, LOG_TRACE, error
    use util, only: imun, pi, mi, qi, c
    use spline, only: spline_coeff, spline_val_0
    use do_magfie_mod, only: do_magfie, s, iota, R0, eps, psi_pr, &
        bphcov, dbthcovds, dbphcovds, q, dqds, sign_theta
    use do_magfie_pert_mod, only: do_magfie_pert_amp
    use neort_profiles, only: vth, Om_tE, dOm_tEds
    use driftorbit, only: etatp, etadt, etamin, etamax, epsmn, mth, mph, m0, mth, &
        init_done, pertfile, magdrift, nonlin, epsst_spl, epssp_spl, epst_spl, epsp_spl, &
        sign_vpar, sign_vpar_htheta

    implicit none

    integer, parameter :: nvar = 7
    real(8) :: th0

    logical :: noshear = .false.      ! neglect magnetic shear


    interface
        subroutine timestep_i(v, eta, neq, t, y, ydot)
            real(8), intent(in) :: v, eta
            integer, intent(in) :: neq
            real(8), intent(in) :: t
            real(8), intent(in) :: y(neq)
            real(8), intent(out) :: ydot(neq)
        end subroutine timestep_i
    end interface

contains

    subroutine dvode_error_context(where, v_in, eta_in, tcur, tout, ist)
        use do_magfie_mod, only: s, iota, R0, q, psi_pr, eps
        use driftorbit, only: etatp, etadt, etamin, etamax, mth, mph, sign_vpar
        use neort_profiles, only: vth, Om_tE
        character(*), intent(in) :: where
        real(8), intent(in) :: v_in, eta_in, tcur, tout
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
        write(0,'(A,1X,I0,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') '  mth=', mth, 'mph=', mph, 'sign_vpar=', dble(sign_vpar)
        write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  eps=', eps, 'psi_pr=', psi_pr
        call error('DVODE MXSTEP')
    end subroutine dvode_error_context

    pure function to_es(x) result(sout)
        real(8), intent(in) :: x
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
        real(8), intent(in) :: v, eta
        real(8), intent(out) :: taub, bounceavg(nvar)
        real(8), optional :: taub_estimate  ! estimated bounce time (user input)
        real(8) :: findroot_res(nvar + 1)
        real(8) :: bmod, htheta
        real(8) :: y0(nvar)

        ! Initialize bounce-averated quantities y0. Their meaning
        ! is defined inside subroutine timestep (thin orbit integration)
        call evaluate_bfield_local(bmod, htheta)
        sign_vpar_htheta = sign(1d0, htheta)*sign_vpar
        y0 = 1d-15
        y0(1) = th0         ! poloidal angle theta
        y0(2) = sign_vpar_htheta*vpar(v, eta, bmod)  ! parallel velocity vpar
        y0(3) = 0d0         ! toroidal velocity v_ph for drift frequency Om_ph
        y0(4) = 0d0         ! perturbed Hamiltonian real part
        y0(5) = 0d0         ! perturbed Hamiltonian imaginary part
        y0(6) = 0d0         ! 1/abs(B)
        ! y0(7) = 0d0       ! abs(B)

        ! If bounce time estimate exists (elliptic integrals),
        ! initialize taub with it, owtherwise estimate here.
        if (present(taub_estimate)) then
            taub = taub_estimate
        else
            taub = 2.0*pi/abs(vperp(v, eta, bmod)*iota/R0*sqrt(eps/2d0))
        end if

        ! Look for exactly one orbit turn via root-finding.
        ! Start by looking for 5 points per turn.
        findroot_res = bounce_integral(v, eta, nvar, y0, taub/5d0, timestep)

        taub = findroot_res(1)
        bounceavg = findroot_res(2:)/taub
    end subroutine bounce

    subroutine evaluate_bfield_local(bmod, htheta)
        real(8), intent(out) :: bmod, htheta
        real(8) :: sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0d0
        x(3) = th0
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        htheta = hctrvr(3)
    end subroutine evaluate_bfield_local

    pure function vpar(v, eta, bmod)
        !   parallel velocity
        real(8) :: vpar
        real(8), intent(in) :: v, eta, bmod
        vpar = v*sqrt(1d0 - eta*bmod)
        if (isnan(vpar)) then
            vpar = 0d0
        end if
    end function vpar

    pure function vperp(v, eta, bmod)
        !   perpendicular velocity
        real(8) :: vperp
        real(8), intent(in) :: v, eta, bmod
        vperp = v*sqrt(eta*bmod)
        if (isnan(vperp)) then
            vperp = 0d0
        end if
    end function vperp

    subroutine bounce_fast(v, eta, taub, bounceavg, ts, istate_out)
        use dvode_f90_m

        real(8), intent(in) :: v, eta, taub
        real(8), intent(out) :: bounceavg(nvar)
        procedure(timestep_i) :: ts
        integer, intent(out), optional :: istate_out

        real(8) :: t1, t2, bmod, htheta
        real(8) :: y(nvar)
        real(8) :: atol(nvar), rtol
        integer :: neq, itask, istate
        type(vode_opts) :: options

        t1 = 0d0
        t2 = taub

        call evaluate_bfield_local(bmod, htheta)
        sign_vpar_htheta = sign(1d0, htheta)*sign_vpar
        y = 1d-15
        y(1) = th0
        y(2) = sign_vpar_htheta*vpar(v, eta, bmod)
        y(3:6) = 0d0

        neq = nvar
        rtol = 1d-9
        atol = 1d-10
        itask = 1
        istate = 1
        options = set_normal_opts(abserr_vector=atol, relerr=rtol)

        call dvode_f90(timestep_wrapper, neq, y, t1, t2, itask, istate, options)
        if (istate == -1) then
            call dvode_error_context('bounce_fast', v, eta, t1, t2, istate)
        end if

        bounceavg = y/taub
        if (present(istate_out)) istate_out = istate

    contains

        subroutine timestep_wrapper(neq_, t_, y_, ydot_)
            ! Wrapper routine for timestep to work with VODE
            integer, intent(in) :: neq_
            real(8), intent(in) :: t_
            real(8), intent(in) :: y_(neq_)
            real(8), intent(out) :: ydot_(neq_)

            call ts(v, eta, neq_, t_, y_, ydot_)
        end subroutine timestep_wrapper
    end subroutine bounce_fast

    function bounce_time(v, eta, taub_estimate) result(taub)
        use dvode_f90_m, only: dvode_f90, set_normal_opts, vode_opts

        real(8), intent(in) :: v, eta
        real(8), intent(in), optional :: taub_estimate
        real(8) :: taub

        integer, parameter :: neq = 2
        real(8) :: y0(neq), roots(neq+1)
        real(8) :: bmod, htheta

        call evaluate_bfield_local(bmod, htheta)
        sign_vpar_htheta = sign(1d0, htheta)*sign_vpar

        y0(1) = th0         ! poloidal angle theta
        y0(2) = sign_vpar_htheta*vpar(v, eta, bmod)  ! parallel velocity vpar

        if (present(taub_estimate)) then
            taub = taub_estimate
        else
            taub = 2.0*pi/abs(vperp(v, eta, bmod)*iota/R0*sqrt(eps/2d0))
        end if

        roots = bounce_integral(v, eta, neq, y0, taub, timestep_poloidal_motion)
        taub = roots(1)

    end function bounce_time

    subroutine timestep_poloidal_motion(v, eta, neq, t, y, ydot)
        real(8), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(8), intent(in) :: t, y(neq)
        real(8), intent(out) :: ydot(neq)

        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)

        x(1) = s
        x(2) = 0d0
        x(3) = y(1)

        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        call poloidal_velocity(v, eta, bmod, hctrvr(3), hder(3), y(2), ydot)
    end subroutine timestep_poloidal_motion

    pure subroutine poloidal_velocity(v, eta, bmod, hthctr, hderth, v_par, ydot)
        real(8), intent(in) :: v, eta, bmod, hthctr, hderth
        real(8), intent(in) :: v_par
        real(8), intent(out) :: ydot(2)

        ydot(1) = v_par*hthctr                                  ! theta
        ydot(2) = -v**2*eta/2d0*hthctr*hderth*bmod              ! v_par
    end subroutine poloidal_velocity

    function bounce_integral(v, eta, neq, y0, dt, ts)
        !
        !  Finds the root of an orbit after the first turn
        !
        use dvode_f90_m

        real(8) :: bounce_integral(neq + 1)
        real(8), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(8), intent(in) :: y0(neq), dt
        procedure(timestep_i) :: ts

        integer :: n

        integer :: k, state, rootstate
        real(8) :: ti, told
        real(8) :: y(neq), yold(neq)

        logical :: passing

        real(8) :: atol(neq), rtol, tout
        integer :: itask, istate
        type(vode_opts) :: options

        rtol = 1d-9
        atol = 1d-10
        itask = 1
        istate = 1
        options = set_normal_opts(abserr_vector=atol, relerr=rtol, nevents=2)

        ! check for passing orbit
        passing = .false.
        if (eta < etatp) then
            passing = .true.
        end if

        n = 500
        rootstate = -1

        y = y0
        yold = y0
        ti = 0d0
        state = 1
        if (get_log_level() >= LOG_TRACE) then
            write(*,'(A,1X,ES12.5,1X,ES12.5,1X,A,1X,A)') &
                '[TRACE] bounce_integral start v,eta=', v, eta, 'pass=', merge('T','F',eta<etatp)
        end if
        do k = 2, n
            yold = y
            told = ti

            tout = ti + dt
            call dvode_f90(timestep_wrapper, neq, y, ti, tout, itask, istate, options, &
                        g_fcn=bounceroots)
            if (istate == -1) then
                call dvode_error_context('bounce_integral', v, eta, ti, tout, istate)
            end if
            if (get_log_level() >= LOG_TRACE) then
                write(*,'(A,I0,A,ES12.5,A,ES12.5,A,I0)') &
                    '[TRACE] step k=', k, ' ti=', ti, ' y1=', y(1), ' istate=', istate
            end if
            if (istate == 3) then
                if (passing .or. (yold(1) - th0) < 0) then
                    exit
                end if
            end if

            istate = 2
        end do
        if (istate /= 3) then
            write(0,'(A)') '[ERROR] bounce_integral: no event after 500 iterations'
            write(0,'(A,1X,A)') '  region =', merge('passing','trapped',eta<etatp)
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  v =', v, 'eta =', eta
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  ti =', ti, 'dt =', dt
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') '  etamin =', etamin, 'etamax =', etamax, 'etatp =', etatp
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5)') '  theta(y1) =', y(1), 'th0 =', th0
            write(0,'(A,1X,I0,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') '  mth =', mth, 'mph =', mph, 'sign_vpar =', dble(sign_vpar)
            write(0,'(A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5,2X,A,1X,ES12.5)') '  s =', s, 'R0 =', R0, 'q =', q, 'iota =', iota
        end if

        bounce_integral(1) = ti
        bounce_integral(2:) = y

    contains

        subroutine timestep_wrapper(neq_, t_, y_, ydot_)
            ! Wrapper routine for timestep to work with VODE
            integer, intent(in) :: neq_
            real(8), intent(in) :: t_
            real(8), intent(in) :: y_(neq_)
            real(8), intent(out) :: ydot_(neq_)

            call ts(v, eta, neq_, t_, y_, ydot_)
        end subroutine timestep_wrapper
    end function bounce_integral

    subroutine bounceroots(NEQ, T, Y, NG, GOUT)
        integer, intent(in) :: NEQ, NG
        real(8), intent(in) :: T, Y(neq)
        real(8), intent(out) :: GOUT(ng)
        associate (dummy => T)
        end associate
        GOUT(1) = sign_vpar_htheta*(Y(1) - th0) ! trapped orbit return to starting point
        GOUT(2) = sign_vpar_htheta*(2d0*pi - (Y(1) - th0))  ! passing orbit return
        return
    end subroutine bounceroots

    subroutine timestep(v, eta, neq, t, y, ydot)
        !
        !  Timestep function for orbit integration.
        !  Includes poloidal angle theta and parallel velocity.
        !  More integrands may be added starting from y(3)
        !

        real(8), intent(in) :: v, eta
        integer, intent(in) :: neq
        real(8), intent(in) :: t
        real(8), intent(in) :: y(neq)
        real(8), intent(out) :: ydot(neq)

        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(8) :: Om_tB_v
        real(8) :: shearterm


        x(1) = s
        x(2) = 0d0
        x(3) = y(1)
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)

        shearterm = Bphcov*dqds
        if (noshear) then
            shearterm = 0
        end if

        Om_tB_v = mi*c*q/(2d0*qi*sign_theta*psi_pr*bmod)*( &      ! Om_tB/v**2
                  -(2d0 - eta*bmod)*bmod*hder(1) &
                  + 2d0*(1d0 - eta*bmod)*hctrvr(3)* &
                  (dBthcovds + q*dBphcovds + shearterm))

        ydot(1) = y(2)*hctrvr(3)                                    ! theta
        ydot(2) = -0.5d0*v**2*eta*hctrvr(3)*hder(3)*bmod            ! v_par
        ydot(3) = Om_tB_v  ! for bounce average of Om_tB/v**2
    end subroutine timestep

end module neort_orbit
