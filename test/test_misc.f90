program test_misc
    use driftorbit
    use neort, only: main, init_plasma, runmode, runname

    implicit none

    call main

    if (runmode == "test_profile") then
        call test_profile
        stop
    elseif (runmode == "test_bounce") then
        call test_bounce
        stop
    elseif (runmode == "test_torfreq") then
        call test_torfreq
        stop
    elseif (runmode == "test_resline") then
        call test_resline
        stop
    elseif (runmode == "test_box") then
        if (orbit_mode_transp <= 0) &
            error stop "need orbit_mode_transp>=0 for test_box"
        call test_box
        stop
    elseif (runmode == "test_torque_integral") then
        call test_torque_integral
        stop
    elseif (runmode == "test_integral") then
        call test_integral
        stop
    end if

contains

    subroutine test_bounce
        use dvode_f90_m2
        use orbit

        real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(8) :: dt
        integer, parameter :: n = 1000

        integer :: k
        real(8) :: t, y0(4), y(4)

        real(8) :: atol(4), rtol, tout
        integer :: neq, itask, istate
        type(vode_opts) :: options

        neq = 4
        t = 0.0d0
        rtol = 1d-11
        atol = 1d-12
        options = set_normal_opts(abserr_vector=atol, relerr=rtol)

        print *, "test_bounce"

        x(1) = s
        x(2) = 0d0
        x(3) = 0d0
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

        !v = 143730130.234083
        !eta = 4.686503826199121E-005

        v = vth
        !eta = (1d0-1.86727d-1)/bmod
        !eta = etatp*(1+1d-2)
        eta = etatp + (etadt - etatp)/1d1
        !eta = etatp*1d-2
        !eta = etatp*0.8
        !eta = etatp*1.1
        !eta = etatp*(1+1d-7)
        !eta = etadt*(1-1d-2)

        !v = 2.2*vth
        !eta = 0.00003

        call bounce
        call disp("test_bounce: Om_tB     = ", v**2*bounceavg(3))
        call disp("test_bounce: taub      = ", taub)
        call disp("test_bounce: taub_est  = ", 2.0*pi/(vperp(bmod)*iota/R0*sqrt(eps/2d0)))
        call disp("test_bounce: bmod      = ", bmod)

        dt = taub/(n - 1)

        ! go until banana tip or inboard side
        t = 0d0
        y(1:3) = x
        y(4) = -vpar(bmod)
        itask = 1
        istate = 1
        if (eta > etatp) then
            tout = -taub/4.0
        else
            tout = -taub/2.0
        end if
        call dvode_f90(step_zeroorder, neq, y, t, tout, itask, istate, options)
        y0 = y

        t = 0d0
        y = y0
        y(2) = 0.0 ! reset toroidal angle
        itask = 1
        istate = 1
        open (unit=9, file=trim(adjustl(runname))//"_bounce_zeroorder.out", recl=1024)
        do k = 2, n
            call do_magfie(y(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            write (9, *) t/taub, y, q, taub, v**2*bounceavg(3), y(1) + c*mi*y(4)*hcovar(2)*q/(qi*psi_pr)
            tout = t + dt
            call dvode_f90(step_zeroorder, neq, y, t, tout, itask, istate, options)
        end do
        call do_magfie(y(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        write (9, *) t/taub, y, q, taub, v**2*bounceavg(3), y(1) + c*mi*y(4)*hcovar(2)*q/(qi*psi_pr)
        close (unit=9)

        t = 0d0
        y = y0
        y(2) = 0.0 ! reset toroidal angle
        y(1) = y(1) + c*mi*y0(4)*hcovar(2)*q/(qi*psi_pr)
        itask = 1
        istate = 1
        open (unit=9, file=trim(adjustl(runname))//"_bounce_full.out", recl=1024)
        do k = 2, n
            call do_magfie(y(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            write (9, *) t/taub, y, q, taub, v**2*bounceavg(3)
            tout = t + dt
            call dvode_f90(step_full, neq, y, t, tout, itask, istate, options)
        end do
        call do_magfie(y(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        write (9, *) t/taub, y, q, taub, v**2*bounceavg(3)
        close (unit=9)

        t = 0d0
        y = y0
        y(2) = 0.0 ! reset toroidal angle
        y(1) = y(1) + c*mi*y0(4)*hcovar(2)*q/(qi*psi_pr)
        itask = 1
        istate = 1
        open (unit=9, file=trim(adjustl(runname))//"_bounce_rela.out", recl=1024)
        do k = 2, n
            call do_magfie(y(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
            write (9, *) t/taub, y, q, taub, v**2*bounceavg(3)
            tout = t + dt
            call dvode_f90(step_rela, neq, y, t, tout, itask, istate, options)
        end do
        call do_magfie(y(1:3), bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        write (9, *) t/taub, y, q, taub, v**2*bounceavg(3)
        close (unit=9)
    end subroutine test_bounce

    subroutine test_box
        ! test box counting
        use dvode_f90_m2
        integer :: k, n

        real(8) :: ti
        real(8) :: y(2), yold(2)

        real(8) :: atol(nvar), rtol, tout, rstats(22)
        integer :: neq, itask, istate, istats(31), numevents
        type(vode_opts) :: options

        real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(8) :: s1old, told
        integer :: sind, sind0 ! s index

        integer :: jroots(2)

        v = 26399452.5418568
        eta = 4.686498216380098e-005

        ! v = vth
        ! eta = etatp*1.01

        x(1) = s
        x(2) = 0d0
        x(3) = 0d0
        call do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)

        call bounce

        neq = 2
        rtol = 1d-12
        atol = 1d-13
        itask = 1
        istate = 1
        numevents = 2
        options = set_normal_opts(abserr_vector=atol, relerr=rtol, nevents=numevents)

        n = 3*size(sbox)
        allocate (taubins(size(sbox) + 1))
        taubins = 0d0

        y(1) = 0d0
        y(2) = sigv*vpar(bmod)
        ti = 0d0

        print *, sbox

        s1 = s + c*mi*y(2)*hcovar(2)*q/(qi*psi_pr)
        sind = size(sbox) + 1
        sprev = -1e5
        snext = 1e5
        do k = 1, size(sbox)
            if (sbox(k) > s1) then
                sind = k
                snext = sbox(k)
                if (k > 1) sprev = sbox(k - 1)
                exit
            end if
        end do
        sind0 = sind

        print *, y(1), sind, s1

        told = 0d0
        do k = 2, n
            yold = y
            s1old = s1
            tout = taub
            call dvode_f90(tsorb, neq, y, ti, tout, itask, istate, options, &
                           g_fcn=sroots)
            if (istate == 2) exit
            if (istate == 3) then
                taubins(sind) = taubins(sind) + ti - told
                told = ti
                call get_stats(rstats, istats, numevents, jroots)
                if (jroots(2) .ne. 0) then
                    sind = sind + 1 ! moving outwards
                    sprev = snext
                    if (sind == size(sbox) + 1) then
                        snext = 1e5
                    else
                        snext = sbox(sind)
                    end if
                end if
                if (jroots(1) .ne. 0) then
                    sind = sind - 1 ! moving inwards
                    snext = sprev
                    if (sind == 1) then
                        sprev = -1e5
                    else
                        sprev = sbox(sind - 1)
                    end if
                end if
                print *, y(1), sind, s1
            end if
        end do
        !if (sind /= sind0) then
        !   write(6,*) s, s1, sind0, sind, v, eta, etatp, etadt
        !   stop("START AND END RADIUS OF ORBIT ARE NOT THE SAME")
        !end if
        print *, ti/taub, sind, s1
        taubins(sind) = taubins(sind) + taub - told

        taubins = taubins/taub
        print *, taubins
        deallocate (taubins)

    end subroutine test_box

    subroutine test_torfreq
        integer, parameter :: n = 100
        integer :: k
        real(8) :: Omph, dOmphdv, dOmphdeta, dOmphds
        real(8) :: Omth, dOmthdv, dOmthdeta, dOmthds
        real(8) :: OmtB, dOmtBdv, dOmtBdeta
        real(8) :: aa, b

        v = vth

        etamin = etatp
        etamax = etatp + (etadt - etatp)*(1d0 - epsst_spl)

        b = log(epst_spl)
        aa = 1d0/(n - 1d0)*(log(etamax/etamin - 1d0) - b)
        eta = etamax

        sigv = 1

        call disp("test_torfreq: vth        = ", vth)
        call disp("test_torfreq: v/vth      = ", v/vth)
        call disp("test_torfreq: mph        = ", 1d0*mph)
        call disp("test_torfreq: mth        = ", 1d0*mth)
        call disp("test_torfreq: Om_tE      = ", Om_tE)
        call disp("test_torfreq: Om_tB_ref  = ", c*mi*vth**2/(2*qi*psi_pr))
        call bounce
        call disp("test_torfreq: Om_tB_ba   = ", vth**2*bounceavg(3))
        call disp("test_torfreq: etamin = ", etamin)
        call disp("test_torfreq: etamax = ", etamax)
        call Om_th(Omth, dOmthdv, dOmthdeta)

        call disp("test_torfreq: Om_th_approx    = ", v/(q*R0*sqrt(2d0/eps)))
        call disp("test_torfreq: Om_th_deeptrap  = ", Omth)

        open (unit=9, file=trim(adjustl(runname))//"_torfreq.out", recl=1024)
        write (9, *) "!1:eta                   "// &
            " 2:etatp                  "// &
            " 3:etadt                  "// &
            " 4:Om_tE                  "// &
            " 5:OmtB                   "// &
            " 6:dOmtbdv                "// &
            " 7:dOmtbdeta              "// &
            " 8:Omth                   "// &
            " 9:dOmthdv                "// &
            "10:dOmthdeta              "// &
            "11:dOmthds                "// &
            "12:Omph                   "// &
            "13:dOmphdv                "// &
            "14:dOmphdeta              "// &
            "15:dOmphds                "
        do k = 0, n - 1
            eta = etamin*(1d0 + exp(aa*k + b))
            call Om_ph(Omph, dOmphdv, dOmphdeta)
            call Om_th(Omth, dOmthdv, dOmthdeta)
            call Om_tB(OmtB, dOmtBdv, dOmtBdeta)
            call d_Om_ds(dOmthds, dOmphds)
            write (9, *) eta, etatp, etadt, &
                Om_tE, OmtB, dOmtbdv, dOmtbdeta, &
                Omth, dOmthdv, dOmthdeta, dOmthds, &
                Omph, dOmphdv, dOmphdeta, dOmphds
        end do
        close (unit=9)
    end subroutine test_torfreq

    subroutine test_resline
        integer, parameter :: n = 500
        integer :: k
        real(8) :: vmin, vmax
        real(8) :: etarest(2), etaresp(2)
        real(8) :: roots(nlev, 3)
        integer :: nroots, kr

        print *, "test_resline"

        vmin = 1d-6*vth
        vmax = 10d0*vth

        etaresp = etatp
        etarest = etatp
        sigv = 1

        open (unit=9, file=trim(adjustl(runname))//"_resline_pco.out", recl=1024)
        open (unit=10, file=trim(adjustl(runname))//"_resline_pct.out", recl=1024)
        open (unit=11, file=trim(adjustl(runname))//"_resline_t.out", recl=1024)
        do k = 0, n - 1
            v = vmin + k/(n - 1d0)*(vmax - vmin)

            if (.not. nopassing) then
                ! resonance (passing)
                sigv = 1
                call driftorbit_coarse(etatp*epsp, etatp*(1 - epsp), roots, nroots)
                do kr = 1, nroots
                    etaresp = driftorbit_root(1d-8*abs(Om_tE), roots(kr, 1), roots(kr, 2))
                    write (9, *) v/vth, kr, etaresp(1), 0d0, etatp
                end do
                sigv = -1
                call driftorbit_coarse(etatp*epsp, etatp*(1 - epsp), roots, nroots)
                do kr = 1, nroots
                    etaresp = driftorbit_root(1d-8*abs(Om_tE), roots(kr, 1), roots(kr, 2))
                    write (10, *) v/vth, kr, etaresp(1), 0d0, etatp
                end do
            end if

            ! resonance (trapped)
            sigv = 1
            call driftorbit_coarse(etatp*(1 + epst), etadt*(1 - epst), roots, nroots)
            do kr = 1, nroots
                etarest = driftorbit_root(1d-8*abs(Om_tE), roots(kr, 1), roots(kr, 2))
                write (11, *) v/vth, kr, etarest(1), etatp, etadt
            end do
        end do
        close (unit=11)
        close (unit=10)
        close (unit=9)
    end subroutine test_resline

    subroutine test_integral
        real(8) :: eta_res(2)
        real(8) :: D11, D12, Dp, vold, dsdreff, D11sum, D12sum
        integer :: ku, nu
        real(8) :: xs, kappa2s
        real(8) :: roots(nlev, 3)
        integer :: nroots, kr
        character(1024) :: fn1, fn2, fn3
        real(8) :: tol0(2)

        tol0 = 1d-9

        call disp("test_integral: Mach num Mt         = ", M_t)
        call disp("test_integral: Poloidal mode mth   = ", 1d0*mth)
        print *, "test_integral: Boxes sbox = ", sbox

        write (fn1, *) trim(adjustl(runname))//"_integral_pco.out"
        write (fn2, *) trim(adjustl(runname))//"_integral_pct.out"
        write (fn3, *) trim(adjustl(runname))//"_integral_t.out"

        nu = 1000

        vmin2 = 1d-6*vth
        vmax2 = 3d0*vth

        print *, "==CO-PASSING=="
        sigv = 1
        etamin = epsp*etatp
        etamax = (1 - epsp)*etatp

        print *, "vrange: ", vmin2/vth, vmax2/vth
        print *, "etarange: ", etamin, etamax
        !if (vsteps>0) print *, "D11/Dp (INT) = ", flux_integral(vmin2, vmax2, tol0)
        !if (vsteps>0) print *, "D11/Dp (MID) = ", flux_integral_mid(vmin2, vmax2)
        !if (orbit_mode_transp>0) print *, fluxint_box(1,:)
        !print *, "D11/Dp (ODE) = ", flux_integral_ode(vmin2, vmax2)
        !vmin2 = vmin2 + 1d-10*(vmax2-vmin2)
        !vmax2 = vmax2 - 1d-10*(vmax2-vmin2)

        D11 = 0d0
        D12 = 0d0
        D11sum = 0d0
        D12sum = 0d0
        dsdreff = 2d0/a*sqrt(s)
        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)

        open (unit=10, file=trim(adjustl(fn1)), recl=1024)
        v = vmax2
        vold = vmax2
        do ku = 0, nu - 1
            vold = v
            v = vmax2 - ku*(vmax2 - vmin2)/(nu - 1)

            call driftorbit_coarse(etamin, etamax, roots, nroots)
            if (nroots == 0) cycle

            do kr = 1, nroots
                eta_res = driftorbit_root(max(1d-9*abs(Om_tE), 1d-12), roots(kr, 1), roots(kr, 2))
                eta = eta_res(1)
                D11 = (vold - v)/vth*D11int_u(v/vth)*dsdreff**(-2)
                D12 = (vold - v)/vth*D12int_u(v/vth)*dsdreff**(-2)
                if ((ku == 1) .or. (ku == nu - 1)) then
                    D11 = D11/2d0
                    D12 = D12/2d0
                end if
                D11sum = D11sum + D11
                D12sum = D12sum + D12
            end do

            xs = (v/vth)**2
            kappa2s = (1d0/eta_res(1) - B0*(1d0 - eps))/(2d0*B0*eps)

            write (10, *) v/vth, eta_res(1), D11/Dp, D12/Dp, &
                (eta_res(1) - etamin)/(etamax - etamin), &
                vth, B0, xs, kappa2s, 0d0, etatp
        end do
        print *, "D11/Dp (SUM) = ", D11sum/Dp, D12sum/Dp
        close (unit=10)

        print *, "==CTR-PASSING=="
        sigv = -1
        etamin = epsp*etatp
        etamax = (1 - epsp)*etatp

        print *, "vrange: ", vmin2/vth, vmax2/vth
        print *, "etarange: ", etamin, etamax
        if (vsteps > 0) print *, "D11/Dp (INT) = ", flux_integral(vmin2, vmax2, tol0)
        if (vsteps > 0) print *, "D11/Dp (MID) = ", flux_integral_mid(vmin2, vmax2)
        if (orbit_mode_transp > 0) print *, fluxint_box(1, :)
        !print *, "D11/Dp (ODE) = ", flux_integral_ode(vmin2, vmax2)
        !vmin2 = vmin2 + 1d-10*(vmax2-vmin2)
        !vmax2 = vmax2 - 1d-10*(vmax2-vmin2)

        D11 = 0d0
        D12 = 0d0
        D11sum = 0d0
        D12sum = 0d0
        dsdreff = 2d0/a*sqrt(s)
        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)

        open (unit=10, file=trim(adjustl(fn2)), recl=1024)
        v = vmax2
        vold = vmax2
        do ku = 0, nu - 1
            vold = v
            v = vmax2 - ku*(vmax2 - vmin2)/(nu - 1)

            call driftorbit_coarse(etamin, etamax, roots, nroots)
            if (nroots == 0) cycle

            do kr = 1, nroots
                eta_res = driftorbit_root(max(1d-9*abs(Om_tE), 1d-12), roots(kr, 1), roots(kr, 2))
                eta = eta_res(1)
                D11 = (vold - v)/vth*D11int_u(v/vth)*dsdreff**(-2)
                D12 = (vold - v)/vth*D12int_u(v/vth)*dsdreff**(-2)
                if ((ku == 1) .or. (ku == nu - 1)) then
                    D11 = D11/2d0
                    D12 = D12/2d0
                end if
                D11sum = D11sum + D11
                D12sum = D12sum + D12
            end do

            xs = (v/vth)**2
            kappa2s = (1d0/eta_res(1) - B0*(1d0 - eps))/(2d0*B0*eps)

            write (10, *) v/vth, eta_res(1), D11/Dp, D12/Dp, &
                (eta_res(1) - etamin)/(etamax - etamin), &
                vth, B0, xs, kappa2s, 0d0, etatp
        end do
        print *, "D11/Dp (SUM) = ", D11sum/Dp, D12sum/Dp
        close (unit=10)

        print *, "==TRAPPED=="
        sigv = 1d0
        etamin = (1 + epst)*etatp
        etamax = (1 - epst)*etadt
        !call find_vlim_t(vmin2, vmax2)
        !print *, "vrange: ", vmin2/vth, vmax2/vth

        print *, "vrange: ", vmin2/vth, vmax2/vth
        if (vsteps > 0) print *, "D11/Dp (INT) = ", flux_integral(vmin2, vmax2, tol0)
        if (vsteps > 0) print *, "D11/Dp (MID) = ", flux_integral_mid(vmin2, vmax2)
        if (orbit_mode_transp > 0) print *, fluxint_box(1, :)

        D11 = 0d0
        D12 = 0d0
        D11sum = 0d0
        D12sum = 0d0
        dsdreff = 2d0/a*sqrt(s)
        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)

        open (unit=10, file=trim(adjustl(fn3)), recl=1024)
        v = vmax2
        vold = vmax2
        do ku = 0, nu - 1
            vold = v
            v = vmax2 - ku*(vmax2 - vmin2)/(nu - 1)

            call driftorbit_coarse(etamin, etamax, roots, nroots)
            if (nroots == 0) cycle

            do kr = 1, nroots
                eta_res = driftorbit_root(max(1d-9*abs(Om_tE), 1d-12), roots(kr, 1), roots(kr, 2))
                eta = eta_res(1)
                D11 = (vold - v)/vth*D11int_u(v/vth)*dsdreff**(-2)
                D12 = (vold - v)/vth*D12int_u(v/vth)*dsdreff**(-2)
                if ((ku == 1) .or. (ku == nu - 1)) then
                    D11 = D11/2d0
                    D12 = D12/2d0
                end if
                D11sum = D11sum + D11
                D12sum = D12sum + D12
            end do

            xs = (v/vth)**2
            kappa2s = (1d0/eta_res(1) - B0*(1d0 - eps))/(2d0*B0*eps)

            write (10, *) v/vth, eta_res(1), D11/Dp, D12/Dp, &
                (eta_res(1) - etamin)/(etamax - etamin), &
                vth, B0, xs, kappa2s, etatp, etadt
        end do
        print *, "D11/Dp (SUM) = ", D11sum/Dp, D12sum/Dp
        close (unit=10)
    end subroutine test_integral

    subroutine test_torque_integral
        real(8) :: eta_res(2)
        real(8) :: Tphi, dTphidu
        integer :: ku, nu
        real(8) :: kappa2s, du, ux
        real(8) :: roots(nlev, 3)
        integer :: nroots, kr
        character(1024) :: fn1, fn2, fn3
        real(8) :: tol0(2)

        tol0 = 1d-9

        call disp("test_torque_integral: Mach num Mt         = ", M_t)
        call disp("test_torque_integral: Poloidal mode mth   = ", 1d0*mth)
        print *, "test_torque_integral: Boxes sbox = ", sbox

        write (fn1, *) trim(adjustl(runname))//"_torque_integral_pco.out"
        write (fn2, *) trim(adjustl(runname))//"_torque_integral_pct.out"
        write (fn3, *) trim(adjustl(runname))//"_torque_integral_t.out"

        nu = 1000

        vmin2 = 0d0
        vmax2 = 3d0*vth

        if (.not. nopassing) then
            print *, "==CO-PASSING=="
            sigv = 1
            etamin = epsp*etatp
            etamax = (1 - epsp)*etatp

            print *, "vrange: ", vmin2/vth, vmax2/vth
            print *, "etarange: ", etamin, etamax

            Tphi = 0d0

            open (unit=10, file=trim(adjustl(fn1)), recl=1024)

            du = (vmax2 - vmin2)/(vth*nu)
            do ku = 0, nu - 1
                ux = vmin2/vth + du/2d0 + ku*du
                v = ux*vth

                call driftorbit_coarse(etamin, etamax, roots, nroots)
                if (nroots == 0) cycle

                open (unit=11, file=trim(adjustl(runname))//"_int_cop.out", recl=8192, position="append")
                do kr = 1, nroots
                    eta_res = driftorbit_root(max(1d-9*abs(Om_tE), 1d-12), roots(kr, 1), roots(kr, 2))
                    eta = eta_res(1)
                    write (11, "(i3, i3, f8.3, f8.3)") mth, kr, ux, eta_res(1)
                    dTphidu = Tphi_int(ux, eta_res)
                    Tphi = Tphi + dTphidu*du
                end do
                close (unit=11)

                kappa2s = (1d0/eta_res(1) - B0*(1d0 - eps))/(2d0*B0*eps)

                write (10, *) ux, eta_res(1), dTphidu, &
                    (eta_res(1) - etamin)/(etamax - etamin), &
                    vth, B0, kappa2s, 0d0, etatp
            end do
            print *, "Tphi (SUM) = ", Tphi/dVds
            close (unit=10)

            print *, "==CTR-PASSING=="
            sigv = -1
            etamin = epsp*etatp
            etamax = (1 - epsp)*etatp

            print *, "vrange: ", vmin2/vth, vmax2/vth
            print *, "etarange: ", etamin, etamax
            !if (vsteps>0) print *, "D11/Dp (INT) = ", flux_integral(vmin2, vmax2, tol0)
            !if (vsteps>0) print *, "D11/Dp (MID) = ", flux_integral_mid(vmin2, vmax2)
            !if(orbit_mode_transp>0) print *, fluxint_box(1,:)
            !print *, "D11/Dp (ODE) = ", flux_integral_ode(vmin2, vmax2)

            Tphi = 0d0

            open (unit=10, file=trim(adjustl(fn2)), recl=1024)

            du = (vmax2 - vmin2)/(vth*nu)
            do ku = 0, nu - 1
                ux = vmin2/vth + du/2d0 + ku*du
                v = ux*vth

                call driftorbit_coarse(etamin, etamax, roots, nroots)
                if (nroots == 0) cycle

                open (unit=11, file=trim(adjustl(runname))//"_int_ctr.out", recl=8192, position="append")
                do kr = 1, nroots
                    eta_res = driftorbit_root(max(1d-9*abs(Om_tE), 1d-12), roots(kr, 1), roots(kr, 2))
                    eta = eta_res(1)
                    write (11, "(i3, i3, f8.3, f8.3)") mth, kr, ux, eta_res(1)
                    dTphidu = Tphi_int(ux, eta_res)
                    Tphi = Tphi + dTphidu*du
                end do
                close (unit=11)

                kappa2s = (1d0/eta_res(1) - B0*(1d0 - eps))/(2d0*B0*eps)

                write (10, *) ux, eta_res(1), dTphidu, &
                    (eta_res(1) - etamin)/(etamax - etamin), &
                    vth, B0, kappa2s, 0d0, etatp
            end do
            print *, "Tphi (SUM) = ", Tphi/dVds
            close (unit=10)
        end if

        print *, "==TRAPPED=="
        sigv = 1d0
        etamin = (1 + epst)*etatp
        etamax = (1 - epst)*etadt

        print *, "vrange: ", vmin2/vth, vmax2/vth

        Tphi = 0d0

        open (unit=10, file=trim(adjustl(fn3)), recl=1024)

        du = (vmax2 - vmin2)/(vth*nu)
        do ku = 0, nu - 1
            ux = vmin2/vth + du/2d0 + ku*du
            v = ux*vth

            call driftorbit_coarse(etamin, etamax, roots, nroots)
            if (nroots == 0) cycle

            open (unit=11, file=trim(adjustl(runname))//"_int_t.out", recl=8192, position="append")
            do kr = 1, nroots
                eta_res = driftorbit_root(max(1d-9*abs(Om_tE), 1d-12), roots(kr, 1), roots(kr, 2))
                eta = eta_res(1)
                write (11, "(i3, i3, f8.3, f8.3)") mth, kr, ux, eta_res(1)
                dTphidu = Tphi_int(ux, eta_res)
                Tphi = Tphi + dTphidu*du
            end do
            close (unit=11)

            kappa2s = (1d0/eta_res(1) - B0*(1d0 - eps))/(2d0*B0*eps)

            write (10, *) ux, eta_res(1), dTphidu, &
                (eta_res(1) - etamin)/(etamax - etamin), &
                vth, B0, kappa2s, 0d0, etatp
        end do
        print *, "Tphi (SUM) = ", Tphi/dVds
        close (unit=10)
    end subroutine test_torque_integral

    subroutine test_profile
        integer :: k
        real(8), allocatable :: data(:, :)

        print *, "test_profile"

        call readdata("profile.in", 3, data)

        open (unit=9, file=trim(adjustl(runname))//"_profile.out", recl=1024)
        write (9, *) "#s deltas=c*mi*vth*hcovar(2)*q/(qi*psi_pr) dVds q psi_pr"
        do k = 1, size(data, 1)
            s = data(k, 1)
            call output_flux_surface_data(9)
        end do
        s = 0.96
        call output_flux_surface_data(9)
        s = 0.97
        call output_flux_surface_data(9)
        s = 0.98
        call output_flux_surface_data(9)
        s = 0.99
        call output_flux_surface_data(9)
        s = 1.0
        call output_flux_surface_data(9)

        close (unit=9)

        deallocate (data)
    end subroutine test_profile

    subroutine output_flux_surface_data(unit)
        integer, intent(in) :: unit
        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        x(1) = s
        x(2) = 0d0
        x(3) = 0d0
        call init_plasma
        call init_fsa
        call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
        write (unit, *) x(1), c*mi*vth*hcovar(2)*q/(qi*psi_pr), dVds, q, psi_pr, vth
    end subroutine output_flux_surface_data
end program test_misc
