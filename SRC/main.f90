program neo_rt
    implicit none

    call main

contains

    subroutine main
        use driftorbit
        use lineint
        use do_magfie_mod, only: s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds, &
                                 q, dqds, iota, R0, a, eps, inp_swi !,B0h, B00
        use do_magfie_pert_mod, only: do_magfie_pert_init, do_magfie_pert_amp, mph
        use polylag_3, only: mp, indef, plag1d
        use util, only: clearfile
        implicit none

        logical :: odeint
        character(len=1024) :: tmp
        character(:), allocatable :: runname
        integer :: tmplen
        logical :: plasmafile, profilefile
        character(len=64) :: runmode = "torque"

        call get_command_argument(1, tmp, tmplen)
        runname = trim(adjustl(tmp))

        call read_control

        call do_magfie_init  ! init axisymmetric part of field from infile
        if (pertfile) then
            ! init non-axisymmetric perturbation of field from infile_pert
            call do_magfie_pert_init
        end if

        ! Init plasma profiles of radial electric field.
        ! Read profile.in in cases where it's needed.
        inquire (file="profile.in", exist=profilefile)
        if (profilefile) then
            call init_profile
        else
            if (orbit_mode_transp > 0) then
                stop "need profile.in for finite orbit width transport"
            elseif (intoutput .or. nonlin) then
                stop "need profile.in for integral output or nonlinear calculation"
            end if
        end if

        inquire (file="plasma.in", exist=plasmafile)
        ! init plasma density and temperature profiles
        if (plasmafile) then
            call init_plasma
        else
            if ((runmode == "torque") .or. nonlin) then
                stop "need plasma.in for nonlinear or torque calculation"
            end if
        end if

        call init_test

        call test_magfie

        if (runmode == "test_profile") then
            if (.not. profilefile) error stop "need profile.in for test_profile"
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
        elseif (runmode == "torque") then
            call compute_torque
            stop
        elseif (runmode == "transport") then
            call test_machrange2
            stop
        end if
    end subroutine main

    subroutine read_control
        character(1) :: dummy
        real(8) :: qs, ms

        namelist /params/ s, M_t, qs, ms, vth, epsmn, m0, mph, mth, supban, &
            magdrift, nopassing, noshear, pertfile, odeint, nonlin, &
            bfac, efac, inp_swi, orbit_mode_avg, orbit_mode_transp, &
            vsteps, intoutput, runmode

        open (unit=9, file=trim(adjustl(runname))//".in", status="old", form="formatted")
        read (9, nml=params)
        close (unit=9)

        M_t = M_t*efac/bfac
        qi = qs*qe
        mi = ms*mu
    end subroutine read_control

    subroutine init_profile
        ! Init s profile for finite orbit width boxes in radial s
        real(8), allocatable :: data(:, :)
        real(8) :: splineval(3)
        integer :: k

        call readdata("profile.in", 3, data)

        allocate (Mt_spl_coeff(size(data(:, 1)), 5))

        Mt_spl_coeff = spline_coeff(data(:, 1), data(:, 2))
        splineval = spline_val_0(Mt_spl_coeff, s)

        M_t = splineval(1)*efac/bfac
        dM_tds = splineval(2)*efac/bfac

        if (orbit_mode_transp > 0) then
            allocate (sbox(size(data, 1) + 1))
            allocate (taubins(size(sbox) + 1))
            sbox(1) = 0.0
            sbox(size(sbox)) = 1.0
            do k = 1, (size(data, 1) - 1)
                sbox(k + 1) = (data(k, 1) + data(k + 1, 1))/2d0
            end do
            allocate (fluxint_box(2, size(sbox) + 1))
            allocate (torque_int_box(size(sbox) + 1))
        end if

        deallocate (data)
    end subroutine init_profile

    subroutine init_plasma
        real(8), parameter :: pmass = 1.6726d-24

        real(8) :: amb, am1, am2, Zb, Z1, Z2, dchichi, slowrate, dchichi_norm, slowrate_norm
        real(8) :: v0, ebeam
        real(8), dimension(:, :), allocatable :: plasma(:, :)
        integer, dimension(mp) :: indu
        real(8), dimension(mp) :: xp, fp
        real(8) :: dxm1
        integer :: nplasma, i

        ! read plasma file
        open (1, file="plasma.in")
        read (1, *)
        read (1, *) nplasma, am1, am2, Z1, Z2
        read (1, *)
        allocate (plasma(nplasma, 6))
        do i = 1, nplasma
            read (1, *) plasma(i, :)
        end do
        dxm1 = 1.d0/(plasma(2, 1) - plasma(1, 1))
        close (1)

        ! interpolate to s value
        call indef(s, plasma(1, 1), dxm1, nplasma, indu)

        xp = plasma(indu, 1)
        fp = plasma(indu, 2)
        call plag1d(s, fp, dxm1, xp, ni1, dni1ds)
        fp = plasma(indu, 3)
        call plag1d(s, fp, dxm1, xp, ni2, dni2ds)
        fp = plasma(indu, 4)
        call plag1d(s, fp, dxm1, xp, Ti1, dTi1ds)
        fp = plasma(indu, 5)
        call plag1d(s, fp, dxm1, xp, Ti2, dTi2ds)
        fp = plasma(indu, 6)
        call plag1d(s, fp, dxm1, xp, Te, dTeds)

        qi = Z1*qe
        mi = am1*mu
        vth = sqrt(2d0*Ti1*ev/mi)
        v0 = vth
        amb = 2d0
        Zb = 1d0
        ebeam = amb*pmass*v0**2/(2d0*ev)

        call loacol_nbi(amb, am1, am2, Zb, Z1, Z2, ni1, ni2, Ti1, Ti2, Te, &
                        ebeam, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)

    end subroutine init_plasma

    subroutine init_test
        call init
        !if (supban) then
        ! set thermal velocity so that ExB = reference toroidal drift
        !print *, "Set vth from: ", vth
        !   vth = abs(2*M_t*qi*psi_pr/(mi*c*R0)) ! thermal velocity
        !   print *, "Set vth to: ", vth
        !end if
        Om_tE = vth*M_t/R0                   ! toroidal ExB drift frequency
        dOm_tEds = vth*dM_tds/R0

        etamin = (1 + epst)*etatp
        etamax = (1 - epst)*etadt
        sigv = 1

        if (runmode == "torque") then
            ! thermodynamic forces
            A1 = dni1ds/ni1 - qi/(Ti1*ev)*psi_pr/(q*c)*Om_tE - 3d0/2d0*dTi1ds/Ti1
            A2 = dTi1ds/Ti1

            print *, "A1 = ", A1, "A2 = ", A2
        end if
    end subroutine init_test

    subroutine test_magfie
!    use do_magfie_neo_mod, only: s_neo => s, do_magfie_neo_init => do_magfie_init,&
!         do_magfie_neo => do_magfie
!    use do_magfie_pert_neo_mod, only: do_magfie_pert_neo_init => do_magfie_pert_init,&
!         do_magfie_pert_neo_amp => do_magfie_pert_amp

        integer, parameter :: nth = 50
        integer :: k
        real(8) :: thmin, thmax
        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(8) :: Dp, Drp
        real(8) :: ux, dpp, dhh, fpeff
        complex(8) :: bn

! comparison with Neo2 magfie
!    s_neo = s
!    call do_magfie_neo_init
! comparison with Neo2 magfie pert
!    if (pertfile) then
!       call do_magfie_pert_neo_init
!    end if

        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
        Drp = 4*mph*q/(eps**2*sqrt(pi))  ! actually this is Drp/Dp

        open (unit=9, file=trim(adjustl(runname))//"_magfie_param.out", recl=1024)

        thmin = -pi
        thmax = pi
        x(1) = s
        x(2) = 0d0
        x(3) = 0d0

        write (9, *) "-------------------------"
        write (9, *) "test_magfie: s         = ", s
        write (9, *) "test_magfie: R0        = ", R0
        write (9, *) "test_magfie: a         = ", a
        write (9, *) "test_magfie: eps       = ", eps
        write (9, *) "test_magfie: A         = ", 1/eps
        write (9, *) "test_magfie: psi_pr    = ", psi_pr
        write (9, *) "test_magfie: B0        = ", B0
        !write(9,*) "test_magfie: B0h       = ", B0h
        !write(9,*) "test_magfie: B00       = ", B00
        write (9, *) "test_magfie: Bthcov    = ", Bthcov
        write (9, *) "test_magfie: Bphcov    = ", Bphcov
        write (9, *) "test_magfie: dBthcovds = ", dBthcovds
        write (9, *) "test_magfie: dBphcovds = ", dBphcovds
        write (9, *) "test_magfie: q         = ", q
        write (9, *) "test_magfie: iota      = ", iota
        write (9, *) "test_magfie: dVds      = ", dVds
        write (9, *) "test_magfie: M_t       = ", M_t
        write (9, *) "test_magfie: Om_tE     = ", Om_tE
        write (9, *) "test_magfie: Om_tBref  = ", c*mi*vth**2/(2d0*qi*psi_pr)
        write (9, *) "test_magfie: vth       = ", vth
        write (9, *) "test_magfie: T [eV]    = ", mi/2d0*vth**2/eV
        write (9, *) "test_magfie: m0        = ", 1d0*m0
        write (9, *) "test_magfie: n0        = ", 1d0*mph
        write (9, *) "test_magfie: Dp        = ", Dp
        write (9, *) "test_magfie: Drp/Dp    = ", Drp
        write (9, *) "test_magfie: etatp     = ", etatp
        write (9, *) "test_magfie: etadt     = ", etadt
        write (9, *) "-------------------------"
        write (9, *) "test_magfie: pertfile  = ", pertfile
        if (nonlin) then
            write (9, *) "-------------------------"
            write (9, *) "test_magfie: dpp       = ", dpp
            write (9, *) "test_magfie: dhh       = ", dhh
            write (9, *) "test_magfie: dfpeff    = ", fpeff
            write (9, *) "-------------------------"
        end if

        close (unit=9)

        open (unit=9, file=trim(adjustl(runname))//"_magfie.out", recl=1024)
!    open(unit=10, file=trim(adjustl(runname))//"_magfie_neo.out", recl=1024)
        do k = 0, nth - 1
            x(3) = thmin + k*(thmax - thmin)/(nth - 1)
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            if (pertfile) then
                call do_magfie_pert_amp(x, bn)
                bn = epsmn*bn/bmod
            else
                bn = epsmn*exp(imun*m0*x(3))
            end if
            write (9, *) x(3), bmod, sqrtg, hder(1), hder(2), hder(3), hcovar(1), &
                hcovar(2), hcovar(3), hctrvr(1), hctrvr(2), hctrvr(3), &
                hcurl(1), hcurl(2), hcurl(3), &
                real(bn), aimag(bn), real(epsmn*exp(imun*m0*x(3))), &
                aimag(epsmn*exp(imun*m0*x(3))) !16
!       call do_magfie_neo( x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl )
!       if (pertfile) then
!          call do_magfie_pert_neo_amp( x, bn )
!          bn = bn/bmod
!       else
!          bn = epsmn*exp(imun*m0*x(3))
!       end if
!      write(10,*) x(3), bmod, sqrtg, hder(1), hder(2), hder(3), hcovar(1),&
!           hcovar(2), hcovar(3), hctrvr(1), hctrvr(2), hctrvr(3),&
!            hcurl(1), hcurl(2), hcurl(3), &
!            real(bn), aimag(bn), real(epsmn*exp(imun*m0*x(3))),&
!            aimag(epsmn*exp(imun*m0*x(3))) !16
        end do
!    close(unit=10)
        close (unit=9)
    end subroutine test_magfie

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

        real(8) :: tol
        integer :: n

        integer :: k, state
        real(8) :: ti
        real(8) :: y(2), ydot(2), yold(2)

        real(8) :: atol(nvar), rtol, tout, rstats(22)
        integer :: neq, itask, istate, istats(31), numevents
        type(vode_opts) :: options

        real(8) :: bmod, sqrtg, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)

        real(8) :: s1old, told, s1dot, sbound
        real(8), allocatable :: taubins(:)
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

    !> Do calculation for a range of mach numbers
    !>
    !> Range of mach numbers depends on input parameters Mtmin, Mtmax and
    !> Mtnum.
    !> Output file names depend on commandline parameter "runname".
    subroutine test_machrange2
        integer :: j, k
        real(8) :: fluxrespco(2), fluxrespctr(2), &
                   fluxrest(2), fluxpco(2), fluxpctr(2), fluxt(2)
        real(8) :: vminp, vmaxp, vmint, vmaxt
        integer :: mthmin, mthmax
        integer :: mth0co, mth0ctr, mth0t, mthco, mthctr, mtht ! mode counters
        character(1024) :: fn1, fn2

        real(8) :: tolpco(2), tolpctr(2), tolt(2), tol0(2)

        logical :: firstloop

        integer :: kb

        firstloop = .true.

        fluxpco = 0d0
        fluxpctr = 0d0
        fluxt = 0d0

        write (fn1, *) trim(adjustl(runname))//".out"
        write (fn2, *) trim(adjustl(runname))//"_integral.out"

        ! absolute tolerances for integrals
        tol0 = 1d-15
        tolpco = tol0
        tolpctr = tol0
        tolt = tol0

        fluxpco = 0d0
        fluxpctr = 0d0
        fluxt = 0d0

        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0

        if (supban) then
            mthmin = 0
            mthmax = 0
        else
            mthmin = -ceiling(2*abs(mph)*q)
            mthmax = ceiling(2*abs(mph)*q)
        end if

        mth0co = 0
        mth0ctr = 0
        mth0t = 0
        mthco = 0
        mthctr = 0
        mtht = 0

        do j = mthmin, mthmax
            mth = j
            fluxrest = 0d0
            fluxrespco = 0d0
            fluxrespctr = 0d0

            vminp = 1d-6*vth
            vmaxp = 3d0*vth
            vmint = 1d-6*vth
            vmaxt = 3d0*vth

            ! superbanana resonance
            if (supban) then
                sigv = 1
                vmint = 0.01*vth
                vmaxt = 5*vth
                etamin = (1 + epst)*etatp
                etamax = (1 - epst)*etadt
                if (odeint) then
                    fluxrest = flux_integral_ode(vmint, vmaxt)
                elseif (vsteps > 0) then
                    fluxrest = flux_integral_mid(vmint, vmaxt)
                else
                    fluxrest = flux_integral(vmint, vmaxt, tol0)
                end if
                fluxt = fluxt + fluxrest
            else
                ! passing resonance (co-passing)
                if (.not. nopassing) then
                    sigv = 1
                    etamin = epsp*etatp
                    etamax = (1 - epsp)*etatp
                    if (odeint) then
                        fluxrespco = flux_integral_ode(vminp, vmaxp)
                    elseif (vsteps > 0) then
                        fluxrespco = flux_integral_mid(vminp, vmaxp)
                    else
                        fluxrespco = flux_integral(vminp, vmaxp, tolpco)
                    end if
                    fluxpco = fluxpco + fluxrespco
                    tolpco(1) = max(1d-6*fluxpco(1), tolpco(1))
                    tolpco(2) = max(1d-6*fluxpco(2), tolpco(2))
                    if (orbit_mode_transp > 0) then
                        open (unit=9, file=trim(adjustl(runname))//"_box_cop.out", recl=8192, position="append")
                        write (9, *) mth, fluxint_box(1, :), fluxint_box(2, :)
                        close (unit=9)
                    end if
                end if

                ! passing resonance (counter-passing)
                if (.not. nopassing) then
                    sigv = -1
                    etamin = epsp*etatp
                    etamax = (1 - epsp)*etatp
                    if (odeint) then
                        fluxrespctr = flux_integral_ode(vminp, vmaxp)
                    elseif (vsteps > 0) then
                        fluxrespctr = flux_integral_mid(vminp, vmaxp)
                    else
                        fluxrespctr = flux_integral(vminp, vmaxp, tolpctr)
                    end if
                    fluxpctr = fluxpctr + fluxrespctr
                    tolpctr(1) = max(1d-6*fluxpctr(1), tolpctr(1))
                    tolpctr(2) = max(1d-6*fluxpctr(2), tolpctr(2))
                    if (orbit_mode_transp > 0) then
                        open (unit=9, file=trim(adjustl(runname))//"_box_ctr.out", recl=8192, position="append")
                        write (9, *) mth, fluxint_box(1, :), fluxint_box(2, :)
                        close (unit=9)
                    end if
                end if

                ! trapped resonance (trapped)
                sigv = 1
                etamin = (1 + epst)*etatp
                etamax = (1 - epst)*etadt
                if (odeint) then
                    fluxrest = flux_integral_ode(vmint, vmaxt)
                elseif (vsteps > 0) then
                    fluxrest = flux_integral_mid(vmint, vmaxt)
                else
                    if (mth == 0) then
                        fluxrest = flux_integral(vmint, vmaxt, tol0)
                    else
                        fluxrest = flux_integral(vmint, vmaxt, tolt)
                    end if
                end if
                fluxt = fluxt + fluxrest
                tolt(1) = max(1d-6*fluxt(1), tolt(1))
                tolt(2) = max(1d-6*fluxt(2), tolt(2))
                if (orbit_mode_transp > 0) then
                    open (unit=9, file=trim(adjustl(runname))//"_box_t.out", recl=8192, position="append")
                    write (9, *) mth, fluxint_box(1, :), fluxint_box(2, :)
                    close (unit=9)
                end if
            end if

            print *, ""
            print *, "test_flux: Mt = ", M_t, ", mth = ", mth
            write (*, "(4ES12.2,2F12.2)") fluxrespco(1), fluxrespctr(1), &
                fluxrest(1), fluxrespco(1) + fluxrespctr(1) + fluxrest(1), &
                vminp/vth, vmint/vth
            write (*, "(4ES12.2,2F12.2)") fluxrespco(2), fluxrespctr(2), &
                fluxrest(2), fluxrespco(2) + fluxrespctr(2) + fluxrest(2), &
                vmaxp/vth, vmaxt/vth

            open (unit=10, file=trim(adjustl(fn2)), recl=1024, position="append")
            write (10, *) M_t, mth, fluxrespco(1), fluxrespctr(1), fluxrest(1), &
                fluxrespco(1) + fluxrespctr(1) + fluxrest(1), &
                fluxrespco(2), fluxrespctr(2), fluxrest(2), &
                fluxrespco(2) + fluxrespctr(2) + fluxrest(2), vminp/vth, vmaxp/vth, &
                vmint/vth, vmaxt/vth
            close (unit=10)

            if (supban) then
                exit
            end if
        end do
        ! Here D11 and D12 are written - just different names have been used.
        ! pco is part from co-passing particles, pctr is part from counter
        ! passing particles, t is part from trapped particles. Sum of all
        ! three terms is also stored.
        ! Index 1 is D11, index 2 D12.
        open (unit=9, file=trim(adjustl(fn1)), recl=1024, position="append")
        write (9, *) M_t, fluxpco(1), fluxpctr(1), fluxt(1), &
            fluxpco(1) + fluxpctr(1) + fluxt(1), &
            fluxpco(2), fluxpctr(2), fluxt(2), &
            fluxpco(2) + fluxpctr(2) + fluxt(2)
        close (unit=9)
    end subroutine test_machrange2

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

    subroutine compute_torque
        ! computes torque per radial distance dTphi/ds at s=sphi

        integer :: j, k
        real(8) :: Tresco, Tresctr, Trest, Tco, Tctr, Tt
        real(8) :: vminp, vmaxp, vmint, vmaxt
        integer :: mthmin, mthmax
        integer :: mth0co, mth0ctr, mth0t, mthco, mthctr, mtht ! mode counters
        character(1024) :: fn1, fn2

        real(8) :: tolpco(2), tolpctr(2), tolt(2), tol0(2)

        logical :: firstloop

        firstloop = .true.

        write (fn1, *) trim(adjustl(runname))//"_torque.out"
        call clearfile(fn1)
        write (fn2, *) trim(adjustl(runname))//"_torque_integral.out"
        call clearfile(fn2)

        call clearfile(trim(adjustl(runname))//"_torque_box_co.out")
        call clearfile(trim(adjustl(runname))//"_torque_box_ctr.out")
        call clearfile(trim(adjustl(runname))//"_torque_box_t.out")

        ! absolute tolerances for integrals
        tol0 = 1d-15
        tolpco = tol0
        tolpctr = tol0
        tolt = tol0

        Tco = 0d0
        Tctr = 0d0
        Tt = 0d0

        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0

        if (supban) then
            mthmin = 0
            mthmax = 0
        else
            mthmin = -ceiling(2*abs(mph)*q)
            mthmax = ceiling(2*abs(mph)*q)
        end if

        mth0co = 0
        mth0ctr = 0
        mth0t = 0
        mthco = 0
        mthctr = 0
        mtht = 0

        if (intoutput) open (unit=11, file=trim(adjustl(runname))//"_intoutput.out", recl=1024)

        do j = mthmin, mthmax
            mth = j

            vminp = 1d-6*vth
            vmaxp = 4d0*vth

            vmint = vminp
            vmaxt = vmaxp

            Tresco = 0d0
            Tresctr = 0d0
            Trest = 0d0

            ! superbanana resonance
            if (supban) then
                sigv = 1
                vmint = 0.01*vth
                vmaxt = 5*vth
                etamin = (1 + epst)*etatp
                etamax = (1 - epst)*etadt
                Trest = torque_integral_mid(vmint, vmaxt)
                Tt = Tt + Trest
            else
                ! passing resonance (co-passing)
                if (.not. nopassing) then
                    sigv = 1
                    etamin = epsp*etatp
                    etamax = (1 - epsp)*etatp
                    Tresco = torque_integral_mid(vminp, vmaxp)
                    Tco = Tco + Tresco
                    if (orbit_mode_transp > 0) then
                        open (unit=9, file=trim(adjustl(runname))//"_torque_box_co.out", recl=8192, &
                              position="append")
                        write (9, *) mth, torque_int_box
                        close (unit=9)
                    end if
                end if

                ! passing resonance (counter-passing)
                if (.not. nopassing) then
                    sigv = -1
                    etamin = epsp*etatp
                    etamax = (1 - epsp)*etatp
                    Tresctr = torque_integral_mid(vminp, vmaxp)
                    Tctr = Tctr + Tresctr
                    if (orbit_mode_transp > 0) then
                        open (unit=9, file=trim(adjustl(runname))//"_torque_box_ctr.out", recl=8192, &
                              position="append")
                        write (9, *) mth, torque_int_box
                        close (unit=9)
                    end if
                end if

                ! trapped resonance (trapped)
                sigv = 1
                etamin = (1 + epst)*etatp
                etamax = (1 - epst)*etadt
                Trest = torque_integral_mid(vmint, vmaxt)
                Tt = Tt + Trest
                if (orbit_mode_transp > 0) then
                    open (unit=9, file=trim(adjustl(runname))//"_torque_box_t.out", recl=8192, &
                          position="append")
                    write (9, *) mth, torque_int_box
                    close (unit=9)
                end if
            end if

            print *, ""
            print *, "compute_torque: Mt = ", M_t, ", mth = ", mth
            write (*, "(4ES12.2,2F12.2)") Tresco, Tresctr, &
                Trest, Tresco + Tresctr + Trest

            open (unit=10, file=trim(adjustl(fn2)), recl=1024, position="append")
            write (10, *) mth, Tresco, Tresctr, Trest
            close (unit=10)

            if (supban) then
                exit
            end if
        end do
        open (unit=9, file=trim(adjustl(fn1)), recl=1024, position="append")
        write (9, *) s, dVds, M_t, Tco, Tctr, Tt
        close (unit=9)
        if (intoutput) close (unit=11)
    end subroutine compute_torque

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

end program neo_rt
