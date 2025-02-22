module neort
    use driftorbit
    implicit none

    character(:), allocatable :: runname
    character(len=64) :: runmode = "torque"

contains

    subroutine main
        use do_magfie_mod, only: do_magfie_init
        use do_magfie_pert_mod, only: do_magfie_pert_init
        use driftorbit, only: pertfile, orbit_mode_transp, intoutput, nonlin

        character(len=1024) :: tmp
        integer :: tmplen
        logical :: use_rotation_profiles, use_thermodynamic_profiles

        call get_command_argument(1, tmp, tmplen)
        runname = trim(adjustl(tmp))

        call read_control

        call do_magfie_init  ! init axisymmetric part of field from infile
        if (pertfile) then
            ! init non-axisymmetric perturbation of field from infile_pert
            call do_magfie_pert_init
        end if

        ! Init plasma profiles of radial electric field.
        ! Read profile.in in cases where it is needed.
        inquire (file="profile.in", exist=use_rotation_profiles)
        if (use_rotation_profiles) then
            call init_profile
        else
            if (orbit_mode_transp > 0) then
                stop "need profile.in for finite orbit width transport"
            elseif (intoutput .or. nonlin) then
                stop "need profile.in for integral output or nonlinear calculation"
            end if
        end if

        inquire (file="plasma.in", exist=use_thermodynamic_profiles)
        ! init plasma density and temperature profiles
        if (use_thermodynamic_profiles) then
            call init_plasma
        else
            if ((runmode == "torque") .or. nonlin) then
                stop "need plasma.in for nonlinear or torque calculation"
            end if
        end if

        call init_run(use_thermodynamic_profiles)
        call check_magfie
        call compute_transport
    end subroutine main

    subroutine read_control
        use driftorbit
        real(8) :: qs, ms

        namelist /params/ runmode, s, M_t, qs, ms, vth, epsmn, m0, mph, supban, &
            magdrift, nopassing, noshear, pertfile, nonlin, bfac, efac, inp_swi, &
            orbit_mode_avg, orbit_mode_transp, vsteps, intoutput

        open (unit=9, file=trim(adjustl(runname))//".in", status="old", form="formatted")
        read (9, nml=params)
        close (unit=9)

        M_t = M_t*efac/bfac
        qi = qs*qe
        mi = ms*mu
    end subroutine read_control

    subroutine init_profile
        ! Init s profile for finite orbit width boxes in radial s
        use driftorbit, only: s, M_t, dM_tds, efac, bfac, orbit_mode_transp
        use spline

        ! For splining electric precession frequency
        real(8), allocatable :: Mt_spl_coeff(:, :)

        real(8), allocatable :: data(:, :)
        real(8) :: splineval(3)

        call readdata("profile.in", 3, data)

        allocate (Mt_spl_coeff(size(data(:, 1)), 5))

        Mt_spl_coeff = spline_coeff(data(:, 1), data(:, 2))
        splineval = spline_val_0(Mt_spl_coeff, s)

        M_t = splineval(1)*efac/bfac
        dM_tds = splineval(2)*efac/bfac

        if (orbit_mode_transp > 0) then
            call init_profile_finite_width(data)
        end if
    end subroutine init_profile

    subroutine init_profile_finite_width(data)
        use driftorbit, only: sbox, taubins, fluxint_box, torque_int_box

        real(8), allocatable, intent(in) :: data(:, :)

        integer :: k

        allocate (sbox(size(data, 1) + 1))
        allocate (taubins(size(sbox) + 1))
        sbox(1) = 0.0
        sbox(size(sbox)) = 1.0
        do k = 1, (size(data, 1) - 1)
            sbox(k + 1) = (data(k, 1) + data(k + 1, 1))/2d0
        end do
        allocate (fluxint_box(2, size(sbox) + 1))
        allocate (torque_int_box(size(sbox) + 1))
    end subroutine init_profile_finite_width

    subroutine init_plasma
        use driftorbit, only: s, ni1, dni1ds, ni2, dni2ds, Ti1, dTi1ds, Ti2, dTi2ds, &
                              Te, dTeds, qi, qe, mi, mu, ev, loacol_nbi, vth
        use polylag_3, only: plag1d, indef, mp

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
        ! dvthds = TODO
        v0 = vth
        amb = 2d0
        Zb = 1d0
        ebeam = amb*pmass*v0**2/(2d0*ev)

        call loacol_nbi(amb, am1, am2, Zb, Z1, Z2, ni1, ni2, Ti1, Ti2, Te, &
                        ebeam, v0, dchichi, slowrate, dchichi_norm, slowrate_norm)

    end subroutine init_plasma

    subroutine init_run(use_thermodynamic_profiles)

        logical, intent(in) :: use_thermodynamic_profiles

        call init

        Om_tE = vth*M_t/R0                   ! toroidal ExB drift frequency
        dOm_tEds = vth*dM_tds/R0 ! TODO + M_t*dvthds/R0

        sigv = 1
        call get_trapped_region(etamin, etamax)

        if (use_thermodynamic_profiles) then
            A1 = dni1ds/ni1 - qi/(Ti1*ev)*psi_pr/(q*c)*Om_tE - 3d0/2d0*dTi1ds/Ti1
            A2 = dTi1ds/Ti1
        end if
    end subroutine init_run

    subroutine check_magfie

        integer, parameter :: nth = 50
        integer :: k
        real(8) :: thmin, thmax
        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(8) :: Dp, Drp
        real(8) :: dpp, dhh, fpeff
        complex(8) :: bn

        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
        Drp = 4*mph*q/(eps**2*sqrt(pi))  ! actually this is Drp/Dp

        open (unit=9, file=trim(adjustl(runname))//"_magfie_param.out", recl=1024)

        thmin = -pi
        thmax = pi
        x(1) = s
        x(2) = 0d0
        x(3) = 0d0

        write (9, *) "-------------------------"
        write (9, *) "check_magfie: s         = ", s
        write (9, *) "check_magfie: R0        = ", R0
        write (9, *) "check_magfie: a         = ", a
        write (9, *) "check_magfie: eps       = ", eps
        write (9, *) "check_magfie: A         = ", 1/eps
        write (9, *) "check_magfie: psi_pr    = ", psi_pr
        write (9, *) "check_magfie: B0        = ", B0
        write (9, *) "check_magfie: Bthcov    = ", Bthcov
        write (9, *) "check_magfie: Bphcov    = ", Bphcov
        write (9, *) "check_magfie: dBthcovds = ", dBthcovds
        write (9, *) "check_magfie: dBphcovds = ", dBphcovds
        write (9, *) "check_magfie: q         = ", q
        write (9, *) "check_magfie: iota      = ", iota
        write (9, *) "check_magfie: dVds      = ", dVds
        write (9, *) "check_magfie: M_t       = ", M_t
        write (9, *) "check_magfie: Om_tE     = ", Om_tE
        write (9, *) "check_magfie: Om_tBref  = ", c*mi*vth**2/(2d0*qi*psi_pr)
        write (9, *) "check_magfie: vth       = ", vth
        write (9, *) "check_magfie: T [eV]    = ", mi/2d0*vth**2/eV
        write (9, *) "check_magfie: m0        = ", 1d0*m0
        write (9, *) "check_magfie: n0        = ", 1d0*mph
        write (9, *) "check_magfie: Dp        = ", Dp
        write (9, *) "check_magfie: Drp/Dp    = ", Drp
        write (9, *) "check_magfie: etatp     = ", etatp
        write (9, *) "check_magfie: etadt     = ", etadt
        write (9, *) "-------------------------"
        write (9, *) "check_magfie: pertfile  = ", pertfile
        if (nonlin) then
            write (9, *) "-------------------------"
            write (9, *) "check_magfie: dpp       = ", dpp
            write (9, *) "check_magfie: dhh       = ", dhh
            write (9, *) "check_magfie: dfpeff    = ", fpeff
            write (9, *) "-------------------------"
        end if

        close (unit=9)

        open (unit=9, file=trim(adjustl(runname))//"_magfie.out", recl=1024)
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
                aimag(epsmn*exp(imun*m0*x(3)))
        end do
        close (unit=9)
    end subroutine check_magfie

    subroutine compute_transport
        integer :: j

        ! Transport coefficients D11, D12 in approximate effective radius r=2d0/a*sqrt(s)
        real(8) :: Dco(2), Dctr(2), Dt(2)

        ! Torque density dTphi_int/ds for integration over normalized toroidal flux s
        real(8) :: Tco, Tctr, Tt
        integer :: mthmin, mthmax

        logical :: firstloop

        firstloop = .true.

        Dco = 0d0
        Dctr = 0d0
        Dt = 0d0
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

        do j = mthmin, mthmax
            call compute_transport_harmonic(j, Dco, Dctr, Dt, Tco, Tctr, Tt)
        end do

        open (unit=9, file=trim(adjustl(runname))//".out", recl=1024)
        write (9, *) "# M_t D11t D11ctr D11co D11 D12t D12ctr D12co D12"
        write (9, *) M_t, Dco(1), Dctr(1), Dt(1), Dco(1) + Dctr(1) + Dt(1), &
            Dco(2), Dctr(2), Dt(2), Dco(2) + Dctr(2) + Dt(2)
        close (unit=9)

        open (unit=9, file=trim(adjustl(runname))//"_torque.out", recl=1024)
        write (9, *) "# s dVds M_t Tco Tctr Tt"
        write (9, *) s, dVds, M_t, Tco, Tctr, Tt
        close (unit=9)
    end subroutine compute_transport

    subroutine compute_transport_harmonic(j, Dco, Dctr, Dt, Tco, Tctr, Tt)
        integer, intent(in) :: j
        real(8), intent(inout) :: Dco(2), Dctr(2), Dt(2), Tco, Tctr, Tt

        real(8) :: Dresco(2), Dresctr(2), Drest(2), Tresco, Tresctr, Trest
        real(8) :: vminp, vmaxp, vmint, vmaxt

        mth = j

        Drest = 0d0
        Dresco = 0d0
        Dresctr = 0d0
        Tresco = 0d0
        Tresctr = 0d0
        Trest = 0d0

        vminp = 1d-6*vth
        vmaxp = 3d0*vth
        vmint = vminp
        vmaxt = vmaxp

        ! Superbanana resonance
        if (supban) then
            vmint = 0.01d0*vth
            vmaxt = 5d0*vth
            sigv = 1
            call get_trapped_region(etamin, etamax)
            call transport_integral_mid(vmint, vmaxt, Drest, Trest)
            Dt = Dt + Drest
            Tt = Tt + Trest
        else
            ! Passing resonance (co-passing)
            if (.not. nopassing) then
                sigv = 1
                call get_passing_region(etamin, etamax)
                call transport_integral_mid(vminp, vmaxp, Dresco, Tresco)
                Dco = Dco + Dresco
                Tco = Tco + Tresco
                if (orbit_mode_transp > 0) then
                    call write_transport_box(trim(adjustl(runname))//"_box_co.out")
                end if
            end if

            ! Passing resonance (counter-passing)
            if (.not. nopassing) then
                sigv = -1
                call get_passing_region(etamin, etamax)
                call transport_integral_mid(vminp, vmaxp, Dresctr, Tresctr)
                Dctr = Dctr + Dresctr
                Tctr = Tctr + Tresctr
                if (orbit_mode_transp > 0) then
                    call write_transport_box(trim(adjustl(runname))//"_box_ctr.out")
                end if
            end if

            ! Trapped resonance (trapped)
            sigv = 1
            call get_trapped_region(etamin, etamax)
            call transport_integral_mid(vmint, vmaxt, Drest, Trest)
            Dt = Dt + Drest
            Tt = Tt + Trest
            if (orbit_mode_transp > 0) then
                call write_transport_box(trim(adjustl(runname))//"_box_t.out")
            end if
        end if

        print *, ""
        print *, "test_flux: Mt = ", M_t, ", mth = ", mth
        write (*, "(4ES12.2,2F12.2)") Dresco(1), Dresctr(1), &
            Drest(1), Dresco(1) + Dresctr(1) + Drest(1), &
            vminp/vth, vmint/vth
        write (*, "(4ES12.2,2F12.2)") Dresco(2), Dresctr(2), &
            Drest(2), Dresco(2) + Dresctr(2) + Drest(2), &
            vmaxp/vth, vmaxt/vth

        open (unit=10, file=trim(adjustl(runname))//"_integral.out", recl=1024, position="append")
        write (10, *) M_t, mth, Dresco(1), Dresctr(1), Drest(1), &
            Dresco(1) + Dresctr(1) + Drest(1), &
            Dresco(2), Dresctr(2), Drest(2), &
            Dresco(2) + Dresctr(2) + Drest(2), vminp/vth, vmaxp/vth, &
            vmint/vth, vmaxt/vth
        close (unit=10)

        print *, ""
        print *, "compute_torque: Mt = ", M_t, ", mth = ", mth
        write (*, "(4ES12.2,2F12.2)") Tresco, Tresctr, &
            Trest, Tresco + Tresctr + Trest

        open (unit=10, file=trim(adjustl(runname))//"_torque_integral.out", recl=1024, position="append")
        write (10, *) mth, Tresco, Tresctr, Trest
        close (unit=10)
    end subroutine compute_transport_harmonic

    subroutine write_transport_box(filename)
        implicit none
        character(len=*), intent(in) :: filename

        open (unit=9, file=filename, recl=8192, position="append")
        write (9, *) mth, fluxint_box(1, :), fluxint_box(2, :)
        close (unit=9)
    end subroutine write_transport_box

    end module neort
