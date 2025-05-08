module neort
    use neort_profiles, only: init_profile_input, init_plasma_input, &
        init_thermodynamic_forces, init_profiles, vth, dvthds, ni1, dni1ds, Ti1, &
        dTi1ds, qi, mi, mu, qe
    use neort_magfie, only: init_fsa, init_misc
    use neort_orbit, only: noshear
    use neort_freq, only: init_Om_spl, init_Om_pass_spl
    use neort_transport, only: compute_transport_integral
    use driftorbit
    implicit none

    character(1024) :: runname

    ! Number of integration steps in v, set 0 for adaptive integration by quadpack
    integer :: vsteps = 256

contains

    subroutine main
        use do_magfie_mod, only: do_magfie_init
        use do_magfie_pert_mod, only: do_magfie_pert_init
        use driftorbit, only: pertfile

        logical :: file_exists

        call get_command_argument(1, runname)
        call read_control
        call do_magfie_init  ! init axisymmetric part of field from infile
        if (pertfile) call do_magfie_pert_init ! else epsmn*exp(imun*(m0*th + mph*ph))
        call init_profiles(R0)

        inquire(file="plasma.in", exist=file_exists)
        if (file_exists) then
            call init_plasma_input(s)
        else if (nonlin .or. comptorque) then
            error stop "plasma.in required for nonlin or comptorque"
        end if

        inquire(file="profile.in", exist=file_exists)
        if (file_exists) then
            call init_profile_input(s, R0, efac, bfac)
        else if (nonlin) then
            error stop "profile.in required for nonlin"
        end if

        call init
        call check_magfie
        call compute_transport
    end subroutine main

    subroutine read_control
        use driftorbit
        real(8) :: qs, ms

        namelist /params/ s, M_t, qs, ms, vth, epsmn, m0, mph, comptorque, magdrift, &
            nopassing, noshear, pertfile, nonlin, bfac, efac, inp_swi, vsteps

        open (unit=9, file=trim(adjustl(runname))//".in", status="old", form="formatted")
        read (9, nml=params)
        close (unit=9)

        M_t = M_t*efac/bfac
        qi = qs*qe
        mi = ms*mu
    end subroutine read_control

    subroutine init
        init_done = .false.
        call init_fsa
        call init_misc
        call init_Om_spl       ! frequencies of trapped orbits
        if (.not. nopassing) call init_Om_pass_spl  ! frequencies of passing orbits
        sigv = 1
        call set_to_trapped_region(etamin, etamax)
        if (comptorque) call init_thermodynamic_forces(psi_pr, q)
        init_done = .true.
    end subroutine init

    pure subroutine set_to_trapped_region(eta_min, eta_max)
        real(8), intent(out) :: eta_min, eta_max
        eta_min = (1 + epst)*etatp
        eta_max = (1 - epst)*etadt
    end subroutine set_to_trapped_region

    pure subroutine set_to_passing_region(eta_min, eta_max)
        real(8), intent(out) :: eta_min, eta_max
        eta_min = epsp*etatp
        eta_max = (1 - epsp)*etatp
    end subroutine set_to_passing_region

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

        Dco = 0d0
        Dctr = 0d0
        Dt = 0d0
        Tco = 0d0
        Tctr = 0d0
        Tt = 0d0

        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0 + M_t*dvthds/R0

        mthmin = -ceiling(2*abs(mph*q))
        mthmax = ceiling(2*abs(mph*q))

        do j = mthmin, mthmax
            call compute_transport_harmonic(j, Dco, Dctr, Dt, Tco, Tctr, Tt)
        end do

        open (unit=9, file=trim(adjustl(runname))//".out", recl=1024)
        write (9, *) "# M_t D11t D11ctr D11co D11 D12t D12ctr D12co D12"
        write (9, *) M_t, Dco(1), Dctr(1), Dt(1), Dco(1) + Dctr(1) + Dt(1), &
            Dco(2), Dctr(2), Dt(2), Dco(2) + Dctr(2) + Dt(2)
        close (unit=9)

        if (comptorque) then
            open (unit=9, file=trim(adjustl(runname))//"_torque.out", recl=1024)
            write (9, *) "# s dVds M_t Tco Tctr Tt"
            write (9, *) s, dVds, M_t, Tco, Tctr, Tt
            close (unit=9)
        end if
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

        ! Passing resonance (co-passing)
        if (.not. nopassing) then
            sigv = 1
            call set_to_passing_region(etamin, etamax)
            call compute_transport_integral(vminp, vmaxp, vsteps, Dresco, Tresco)
            Dco = Dco + Dresco
            Tco = Tco + Tresco
        end if

        ! Passing resonance (counter-passing)
        if (.not. nopassing) then
            sigv = -1
            call set_to_passing_region(etamin, etamax)
            call compute_transport_integral(vminp, vmaxp, vsteps, Dresctr, Tresctr)
            Dctr = Dctr + Dresctr
            Tctr = Tctr + Tresctr
        end if

        ! Trapped resonance (trapped)
        sigv = 1
        call set_to_trapped_region(etamin, etamax)
        call compute_transport_integral(vmint, vmaxt, vsteps, Drest, Trest)
        Dt = Dt + Drest
        Tt = Tt + Trest

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

    end module neort
