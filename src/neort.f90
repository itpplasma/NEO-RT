module neort
    use logger, only: debug, set_log_level, get_log_level, log_result, LOG_INFO
    use neort_datatypes, only: magfie_data_t, transport_data_t, transport_harmonic_t
    use neort_profiles, only: read_and_init_profile_input, read_and_init_plasma_input, &
                              init_thermodynamic_forces, init_profiles, vth, dvthds, ni1, &
                              dni1ds, Ti1, dTi1ds, qi, mi, mu, qe
    use neort_magfie, only: init_flux_surface_average
    use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline
    use neort_transport, only: compute_transport_integral
    use do_magfie_mod, only: s

    implicit none

    ! Number of integration steps in v, set 0 for adaptive integration by quadpack
    integer :: vsteps = 256

contains

    subroutine read_and_set_control(base_path)
        ! set global control parameters directly from file
        ! TODO: warn if unused values are non-zero
        use do_magfie_mod, only: s, bfac, inp_swi
        use do_magfie_pert_mod, only: mph, set_mph
        use driftorbit, only: epsmn, m0, comptorque, magdrift, nopassing, pertfile, nonlin, efac
        use logger, only: set_log_level
        use neort_orbit, only: noshear
        use neort_profiles, only: M_t, vth

        character(len=*), intent(in) :: base_path
        real(8) :: qs, ms
        integer :: log_level

        namelist /params/ s, M_t, qs, ms, vth, epsmn, m0, mph, comptorque, magdrift, &
            nopassing, noshear, pertfile, nonlin, bfac, efac, inp_swi, vsteps, log_level

        open (unit=9, file=trim(adjustl(base_path))//".in", status="old", form="formatted")
        read (9, nml=params)
        close (unit=9)

        M_t = M_t * efac / bfac
        qi = qs * qe
        mi = ms * mu
        call set_mph(mph)
        call set_log_level(log_level)
    end subroutine read_and_set_control

    subroutine init
        use do_magfie_mod, only: psi_pr, q
        use driftorbit, only: nopassing, sign_vpar, etamin, etamax, comptorque

        call debug('init')
        call init_flux_surface_average(s)
        call init_canon_freq_trapped_spline
        if (.not. nopassing) call init_canon_freq_passing_spline
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax)
        if (comptorque) call init_thermodynamic_forces(psi_pr, q)
        call debug('init complete')
    end subroutine init

    pure subroutine set_to_trapped_region(eta_min, eta_max)
        use driftorbit, only: epst, etatp, etadt

        real(8), intent(out) :: eta_min, eta_max

        eta_min = (1 + epst)*etatp
        eta_max = (1 - epst)*etadt
    end subroutine set_to_trapped_region

    pure subroutine set_to_passing_region(eta_min, eta_max)
        use driftorbit, only: epsp, etatp

        real(8), intent(out) :: eta_min, eta_max

        eta_min = epsp*etatp
        eta_max = (1 - epsp)*etatp
    end subroutine set_to_passing_region

    subroutine check_magfie(data)
        use do_magfie_mod, only: do_magfie, sign_theta
        use do_magfie_pert_mod, only: do_magfie_pert_amp
        use driftorbit

        type(magfie_data_t), intent(out) :: data

        integer, parameter :: nth = 50
        integer :: k
        real(8) :: thmin, thmax
        real(8) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(8) :: Dp, Drp
        real(8) :: dpp, dhh, fpeff
        complex(8) :: bn

        Dp = pi*vth**3/(16d0*R0*iota*(qi*B0/(mi*c))**2)
        Drp = 4*mph*q/(eps**2*sqrt(pi))  ! actually this is Drp/Dp

        thmin = -pi
        thmax = pi
        x(1) = s
        x(2) = 0d0
        x(3) = 0d0

        data%params%s = s
        data%params%R0 = R0
        data%params%a = a
        data%params%eps = eps
        data%params%psi_pr = psi_pr
        data%params%B0 = B0
        data%params%Bthcov = Bthcov
        data%params%Bphcov = Bphcov
        data%params%dBthcovds = dBthcovds
        data%params%dBphcovds = dBphcovds
        data%params%q = q
        data%params%iota = iota
        data%params%dVds = dVds
        data%params%M_t = M_t
        data%params%Om_tE = Om_tE
        data%params%Om_tBref = c*mi*vth**2/(2d0*qi*sign_theta*psi_pr)
        data%params%vth = vth
        data%params%T_in_eV = mi/2d0*vth**2/eV
        data%params%m0 = 1d0*m0
        data%params%n0 = 1d0*mph
        data%params%Dp = Dp
        data%params%Drp_over_Dp = Drp
        data%params%etatp = etatp
        data%params%etadt = etadt
        data%params%pertfile = pertfile
        data%params%nonlin = nonlin
        dpp = 0d0
        dhh = 0d0
        fpeff = 0d0
        data%params%dpp = dpp
        data%params%dhh = dhh
        data%params%fpeff = fpeff

        data%n_points = nth
        if (allocated(data%tensors)) deallocate(data%tensors)
        allocate(data%tensors(nth))

        do k = 1, nth
            x(3) = thmin + (k - 1)*(thmax - thmin)/(nth - 1)
            call do_magfie(x, bmod, sqrtg, hder, hcovar, hctrvr, hcurl)
            if (pertfile) then
                call do_magfie_pert_amp(x, bn)
                bn = epsmn*bn/bmod
            else
                bn = epsmn*exp(imun*m0*x(3))
            end if
            data%tensors(k)%theta = x(3)
            data%tensors(k)%bmod = bmod
            data%tensors(k)%sqrtg = sqrtg
            data%tensors(k)%hder = hder
            data%tensors(k)%hcovar = hcovar
            data%tensors(k)%hctrvr = hctrvr
            data%tensors(k)%hcurl = hcurl
            data%tensors(k)%bn = bn
            data%tensors(k)%eps_exp = epsmn*exp(imun*m0*x(3))
        end do
    end subroutine check_magfie

    subroutine compute_transport(result_)
        use driftorbit

        type(transport_data_t), intent(out) :: result_

        integer :: j, idx

        ! Transport coefficients D11, D12 in approximate effective radius r=2d0/a*sqrt(s)
        real(8) :: Dco(2), Dctr(2), Dt(2)

        ! Torque density dTphi_int/ds for integration over normalized toroidal flux s
        real(8) :: Tco, Tctr, Tt
        integer :: mthmin, mthmax

        call debug('compute_transport')

        if (allocated(result_%harmonics)) deallocate(result_%harmonics)

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

        if (mthmax < mthmin) then
            ! Edge case: no valid harmonics (upper bound below lower bound).
            ! In this situation we intentionally allocate an empty array.
            allocate(result_%harmonics(0))
        else
            allocate(result_%harmonics(mthmax - mthmin + 1))
        end if

        idx = 0
        do j = mthmin, mthmax
            idx = idx + 1
            call compute_transport_harmonic(j, Dco, Dctr, Dt, Tco, Tctr, Tt, result_%harmonics(idx))
        end do

        result_%summary%M_t = M_t
        result_%summary%Dco = Dco
        result_%summary%Dctr = Dctr
        result_%summary%Dt = Dt

        result_%torque%has_torque = comptorque
        if (comptorque) then
            result_%torque%s = s
            result_%torque%dVds = dVds
            result_%torque%M_t = M_t
            result_%torque%Tco = Tco
            result_%torque%Tctr = Tctr
            result_%torque%Tt = Tt
        else
            result_%torque%s = 0d0
            result_%torque%dVds = 0d0
            result_%torque%M_t = 0d0
            result_%torque%Tco = 0d0
            result_%torque%Tctr = 0d0
            result_%torque%Tt = 0d0
        end if
        call debug('compute_transport complete')
    end subroutine compute_transport

    subroutine compute_transport_harmonic(j, Dco, Dctr, Dt, Tco, Tctr, Tt, harmonic)
        use driftorbit

        integer, intent(in) :: j
        real(8), intent(inout) :: Dco(2), Dctr(2), Dt(2), Tco, Tctr, Tt
        type(transport_harmonic_t), intent(out) :: harmonic

        real(8) :: Dresco(2), Dresctr(2), Drest(2), Tresco, Tresctr, Trest
        real(8) :: vminp, vmaxp, vmint, vmaxt
        character(len=256) :: buffer

        write(buffer, '(A,ES12.5,A,I0)') "compute_transport_harmonic: M_t = ", M_t, ", mth = ", j
        call debug(buffer)

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
            sign_vpar = 1
            call set_to_passing_region(etamin, etamax)
            call compute_transport_integral(vminp, vmaxp, vsteps, Dresco, Tresco)
            Dco = Dco + Dresco
            Tco = Tco + Tresco
        end if

        ! Passing resonance (counter-passing)
        if (.not. nopassing) then
            sign_vpar = -1
            call set_to_passing_region(etamin, etamax)
            call compute_transport_integral(vminp, vmaxp, vsteps, Dresctr, Tresctr)
            Dctr = Dctr + Dresctr
            Tctr = Tctr + Tresctr
        end if

        ! Trapped resonance (trapped)
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax)
        call compute_transport_integral(vmint, vmaxt, vsteps, Drest, Trest)
        Dt = Dt + Drest
        Tt = Tt + Trest

        call log_result("")
        write(buffer, '(A,ES12.5,A,I0)') "test_flux: Mt = ", M_t, ", mth = ", mth
        call log_result(buffer)
        write(buffer, '(4ES12.2,2F12.2)') Dresco(1), Dresctr(1), &
            Drest(1), Dresco(1) + Dresctr(1) + Drest(1), &
            vminp/vth, vmint/vth
        call log_result(buffer)
        write(buffer, '(4ES12.2,2F12.2)') Dresco(2), Dresctr(2), &
            Drest(2), Dresco(2) + Dresctr(2) + Drest(2), &
            vmaxp/vth, vmaxt/vth
        call log_result(buffer)

        call log_result("")
        write(buffer, '(A,ES12.5,A,I0)') "compute_torque: Mt = ", M_t, ", mth = ", mth
        call log_result(buffer)
        write(buffer, '(4ES12.2)') Tresco, Tresctr, Trest, Tresco + Tresctr + Trest
        call log_result(buffer)

        harmonic%mth = mth
        harmonic%Dresco = Dresco
        harmonic%Dresctr = Dresctr
        harmonic%Drest = Drest
        harmonic%Tresco = Tresco
        harmonic%Tresctr = Tresctr
        harmonic%Trest = Trest
        harmonic%vminp_over_vth = vminp/vth
        harmonic%vmaxp_over_vth = vmaxp/vth
        harmonic%vmint_over_vth = vmint/vth
        harmonic%vmaxt_over_vth = vmaxt/vth

        call debug('compute_transport_harmonic complete')
    end subroutine compute_transport_harmonic

    subroutine write_magfie_data_to_files(data, base_path)
        type(magfie_data_t), intent(in) :: data
        character(len=*), intent(in) :: base_path

        integer :: k
        integer, parameter :: unit = 9

        open (unit=unit, file=trim(adjustl(base_path))//"_magfie_param.out", recl=1024)
        write (unit, *) "-------------------------"
        write (unit, *) "check_magfie: s         = ", data%params%s
        write (unit, *) "check_magfie: R0        = ", data%params%R0
        write (unit, *) "check_magfie: a         = ", data%params%a
        write (unit, *) "check_magfie: eps       = ", data%params%eps
        write (unit, *) "check_magfie: A         = ", 1.0d0 / data%params%eps
        write (unit, *) "check_magfie: psi_pr    = ", data%params%psi_pr
        write (unit, *) "check_magfie: B0        = ", data%params%B0
        write (unit, *) "check_magfie: Bthcov    = ", data%params%Bthcov
        write (unit, *) "check_magfie: Bphcov    = ", data%params%Bphcov
        write (unit, *) "check_magfie: dBthcovds = ", data%params%dBthcovds
        write (unit, *) "check_magfie: dBphcovds = ", data%params%dBphcovds
        write (unit, *) "check_magfie: q         = ", data%params%q
        write (unit, *) "check_magfie: iota      = ", data%params%iota
        write (unit, *) "check_magfie: dVds      = ", data%params%dVds
        write (unit, *) "check_magfie: M_t       = ", data%params%M_t
        write (unit, *) "check_magfie: Om_tE     = ", data%params%Om_tE
        write (unit, *) "check_magfie: Om_tBref  = ", data%params%Om_tBref
        write (unit, *) "check_magfie: vth       = ", data%params%vth
        write (unit, *) "check_magfie: T [eV]    = ", data%params%T_in_eV
        write (unit, *) "check_magfie: m0        = ", data%params%m0
        write (unit, *) "check_magfie: n0        = ", data%params%n0
        write (unit, *) "check_magfie: Dp        = ", data%params%Dp
        write (unit, *) "check_magfie: Drp/Dp    = ", data%params%Drp_over_Dp
        write (unit, *) "check_magfie: etatp     = ", data%params%etatp
        write (unit, *) "check_magfie: etadt     = ", data%params%etadt
        write (unit, *) "-------------------------"
        write (unit, *) "check_magfie: pertfile  = ", data%params%pertfile
        if (data%params%nonlin) then
            write (unit, *) "-------------------------"
            write (unit, *) "check_magfie: dpp       = ", data%params%dpp
            write (unit, *) "check_magfie: dhh       = ", data%params%dhh
            write (unit, *) "check_magfie: dfpeff    = ", data%params%fpeff
            write (unit, *) "-------------------------"
        end if
        close (unit=unit)

        open (unit=unit, file=trim(adjustl(base_path))//"_magfie.out", recl=1024)
        do k = 1, data%n_points
            write (unit, *) data%tensors(k)%theta, data%tensors(k)%bmod, data%tensors(k)%sqrtg, &
                            data%tensors(k)%hder(1), data%tensors(k)%hder(2), &
                            data%tensors(k)%hder(3), data%tensors(k)%hcovar(1), &
                            data%tensors(k)%hcovar(2), data%tensors(k)%hcovar(3), &
                            data%tensors(k)%hctrvr(1), data%tensors(k)%hctrvr(2), &
                            data%tensors(k)%hctrvr(3), data%tensors(k)%hcurl(1), &
                            data%tensors(k)%hcurl(2), data%tensors(k)%hcurl(3), &
                            real(data%tensors(k)%bn), aimag(data%tensors(k)%bn), &
                            real(data%tensors(k)%eps_exp), aimag(data%tensors(k)%eps_exp)
        end do
        close (unit=unit)
    end subroutine write_magfie_data_to_files

    subroutine write_transport_data_to_files(data, base_path)
        type(transport_data_t), intent(in) :: data
        character(len=*), intent(in) :: base_path

        integer :: k
        real(8) :: total_D1, total_D2
        integer, parameter :: unit1 = 9
        integer, parameter :: unit2 = 10

        open (unit=unit1, file=trim(adjustl(base_path))//".out", recl=1024)
        write (unit1, *) "# M_t D11co D11ctr D11t D11 D12co D12ctr D12t D12"
        total_D1 = data%summary%Dco(1) + data%summary%Dctr(1) + data%summary%Dt(1)
        total_D2 = data%summary%Dco(2) + data%summary%Dctr(2) + data%summary%Dt(2)
        write (unit1, *) data%summary%M_t, data%summary%Dco(1), data%summary%Dctr(1), &
                        data%summary%Dt(1), total_D1, data%summary%Dco(2), data%summary%Dctr(2), &
                        data%summary%Dt(2), total_D2
        close (unit=unit1)

        if (data%torque%has_torque) then
            open (unit=unit1, file=trim(adjustl(base_path))//"_torque.out", recl=1024)
            write (unit1, *) "# s dVds M_t Tco Tctr Tt"
            write (unit1, *) data%torque%s, data%torque%dVds, data%torque%M_t, data%torque%Tco, &
                             data%torque%Tctr, data%torque%Tt
            close (unit=unit1)
        end if

        open (unit=unit1, file=trim(adjustl(base_path))//"_integral.out", recl=1024)
        open (unit=unit2, file=trim(adjustl(base_path))//"_torque_integral.out", recl=1024)
        do k = 1, size(data%harmonics)
            total_D1 = data%harmonics(k)%Dresco(1) + data%harmonics(k)%Dresctr(1) + &
                       data%harmonics(k)%Drest(1)
            total_D2 = data%harmonics(k)%Dresco(2) + data%harmonics(k)%Dresctr(2) + &
                       data%harmonics(k)%Drest(2)
            write (unit1, *) data%summary%M_t, data%harmonics(k)%mth, data%harmonics(k)%Dresco(1), &
                             data%harmonics(k)%Dresctr(1), data%harmonics(k)%Drest(1), &
                             total_D1, data%harmonics(k)%Dresco(2), data%harmonics(k)%Dresctr(2), &
                             data%harmonics(k)%Drest(2), total_D2, &
                             data%harmonics(k)%vminp_over_vth, data%harmonics(k)%vmaxp_over_vth, &
                             data%harmonics(k)%vmint_over_vth, data%harmonics(k)%vmaxt_over_vth

            write (unit2, *) data%harmonics(k)%mth, data%harmonics(k)%Tresco, &
                             data%harmonics(k)%Tresctr, data%harmonics(k)%Trest
        end do
        close (unit=unit1)
        close (unit=unit2)
    end subroutine write_transport_data_to_files

end module neort
