module neort
    use iso_fortran_env, only: dp => real64
    use logger, only: debug, set_log_level, get_log_level, log_result, LOG_INFO
    use neort_datatypes, only: magfie_data_t, transport_data_t, transport_harmonic_t
    use neort_profiles, only: init_thermodynamic_forces, init_profiles, vth, dvthds, ni1, &
                              dni1ds, Ti1, dTi1ds, qi, mi, mu, qe
    use neort_magfie, only: init_flux_surface_average
    use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline
    use neort_transport, only: compute_transport_integral
    use do_magfie_mod, only: s

    implicit none

    ! Number of integration steps in v, set 0 for adaptive integration by quadpack
    integer :: vsteps = 256

    ! Maximum absolute poloidal/orbit resonance harmonic. A negative value
    ! preserves the historical q-dependent automatic range.
    integer :: mth_max_abs = -1

    ! Upper velocity-space cutoff in units of the thermal velocity. Default 4.0
    ! captures the far-tail resonance; the historical hard-coded 3.0 bound
    ! truncated it. Set vmax_over_vth = 3.0 to reproduce pre-2026-07-20 results.
    real(dp) :: vmax_over_vth = 4.0_dp

    ! Resonance components of a harmonic: co-passing, counter-passing, trapped
    integer, parameter :: component_co = 1, component_ctr = 2, component_trapped = 3
    integer, parameter :: ncomponent = 3

contains

    pure subroutine harmonic_bounds(mph_value, q_value, max_abs, mth_min, mth_max)
        integer, intent(in) :: mph_value, max_abs
        real(dp), intent(in) :: q_value
        integer, intent(out) :: mth_min, mth_max

        if (max_abs >= 0) then
            mth_min = -max_abs
            mth_max = max_abs
        else
            mth_min = -ceiling(2*abs(mph_value*q_value))
            mth_max = ceiling(2*abs(mph_value*q_value))
        end if
    end subroutine harmonic_bounds

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

        real(dp), intent(out) :: eta_min, eta_max

        eta_min = (1 + epst)*etatp
        eta_max = (1 - epst)*etadt
    end subroutine set_to_trapped_region

    pure subroutine set_to_passing_region(eta_min, eta_max)
        use driftorbit, only: epsp, etatp

        real(dp), intent(out) :: eta_min, eta_max

        eta_min = epsp*etatp
        eta_max = (1 - epsp)*etatp
    end subroutine set_to_passing_region

    subroutine check_magfie(data)
        use util, only: pi, c, imun, eV
        use do_magfie_mod, only: do_magfie, sign_theta, R0, a, psi_pr, iota, q, Bthcov, Bphcov, &
            dBthcovds, dBphcovds, eps
        use do_magfie_pert_mod, only: do_magfie_pert_amp, mph
        use driftorbit, only: B0, etatp, etadt, M_t, Om_tE, m0, pertfile, nonlin, &
                              epsmn, pertfile_scale, dVds

        type(magfie_data_t), intent(out) :: data

        integer, parameter :: nth = 50
        integer :: k
        real(dp) :: thmin, thmax
        real(dp) :: bmod, sqrtg, x(3), hder(3), hcovar(3), hctrvr(3), hcurl(3)
        real(dp) :: D_plateau, Drp
        real(dp) :: dpp, dhh, fpeff
        complex(dp) :: bn

        D_plateau = pi * vth**3 / (16.0_dp * R0 * iota * (qi * B0 / (mi * c))**2)
        Drp = 4 * mph * q / (eps**2 * sqrt(pi))  ! actually this is Drp/Dp

        thmin = -pi
        thmax = pi
        x(1) = s
        x(2) = 0.0_dp
        x(3) = 0.0_dp

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
        data%params%Om_tBref = c * mi * vth**2 / (2.0_dp * qi * sign_theta * psi_pr)
        data%params%vth = vth
        data%params%T_in_eV = mi / 2.0_dp * vth**2 / eV
        data%params%m0 = 1.0_dp * m0
        data%params%n0 = 1.0_dp * mph
        data%params%Dp = D_plateau
        data%params%Drp_over_Dp = Drp
        data%params%etatp = etatp
        data%params%etadt = etadt
        data%params%pertfile = pertfile
        data%params%nonlin = nonlin
        dpp = 0.0_dp
        dhh = 0.0_dp
        fpeff = 0.0_dp
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
                bn = pertfile_scale*bn/bmod
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

    subroutine prepare_transport(mthmin, mthmax)
        ! Per-species state shared by all harmonics and the harmonic range.
        ! Every thread computing a resonance component of this species must
        ! hold it (see compute_transport_component).
        use driftorbit, only: M_t, R0, mph, dOm_tEds, dM_tds, supban
        use neort_profiles, only: Om_tE
        use do_magfie_mod, only: q

        integer, intent(out) :: mthmin, mthmax

        Om_tE = vth*M_t/R0
        dOm_tEds = vth*dM_tds/R0 + M_t*dvthds/R0

        call harmonic_bounds(mph, q, mth_max_abs, mthmin, mthmax)

        if (supban) then
            ! Superbanana plateau is the trapped ell=0 resonance only.
            mthmin = 0
            mthmax = 0
        end if
    end subroutine prepare_transport

    logical function component_active(component)
        ! Superbanana plateau and nopassing runs only have the trapped component.
        use driftorbit, only: nopassing, supban
        integer, intent(in) :: component
        if (component == component_trapped) then
            component_active = .true.
        else
            component_active = .not. (supban .or. nopassing)
        end if
    end function component_active

    subroutine compute_transport(result_)
        type(transport_data_t), intent(out) :: result_

        integer :: j, component, mthmin, mthmax
        real(dp), allocatable :: Dres(:, :, :), Tres(:, :)

        call debug('compute_transport')

        call prepare_transport(mthmin, mthmax)
        allocate(Dres(2, ncomponent, mthmin:mthmax), Tres(ncomponent, mthmin:mthmax))
        do j = mthmin, mthmax
            do component = 1, ncomponent
                call compute_transport_component(j, component, Dres(:, component, j), &
                    Tres(component, j))
            end do
        end do
        call assemble_transport(mthmin, mthmax, Dres, Tres, result_)

        call debug('compute_transport complete')
    end subroutine compute_transport

    subroutine compute_transport_harmonic(j, Dco, Dctr, Dt, Tco, Tctr, Tt, harmonic)
        ! All components of harmonic j, added to the running sums.
        integer, intent(in) :: j
        real(dp), intent(inout) :: Dco(2), Dctr(2), Dt(2), Tco, Tctr, Tt
        type(transport_harmonic_t), intent(out) :: harmonic

        real(dp) :: Dres(2, ncomponent), Tres(ncomponent)
        integer :: component

        do component = 1, ncomponent
            call compute_transport_component(j, component, Dres(:, component), &
                Tres(component))
        end do
        call assemble_harmonic(j, Dres, Tres, Dco, Dctr, Dt, Tco, Tctr, Tt, harmonic)
    end subroutine compute_transport_harmonic

    subroutine compute_transport_component(j, component, D, T)
        ! Resonance integral of harmonic j for one orbit class. Depends only on
        ! the per-species state of prepare_transport, so the components of all
        ! harmonics may be computed in any order or on any thread holding it.
        use driftorbit, only: mth, etamin, etamax, sign_vpar

        integer, intent(in) :: j, component
        real(dp), intent(out) :: D(2), T

        real(dp) :: vmin, vmax
        character(len=256) :: buffer

        D = 0.0_dp
        T = 0.0_dp
        if (.not. component_active(component)) return

        write(buffer, '(A,I0,A,I0)') "compute_transport_component: mth = ", j, &
            ", component = ", component
        call debug(buffer)

        mth = j
        call velocity_bounds(component, vmin, vmax)
        select case (component)
        case (component_co)
            sign_vpar = 1
            call set_to_passing_region(etamin, etamax)
        case (component_ctr)
            sign_vpar = -1
            call set_to_passing_region(etamin, etamax)
        case default
            sign_vpar = 1
            call set_to_trapped_region(etamin, etamax)
        end select
        call compute_transport_integral(vmin, vmax, vsteps, D, T)
    end subroutine compute_transport_component

    subroutine velocity_bounds(component, vmin, vmax)
        use driftorbit, only: supban

        integer, intent(in) :: component
        real(dp), intent(out) :: vmin, vmax

        vmin = 1.0e-6_dp * vth
        vmax = vmax_over_vth * vth
        if (supban .and. component == component_trapped) then
            ! Shaing superbanana plateau: trapped ell=0 resonance only, with a
            ! wide velocity bracket so the analytic Omph = 0 resonance is
            ! captured across the thermal range.
            vmin = 0.01_dp * vth
            vmax = 5.0_dp * vth
        end if
    end subroutine velocity_bounds

    subroutine assemble_transport(mthmin, mthmax, Dres, Tres, result_)
        ! Sums the component integrals of all harmonics in harmonic order and
        ! leaves the state where the harmonic loop historically ended.
        use driftorbit, only: mth, etamin, etamax, sign_vpar, comptorque, dVds, M_t

        integer, intent(in) :: mthmin, mthmax
        real(dp), intent(in) :: Dres(:, :, mthmin:), Tres(:, mthmin:)
        type(transport_data_t), intent(out) :: result_

        integer :: j, idx

        ! Transport coefficients D11, D12 in approximate effective radius r=2.0_dp/a*sqrt(s)
        real(dp) :: Dco(2), Dctr(2), Dt(2)

        ! Torque density dTphi_int/ds for integration over normalized toroidal flux s
        real(dp) :: Tco, Tctr, Tt

        if (allocated(result_%harmonics)) deallocate(result_%harmonics)

        Dco = 0.0_dp
        Dctr = 0.0_dp
        Dt = 0.0_dp
        Tco = 0.0_dp
        Tctr = 0.0_dp
        Tt = 0.0_dp

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
            call assemble_harmonic(j, Dres(:, :, j), Tres(:, j), Dco, Dctr, Dt, &
                Tco, Tctr, Tt, result_%harmonics(idx))
        end do
        if (mthmax >= mthmin) then
            mth = mthmax
            sign_vpar = 1
            call set_to_trapped_region(etamin, etamax)
        end if

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
            result_%torque%s = 0.0_dp
            result_%torque%dVds = 0.0_dp
            result_%torque%M_t = 0.0_dp
            result_%torque%Tco = 0.0_dp
            result_%torque%Tctr = 0.0_dp
            result_%torque%Tt = 0.0_dp
        end if
    end subroutine assemble_transport

    subroutine assemble_harmonic(j, Dres, Tres, Dco, Dctr, Dt, Tco, Tctr, Tt, harmonic)
        use driftorbit, only: M_t, supban

        integer, intent(in) :: j
        real(dp), intent(in) :: Dres(2, ncomponent), Tres(ncomponent)
        real(dp), intent(inout) :: Dco(2), Dctr(2), Dt(2), Tco, Tctr, Tt
        type(transport_harmonic_t), intent(out) :: harmonic

        real(dp) :: vminp, vmaxp, vmint, vmaxt
        character(len=256) :: buffer

        call velocity_bounds(component_co, vminp, vmaxp)
        call velocity_bounds(component_trapped, vmint, vmaxt)

        if (component_active(component_co)) then
            Dco = Dco + Dres(:, component_co)
            Tco = Tco + Tres(component_co)
        end if
        if (component_active(component_ctr)) then
            Dctr = Dctr + Dres(:, component_ctr)
            Tctr = Tctr + Tres(component_ctr)
        end if
        Dt = Dt + Dres(:, component_trapped)
        Tt = Tt + Tres(component_trapped)

        harmonic%mth = j
        harmonic%Dresco = Dres(:, component_co)
        harmonic%Dresctr = Dres(:, component_ctr)
        harmonic%Drest = Dres(:, component_trapped)
        harmonic%Tresco = Tres(component_co)
        harmonic%Tresctr = Tres(component_ctr)
        harmonic%Trest = Tres(component_trapped)
        harmonic%vminp_over_vth = vminp/vth
        harmonic%vmaxp_over_vth = vmaxp/vth
        harmonic%vmint_over_vth = vmint/vth
        harmonic%vmaxt_over_vth = vmaxt/vth

        if (supban) return

        call log_result("")
        write(buffer, '(A,ES12.5,A,I0)') "test_flux: Mt = ", M_t, ", mth = ", j
        call log_result(buffer)
        write(buffer, '(4ES12.2,2F12.2)') Dres(1, :), &
            Dres(1, component_co) + Dres(1, component_ctr) &
            + Dres(1, component_trapped), vminp/vth, vmint/vth
        call log_result(buffer)
        write(buffer, '(4ES12.2,2F12.2)') Dres(2, :), &
            Dres(2, component_co) + Dres(2, component_ctr) &
            + Dres(2, component_trapped), vmaxp/vth, vmaxt/vth
        call log_result(buffer)

        call log_result("")
        write(buffer, '(A,ES12.5,A,I0)') "compute_torque: Mt = ", M_t, ", mth = ", j
        call log_result(buffer)
        write(buffer, '(4ES12.2)') Tres(component_co), Tres(component_ctr), &
            Tres(component_trapped), &
            Tres(component_co) + Tres(component_ctr) + Tres(component_trapped)
        call log_result(buffer)
    end subroutine assemble_harmonic

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
        write (unit, *) "check_magfie: A         = ", 1.0_dp / data%params%eps
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
        real(dp) :: total_D1, total_D2
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
