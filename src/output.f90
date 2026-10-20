module neort_output
    ! Consolidated run output. The HDF5 backend replaces the six plain-text
    ! files with a single self-describing <runname>.h5 carrying units, the
    ! resolved run configuration, and the harmonic- and theta-resolved data
    ! that the text pipeline used to drop on the floor.
    use iso_fortran_env, only: dp => real64
    use neort_datatypes, only: magfie_data_t, transport_data_t, mixture_result_t, species_result_t
    use neort_species, only: species_state_t

    implicit none

    private
    public :: write_output, write_mixture_output, FORMAT_VERSION

    ! Bump when the group/dataset layout changes incompatibly.
    integer, parameter :: FORMAT_VERSION = 1

contains

    subroutine write_output(magfie_data, transport_data, base_path, output_format)
        use neort, only: write_magfie_data_to_files, write_transport_data_to_files

        type(magfie_data_t), intent(in) :: magfie_data
        type(transport_data_t), intent(in) :: transport_data
        character(len=*), intent(in) :: base_path, output_format

        select case (trim(output_format))
        case ("hdf5")
            call write_hdf5(magfie_data, transport_data, base_path)
        case ("text")
            call write_magfie_data_to_files(magfie_data, base_path)
            call write_transport_data_to_files(transport_data, base_path)
        case ("both")
            call write_hdf5(magfie_data, transport_data, base_path)
            call write_magfie_data_to_files(magfie_data, base_path)
            call write_transport_data_to_files(transport_data, base_path)
        case default
            error stop "output_format must be 'hdf5', 'text' or 'both'"
        end select
    end subroutine write_output

    subroutine write_hdf5(magfie_data, transport_data, base_path)
        use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5_create, h5_close, &
                              h5_define_group, h5_close_group, h5_add

        type(magfie_data_t), intent(in) :: magfie_data
        type(transport_data_t), intent(in) :: transport_data
        character(len=*), intent(in) :: base_path

        integer(HID_T) :: file_id, group_id

        call h5_init()
        call h5_create(trim(adjustl(base_path))//".h5", file_id)

        call h5_add(file_id, "format_version", FORMAT_VERSION, &
                    comment="NEO-RT output layout version")
        call h5_add(file_id, "created", timestamp(), comment="ISO 8601 local time")

        call h5_define_group(file_id, "config", group_id)
        call write_config_group(group_id)
        call h5_close_group(group_id)

        call h5_define_group(file_id, "magfie", group_id)
        call write_magfie_group(group_id, magfie_data)
        call h5_close_group(group_id)

        call h5_define_group(file_id, "transport", group_id)
        call write_transport_group(group_id, transport_data)
        call h5_close_group(group_id)

        call h5_define_group(file_id, "torque", group_id)
        call write_torque_group(group_id, transport_data)
        call h5_close_group(group_id)

        call h5_define_group(file_id, "harmonics", group_id)
        call write_harmonics_group(group_id, transport_data)
        call h5_close_group(group_id)

        call h5_close(file_id)
        call h5_deinit()
    end subroutine write_hdf5

    subroutine write_config_group(group_id, species)
        ! Resolved run settings, i.e. what the solver actually used after
        ! defaults and efac/bfac scaling were applied.
        use hdf5_tools, only: HID_T, h5_add
        use do_magfie_mod, only: bfac, inp_swi
        use do_magfie_pert_mod, only: mph, inp_swi_pert
        use driftorbit, only: epsmn, m0, comptorque, magdrift, magdrift_passing, &
                              nopassing, pertfile, nonlin, efac, supban, drift_phase, &
                              fow_radial, pertfile_curvature
        use neort, only: vsteps, mth_max_abs, vmax_over_vth, monoenergetic_x
        use neort_orbit, only: noshear
        use util, only: qe, mu, qi, mi

        integer(HID_T), intent(in) :: group_id
        type(species_state_t), intent(in), optional :: species

        call h5_add(group_id, "epsmn", epsmn, comment="perturbation amplitude B1/B0")
        call h5_add(group_id, "m0", m0, comment="poloidal perturbation mode")
        call h5_add(group_id, "mph", mph, comment="toroidal perturbation mode")
        if (present(species)) then
            call h5_add(group_id, "qs", species%charge_e, unit="1", &
                comment="charge / elementary charge")
            call h5_add(group_id, "ms", species%mass_kg*1.0e3_dp/mu, unit="1", &
                comment="mass / atomic mass unit")
        else
            call h5_add(group_id, "qs", qi/qe, unit="1", comment="charge / elementary charge")
            call h5_add(group_id, "ms", mi/mu, unit="1", comment="mass / atomic mass unit")
        end if
        call h5_add(group_id, "comptorque", comptorque)
        call h5_add(group_id, "supban", supban, comment="Shaing superbanana plateau only")
        call h5_add(group_id, "drift_phase", drift_phase, &
                    comment="oscillating drift phase kept in orbit integral")
        call h5_add(group_id, "fow_radial", fow_radial, &
                    comment="perturbation sampled at first-order orbit radius")
        call h5_add(group_id, "magdrift", magdrift)
        call h5_add(group_id, "magdrift_passing", magdrift_passing)
        call h5_add(group_id, "nopassing", nopassing)
        call h5_add(group_id, "noshear", noshear)
        call h5_add(group_id, "pertfile", pertfile)
        call h5_add(group_id, "pertfile_curvature", pertfile_curvature, &
                    comment="general ideal drive with kappa.xi_perp table")
        call h5_add(group_id, "nonlin", nonlin)
        call h5_add(group_id, "bfac", bfac, comment="B field scaling factor")
        call h5_add(group_id, "efac", efac, comment="E field scaling factor")
        call h5_add(group_id, "inp_swi", inp_swi, comment="Boozer input switch")
        call h5_add(group_id, "inp_swi_pert", inp_swi_pert, &
                    comment="Boozer input switch for perturbation")
        call h5_add(group_id, "vsteps", vsteps, comment="velocity integration steps")
        call h5_add(group_id, "mth_max_abs", mth_max_abs, &
                    comment="max |mth|; negative means q-dependent range")
        call h5_add(group_id, "vmax_over_vth", vmax_over_vth, unit="1", &
                    comment="upper velocity cutoff / thermal velocity")
        call write_krook_config_group(group_id)
        call h5_add(group_id, "monoenergetic_x", monoenergetic_x, unit="1", &
                    comment="E/T; -1 selects Maxwellian integration")
    end subroutine write_config_group

    subroutine write_krook_config_group(group_id)
        use neort_mixture_collisions, only: krook_rate_source, krook_rate_scale
        use hdf5_tools, only: HID_T, h5_add
        use neort_profiles, only: profile_species
        use neort_krook, only: krook_collisions, krook_nu, krook_model, krook_harmonic, &
            krook_geometry_trapped, krook_geometry_passing, krook_eta_steps, &
            krook_max_depth, krook_rtol, krook_atol
        integer(HID_T), intent(in) :: group_id
        call h5_add(group_id, "profile_species", profile_species)
        call h5_add(group_id, "krook_collisions", krook_collisions)
        call h5_add(group_id, "krook_relaxed_variable", 'nonadiabatic_g', &
            comment='g=delta_f-delta_H*dF0/dH; distinct from plain delta_f damping')
        call h5_add(group_id, "krook_nu", krook_nu, unit="s^-1", &
            comment="supplied Krook damping rate; no 2pi conversion")
        call h5_add(group_id, "krook_rate_source", krook_rate_source)
        call h5_add(group_id, "krook_rate_scale", krook_rate_scale)
        call h5_add(group_id, "krook_model", krook_model)
        call h5_add(group_id, "krook_harmonic", krook_harmonic)
        call h5_add(group_id, "krook_geometry_trapped", krook_geometry_trapped)
        call h5_add(group_id, "krook_geometry_passing", krook_geometry_passing)
        call h5_add(group_id, "krook_eta_steps", krook_eta_steps)
        call h5_add(group_id, "krook_max_depth", krook_max_depth)
        call h5_add(group_id, "krook_rtol", krook_rtol)
        call h5_add(group_id, "krook_atol", krook_atol)
    end subroutine write_krook_config_group

    subroutine write_magfie_group(group_id, data)
        use hdf5_tools, only: HID_T, h5_add

        integer(HID_T), intent(in) :: group_id
        type(magfie_data_t), intent(in) :: data

        integer :: n, k
        real(dp), allocatable :: theta(:), bmod(:), sqrtg(:)
        real(dp), allocatable :: hder(:, :), hcovar(:, :), hctrvr(:, :), hcurl(:, :)
        complex(dp), allocatable :: bn(:), eps_exp(:)

        call h5_add(group_id, "s", data%params%s, unit="1", comment="normalised toroidal flux")
        call h5_add(group_id, "R0", data%params%R0, unit="cm", comment="major radius")
        call h5_add(group_id, "a", data%params%a, unit="cm", comment="minor radius")
        call h5_add(group_id, "eps", data%params%eps, unit="1", comment="inverse aspect ratio")
        call h5_add(group_id, "psi_pr", data%params%psi_pr, &
                    comment="radial derivative of poloidal flux")
        call h5_add(group_id, "B0", data%params%B0, unit="G", comment="reference |B|")
        call h5_add(group_id, "Bthcov", data%params%Bthcov, comment="covariant poloidal B")
        call h5_add(group_id, "Bphcov", data%params%Bphcov, comment="covariant toroidal B")
        call h5_add(group_id, "dBthcovds", data%params%dBthcovds)
        call h5_add(group_id, "dBphcovds", data%params%dBphcovds)
        call h5_add(group_id, "q", data%params%q, unit="1", comment="safety factor")
        call h5_add(group_id, "iota", data%params%iota, unit="1", comment="rotational transform")
        call h5_add(group_id, "dVds", data%params%dVds, unit="cm^3", &
                    comment="flux surface volume derivative")
        call h5_add(group_id, "M_t", data%params%M_t, unit="1", comment="toroidal Mach number")
        call h5_add(group_id, "Om_tE", data%params%Om_tE, unit="1/s", &
                    comment="toroidal rotation frequency from E field")
        call h5_add(group_id, "Om_tBref", data%params%Om_tBref, unit="1/s", &
                    comment="reference magnetic drift frequency")
        call h5_add(group_id, "vth", data%params%vth, unit="cm/s", comment="thermal velocity")
        call h5_add(group_id, "T", data%params%T_in_eV, unit="eV", comment="temperature")
        call h5_add(group_id, "m0", data%params%m0, comment="poloidal perturbation mode")
        call h5_add(group_id, "n0", data%params%n0, comment="toroidal perturbation mode")
        call h5_add(group_id, "Dp", data%params%Dp, comment="plateau diffusion coefficient")
        call h5_add(group_id, "Drp_over_Dp", data%params%Drp_over_Dp, unit="1")
        call h5_add(group_id, "etatp", data%params%etatp, &
                    comment="eta at trapped-passing boundary")
        call h5_add(group_id, "etadt", data%params%etadt, &
                    comment="eta at deeply trapped boundary")
        call h5_add(group_id, "pertfile", data%params%pertfile)
        call h5_add(group_id, "nonlin", data%params%nonlin)

        if (data%params%nonlin) then
            call h5_add(group_id, "dpp", data%params%dpp)
            call h5_add(group_id, "dhh", data%params%dhh)
            call h5_add(group_id, "fpeff", data%params%fpeff)
        end if

        n = data%n_points
        allocate (theta(n), bmod(n), sqrtg(n))
        allocate (hder(3, n), hcovar(3, n), hctrvr(3, n), hcurl(3, n))
        allocate (bn(n), eps_exp(n))

        do k = 1, n
            theta(k) = data%tensors(k)%theta
            bmod(k) = data%tensors(k)%bmod
            sqrtg(k) = data%tensors(k)%sqrtg
            hder(:, k) = data%tensors(k)%hder
            hcovar(:, k) = data%tensors(k)%hcovar
            hctrvr(:, k) = data%tensors(k)%hctrvr
            hcurl(:, k) = data%tensors(k)%hcurl
            bn(k) = data%tensors(k)%bn
            eps_exp(k) = data%tensors(k)%eps_exp
        end do

        call h5_add(group_id, "theta", theta, [1], [n], unit="rad", &
                    comment="Boozer poloidal angle")
        call h5_add(group_id, "bmod", bmod, [1], [n], unit="G", comment="|B|")
        call h5_add(group_id, "sqrtg", sqrtg, [1], [n], comment="Jacobian determinant")
        call h5_add(group_id, "hder", hder, [1, 1], [3, n], comment="derivatives of |B| direction")
        call h5_add(group_id, "hcovar", hcovar, [1, 1], [3, n], comment="covariant b components")
        call h5_add(group_id, "hctrvr", hctrvr, [1, 1], [3, n], &
                    comment="contravariant b components")
        call h5_add(group_id, "hcurl", hcurl, [1, 1], [3, n], comment="curl of b")
        call h5_add(group_id, "bn", bn, [1], [n], comment="perturbation amplitude / |B|")
        call h5_add(group_id, "eps_exp", eps_exp, [1], [n], &
                    comment="reference analytic epsmn*exp(i*m0*theta)")
    end subroutine write_magfie_group

    subroutine write_transport_group(group_id, data)
        use hdf5_tools, only: HID_T, h5_add

        integer(HID_T), intent(in) :: group_id
        type(transport_data_t), intent(in) :: data

        call h5_add(group_id, "M_t", data%summary%M_t, unit="1", comment="toroidal Mach number")
        call h5_add(group_id, "D11co", data%summary%Dco(1), comment="co-passing D11")
        call h5_add(group_id, "D11ctr", data%summary%Dctr(1), comment="counter-passing D11")
        call h5_add(group_id, "D11t", data%summary%Dt(1), comment="trapped D11")
        call h5_add(group_id, "D11", &
                    data%summary%Dco(1) + data%summary%Dctr(1) + data%summary%Dt(1), &
                    comment="particle flux coefficient, co+ctr+trapped")
        call h5_add(group_id, "D12co", data%summary%Dco(2), comment="co-passing D12")
        call h5_add(group_id, "D12ctr", data%summary%Dctr(2), comment="counter-passing D12")
        call h5_add(group_id, "D12t", data%summary%Dt(2), comment="trapped D12")
        call h5_add(group_id, "D12", &
                    data%summary%Dco(2) + data%summary%Dctr(2) + data%summary%Dt(2), &
                    comment="momentum flux coefficient, co+ctr+trapped")
    end subroutine write_transport_group

    subroutine write_torque_group(group_id, data)
        use hdf5_tools, only: HID_T, h5_add

        integer(HID_T), intent(in) :: group_id
        type(transport_data_t), intent(in) :: data

        call h5_add(group_id, "has_torque", data%torque%has_torque, &
                    comment="false means the remaining torque values are unset")
        call h5_add(group_id, "s", data%torque%s, unit="1")
        call h5_add(group_id, "dVds", data%torque%dVds, unit="cm^3")
        call h5_add(group_id, "M_t", data%torque%M_t, unit="1")
        call h5_add(group_id, "Tco", data%torque%Tco, unit="erg", &
                    comment="co-passing torque density dTphi_int/ds")
        call h5_add(group_id, "Tctr", data%torque%Tctr, unit="erg", &
                    comment="counter-passing torque density dTphi_int/ds")
        call h5_add(group_id, "Tt", data%torque%Tt, unit="erg", &
                    comment="trapped torque density dTphi_int/ds")
    end subroutine write_torque_group

    subroutine write_harmonics_group(group_id, data)
        ! Per-harmonic transport and torque. The text pipeline split these
        ! across _integral.out and _torque_integral.out and the collector
        ! script read neither.
        use hdf5_tools, only: HID_T, h5_add

        integer(HID_T), intent(in) :: group_id
        type(transport_data_t), intent(in) :: data

        integer :: n, k
        integer, allocatable :: mth(:)
        real(dp), allocatable :: D11co(:), D11ctr(:), D11t(:)
        real(dp), allocatable :: D12co(:), D12ctr(:), D12t(:)
        real(dp), allocatable :: Tco(:), Tctr(:), Tt(:)
        real(dp), allocatable :: vminp(:), vmaxp(:), vmint(:), vmaxt(:)

        n = size(data%harmonics)
        call h5_add(group_id, "n_harmonics", n)
        if (n == 0) return

        allocate (mth(n))
        allocate (D11co(n), D11ctr(n), D11t(n), D12co(n), D12ctr(n), D12t(n))
        allocate (Tco(n), Tctr(n), Tt(n))
        allocate (vminp(n), vmaxp(n), vmint(n), vmaxt(n))

        do k = 1, n
            mth(k) = data%harmonics(k)%mth
            D11co(k) = data%harmonics(k)%Dresco(1)
            D11ctr(k) = data%harmonics(k)%Dresctr(1)
            D11t(k) = data%harmonics(k)%Drest(1)
            D12co(k) = data%harmonics(k)%Dresco(2)
            D12ctr(k) = data%harmonics(k)%Dresctr(2)
            D12t(k) = data%harmonics(k)%Drest(2)
            Tco(k) = data%harmonics(k)%Tresco
            Tctr(k) = data%harmonics(k)%Tresctr
            Tt(k) = data%harmonics(k)%Trest
            vminp(k) = data%harmonics(k)%vminp_over_vth
            vmaxp(k) = data%harmonics(k)%vmaxp_over_vth
            vmint(k) = data%harmonics(k)%vmint_over_vth
            vmaxt(k) = data%harmonics(k)%vmaxt_over_vth
        end do

        call h5_add(group_id, "mth", mth, [1], [n], comment="poloidal resonance harmonic")
        call h5_add(group_id, "D11co", D11co, [1], [n])
        call h5_add(group_id, "D11ctr", D11ctr, [1], [n])
        call h5_add(group_id, "D11t", D11t, [1], [n])
        call h5_add(group_id, "D11", D11co + D11ctr + D11t, [1], [n], &
                    comment="co+ctr+trapped per harmonic")
        call h5_add(group_id, "D12co", D12co, [1], [n])
        call h5_add(group_id, "D12ctr", D12ctr, [1], [n])
        call h5_add(group_id, "D12t", D12t, [1], [n])
        call h5_add(group_id, "D12", D12co + D12ctr + D12t, [1], [n], &
                    comment="co+ctr+trapped per harmonic")
        call h5_add(group_id, "Tco", Tco, [1], [n], unit="erg")
        call h5_add(group_id, "Tctr", Tctr, [1], [n], unit="erg")
        call h5_add(group_id, "Tt", Tt, [1], [n], unit="erg")
        call h5_add(group_id, "vminp_over_vth", vminp, [1], [n], unit="1")
        call h5_add(group_id, "vmaxp_over_vth", vmaxp, [1], [n], unit="1")
        call h5_add(group_id, "vmint_over_vth", vmint, [1], [n], unit="1")
        call h5_add(group_id, "vmaxt_over_vth", vmaxt, [1], [n], unit="1")
    end subroutine write_harmonics_group

    subroutine write_mixture_output(result_, base_path, output_format)
        type(mixture_result_t), intent(in) :: result_
        character(len=*), intent(in) :: base_path, output_format
        integer :: k
        select case (trim(output_format))
        case ('hdf5', 'both')
            call write_mixture_hdf5(result_, base_path)
        case ('text')
        case default
            error stop "output_format must be 'hdf5', 'text' or 'both'"
        end select
        if (trim(output_format) == 'hdf5') return
        do k = 1, size(result_%species)
            call write_output(result_%species(k)%magfie, result_%species(k)%transport, &
                trim(base_path)//'_'//trim(result_%species(k)%species%name), 'text')
        end do
        call write_mixture_text(result_, base_path)
    end subroutine write_mixture_output

    subroutine write_mixture_hdf5(result_, base_path)
        use hdf5_tools, only: HID_T, h5_init, h5_deinit, h5_create, h5_close, &
            h5_add, h5_define_group, h5_close_group
        type(mixture_result_t), intent(in) :: result_
        character(len=*), intent(in) :: base_path
        integer(HID_T) :: file_id, group_id, species_id
        integer :: k
        call h5_init()
        call h5_create(trim(adjustl(base_path))//'.h5', file_id)
        call h5_add(file_id, 'format_version', 2, comment='NEO-RT multispecies output layout')
        call h5_add(file_id, 'created', timestamp(), comment='ISO 8601 local time')
        call h5_add(file_id, 'n_species', size(result_%species))
        call write_mixture_metadata(file_id, result_)
        call h5_define_group(file_id, 'species', group_id)
        do k = 1, size(result_%species)
            call h5_define_group(group_id, trim(result_%species(k)%species%name), species_id)
            call write_species_result_group(species_id, result_%species(k))
            call h5_close_group(species_id)
        end do
        call h5_close_group(group_id)
        call h5_close(file_id)
        call h5_deinit()
    end subroutine write_mixture_hdf5

    subroutine write_mixture_metadata(file_id, result_)
        use hdf5_tools, only: HID_T, h5_add, h5_define_group, h5_close_group
        use neort_species, only: species_mixture
        integer(HID_T), intent(in) :: file_id
        type(mixture_result_t), intent(in) :: result_
        integer(HID_T) :: group_id
        call h5_define_group(file_id, 'mixture', group_id)
        call h5_add(group_id, 'reference_s', species_mixture%reference_s, unit='1')
        call h5_add(group_id, 'reference_omega_e', species_mixture%omega_e, unit='s^-1')
        call h5_add(group_id, 'reference_domega_e_ds', species_mixture%domega_e_ds, unit='s^-1')
        call h5_add(group_id, 'torque_observable', 'electromagnetic_perturbation_work', &
            comment='nonadiabatic-g Krook; evaluated by legacy flux-force-equivalent formula')
        call h5_add(group_id, 'coulomb_log', species_mixture%coulomb_log, unit='1', &
            comment='explicit common Coulomb logarithm for pairwise rate models')
        call h5_add(group_id, 'profile_model', 'affine local expansion in normalized toroidal flux')
        call h5_add(group_id, 'torque_total', sum(result_%torque_total), unit='erg', &
            comment='sum of species co+counter+trapped dTphi_int/ds')
        call h5_add(group_id, 'torque_by_population', result_%torque_total, [1], [3], unit='erg', &
            comment='species sum: co-passing, counter-passing, trapped dTphi_int/ds')
        call h5_close_group(group_id)
    end subroutine write_mixture_metadata

    subroutine write_species_result_group(species_id, result_)
        use hdf5_tools, only: HID_T, h5_define_group, h5_close_group
        integer(HID_T), intent(in) :: species_id
        type(species_result_t), intent(in) :: result_
        integer(HID_T) :: group_id
        call write_species_metadata(species_id, result_)
        call h5_define_group(species_id, 'config', group_id)
        call write_config_group(group_id, result_%species)
        call h5_close_group(group_id)
        call h5_define_group(species_id, 'magfie', group_id)
        call write_magfie_group(group_id, result_%magfie)
        call h5_close_group(group_id)
        call h5_define_group(species_id, 'transport', group_id)
        call write_transport_group(group_id, result_%transport)
        call h5_close_group(group_id)
        call h5_define_group(species_id, 'torque', group_id)
        call write_torque_group(group_id, result_%transport)
        call h5_close_group(group_id)
        call h5_define_group(species_id, 'harmonics', group_id)
        call write_harmonics_group(group_id, result_%transport)
        call h5_close_group(group_id)
    end subroutine write_species_result_group

    subroutine write_species_metadata(species_id, result_)
        use hdf5_tools, only: HID_T, h5_add
        integer(HID_T), intent(in) :: species_id
        type(species_result_t), intent(in) :: result_
        call h5_add(species_id, 'thermal_deflection_rate', &
            result_%collision_thermal_deflection_rate, unit='s^-1', &
            comment='sum_b nu_D_ab(v_th_a), before scalar rate scaling/geometry')
        call h5_add(species_id, 'collision_reference_frequency', &
            result_%collision_reference_frequency, unit='s^-1', &
            comment='Gamma_aa/v_th_a^3, distinct from summed deflection rate')
        call h5_add(species_id, 'mass_kg', result_%species%mass_kg, unit='kg')
        call h5_add(species_id, 'charge_e', result_%species%charge_e, unit='1')
        call h5_add(species_id, 'density_m3', result_%species%density_m3, unit='m^-3')
        call h5_add(species_id, 'temperature_ev', result_%species%temperature_ev, unit='eV')
        call h5_add(species_id, 'density_gradient_m3', &
            result_%species%density_gradient_m3, unit='m^-3', comment='dn/ds')
        call h5_add(species_id, 'temperature_gradient_ev', &
            result_%species%temperature_gradient_ev, unit='eV', comment='dT/ds')
        call h5_add(species_id, 'domega_e_ds', result_%domega_e_ds, unit='s^-1')
        call h5_add(species_id, 'thermodynamic_forces', result_%force, [1], [2], unit='1')
    end subroutine write_species_metadata

    subroutine write_mixture_text(result_, base_path)
        type(mixture_result_t), intent(in) :: result_
        character(len=*), intent(in) :: base_path
        integer :: unit, k
        open(newunit=unit, file=trim(base_path)//'_species.out', status='replace')
        write(unit, '(a)') '# species mass[kg] Z n[m^-3] T[eV] dn/ds[m^-3] dT/ds[eV]'// &
            ' OmegaE[s^-1] dOmegaE/ds[s^-1] Mach Tco[erg] Tctr[erg] Tt[erg] Ttotal[erg]'
        do k = 1, size(result_%species)
            call write_species_text_row(unit, result_%species(k))
        end do
        close(unit)
        open(newunit=unit, file=trim(base_path)//'_mixture_torque.out', status='replace')
        write(unit, '(a)') '# s Tco[erg] Tctr[erg] Tt[erg] Ttotal[erg]; sum of species dTphi/ds'
        write(unit, '(5(es24.16e3,1x))') result_%species(1)%magfie%params%s, &
            result_%torque_total, sum(result_%torque_total)
        close(unit)
    end subroutine write_mixture_text

    subroutine write_species_text_row(unit, result_)
        integer, intent(in) :: unit
        type(species_result_t), intent(in) :: result_
        write(unit, '(a,1x,13(es24.16e3,1x))') trim(result_%species%name), &
            result_%species%mass_kg, result_%species%charge_e, result_%species%density_m3, &
            result_%species%temperature_ev, result_%species%density_gradient_m3, &
            result_%species%temperature_gradient_ev, result_%magfie%params%Om_tE, &
            result_%domega_e_ds, result_%magfie%params%M_t, result_%transport%torque%Tco, &
            result_%transport%torque%Tctr, result_%transport%torque%Tt, &
            result_%transport%torque%Tco + result_%transport%torque%Tctr + &
            result_%transport%torque%Tt
    end subroutine write_species_text_row

    function timestamp() result(stamp)
        character(len=19) :: stamp
        character(len=8) :: date
        character(len=10) :: time

        call date_and_time(date=date, time=time)
        stamp = date(1:4)//"-"//date(5:6)//"-"//date(7:8)//"T"// &
                time(1:2)//":"//time(3:4)//":"//time(5:6)
    end function timestamp

end module neort_output
