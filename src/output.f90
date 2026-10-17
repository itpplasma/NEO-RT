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
        use neort_lorentz_controls, only: lorentz_pitch_collisions
        use neort_fp_mode, only: tp_fp_collisions

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
        if (lorentz_pitch_collisions .and. trim(output_format) /= 'hdf5') &
            call write_lorentz_text_metadata(transport_data, base_path)
        if (tp_fp_collisions .and. trim(output_format) /= 'hdf5') &
            call write_fp_text_metadata(transport_data, base_path)
        if (transport_data%conserving_active .and. trim(output_format) /= 'hdf5') &
            call write_conserving_text_metadata(transport_data, base_path)
    end subroutine write_output

    subroutine write_conserving_text_metadata(data, base_path)
        type(transport_data_t), intent(in) :: data
        character(*), intent(in) :: base_path
        integer :: unit, driver
        open(newunit=unit, file=trim(base_path)//'_collision_model.out', &
            status='replace')
        write(unit, '(a)') '# source_feature=all-pair-linearized-Landau-commonT-v1'
        write(unit, '(a)') '# susceptibility_valid=F; all D columns are undefined placeholders'
        write(unit, '(a)') '# joint common-T enriched response; direct electromagnetic work'
        write(unit, '(a)') '# work branch order: trapped, co-passing, counter-passing'
        write(unit, '(a)') '# velocity domain: 0<=u<infinity; '// &
            'graph constant beyond outer knot'
        write(unit, '(a,i0)') '# conserving_rank=', data%conserving_rank
        write(unit, '(a,i0)') '# conserving_nullity=', data%conserving_nullity
        write(unit, '(a,es24.16)') '# conserving_singular_cutoff=', &
            data%conserving_singular_cutoff
        write(unit, '(a,es24.16)') '# conserving_residual=', data%conserving_residual
        write(unit, '(a,es24.16)') '# conserving_compatibility=', &
            data%conserving_compatibility
        write(unit, '(a,es24.16)') '# conserving_mass_condition=', &
            data%conserving_mass_condition
        write(unit, '(a,es24.16)') '# conserving_entropy=', data%conserving_entropy
        write(unit, '(a,6es24.16)') '# conserving_work=', data%conserving_work
        if (allocated(data%conserving_cross_work)) then
            do driver = 1, size(data%conserving_cross_work, 2)
                write(unit, '(a,i0,6es24.16)') '# driver_work ', driver, &
                    data%conserving_cross_work(:, driver)
            end do
        end if
        if (allocated(data%conserving_cross_work_qualified)) then
            do driver = 1, size(data%conserving_cross_work_qualified)
                write(unit, '(a,i0,1x,l1)') '# driver_work_qualified ', driver, &
                    data%conserving_cross_work_qualified(driver)
            end do
        end if
        close(unit)
    end subroutine write_conserving_text_metadata

    subroutine write_fp_text_metadata(data, base_path)
        type(transport_data_t), intent(in) :: data
        character(*), intent(in) :: base_path
        integer :: unit
        open(newunit=unit, file=trim(base_path)//'_collision_model.out', &
            status='replace')
        write(unit, '(a)') '# model=equal_temperature_test_particle_fp'
        write(unit, '(a)') '# susceptibility_valid=F; '// &
            'all D columns are undefined placeholders'
        write(unit, '(a,l1)') '# fp_auxiliary_valid=', data%fp_auxiliary_valid
        write(unit, '(a,l1)') '# fp_physical_response=', data%fp_physical_response
        write(unit, '(a)') '# finite speed endpoint reflects; '// &
            'common physical origin u=0'
        write(unit, '(a)') '# torque is direct electromagnetic perturbation work'
        close(unit)
    end subroutine write_fp_text_metadata

    subroutine write_lorentz_text_metadata(data, base_path)
        type(transport_data_t), intent(in) :: data
        character(*), intent(in) :: base_path
        integer :: unit
        open(newunit=unit, file=trim(base_path)//'_collision_model.out', &
            status='replace')
        write(unit, '(a)') '# model=common_phase_bounce_averaged_lorentz'
        write(unit, '(a,l1)') '# susceptibility_valid=', data%susceptibility_valid
        write(unit, '(a)') '# D columns are susceptibility weights, '// &
            'not direct flux moments'
        if (.not. data%susceptibility_valid) &
            write(unit, '(a)') '# D columns are undefined placeholders; '// &
                'physical zero-drive torque=0'
        close(unit)
    end subroutine write_lorentz_text_metadata

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
                              fow_radial
        use neort, only: vsteps, mth_max_abs, vmax_over_vth
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
        call write_lorentz_config_group(group_id)
        call write_fp_config_group(group_id)
        call write_conserving_config_group(group_id)
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

    subroutine write_lorentz_config_group(group_id)
        use hdf5_tools, only: HID_T, h5_add
        use neort_lorentz_controls, only: lorentz_pitch_collisions, lorentz_nu, &
            lorentz_rate_source, lorentz_rate_scale, lorentz_pitch_steps
        integer(HID_T), intent(in) :: group_id
        call h5_add(group_id, 'lorentz_pitch_collisions', lorentz_pitch_collisions)
        call h5_add(group_id, 'lorentz_nu', lorentz_nu, unit='s^-1', &
            comment='supplied physical deflection frequency; '// &
                'C_L=(nu_D/2)d_xi[(1-xi2)d_xi]')
        call h5_add(group_id, 'lorentz_rate_source', lorentz_rate_source, &
            comment='0 supplied constant nu_D; '// &
                '1 all-background nu_D at physical test speed')
        call h5_add(group_id, 'lorentz_rate_scale', lorentz_rate_scale)
        call h5_add(group_id, 'lorentz_pitch_steps', lorentz_pitch_steps, &
            comment='P1 cells per branch, physical endpoints and one shared TP value')
        if (.not. lorentz_pitch_collisions) return
        call h5_add(group_id, 'lorentz_response_model', 'common_phase_bounce_averaged')
        call h5_add(group_id, 'lorentz_relaxed_variable', 'nonadiabatic_g', &
            comment='coordinate-consistent delta C F0=-C0[delta H*dF0/dH]')
        call h5_add(group_id, 'lorentz_harmonic_policy', 'trapped_0_passing_minus_nq', &
            comment='fixed common physical phase; '// &
                'mth_max_abs controls only legacy mode')
        call h5_add(group_id, 'lorentz_ordering', 'low_collision_thin_orbit', &
            comment='not a full arbitrary-collisionality DKE or Fokker-Planck operator')
    end subroutine write_lorentz_config_group

    subroutine write_fp_config_group(group_id)
        use hdf5_tools, only: HID_T, h5_add
        use neort_fp_controls, only: tp_fp_collisions, fp_rate_scale, &
            fp_speed_steps, fp_pitch_steps, fp_umax
        integer(HID_T), intent(in) :: group_id
        if (.not. tp_fp_collisions) return
        call h5_add(group_id, 'tp_fp_collisions', tp_fp_collisions)
        call h5_add(group_id, 'fp_rate_scale', fp_rate_scale)
        call h5_add(group_id, 'fp_speed_steps', fp_speed_steps)
        call h5_add(group_id, 'fp_pitch_steps', fp_pitch_steps)
        call h5_add(group_id, 'fp_umax', fp_umax)
        call h5_add(group_id, 'fp_operator', 'equal_T_Maxwellian_test_particle')
        call h5_add(group_id, 'fp_relaxed_variable', 'nonadiabatic_g')
        call h5_add(group_id, 'fp_rate_source', 'all_background_physical_Coulomb')
        call h5_add(group_id, 'fp_speed_boundary', 'common_origin_reflecting_umax')
        call h5_add(group_id, 'fp_ordering', 'low_collision_common_phase_thin_orbit')
        call h5_add(group_id, 'fp_conservation', 'test_particle_number_only')
    end subroutine write_fp_config_group

    subroutine write_conserving_config_group(group_id)
        use hdf5_tools, only: HID_T, h5_add
        use neort_fp_controls, only: conserving_fp_collisions, &
            conserving_speed_steps, conserving_pitch_steps, conserving_umax, &
            conserving_rate_scale, &
            conserving_theta_order, conserving_velocity_order, conserving_cm_order, &
            conserving_rel_order, conserving_polar_order, conserving_azimuth_order, &
            conserving_rank_rtol, conserving_compatibility_rtol
        integer(HID_T), intent(in) :: group_id
        if (.not. conserving_fp_collisions) return
        call h5_add(group_id, 'conserving_fp_collisions', conserving_fp_collisions)
        call h5_add(group_id, 'conserving_speed_steps', conserving_speed_steps)
        call h5_add(group_id, 'conserving_pitch_steps', conserving_pitch_steps)
        call h5_add(group_id, 'conserving_umax', conserving_umax)
        call h5_add(group_id, 'conserving_rate_scale', conserving_rate_scale, &
            comment='one positive global multiplier on the entire all-pair operator')
        call h5_add(group_id, 'conserving_theta_order', conserving_theta_order)
        call h5_add(group_id, 'conserving_velocity_order', conserving_velocity_order)
        call h5_add(group_id, 'conserving_cm_order', conserving_cm_order)
        call h5_add(group_id, 'conserving_rel_order', conserving_rel_order)
        call h5_add(group_id, 'conserving_polar_order', conserving_polar_order)
        call h5_add(group_id, 'conserving_azimuth_order', conserving_azimuth_order)
        call h5_add(group_id, 'conserving_rank_rtol', conserving_rank_rtol)
        call h5_add(group_id, 'conserving_compatibility_rtol', &
            conserving_compatibility_rtol)
        call h5_add(group_id, 'conserving_basis', 'joined_pitch_speed_plus_flow_energy')
        call h5_add(group_id, 'conserving_speed_tail', &
            'constant_graph_analytic_moments')
        call h5_add(group_id, 'conserving_ordering', 'common_phase_thin_orbit_common_T')
        call h5_add(group_id, 'conserving_conservation', &
            'species_number_total_parallel_momentum_total_energy')
    end subroutine write_conserving_config_group

    subroutine write_conserving_transport(group_id, data)
        use hdf5_tools, only: HID_T, h5_add
        integer(HID_T), intent(in) :: group_id
        type(transport_data_t), intent(in) :: data
        integer :: drivers
        call h5_add(group_id, 'conserving_active', .true.)
        call h5_add(group_id, 'source_feature', &
            'all-pair-linearized-Landau-commonT-v1')
        call h5_add(group_id, 'conserving_velocity_domain', &
            '0<=u<infinity; graph constant beyond outer knot; analytic flow/energy tails')
        call h5_add(group_id, 'conserving_projection', &
            'common spatial phase; joined graph plus exact parallel flow and energy')
        call h5_add(group_id, 'susceptibility_valid', .false., &
            comment='all legacy D weights are undefined placeholders')
        call h5_add(group_id, 'conserving_work_real', &
            real(data%conserving_work, dp), [1], [3], unit='erg', &
            comment='integral J F0 conj(Herg) g dtheta d3v; trapped,co,counter')
        call h5_add(group_id, 'conserving_work_imag', &
            aimag(data%conserving_work), [1], [3], unit='erg', &
            comment='physical torque by branch is -n*pi times this integral')
        if (allocated(data%conserving_cross_work)) then
            drivers = size(data%conserving_cross_work, 2)
            call h5_add(group_id, 'conserving_cross_work_real', &
                real(data%conserving_cross_work, dp), [1, 1], [3, drivers], &
                unit='erg', &
                comment='receiving branch by physical source species')
            call h5_add(group_id, 'conserving_cross_work_imag', &
                aimag(data%conserving_cross_work), [1, 1], [3, drivers], unit='erg', &
                comment='unqualified driver columns are undefined placeholders')
        end if
        if (allocated(data%conserving_cross_work_qualified)) then
            drivers = size(data%conserving_cross_work_qualified)
            call h5_add(group_id, 'conserving_cross_work_qualified', &
                merge(1, 0, data%conserving_cross_work_qualified), [1], [drivers], &
                comment='1: independently compatible physical source-species response')
        end if
        call h5_add(group_id, 'conserving_rank', data%conserving_rank)
        call h5_add(group_id, 'conserving_nullity', data%conserving_nullity)
        call h5_add(group_id, 'conserving_singular_cutoff', &
            data%conserving_singular_cutoff, unit='1/s')
        call h5_add(group_id, 'conserving_residual', data%conserving_residual)
        call h5_add(group_id, 'conserving_compatibility', data%conserving_compatibility)
        call h5_add(group_id, 'conserving_mass_condition', &
            data%conserving_mass_condition)
        call h5_add(group_id, 'conserving_entropy', data%conserving_entropy, &
            unit='1/s', comment='joint all-species number-weighted collision '// &
                'quadratic form; repeated for each receiving species')
    end subroutine write_conserving_transport

    subroutine write_magfie_group(group_id, data)
        use neort_fp_mode, only: fp_physical_mode
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
        if (fp_physical_mode()) then
            call h5_add(group_id, 'fp_thermal_energy_J', &
                data%params%fp_thermal_energy_j, unit='J')
            call h5_add(group_id, 'fp_vth_m_s', data%params%vth/100.0_dp, unit='m/s')
            call h5_add(group_id, 'fp_charge_over_c', data%params%fp_charge_over_c, &
                comment='native EM coefficient qi/c = physical SI charge/10')
        end if
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
        use neort_lorentz_controls, only: lorentz_pitch_collisions
        use neort_fp_mode, only: tp_fp_collisions, fp_physical_mode

        integer(HID_T), intent(in) :: group_id
        type(transport_data_t), intent(in) :: data

        if (lorentz_pitch_collisions) call h5_add(group_id, 'susceptibility_valid', &
            data%susceptibility_valid, &
            comment='false: D weights are undefined placeholders; '// &
                'physical zero-drive torque=0')
        if (tp_fp_collisions) then
            call h5_add(group_id, 'susceptibility_valid', .false., &
                comment='all legacy D weights are undefined placeholders')
            call h5_add(group_id, 'fp_auxiliary_valid', data%fp_auxiliary_valid, &
                comment='true: work columns solve original H and u^2 H sources')
            call h5_add(group_id, 'fp_physical_response', data%fp_physical_response, &
                comment='true: column 1 contains physical g; column 2 is undefined')
            call h5_add(group_id, 'fp_source_work_real', &
                real(data%fp_source_work, dp), [1, 1], [2, 3], &
                comment='Re integral W conj(H) response; branches trapped,co,counter')
            call h5_add(group_id, 'fp_source_work_imag', &
                aimag(data%fp_source_work), [1, 1], [2, 3], &
                comment='Im integral W conj(H) response; original H receiving dual')
        end if
        if (data%conserving_active) call write_conserving_transport(group_id, data)
        call h5_add(group_id, "M_t", data%summary%M_t, unit="1", comment="toroidal Mach number")
        call h5_add(group_id, "D11co", data%summary%Dco(1), comment="co-passing D11")
        call h5_add(group_id, "D11ctr", data%summary%Dctr(1), comment="counter-passing D11")
        call h5_add(group_id, "D11t", data%summary%Dt(1), comment="trapped D11")
        if (fp_physical_mode()) then
            call h5_add(group_id, 'D11', 0.0_dp, &
                comment='undefined placeholder; use direct electromagnetic torque')
        else if (lorentz_pitch_collisions) then
            call h5_add(group_id, "D11", &
                data%summary%Dco(1) + data%summary%Dctr(1) + data%summary%Dt(1), &
                comment="susceptibility-derived legacy D11 weight; "// &
                    "not a direct flux moment")
        else
            call h5_add(group_id, "D11", &
                data%summary%Dco(1) + data%summary%Dctr(1) + data%summary%Dt(1), &
                comment="particle flux coefficient, co+ctr+trapped")
        end if
        call h5_add(group_id, "D12co", data%summary%Dco(2), comment="co-passing D12")
        call h5_add(group_id, "D12ctr", data%summary%Dctr(2), comment="counter-passing D12")
        call h5_add(group_id, "D12t", data%summary%Dt(2), comment="trapped D12")
        if (fp_physical_mode()) then
            call h5_add(group_id, 'D12', 0.0_dp, &
                comment='undefined placeholder; use direct electromagnetic torque')
        else if (lorentz_pitch_collisions) then
            call h5_add(group_id, "D12", &
                data%summary%Dco(2) + data%summary%Dctr(2) + data%summary%Dt(2), &
                comment="susceptibility-derived legacy D12 weight; "// &
                    "not a direct flux moment")
        else
            call h5_add(group_id, "D12", &
                data%summary%Dco(2) + data%summary%Dctr(2) + data%summary%Dt(2), &
                comment="momentum flux coefficient, co+ctr+trapped")
        end if
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
        use neort_lorentz_controls, only: lorentz_pitch_collisions
        use neort_fp_mode, only: tp_fp_collisions, conserving_fp_collisions
        integer(HID_T), intent(in) :: file_id
        type(mixture_result_t), intent(in) :: result_
        integer(HID_T) :: group_id
        call h5_define_group(file_id, 'mixture', group_id)
        call h5_add(group_id, 'reference_s', species_mixture%reference_s, unit='1')
        call h5_add(group_id, 'reference_omega_e', species_mixture%omega_e, unit='s^-1')
        call h5_add(group_id, 'reference_domega_e_ds', species_mixture%domega_e_ds, unit='s^-1')
        if (conserving_fp_collisions) then
            call h5_add(group_id, 'torque_observable', &
                'electromagnetic_perturbation_work', &
                comment='joint all-pair Landau projected common-phase response')
        else if (tp_fp_collisions) then
            call h5_add(group_id, 'torque_observable', &
                'electromagnetic_perturbation_work', &
                comment='nonadiabatic-g two-dimensional test-particle FP response')
        else if (lorentz_pitch_collisions) then
            call h5_add(group_id, 'torque_observable', &
                'electromagnetic_perturbation_work', &
                comment='nonadiabatic-g Lorentz common-phase complex response work')
        else
            call h5_add(group_id, 'torque_observable', &
                'electromagnetic_perturbation_work', &
                comment='nonadiabatic-g Krook; '// &
                    'evaluated by legacy flux-force-equivalent formula')
        end if
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
