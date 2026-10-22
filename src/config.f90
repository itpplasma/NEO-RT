module neort_config
    use iso_fortran_env, only: dp => real64
    use neort_phi1, only: phi1, phi1file, phi1_scale, phi_mn, m0_phi, reset_phi1
    use neort_krook, only: krook_collisions, krook_nu, krook_model, krook_harmonic, &
        krook_geometry_trapped, krook_geometry_passing, krook_eta_steps, &
        krook_max_depth, krook_rtol, krook_atol, validate_krook, reset_krook

    use neort_mixture_collisions, only: krook_rate_source, krook_rate_scale, &
        validate_collision_background, finite_krook_active
    use neort_resonance_cutout, only: resonance_cutout, resonance_cutout_width, &
        validate_resonance_cutout, reset_resonance_cutout, CUTOUT_OFF

    implicit none

    type :: config_t
        real(dp) :: s = 0.0_dp  ! radial coordinate (flux surface)
        real(dp) :: M_t = 0.0_dp  ! Mach number (for single Mach no. run) !*!
        real(dp) :: qs = 0.0_dp  ! particle charge / elementary charge !*!
        real(dp) :: ms = 0.0_dp  ! particle mass / u !*!
        real(dp) :: vth = 0.0_dp  ! thermal velocity / cm/s !*!
        real(dp) :: epsmn = 1.0_dp  ! perturbation amplitude B1/B0 (if pertfile==F)
        real(dp) :: pertfile_scale = 1.0_dp  ! scale factor for perturbations read from file
        integer :: m0 = 0  ! poloidal perturbation mode (if pertfile==F)
        integer :: mph = 0  ! toroidal perturbation mode (if pertfile==F, n>0!)
        logical :: phi1 = .false.  ! prescribed electrostatic perturbation
        character(1024) :: phi1file = ''  ! empty selects analytic phi_mn
        real(dp) :: phi1_scale = 1.0_dp
        complex(dp) :: phi_mn = (0.0_dp, 0.0_dp)  ! statvolt, exp(+i*m0_phi*theta)
        integer :: m0_phi = 0
        logical :: comptorque = .false.  ! compute torque
        logical :: collisional_layer = .false.  ! optional boundary-layer closure
        integer :: resonance_cutout = 0  ! 0 off, 1 zero, 2 boundary fill in the fold layer
        real(dp) :: resonance_cutout_width = 1.0_dp  ! fold-layer width factor >= 1
        logical :: supban = .false.  ! Shaing superbanana-plateau (trapped ell=0) only
        logical :: drift_phase = .false.  ! keep oscillating drift phase in orbit integral
        logical :: fow_radial = .false.  ! first-order radial orbit width in amplitude
        logical :: magdrift = .false.  ! consider magnetic drift
        ! Negative means "follow magdrift", preserving every existing deck.
        integer :: magdrift_passing = -1
        logical :: nopassing = .false.  ! neglect passing particles
        logical :: noshear = .false.  ! neglect magnetic shear term with dqds
        logical :: pertfile = .false.  ! read perturbation from file
        logical :: pertfile_curvature = .false.  ! read kappa.xi_perp table (general drive)
        logical :: nonlin = .false.  ! do nonlinear calculation
        integer :: profile_species = 1  ! 0 electrons, 1 first ion, 2 second ion
        character(len=1024) :: species_file = ''  ! SI mixture at a reference surface
        logical :: krook_collisions = .false.
        real(dp) :: krook_nu = 0.0_dp  ! supplied Krook damping rate [s^-1]
        integer :: krook_rate_source = 0
        real(dp) :: krook_rate_scale = 1.0_dp
        integer :: krook_model = 0  ! 0 constant, 1 energy+harmonic, 2 energy
        integer :: krook_harmonic = 0  ! 0 integer mth, 1 canonical passing mth+nq
        real(dp) :: krook_geometry_trapped = 1.0_dp
        real(dp) :: krook_geometry_passing = 1.0_dp
        integer :: krook_eta_steps = 32
        integer :: krook_max_depth = 40
        real(dp) :: krook_rtol = 1.0e-5_dp
        real(dp) :: krook_atol = 0.0_dp
        real(dp) :: bfac = 1.0_dp  ! scale B field by factor
        real(dp) :: efac = 1.0_dp  ! scale E field by factor
        integer :: inp_swi = 0  ! input switch for Boozer file
        integer :: inp_swi_pert = -1  ! negative inherits inp_swi
        integer :: vsteps = 0  ! integration steps in velocity space
        integer :: mth_max_abs = -1 ! negative: historical q-dependent range
        real(dp) :: vmax_over_vth = 4.0_dp  ! upper velocity cutoff / vth
        real(dp) :: monoenergetic_x = -1.0_dp  ! -1: thermal; positive: E/T
        integer :: log_level = 0  ! how much to log
        character(len=8) :: output_format = "both"  ! 'hdf5', 'text' or 'both'
        !*! will be overwritten if using splines from plasma.in and profile.in files
    end type config_t

    ! Resolved output format, read by neort_output at write time.
    character(len=8) :: output_format = "hdf5"
    character(len=1024) :: species_file = ''

contains

    subroutine set_config(config)
        ! Set global control parameters via config struct
        use do_magfie_mod, only: s, bfac, inp_swi
        use do_magfie_pert_mod, only: mph, set_mph, perturbation_switch => inp_swi_pert
        use driftorbit, only: epsmn, pertfile_scale, m0, comptorque, magdrift, &
            magdrift_passing, nopassing, pertfile, pertfile_curvature, &
            nonlin, efac, supban, drift_phase, fow_radial
        use neort_collisional_layer, only: collisional_layer
        use logger, only: set_log_level
        use neort, only: vsteps, mth_max_abs, vmax_over_vth, monoenergetic_x
        use kinetic_energy_selection, only: validate_energy_ratio
        use neort_orbit, only: noshear
        use neort_profiles, only: M_t, vth, profile_species
        use util, only: qe, mu, qi, mi
        use neort_species, only: reset_mixture, read_species_input

        type(config_t), intent(in) :: config

        call reset_phi1()
        phi1 = config%phi1
        phi1file = config%phi1file
        phi1_scale = config%phi1_scale
        phi_mn = config%phi_mn
        m0_phi = config%m0_phi
        call reset_mixture()
        species_file = config%species_file
        if (len_trim(species_file) > 0) call read_species_input(trim(species_file))
        call validate_mixture_mode(config%nonlin)
        s = config%s
        M_t = config%M_t * config%efac / config%bfac
        vth = config%vth
        epsmn = config%epsmn
        pertfile_scale = config%pertfile_scale
        m0 = config%m0
        mph = config%mph
        comptorque = config%comptorque
        collisional_layer = config%collisional_layer
        resonance_cutout = config%resonance_cutout
        resonance_cutout_width = config%resonance_cutout_width
        supban = config%supban
        drift_phase = config%drift_phase
        fow_radial = config%fow_radial
        magdrift = config%magdrift
        magdrift_passing = config%magdrift_passing
        if (magdrift_passing < 0) magdrift_passing = merge(1, 0, magdrift)
        nopassing = config%nopassing
        noshear = config%noshear
        pertfile = config%pertfile
        pertfile_curvature = config%pertfile_curvature
        call validate_curvature_drive(pertfile, pertfile_curvature)
        nonlin = config%nonlin
        profile_species = config%profile_species
        call set_krook_config(config)
        call validate_resonance_cutout()
        call validate_cutout_combination()
        bfac = config%bfac
        efac = config%efac
        inp_swi = config%inp_swi
        perturbation_switch = resolve_perturbation_switch(config%inp_swi, config%inp_swi_pert)
        call validate_perturbation_switch(pertfile, perturbation_switch)
        vsteps = config%vsteps
        if (config%mth_max_abs < -1) error stop "mth_max_abs must be -1 or nonnegative"
        mth_max_abs = config%mth_max_abs
        if (config%vmax_over_vth <= 0.0_dp) error stop "vmax_over_vth must be positive"
        vmax_over_vth = config%vmax_over_vth
        call validate_energy_ratio(config%monoenergetic_x)
        monoenergetic_x = config%monoenergetic_x
        output_format = config%output_format
        call validate_output_format(output_format)

        qi = config%qs * qe
        mi = config%ms * mu
        call set_mph(config%mph)
        call set_log_level(config%log_level)
    end subroutine set_config

    subroutine read_and_set_config(config_file)
        ! Set global control parameters directly from a file
        use do_magfie_mod, only: s, bfac, inp_swi
        use do_magfie_pert_mod, only: mph, set_mph, inp_swi_pert
        use driftorbit, only: epsmn, pertfile_scale, m0, comptorque, magdrift, &
            magdrift_passing, nopassing, pertfile, pertfile_curvature, &
            nonlin, efac, supban, drift_phase, fow_radial
        use neort_collisional_layer, only: collisional_layer
        use logger, only: set_log_level
        use neort, only: vsteps, mth_max_abs, vmax_over_vth, monoenergetic_x
        use kinetic_energy_selection, only: validate_energy_ratio
        use neort_orbit, only: noshear
        use neort_profiles, only: M_t, vth, profile_species
        use util, only: qe, mu, qi, mi
        use neort_species, only: reset_mixture, read_species_input

        character(len=*), intent(in) :: config_file
        real(dp) :: qs, ms
        integer :: log_level = 0

        namelist /params/ phi1, phi1file, phi1_scale, phi_mn, m0_phi, &
            s, M_t, qs, ms, vth, epsmn, pertfile_scale, m0, mph, comptorque, &
            supban, collisional_layer, drift_phase, fow_radial, &
            magdrift, magdrift_passing, nopassing, noshear, pertfile, pertfile_curvature, nonlin, bfac, efac, inp_swi, &
            inp_swi_pert, vsteps, mth_max_abs, vmax_over_vth, monoenergetic_x, log_level, output_format, &
            profile_species, species_file, krook_collisions, krook_nu, krook_model, krook_harmonic, &
            krook_geometry_trapped, krook_geometry_passing, krook_eta_steps, &
            krook_max_depth, krook_rtol, krook_atol, krook_rate_source, krook_rate_scale, &
            resonance_cutout, resonance_cutout_width

        call reset_phi1()
        collisional_layer = .false.
        call reset_resonance_cutout()
        mth_max_abs = -1
        vmax_over_vth = 4.0_dp
        monoenergetic_x = -1.0_dp
        inp_swi_pert = -1
        pertfile_curvature = .false.
        output_format = "both"
        profile_species = 1
        species_file = ''
        qs = 0.0_dp
        ms = 0.0_dp
        call reset_mixture()
        call reset_krook()
        krook_rate_source = 0
        krook_rate_scale = 1.0_dp
        open (unit=9, file=config_file, status="old", form="formatted")
        read (9, nml=params)
        close (unit=9)

        if (len_trim(species_file) > 0) call read_species_input(trim(species_file))
        call validate_mixture_mode(nonlin)
        call validate_krook()
        call validate_collision_background(.true.)
        call validate_resonance_cutout()
        call validate_cutout_combination()
        call validate_selected_physics(profile_species, nonlin)

        if (magdrift_passing < 0) magdrift_passing = merge(1, 0, magdrift)
        if (mth_max_abs < -1) error stop "mth_max_abs must be -1 or nonnegative"
        if (vmax_over_vth <= 0.0_dp) error stop "vmax_over_vth must be positive"
        call validate_energy_ratio(monoenergetic_x)
        call validate_output_format(output_format)
        inp_swi_pert = resolve_perturbation_switch(inp_swi, inp_swi_pert)
        call validate_perturbation_switch(pertfile, inp_swi_pert)
        call validate_curvature_drive(pertfile, pertfile_curvature)

        M_t = M_t * efac / bfac
        qi = qs * qe
        mi = ms * mu
        call set_mph(mph)
        call set_log_level(log_level)
    end subroutine read_and_set_config

    subroutine validate_cutout_combination()
        ! The cutout regularizes the collisionless delta-function weight only,
        ! and its layer overlaps the boundary-layer closure next to eta_tp.
        use neort_collisional_layer, only: collisional_layer
        if (resonance_cutout /= CUTOUT_OFF .and. finite_krook_active()) &
            error stop 'resonance_cutout and finite Krook broadening cannot be combined'
        if (resonance_cutout /= CUTOUT_OFF .and. collisional_layer) &
            error stop 'resonance_cutout and collisional_layer cannot be combined'
    end subroutine validate_cutout_combination

    subroutine validate_mixture_mode(nonlinear)
        use neort_species, only: mixture_enabled
        logical, intent(in) :: nonlinear
        if (mixture_enabled .and. nonlinear) &
            error stop 'general mixtures require linear transport; legacy nonlin is ion-only'
    end subroutine validate_mixture_mode

    subroutine set_krook_config(config)
        type(config_t), intent(in) :: config
        krook_collisions = config%krook_collisions
        krook_nu = config%krook_nu
        krook_rate_source = config%krook_rate_source
        krook_rate_scale = config%krook_rate_scale
        krook_model = config%krook_model
        krook_harmonic = config%krook_harmonic
        krook_geometry_trapped = config%krook_geometry_trapped
        krook_geometry_passing = config%krook_geometry_passing
        krook_eta_steps = config%krook_eta_steps
        krook_max_depth = config%krook_max_depth
        krook_rtol = config%krook_rtol
        krook_atol = config%krook_atol
        call validate_krook()
        call validate_collision_background()
        call validate_selected_physics(config%profile_species, config%nonlin)
    end subroutine set_krook_config

    subroutine validate_selected_physics(selected, nonlinear)
        integer, intent(in) :: selected
        logical, intent(in) :: nonlinear
        if (selected < 0 .or. selected > 2) error stop 'profile_species must be 0, 1 or 2'
        if (selected == 0 .and. nonlinear) &
            error stop 'nonlinear electron collision model is not implemented'
        if (finite_krook_active() .and. nonlinear) &
            error stop 'finite-width Krook transport requires nonlin=.false.'
    end subroutine validate_selected_physics

    pure integer function resolve_perturbation_switch(axisymmetric_switch, perturbation_switch)
        integer, intent(in) :: axisymmetric_switch, perturbation_switch

        if (perturbation_switch < 0) then
            resolve_perturbation_switch = axisymmetric_switch
        else
            resolve_perturbation_switch = perturbation_switch
        end if
    end function resolve_perturbation_switch

    subroutine validate_output_format(format)
        character(len=*), intent(in) :: format

        if (trim(format) /= "hdf5" .and. trim(format) /= "text" .and. &
            trim(format) /= "both") then
            error stop "output_format must be 'hdf5', 'text' or 'both'"
        end if
    end subroutine validate_output_format

    subroutine validate_curvature_drive(has_perturbation_file, has_curvature_file)
        logical, intent(in) :: has_perturbation_file, has_curvature_file

        if (has_curvature_file .and. .not. has_perturbation_file) &
            error stop "pertfile_curvature requires pertfile = .true."
    end subroutine validate_curvature_drive

    subroutine validate_perturbation_switch(has_perturbation_file, perturbation_switch)
        logical, intent(in) :: has_perturbation_file
        integer, intent(in) :: perturbation_switch

        if (has_perturbation_file .and. perturbation_switch /= 8 .and. &
            perturbation_switch /= 9) then
            error stop "inp_swi_pert must be 8 or 9 for a perturbation .bc file"
        end if
    end subroutine validate_perturbation_switch

end module neort_config
