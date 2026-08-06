module neort_config
    use iso_fortran_env, only: dp => real64

    implicit none

    type :: config_t
        real(dp) :: s = 0.0_dp  ! radial coordinate (flux surface)
        real(dp) :: M_t = 0.0_dp  ! Mach number (for single Mach no. run) !*!
        real(dp) :: qs = 0.0_dp  ! particle charge / elementary charge !*!
        real(dp) :: ms = 0.0_dp  ! particle mass / u !*!
        real(dp) :: vth = 0.0_dp  ! thermal velocity / cm/s !*!
        real(dp) :: epsmn = 0.0_dp  ! perturbation amplitude B1/B0 (if pertfile==F)
        integer :: m0 = 0  ! poloidal perturbation mode (if pertfile==F)
        integer :: mph = 0  ! toroidal perturbation mode (if pertfile==F, n>0!)
        logical :: comptorque = .false.  ! compute torque
        logical :: supban = .false.  ! Shaing superbanana-plateau (trapped ell=0) only
        logical :: magdrift = .false.  ! consider magnetic drift
        ! 0: historical local bounce/drift; 1: real-space GC thin-orbit limit.
        integer :: frequency_model = 0
        ! Negative means "follow magdrift", preserving every existing deck.
        integer :: magdrift_passing = -1
        logical :: nopassing = .false.  ! neglect passing particles
        logical :: noshear = .false.  ! neglect magnetic shear term with dqds
        logical :: pertfile = .false.  ! read perturbation from file
        logical :: nonlin = .false.  ! do nonlinear calculation
        real(dp) :: bfac = 1.0_dp  ! scale B field by factor
        real(dp) :: efac = 1.0_dp  ! scale E field by factor
        integer :: inp_swi = 0  ! input switch for Boozer file
        integer :: inp_swi_pert = -1  ! negative inherits inp_swi
        ! Axisymmetric Boozer file supplying theta_B(s, theta_geo) and the
        ! stream function.  Only meaningful for inp_swi=11 with a Boozer
        ! perturbation, where the orbit's angles are not the ones the
        ! perturbation is indexed in.
        character(len=1024) :: pert_angle_map = ''
        integer :: vsteps = 0  ! integration steps in velocity space
        integer :: mth_max_abs = -1 ! negative: historical q-dependent range
        real(dp) :: vmax_over_vth = 4.0_dp ! upper velocity cutoff / vth
        integer :: log_level = 0 ! how much to log
        !*! will be overwritten if using splines from plasma.in and profile.in files
    end type config_t

contains

    pure logical function perturbation_chart_is_compatible(axis_switch, &
            perturbation_switch, have_angle_map)
        !! A Boozer-Fourier perturbation may only be combined with the direct
        !! GEQDSK chart when the theta_B(s, theta_geo) map is supplied.  Without
        !! it the Fourier sum would be taken against the geometric angle, which
        !! returns a smooth and entirely wrong answer rather than failing.
        integer, intent(in) :: axis_switch, perturbation_switch
        logical, intent(in), optional :: have_angle_map

        logical :: mapped

        mapped = .false.
        if (present(have_angle_map)) mapped = have_angle_map

        if (axis_switch == 11 .and. perturbation_switch == 9) then
            perturbation_chart_is_compatible = mapped
            return
        end if
        perturbation_chart_is_compatible = &
            (axis_switch == 11) .eqv. (perturbation_switch == 11)
    end function perturbation_chart_is_compatible

    subroutine set_config(config)
        ! Set global control parameters via config struct
        use do_magfie_mod, only: s, bfac, inp_swi, has_direct_eqdsk_gc
        use do_magfie_pert_mod, only: mph, set_mph, set_pert_angle_map_path, &
            perturbation_switch => inp_swi_pert
        use driftorbit, only: epsmn, m0, comptorque, magdrift, &
            magdrift_passing, frequency_model, nopassing, pertfile, &
            nonlin, efac, supban
        use logger, only: set_log_level
        use neort, only: vsteps, mth_max_abs, vmax_over_vth
        use neort_orbit, only: noshear
        use neort_profiles, only: M_t, vth
        use util, only: qe, mu, qi, mi

        type(config_t), intent(in) :: config

        s = config%s
        M_t = config%M_t * config%efac / config%bfac
        vth = config%vth
        epsmn = config%epsmn
        m0 = config%m0
        mph = config%mph
        comptorque = config%comptorque
        supban = config%supban
        magdrift = config%magdrift
        frequency_model = config%frequency_model
        magdrift_passing = config%magdrift_passing
        if (magdrift_passing < 0) magdrift_passing = merge(1, 0, magdrift)
        nopassing = config%nopassing
        noshear = config%noshear
        pertfile = config%pertfile
        nonlin = config%nonlin
        bfac = config%bfac
        efac = config%efac
        inp_swi = config%inp_swi
        if (frequency_model < 0 .or. frequency_model > 1) then
            error stop "frequency_model must be 0 (legacy) or 1 (GC thin limit)"
        end if
        if (frequency_model == 1 .and. inp_swi /= 11) then
            error stop "frequency_model=1 requires direct GEQDSK input (inp_swi=11)"
        end if
        if (frequency_model == 1 .and. .not. has_direct_eqdsk_gc) then
            error stop "frequency_model=1 requires the standalone direct-GEQDSK build"
        end if
        if (frequency_model == 1 .and. supban) then
            error stop "frequency_model=1 and supban are mutually exclusive"
        end if
        if (config%inp_swi_pert < 0) then
            perturbation_switch = config%inp_swi
        else
            perturbation_switch = config%inp_swi_pert
        end if
        if (pertfile .and. perturbation_switch /= 8 .and. &
            perturbation_switch /= 9 .and. perturbation_switch /= 11) then
            error stop "inp_swi_pert must be 8/9 (.bc) or 11 (POTATO R-Z grid)"
        end if
        if (pertfile) then
            if (.not. perturbation_chart_is_compatible(inp_swi, &
                    perturbation_switch, len_trim(config%pert_angle_map) > 0)) then
                error stop "direct GEQDSK needs an R-Z perturbation or a Boozer angle map"
            end if
            call set_pert_angle_map_path(config%pert_angle_map)
        end if
        vsteps = config%vsteps
        if (config%mth_max_abs < -1) error stop "mth_max_abs must be -1 or nonnegative"
        mth_max_abs = config%mth_max_abs
        if (config%vmax_over_vth <= 0.0_dp) error stop "vmax_over_vth must be positive"
        vmax_over_vth = config%vmax_over_vth

        qi = config%qs * qe
        mi = config%ms * mu
        call set_mph(config%mph)
        call set_log_level(config%log_level)
    end subroutine set_config

    subroutine read_and_set_config(config_file)
        ! Set global control parameters directly from a file
        use do_magfie_mod, only: s, bfac, inp_swi, has_direct_eqdsk_gc
        use do_magfie_pert_mod, only: mph, set_mph, inp_swi_pert, &
            set_pert_angle_map_path
        use driftorbit, only: epsmn, m0, comptorque, magdrift, &
            magdrift_passing, frequency_model, nopassing, pertfile, &
            nonlin, efac, supban
        use logger, only: set_log_level
        use neort, only: vsteps, mth_max_abs, vmax_over_vth
        use neort_orbit, only: noshear
        use neort_profiles, only: M_t, vth
        use util, only: qe, mu, qi, mi

        character(len=*), intent(in) :: config_file
        real(dp) :: qs, ms
        integer :: log_level = 0

        character(len=1024) :: pert_angle_map

        namelist /params/ s, M_t, qs, ms, vth, epsmn, m0, mph, comptorque, supban, &
            magdrift, magdrift_passing, frequency_model, nopassing, noshear, pertfile, nonlin, bfac, efac, inp_swi, &
            inp_swi_pert, vsteps, mth_max_abs, vmax_over_vth, log_level, pert_angle_map

        mth_max_abs = -1
        vmax_over_vth = 4.0_dp
        inp_swi_pert = -1
        frequency_model = 0
        pert_angle_map = ''
        open (unit=9, file=config_file, status="old", form="formatted")
        read (9, nml=params)
        close (unit=9)

        if (magdrift_passing < 0) magdrift_passing = merge(1, 0, magdrift)
        if (mth_max_abs < -1) error stop "mth_max_abs must be -1 or nonnegative"
        if (vmax_over_vth <= 0.0_dp) error stop "vmax_over_vth must be positive"
        if (inp_swi_pert < 0) inp_swi_pert = inp_swi
        if (frequency_model < 0 .or. frequency_model > 1) then
            error stop "frequency_model must be 0 (legacy) or 1 (GC thin limit)"
        end if
        if (frequency_model == 1 .and. inp_swi /= 11) then
            error stop "frequency_model=1 requires direct GEQDSK input (inp_swi=11)"
        end if
        if (frequency_model == 1 .and. .not. has_direct_eqdsk_gc) then
            error stop "frequency_model=1 requires the standalone direct-GEQDSK build"
        end if
        if (frequency_model == 1 .and. supban) then
            error stop "frequency_model=1 and supban are mutually exclusive"
        end if
        if (pertfile .and. inp_swi_pert /= 8 .and. inp_swi_pert /= 9 .and. &
            inp_swi_pert /= 11) then
            error stop "inp_swi_pert must be 8/9 (.bc) or 11 (POTATO R-Z grid)"
        end if
        if (pertfile) then
            if (.not. perturbation_chart_is_compatible(inp_swi, inp_swi_pert, &
                    len_trim(pert_angle_map) > 0)) then
                error stop "direct GEQDSK needs an R-Z perturbation or a Boozer angle map"
            end if
            call set_pert_angle_map_path(pert_angle_map)
        end if

        M_t = M_t * efac / bfac
        qi = qs * qe
        mi = ms * mu
        call set_mph(mph)
        call set_log_level(log_level)
    end subroutine read_and_set_config

end module neort_config
