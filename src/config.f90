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
        logical :: magdrift = .false.  ! consider magnetic drift
        logical :: nopassing = .false.  ! neglect passing particles
        logical :: noshear = .false.  ! neglect magnetic shear term with dqds
        logical :: pertfile = .false.  ! read perturbation from file
        logical :: nonlin = .false.  ! do nonlinear calculation
        real(dp) :: bfac = 1.0_dp  ! scale B field by factor
        real(dp) :: efac = 1.0_dp  ! scale E field by factor
        integer :: inp_swi = 0  ! input switch for Boozer file
        integer :: vsteps = 0  ! integration steps in velocity space
        integer :: log_level = 0  ! how much to log
        !*! will be overwritten if using splines from plasma.in and profile.in files
    end type config_t

contains

    subroutine set_config(config)
        ! Set global control parameters via config struct
        use do_magfie_mod, only: s, bfac, inp_swi
        use do_magfie_pert_mod, only: mph, set_mph
        use driftorbit, only: epsmn, m0, comptorque, magdrift, nopassing, pertfile, nonlin, efac
        use logger, only: set_log_level
        use neort, only: vsteps
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
        magdrift = config%magdrift
        nopassing = config%nopassing
        noshear = config%noshear
        pertfile = config%pertfile
        nonlin = config%nonlin
        bfac = config%bfac
        efac = config%efac
        inp_swi = config%inp_swi
        vsteps = config%vsteps

        qi = config%qs * qe
        mi = config%ms * mu
        call set_mph(config%mph)
        call set_log_level(config%log_level)
    end subroutine set_config

    subroutine read_and_set_config(config_file)
        ! Set global control parameters directly from a file
        use do_magfie_mod, only: s, bfac, inp_swi
        use do_magfie_pert_mod, only: mph, set_mph
        use driftorbit, only: epsmn, m0, comptorque, magdrift, nopassing, pertfile, nonlin, efac
        use logger, only: set_log_level
        use neort, only: vsteps
        use neort_orbit, only: noshear
        use neort_profiles, only: M_t, vth
        use util, only: qe, mu, qi, mi

        character(len=*), intent(in) :: config_file
        real(dp) :: qs, ms
        integer :: log_level

        namelist /params/ s, M_t, qs, ms, vth, epsmn, m0, mph, comptorque, magdrift, &
            nopassing, noshear, pertfile, nonlin, bfac, efac, inp_swi, vsteps, log_level

        open (unit=9, file=config_file, status="old", form="formatted")
        read (9, nml=params)
        close (unit=9)

        M_t = M_t * efac / bfac
        qi = qs * qe
        mi = ms * mu
        call set_mph(mph)
        call set_log_level(log_level)
    end subroutine read_and_set_config

end module neort_config
