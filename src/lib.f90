module neort_lib
    use iso_fortran_env, only: dp => real64
    use neort_config, only: config_t
    use neort_datatypes, only: transport_data_t, mixture_result_t
    use neort_species, only: species_state_t, mixture_t

    implicit none

    private

    public :: config_t
    public :: transport_data_t
    public :: neort_init
    public :: neort_prepare_splines
    public :: neort_setup_at_s
    public :: neort_compute_at_s
    public :: neort_compute_no_splines
    public :: species_state_t, mixture_t, mixture_result_t
    public :: neort_prepare_mixture, neort_compute_mixture_at_s

    interface neort_init
        module procedure init_from_config_file
        module procedure init_from_config_struct
    end interface neort_init

    interface neort_prepare_splines
        module procedure prepare_splines_from_files
        module procedure prepare_splines_from_memory
    end interface neort_prepare_splines

contains

    subroutine init_from_config_file(config_file, boozer_file, boozer_pert_file)
        use do_magfie_mod, only: read_boozer_file
        use do_magfie_pert_mod, only: read_boozer_pert_file
        use driftorbit, only: pertfile, nonlin
        use neort_config, only: read_and_set_config
        use thetadata_mod, only: init_attenuation_data
        use util, only: check_file

        character(len=*), intent(in) :: config_file
        character(len=*), intent(in) :: boozer_file
        character(len=*), intent(in), optional :: boozer_pert_file

        call check_file(config_file, .true., " could not be found")
        call check_file(boozer_file, .true., " could not be found")

        call read_and_set_config(config_file)
        call read_boozer_file(boozer_file)

        if (present(boozer_pert_file) .and. pertfile) then
            call check_file(boozer_pert_file, .true., " required if pertfile = .true.")
            call read_boozer_pert_file(boozer_pert_file)
        end if

        if (nonlin) call init_attenuation_data()
    end subroutine init_from_config_file

    subroutine init_from_config_struct(config, boozer_file, boozer_pert_file)
        use do_magfie_mod, only: read_boozer_file
        use do_magfie_pert_mod, only: read_boozer_pert_file
        use driftorbit, only: pertfile, nonlin
        use neort_config, only: set_config
        use thetadata_mod, only: init_attenuation_data
        use util, only: check_file

        type(config_t), intent(in) :: config
        character(len=*), intent(in) :: boozer_file
        character(len=*), intent(in), optional :: boozer_pert_file

        call set_config(config)

        call check_file(boozer_file, .true., " could not be found")
        call read_boozer_file(boozer_file)

        if (present(boozer_pert_file) .and. pertfile) then
            call check_file(boozer_pert_file, .true., " required if pertfile = .true.")
            call read_boozer_pert_file(boozer_pert_file)
        end if

        if (nonlin) call init_attenuation_data()
    end subroutine init_from_config_struct

    subroutine prepare_splines_from_files(plasma_file, profile_file)
        use driftorbit, only: nonlin, comptorque
        use neort_profiles, only: read_plasma_input, prepare_plasma_splines, prepare_profile_splines
        use util, only: check_file, readdata

        character(len=*), intent(in) :: plasma_file
        character(len=*), intent(in) :: profile_file

        integer :: nplasma
        real(dp) :: am1, am2, Z1, Z2
        real(dp), allocatable :: plasma_data(:, :)
        real(dp), allocatable :: profile_data(:, :)

        call check_file(plasma_file, nonlin .or. comptorque, &
            " required in spline mode (nonlin = .true. or comptorque = .true.)")
        call check_file(profile_file, nonlin, " required in spline mode (nonlin = .true.)")

        call read_plasma_input(plasma_file, nplasma, am1, am2, Z1, Z2, plasma_data)
        call prepare_plasma_splines(nplasma, am1, am2, Z1, Z2, plasma_data)
        deallocate (plasma_data)

        call readdata(profile_file, 2, profile_data)
        call prepare_profile_splines(profile_data)
        deallocate (profile_data)
    end subroutine prepare_splines_from_files

    subroutine prepare_splines_from_memory(nplasma, am1, am2, Z1, Z2, plasma_data, profile_data)
        use neort_profiles, only: prepare_plasma_splines, prepare_profile_splines

        integer, intent(in) :: nplasma
        real(dp), intent(in) :: am1, am2, Z1, Z2
        real(dp), intent(in) :: plasma_data(:, :)
        real(dp), intent(in) :: profile_data(:, :)

        call prepare_plasma_splines(nplasma, am1, am2, Z1, Z2, plasma_data)
        call prepare_profile_splines(profile_data)
    end subroutine prepare_splines_from_memory

    subroutine neort_setup_at_s(s_val)
        use neort_fp_controls, only: tp_fp_collisions, validate_fp_controls, &
            validate_fp_background
        use do_magfie_mod, only: set_s, init_magfie_at_s, R0, psi_pr, q
        use do_magfie_pert_mod, only: init_magfie_pert_at_s, init_mph_from_shared
        use driftorbit, only: pertfile, nopassing, sign_vpar, etamin, etamax, efac, bfac, comptorque
        use neort, only: set_to_trapped_region
        use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline
        use neort_magfie, only: init_flux_surface_average
        use neort_profiles, only: init_plasma_at_s, init_profile_at_s, init_thermodynamic_forces

        real(dp), intent(in) :: s_val

        if (tp_fp_collisions) call validate_fp_controls()
        call set_s(s_val)

        call init_magfie_at_s()
        if (pertfile) call init_magfie_pert_at_s() ! else epsmn*exp(imun*(m0*th + mph*ph))
        call init_mph_from_shared()

        call init_plasma_at_s() ! overwrites qi, mi and vth from config
        if (tp_fp_collisions) call validate_fp_background()
        call init_profile_at_s(R0, efac, bfac) ! overwrites M_t

        call init_flux_surface_average(s_val)
        call init_canon_freq_trapped_spline() ! sets etamin and etamax
        if (.not. nopassing) call init_canon_freq_passing_spline()
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax) ! overwrites etamin and etamax
        ! psi_pr is the torodial flux at plasma boundary, fixed for all s
        ! psi_pr=const., q is set in do_magfie
        if (comptorque .or. tp_fp_collisions) &
            call init_thermodynamic_forces(psi_pr, q)
    end subroutine neort_setup_at_s

    subroutine neort_compute_at_s(s_val, transport_data_out)
        use neort, only: compute_transport

        real(dp), intent(in) :: s_val
        type(transport_data_t), intent(out) :: transport_data_out

        call neort_setup_at_s(s_val)
        call compute_transport(transport_data_out)
    end subroutine neort_compute_at_s

    subroutine neort_prepare_mixture(mixture)
        use neort_species, only: configure_mixture
        use neort_config, only: validate_mixture_mode
        use driftorbit, only: nonlin
        type(mixture_t), intent(in) :: mixture
        call configure_mixture(mixture)
        call validate_mixture_mode(nonlin)
    end subroutine neort_prepare_mixture

    subroutine neort_compute_mixture_at_s(s_val, result_)
        use neort_species, only: mixture_enabled, species_mixture
        real(dp), intent(in) :: s_val
        type(mixture_result_t), intent(out) :: result_
        integer :: k
        if (.not. mixture_enabled) error stop 'prepare a mixture before computing species'
        allocate(result_%species(size(species_mixture%species)))
        !$omp parallel do default(shared) private(k) schedule(dynamic, 1)
        do k = 1, size(result_%species)
            call compute_mixture_species_at_s(s_val, k, result_%species(k))
        end do
        !$omp end parallel do
        result_%torque_total = 0.0_dp
        do k = 1, size(result_%species)
            result_%torque_total = result_%torque_total + &
                [result_%species(k)%transport%torque%Tco, &
                 result_%species(k)%transport%torque%Tctr, &
                 result_%species(k)%transport%torque%Tt]
        end do
    end subroutine neort_compute_mixture_at_s

    subroutine compute_mixture_species_at_s(s_val, species_index, result_)
        use neort_species, only: selected_species, mixture_at_s
        use neort_profiles, only: A1, A2, dOm_tEds
        use neort_mixture_collisions, only: selected_thermal_deflection_rate, &
            selected_reference_frequency
        use neort_datatypes, only: species_result_t
        use neort, only: check_magfie, compute_transport
        real(dp), intent(in) :: s_val
        integer, intent(in) :: species_index
        type(species_result_t), intent(out) :: result_
        selected_species = species_index
        call neort_setup_at_s(s_val)
        call compute_transport(result_%transport)
        call check_magfie(result_%magfie)
        result_%species = mixture_at_s(species_index)
        result_%collision_thermal_deflection_rate = selected_thermal_deflection_rate()
        result_%collision_reference_frequency = selected_reference_frequency()
        result_%domega_e_ds = dOm_tEds
        if (result_%transport%torque%has_torque) result_%force = [A1, A2]
    end subroutine compute_mixture_species_at_s

    subroutine neort_compute_no_splines(transport_data_out)
        use neort_fp_mode, only: tp_fp_collisions
        use do_magfie_mod, only: init_magfie_at_s, R0, psi_pr, q, s
        use do_magfie_pert_mod, only: init_magfie_pert_at_s, init_mph_from_shared
        use driftorbit, only: pertfile, nopassing, sign_vpar, etamin, etamax, comptorque
        use neort, only: compute_transport, set_to_trapped_region
        use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline
        use neort_magfie, only: init_flux_surface_average
        use neort_profiles, only: init_profiles, init_thermodynamic_forces

        type(transport_data_t), intent(out) :: transport_data_out

        if (tp_fp_collisions) &
            error stop 'TP FP requires explicit mixture profiles; '// &
                'use neort_compute_at_s'
        call init_magfie_at_s()
        if (pertfile) call init_magfie_pert_at_s()
        call init_mph_from_shared()

        call init_profiles(R0)

        call init_flux_surface_average(s)
        call init_canon_freq_trapped_spline()
        if (.not. nopassing) call init_canon_freq_passing_spline()
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax)
        if (comptorque) call init_thermodynamic_forces(psi_pr, q)

        call compute_transport(transport_data_out)
    end subroutine neort_compute_no_splines

end module neort_lib
