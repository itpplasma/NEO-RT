module neort_lib
    use iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: neort_init
    public :: neort_prepare_splines
    public :: neort_prepare_splines_from_files
    public :: neort_compute_at_s
    public :: neort_compute_no_splines
    public :: neort_deinit

contains

    subroutine neort_init(config_file_base, boozer_file, boozer_pert_file)
        use do_magfie_mod, only: read_boozer_file
        use do_magfie_pert_mod, only: read_boozer_pert_file
        use driftorbit, only: pertfile, nonlin
        use neort, only: read_and_set_control
        use thetadata_mod, only: init_attenuation_data
        use util, only: check_file

        character(len=*), intent(in) :: config_file_base
        character(len=*), intent(in) :: boozer_file
        character(len=*), intent(in), optional :: boozer_pert_file

        call check_file(trim(config_file_base) // ".in", .true., ".in could not be found")
        call check_file(boozer_file, .true., " could not be found")

        call read_and_set_control(config_file_base)
        call read_boozer_file(boozer_file)

        if (present(boozer_pert_file) .and. pertfile) then
            call check_file(boozer_pert_file, .true., " required if pertfile = .true.")
            call read_boozer_pert_file(boozer_pert_file)
        end if

        if (nonlin) call init_attenuation_data()
    end subroutine neort_init

    subroutine neort_prepare_splines(nplasma, am1, am2, Z1, Z2, plasma_data, profile_data)
        use neort_profiles, only: prepare_plasma_splines, prepare_profile_splines

        integer, intent(in) :: nplasma
        real(dp), intent(in) :: am1, am2, Z1, Z2
        real(dp), intent(in) :: plasma_data(:, :)
        real(dp), intent(in) :: profile_data(:, :)

        call prepare_plasma_splines(nplasma, am1, am2, Z1, Z2, plasma_data)
        call prepare_profile_splines(profile_data)
    end subroutine neort_prepare_splines

    subroutine neort_prepare_splines_from_files(plasma_file, profile_file)
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
    end subroutine neort_prepare_splines_from_files

    subroutine neort_compute_at_s(s_val, transport_data_out)
        use do_magfie_mod, only: set_s, init_magfie_at_s, R0, psi_pr, q
        use do_magfie_pert_mod, only: init_magfie_pert_at_s, init_mph_from_shared
        use driftorbit, only: pertfile, nopassing, sign_vpar, etamin, etamax, efac, bfac, comptorque
        use neort, only: compute_transport, set_to_trapped_region
        use neort_datatypes, only: transport_data_t
        use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline
        use neort_magfie, only: init_flux_surface_average
        use neort_profiles, only: init_plasma_at_s, init_profile_at_s, init_thermodynamic_forces

        real(dp), intent(in) :: s_val
        type(transport_data_t), intent(out) :: transport_data_out

        call set_s(s_val)

        call init_magfie_at_s()
        if (pertfile) call init_magfie_pert_at_s()
        call init_mph_from_shared()

        call init_plasma_at_s()
        call init_profile_at_s(R0, efac, bfac)

        call init_flux_surface_average(s_val)
        call init_canon_freq_trapped_spline()
        if (.not. nopassing) call init_canon_freq_passing_spline()
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax)
        if (comptorque) call init_thermodynamic_forces(psi_pr, q)

        call compute_transport(transport_data_out)
    end subroutine neort_compute_at_s

    subroutine neort_compute_no_splines(transport_data_out)
        use do_magfie_mod, only: init_magfie_at_s, R0, psi_pr, q, s
        use do_magfie_pert_mod, only: init_magfie_pert_at_s, init_mph_from_shared
        use driftorbit, only: pertfile, nopassing, sign_vpar, etamin, etamax, comptorque
        use neort, only: compute_transport, set_to_trapped_region
        use neort_datatypes, only: transport_data_t
        use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline
        use neort_magfie, only: init_flux_surface_average
        use neort_profiles, only: init_profiles, init_thermodynamic_forces

        type(transport_data_t), intent(out) :: transport_data_out

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

    subroutine neort_deinit()
        use neort_profiles, only: plasma_spl_coeff, Mt_spl_coeff

        if (allocated(plasma_spl_coeff)) deallocate (plasma_spl_coeff)
        if (allocated(Mt_spl_coeff)) deallocate (Mt_spl_coeff)
    end subroutine neort_deinit

end module neort_lib
