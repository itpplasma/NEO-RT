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
            validate_fp_background, conserving_fp_collisions, &
            validate_conserving_controls
        use do_magfie_mod, only: set_s, init_magfie_at_s, R0, psi_pr, q
        use do_magfie_pert_mod, only: init_magfie_pert_at_s, init_mph_from_shared
        use driftorbit, only: pertfile, nopassing, sign_vpar, etamin, etamax, efac, bfac, comptorque
        use neort, only: set_to_trapped_region
        use neort_freq, only: init_canon_freq_trapped_spline, init_canon_freq_passing_spline
        use neort_magfie, only: init_flux_surface_average
        use neort_profiles, only: init_plasma_at_s, init_profile_at_s, init_thermodynamic_forces

        real(dp), intent(in) :: s_val

        if (tp_fp_collisions) call validate_fp_controls()
        if (conserving_fp_collisions) call validate_conserving_controls()
        call set_s(s_val)

        call init_magfie_at_s()
        if (pertfile) call init_magfie_pert_at_s() ! else epsmn*exp(imun*(m0*th + mph*ph))
        call init_mph_from_shared()

        call init_plasma_at_s() ! overwrites qi, mi and vth from config
        if (tp_fp_collisions .or. conserving_fp_collisions) &
            call validate_fp_background()
        call init_profile_at_s(R0, efac, bfac) ! overwrites M_t

        call init_flux_surface_average(s_val)
        call init_canon_freq_trapped_spline() ! sets etamin and etamax
        if (.not. nopassing) call init_canon_freq_passing_spline()
        sign_vpar = 1
        call set_to_trapped_region(etamin, etamax) ! overwrites etamin and etamax
        ! psi_pr is the torodial flux at plasma boundary, fixed for all s
        ! psi_pr=const., q is set in do_magfie
        if (comptorque .or. tp_fp_collisions .or. conserving_fp_collisions) &
            call init_thermodynamic_forces(psi_pr, q)
    end subroutine neort_setup_at_s

    subroutine neort_compute_at_s(s_val, transport_data_out)
        use neort, only: compute_transport
        use neort_fp_mode, only: conserving_fp_collisions

        real(dp), intent(in) :: s_val
        type(transport_data_t), intent(out) :: transport_data_out

        if (conserving_fp_collisions) &
            error stop 'conserving FP requires the joint mixture API'
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
        use neort_fp_mode, only: conserving_fp_collisions
        real(dp), intent(in) :: s_val
        type(mixture_result_t), intent(out) :: result_
        integer :: k
        if (.not. mixture_enabled) error stop 'prepare a mixture before computing species'
        if (conserving_fp_collisions) then
            call compute_conserving_mixture_at_s(s_val, result_)
            return
        end if
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

    subroutine compute_conserving_mixture_at_s(s_val, result_)
        use neort_species, only: species_mixture
        use neort_conserving_native, only: conserving_species_t, &
            prepare_conserving_species, assemble_conserving_collisions
        use neort_conserving_solver, only: solve_conserving_system, &
            conserving_solver_diagnostics_t
        use neort_landau_pair_quadrature, only: landau_quadrature_t, &
            initialize_landau_quadrature
        use neort_fp_controls, only: conserving_speed_steps, &
            conserving_pitch_steps, conserving_umax, conserving_theta_order, &
            conserving_velocity_order, conserving_cm_order, &
            conserving_rel_order, conserving_polar_order, &
            conserving_azimuth_order, conserving_rank_rtol, conserving_rate_scale, &
            conserving_compatibility_rtol, validate_conserving_controls
        use driftorbit, only: mph
        use util, only: pi
        real(dp), intent(in) :: s_val
        type(mixture_result_t), intent(out) :: result_
        type(conserving_species_t), allocatable :: cache(:)
        type(conserving_solver_diagnostics_t) :: diagnostic
        type(landau_quadrature_t) :: quadrature
        real(dp), allocatable :: mass(:, :), streaming(:, :), precession(:, :)
        complex(dp), allocatable :: collision(:, :), source(:, :), coefficient(:)
        complex(dp), allocatable :: cross_coefficient(:, :)
        integer, allocatable :: offsets(:), status(:)
        integer :: species, receivers, first, last, ndof, info, branch, driver
        complex(dp) :: work(3)

        call validate_conserving_controls()
        receivers = size(species_mixture%species)
        allocate(result_%species(receivers), cache(receivers), status(receivers))
        ! Every native callback runs only after complete per-thread geometry,
        ! profile and signed species setup. Pair assembly subsequently reads
        ! these immutable caches without touching native global variables.
        !$omp parallel do default(shared) private(species) schedule(dynamic, 1)
        do species = 1, receivers
            call prepare_species(species)
        end do
        !$omp end parallel do
        do species = 1, receivers
            if (status(species) /= 0) then
                write(*, '(a,i0,a,i0)') 'Conserving FP native species ', &
                    species, ' rejected; status=', status(species)
                error stop 'see doc/conserving_fp.md for the supported projection'
            end if
        end do
        allocate(offsets(receivers + 1))
        offsets(1) = 1
        do species = 1, receivers
            offsets(species + 1) = offsets(species) + cache(species)%basis%ndof
        end do
        ndof = offsets(receivers + 1) - 1
        allocate(mass(ndof, ndof), streaming(ndof, ndof), precession(ndof, ndof))
        allocate(collision(ndof, ndof), source(ndof, receivers))
        mass = 0.0_dp
        streaming = 0.0_dp
        precession = 0.0_dp
        source = (0.0_dp, 0.0_dp)
        do species = 1, receivers
            first = offsets(species)
            last = offsets(species + 1) - 1
            mass(first:last, first:last) = cache(species)%mass
            streaming(first:last, first:last) = cache(species)%streaming
            precession(first:last, first:last) = cache(species)%precession
            source(first:last, species) = cmplx(0.0_dp, &
                -cache(species)%beta, dp)*(cache(species)%A1* &
                sum(cache(species)%hamiltonian, dim=2) + cache(species)%A2* &
                sum(cache(species)%hamiltonian2, dim=2))
        end do
        call initialize_landau_quadrature(conserving_cm_order, &
            conserving_rel_order, conserving_polar_order, &
            conserving_azimuth_order, quadrature)
        call assemble_conserving_collisions(cache, conserving_theta_order, &
            quadrature, collision, info)
        if (info /= 0) then
            write(*, '(a,i0)') 'Conserving FP pair assembly failed; status=', info
            error stop 'invalid common-temperature Landau assembly'
        end if
        collision = conserving_rate_scale*collision
        call solve_conserving_system(mass, collision, streaming, precession, &
            source, coefficient, cross_coefficient, diagnostic, info, &
            conserving_rank_rtol, conserving_compatibility_rtol)
        if (info /= 0) then
            write(*, '(a,i0,a,es14.6)') 'Conserving FP joint solve failed; status=', &
                info, ' physical compatibility=', diagnostic%compatibility
            error stop 'incompatible or unresolved conserving FP response'
        end if
        result_%torque_total = 0.0_dp
        do species = 1, receivers
            first = offsets(species)
            last = offsets(species + 1) - 1
            do branch = 1, 3
                work(branch) = dot_product(cache(species)%hamiltonian(:, branch), &
                    coefficient(first:last))
            end do
            call store_response(species, work)
            result_%torque_total = result_%torque_total + &
                [result_%species(species)%transport%torque%Tco, &
                 result_%species(species)%transport%torque%Tctr, &
                 result_%species(species)%transport%torque%Tt]
        end do

    contains

        subroutine prepare_species(index)
            use neort_species, only: selected_species, mixture_at_s
            use neort_profiles, only: A1, A2, dOm_tEds
            use neort_mixture_collisions, only: selected_thermal_deflection_rate, &
                selected_reference_frequency
            use neort, only: check_magfie
            integer, intent(in) :: index
            selected_species = index
            call neort_setup_at_s(s_val)
            call prepare_conserving_species(conserving_speed_steps, &
                conserving_pitch_steps, conserving_umax, conserving_theta_order, &
                conserving_velocity_order, cache(index), status(index))
            if (status(index) /= 0) return
            call check_magfie(result_%species(index)%magfie)
            result_%species(index)%species = mixture_at_s(index)
            result_%species(index)%collision_thermal_deflection_rate = &
                selected_thermal_deflection_rate()
            result_%species(index)%collision_reference_frequency = &
                selected_reference_frequency()
            result_%species(index)%domega_e_ds = dOm_tEds
            result_%species(index)%force = [A1, A2]
        end subroutine prepare_species

        subroutine store_response(index, physical_work)
            use driftorbit, only: comptorque
            type(transport_data_t) :: data
            integer, intent(in) :: index
            complex(dp), intent(in) :: physical_work(3)
            real(dp) :: torque(3)
            integer :: k
            torque = -real(mph, dp)*pi*aimag(physical_work)
            if (.not. comptorque) torque = 0.0_dp
            data%susceptibility_valid = .false.
            data%conserving_active = .true.
            data%conserving_work = physical_work
            allocate(data%conserving_cross_work(3, receivers))
            allocate(data%conserving_cross_work_qualified(receivers))
            data%conserving_cross_work_qualified = diagnostic%cross_qualified
            data%conserving_cross_work = (0.0_dp, 0.0_dp)
            do driver = 1, receivers
                if (.not. diagnostic%cross_qualified(driver)) cycle
                do k = 1, 3
                    data%conserving_cross_work(k, driver) = dot_product( &
                        cache(index)%hamiltonian(:, k), &
                        cross_coefficient(first:last, driver))
                end do
            end do
            data%conserving_rank = diagnostic%rank
            data%conserving_nullity = diagnostic%nullity
            data%conserving_singular_cutoff = diagnostic%singular_cutoff
            data%conserving_residual = diagnostic%residual
            data%conserving_compatibility = diagnostic%compatibility
            data%conserving_mass_condition = diagnostic%mass_condition
            data%conserving_entropy = diagnostic%entropy
            data%summary%M_t = result_%species(index)%magfie%params%M_t
            data%summary%Dco = 0.0_dp
            data%summary%Dctr = 0.0_dp
            data%summary%Dt = 0.0_dp
            data%torque%has_torque = comptorque
            data%torque%s = s_val
            data%torque%dVds = result_%species(index)%magfie%params%dVds
            data%torque%M_t = data%summary%M_t
            data%torque%Tco = torque(2)
            data%torque%Tctr = torque(3)
            data%torque%Tt = torque(1)
            allocate(data%harmonics(2))
            do k = 1, 2
                data%harmonics(k)%mth = 0
                data%harmonics(k)%Dresco = 0.0_dp
                data%harmonics(k)%Dresctr = 0.0_dp
                data%harmonics(k)%Drest = 0.0_dp
                data%harmonics(k)%Tresco = 0.0_dp
                data%harmonics(k)%Tresctr = 0.0_dp
                data%harmonics(k)%Trest = 0.0_dp
                data%harmonics(k)%vminp_over_vth = 0.0_dp
                data%harmonics(k)%vmint_over_vth = 0.0_dp
                ! No finite integration cutoff: these legacy velocity-range
                ! columns are undefined, like the auxiliary D columns.
                data%harmonics(k)%vmaxp_over_vth = 0.0_dp
                data%harmonics(k)%vmaxt_over_vth = 0.0_dp
            end do
            data%harmonics(1)%Trest = torque(1)
            data%harmonics(2)%mth = -nint( &
                mph*result_%species(index)%magfie%params%q)
            data%harmonics(2)%Tresco = torque(2)
            data%harmonics(2)%Tresctr = torque(3)
            result_%species(index)%transport = data
        end subroutine store_response
    end subroutine compute_conserving_mixture_at_s

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
        use neort_fp_mode, only: tp_fp_collisions, conserving_fp_collisions
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
        if (conserving_fp_collisions) &
            error stop 'conserving FP requires the joint mixture API'
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
