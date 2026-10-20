module neort_species
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_is_finite, ieee_value, ieee_quiet_nan
    implicit none
    private
    public :: species_state_t, mixture_t, configure_mixture, reset_mixture
    public :: read_species_input, evaluate_mixture_at_s, mixture_rotation
    public :: mixture_enabled, species_mixture, mixture_at_s, selected_species
    public :: common_omega_e, common_domega_e_ds

    type :: species_state_t
        character(len=32) :: name = ''
        real(dp) :: mass_kg = 0.0_dp
        real(dp) :: charge_e = 0.0_dp
        real(dp) :: density_m3 = 0.0_dp
        real(dp) :: temperature_ev = 0.0_dp
        real(dp) :: density_gradient_m3 = 0.0_dp
        real(dp) :: temperature_gradient_ev = 0.0_dp
    end type species_state_t

    type :: mixture_t
        real(dp) :: reference_s = 0.0_dp
        real(dp) :: omega_e = 0.0_dp
        real(dp) :: domega_e_ds = 0.0_dp
        real(dp) :: coulomb_log = 0.0_dp
        type(species_state_t), allocatable :: species(:)
    end type mixture_t

    logical :: mixture_enabled = .false.
    type(mixture_t) :: species_mixture
    type(species_state_t), allocatable :: mixture_at_s(:)
    integer :: selected_species = 1
    real(dp) :: common_omega_e = 0.0_dp, common_domega_e_ds = 0.0_dp
    real(dp), parameter :: neutrality_rtol = 1.0e-12_dp
    !$omp threadprivate(mixture_at_s, selected_species, common_omega_e, common_domega_e_ds)

contains

    subroutine configure_mixture(mixture)
        type(mixture_t), intent(in) :: mixture
        if (.not. allocated(mixture%species)) error stop 'mixture species are required'
        if (size(mixture%species) < 1) error stop 'mixture species count must be positive'
        if (.not. all(ieee_is_finite([mixture%reference_s, mixture%omega_e, &
                                     mixture%domega_e_ds, mixture%coulomb_log]))) &
            error stop 'mixture radial coordinate and physical rotation must be finite'
        if (mixture%coulomb_log <= 0.0_dp) error stop 'mixture Coulomb logarithm must be positive'
        call validate_species(mixture%species)
        species_mixture = mixture
        mixture_enabled = .true.
        selected_species = 1
    end subroutine configure_mixture

    subroutine reset_mixture()
        if (allocated(species_mixture%species)) deallocate(species_mixture%species)
        if (allocated(mixture_at_s)) deallocate(mixture_at_s)
        mixture_enabled = .false.
        selected_species = 1
    end subroutine reset_mixture

    subroutine read_species_input(path)
        character(len=*), intent(in) :: path
        integer :: n_species, unit, k
        real(dp) :: reference_s, omega_e, domega_e_ds, coulomb_log
        character(len=32), allocatable :: species_name(:)
        real(dp), allocatable :: mass_kg(:), charge_e(:), density_m3(:), temperature_ev(:)
        real(dp), allocatable :: density_gradient_m3(:), temperature_gradient_ev(:)
        type(mixture_t) :: mixture
        namelist /species_count/ n_species
        namelist /species/ reference_s, omega_e, domega_e_ds, coulomb_log, species_name, mass_kg, &
            charge_e, density_m3, temperature_ev, density_gradient_m3, temperature_gradient_ev
        n_species = 0
        open(newunit=unit, file=path, status='old', action='read')
        read(unit, nml=species_count)
        if (n_species < 1) error stop 'n_species must be positive'
        allocate(species_name(n_species), mass_kg(n_species), charge_e(n_species), &
            density_m3(n_species), temperature_ev(n_species), &
            density_gradient_m3(n_species), temperature_gradient_ev(n_species))
        call initialize_input_arrays()
        read(unit, nml=species)
        close(unit)
        mixture%reference_s = reference_s
        mixture%omega_e = omega_e
        mixture%domega_e_ds = domega_e_ds
        mixture%coulomb_log = coulomb_log
        allocate(mixture%species(n_species))
        do k = 1, n_species
            mixture%species(k) = species_state_t(species_name(k), mass_kg(k), charge_e(k), &
                density_m3(k), temperature_ev(k), density_gradient_m3(k), &
                temperature_gradient_ev(k))
        end do
        call configure_mixture(mixture)
    contains
        subroutine initialize_input_arrays()
            real(dp) :: missing
            missing = ieee_value(0.0_dp, ieee_quiet_nan)
            reference_s = missing
            omega_e = missing
            domega_e_ds = missing
            coulomb_log = missing
            species_name = ''
            mass_kg = missing
            charge_e = missing
            density_m3 = missing
            temperature_ev = missing
            density_gradient_m3 = missing
            temperature_gradient_ev = missing
        end subroutine initialize_input_arrays
    end subroutine read_species_input

    subroutine evaluate_mixture_at_s(s_val)
        real(dp), intent(in) :: s_val
        real(dp) :: delta_s
        integer :: k
        if (.not. mixture_enabled) error stop 'mixture is not configured'
        if (.not. ieee_is_finite(s_val)) error stop 'mixture surface must be finite'
        delta_s = s_val - species_mixture%reference_s
        mixture_at_s = species_mixture%species
        do k = 1, size(mixture_at_s)
            call species_affine_leaf(species_mixture%species(k)%density_m3, &
                species_mixture%species(k)%density_gradient_m3, delta_s, &
                mixture_at_s(k)%density_m3)
            call species_affine_leaf(species_mixture%species(k)%temperature_ev, &
                species_mixture%species(k)%temperature_gradient_ev, delta_s, &
                mixture_at_s(k)%temperature_ev)
        end do
        call species_affine_leaf(species_mixture%omega_e, species_mixture%domega_e_ds, &
            delta_s, common_omega_e)
        common_domega_e_ds = species_mixture%domega_e_ds
        if (.not. ieee_is_finite(common_omega_e)) error stop 'evaluated rotation must be finite'
        call validate_species(mixture_at_s)
    end subroutine evaluate_mixture_at_s

    subroutine mixture_rotation(radius, speed, dspeed_ds, efac, bfac, mach, dmach_ds)
        real(dp), intent(in) :: radius, speed, dspeed_ds, efac, bfac
        real(dp), intent(out) :: mach, dmach_ds
        call species_rotation_leaf(radius, speed, dspeed_ds, &
            common_omega_e*efac/bfac, common_domega_e_ds*efac/bfac, mach, dmach_ds)
    end subroutine mixture_rotation

    subroutine validate_species(species)
        type(species_state_t), intent(in) :: species(:)
        integer :: k, j
        do k = 1, size(species)
            call validate_particle(species(k))
            do j = 1, k - 1
                if (species(k)%name == species(j)%name) &
                    error stop 'mixture species names must be unique'
            end do
        end do
        call validate_neutrality(species)
    end subroutine validate_species

    subroutine validate_particle(species)
        type(species_state_t), intent(in) :: species
        if (len_trim(species%name) == 0) error stop 'mixture species name is required'
        if (scan(trim(species%name), '/ ') /= 0) &
            error stop 'mixture species name cannot contain slash or spaces'
        if (.not. all(ieee_is_finite([species%mass_kg, species%charge_e, &
            species%density_m3, species%temperature_ev, species%density_gradient_m3, &
            species%temperature_gradient_ev]))) error stop 'mixture species must be finite'
        if (species%mass_kg <= 0.0_dp) error stop 'mixture mass must be positive'
        if (species%charge_e == 0.0_dp) error stop 'mixture charge must be nonzero'
        if (species%density_m3 <= 0.0_dp) error stop 'mixture density must be positive'
        if (species%temperature_ev <= 0.0_dp) error stop 'mixture temperature must be positive'
    end subroutine validate_particle

    subroutine validate_neutrality(species)
        type(species_state_t), intent(in) :: species(:)
        real(dp) :: charge_density(size(species)), charge_gradient(size(species))
        charge_density = species%charge_e*species%density_m3
        charge_gradient = species%charge_e*species%density_gradient_m3
        if (abs(sum(charge_density)) > neutrality_rtol*sum(abs(charge_density))) &
            error stop 'mixture density must be charge neutral'
        if (abs(sum(charge_gradient)) > neutrality_rtol*sum(abs(charge_gradient))) &
            error stop 'mixture density gradient must be charge neutral'
    end subroutine validate_neutrality

    include 'generated/species_leaf.inc'
end module neort_species
