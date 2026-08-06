module neort_gc_models
    !! Field, potential, and invariant interfaces for guiding-center orbits.
    !! Concrete equilibrium adapters live outside this physics layer.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_dynamics, only: gc_field_sample_t
    use util, only: c

    implicit none
    private

    integer, parameter, public :: GC_MODEL_SUCCESS = 0
    integer, parameter, public :: GC_MODEL_INVALID_STATE = 1

    type, abstract, public :: gc_field_t
    contains
        procedure(gc_field_evaluate_i), deferred :: evaluate
    end type gc_field_t

    type, abstract, public :: gc_potential_t
    contains
        procedure(gc_potential_evaluate_i), deferred :: evaluate
    end type gc_potential_t

    type, extends(gc_potential_t), public :: gc_zero_potential_t
    contains
        procedure :: evaluate => evaluate_zero_potential
    end type gc_zero_potential_t

    type, extends(gc_potential_t), public :: gc_linear_flux_potential_t
        !! Dimensionless potential q_s Phi/E_ref, E_ref=m v_ref**2/2.
        !! The source convention is dPhi/dpsi_pol=+Omega_E/c.
        real(dp) :: coefficient = 0.0_dp
        real(dp) :: psi_reference = 0.0_dp
    contains
        procedure :: evaluate => evaluate_linear_flux_potential
    end type gc_linear_flux_potential_t

    type, public :: gc_invariants_t
        real(dp) :: energy = 0.0_dp
        real(dp) :: magnetic_moment = 0.0_dp
        real(dp) :: canonical_flux = 0.0_dp
    end type gc_invariants_t

    abstract interface
        subroutine gc_field_evaluate_i(self, position, sample, status)
            import :: dp, gc_field_t, gc_field_sample_t
            class(gc_field_t), intent(in) :: self
            real(dp), intent(in) :: position(3)
            type(gc_field_sample_t), intent(out) :: sample
            integer, intent(out) :: status
        end subroutine gc_field_evaluate_i

        pure subroutine gc_potential_evaluate_i(self, position, sample, &
                potential, gradient, status)
            import :: dp, gc_potential_t, gc_field_sample_t
            class(gc_potential_t), intent(in) :: self
            real(dp), intent(in) :: position(3)
            type(gc_field_sample_t), intent(in) :: sample
            real(dp), intent(out) :: potential, gradient(3)
            integer, intent(out) :: status
        end subroutine gc_potential_evaluate_i
    end interface

    public :: make_linear_flux_potential
    public :: invariants_from_state, state_from_invariants
    public :: canonical_flux_from_state

contains

    pure function make_linear_flux_potential(omega_e, particle_charge, &
            particle_mass, reference_velocity, psi_reference, amplitude) &
            result(model)
        real(dp), intent(in) :: omega_e, particle_charge, particle_mass
        real(dp), intent(in) :: reference_velocity, psi_reference
        real(dp), intent(in), optional :: amplitude
        type(gc_linear_flux_potential_t) :: model
        real(dp) :: scale

        scale = 1.0_dp
        if (present(amplitude)) scale = amplitude
        model%coefficient = scale*2.0_dp*particle_charge*omega_e &
            /(particle_mass*reference_velocity**2*c)
        model%psi_reference = psi_reference
    end function make_linear_flux_potential

    pure subroutine evaluate_zero_potential(self, position, sample, &
            potential, gradient, status)
        class(gc_zero_potential_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_field_sample_t), intent(in) :: sample
        real(dp), intent(out) :: potential, gradient(3)
        integer, intent(out) :: status

        associate (unused_self => self, unused_position => position, &
                unused_sample => sample)
        end associate
        potential = 0.0_dp
        gradient = 0.0_dp
        status = GC_MODEL_SUCCESS
    end subroutine evaluate_zero_potential

    pure subroutine evaluate_linear_flux_potential(self, position, sample, &
            potential, gradient, status)
        class(gc_linear_flux_potential_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_field_sample_t), intent(in) :: sample
        real(dp), intent(out) :: potential, gradient(3)
        integer, intent(out) :: status

        associate (unused_position => position)
        end associate
        potential = self%coefficient*(sample%psi - self%psi_reference)
        gradient = self%coefficient*sample%grad_psi
        status = GC_MODEL_SUCCESS
    end subroutine evaluate_linear_flux_potential

    pure subroutine invariants_from_state(sample, potential, rho0, &
            orbit_width_scale, p, xi, invariants, status)
        type(gc_field_sample_t), intent(in) :: sample
        real(dp), intent(in) :: potential, rho0, orbit_width_scale, p, xi
        type(gc_invariants_t), intent(out) :: invariants
        integer, intent(out) :: status

        invariants = gc_invariants_t()
        status = GC_MODEL_INVALID_STATE
        if (sample%bmod <= 0.0_dp .or. p <= 0.0_dp .or. abs(xi) > 1.0_dp) return
        invariants%energy = p**2 + potential
        invariants%magnetic_moment = p**2*(1.0_dp - xi**2)/sample%bmod
        invariants%canonical_flux = canonical_flux_from_state(sample, rho0, &
            orbit_width_scale, p, xi)
        status = GC_MODEL_SUCCESS
    end subroutine invariants_from_state

    pure subroutine state_from_invariants(sample, potential, invariants, &
            parallel_sign, p, xi, status)
        type(gc_field_sample_t), intent(in) :: sample
        real(dp), intent(in) :: potential
        type(gc_invariants_t), intent(in) :: invariants
        integer, intent(in) :: parallel_sign
        real(dp), intent(out) :: p, xi
        integer, intent(out) :: status

        real(dp) :: p_squared, xi_squared

        p = 0.0_dp
        xi = 0.0_dp
        status = GC_MODEL_INVALID_STATE
        if (sample%bmod <= 0.0_dp .or. parallel_sign == 0) return
        p_squared = invariants%energy - potential
        if (p_squared <= 0.0_dp) return
        xi_squared = 1.0_dp &
            - invariants%magnetic_moment*sample%bmod/p_squared
        if (xi_squared < -100.0_dp*epsilon(xi_squared)) return
        p = sqrt(p_squared)
        xi = sign(sqrt(max(0.0_dp, xi_squared)), real(parallel_sign, dp))
        status = GC_MODEL_SUCCESS
    end subroutine state_from_invariants

    pure function canonical_flux_from_state(sample, rho0, orbit_width_scale, &
            p, xi) result(canonical_flux)
        type(gc_field_sample_t), intent(in) :: sample
        real(dp), intent(in) :: rho0, orbit_width_scale, p, xi
        real(dp) :: canonical_flux

        ! POTATO/libneo source convention: psi_star=psi+rho p_parallel h_phi.
        canonical_flux = sample%psi &
            + orbit_width_scale*rho0*p*xi*sample%hcov(2)
    end function canonical_flux_from_state

end module neort_gc_models
