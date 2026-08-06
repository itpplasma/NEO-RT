module manufactured_gc_field
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use manufactured_torus_kernel, only: evaluate_manufactured_torus
    use neort_gc_dynamics, only: gc_field_sample_t
    use neort_gc_models, only: GC_MODEL_SUCCESS, GC_MODEL_INVALID_STATE, gc_field_t

    implicit none
    private

    type, extends(gc_field_t), public :: manufactured_torus_field_t
        real(dp) :: major_radius = 0.0_dp
        real(dp) :: toroidal_flux = 0.0_dp
        real(dp) :: q_axis = 0.0_dp
        real(dp) :: shear_shape = 0.0_dp
    contains
        procedure :: evaluate => evaluate_field
        procedure :: radial_profiles
    end type manufactured_torus_field_t

contains

    subroutine evaluate_field(self, position, sample, status)
        class(manufactured_torus_field_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_field_sample_t), intent(out) :: sample
        integer, intent(out) :: status
        real(dp) :: q_unused, shear_unused

        sample = gc_field_sample_t()
        status = GC_MODEL_INVALID_STATE
        if (position(1) <= 0.0_dp .or. position(1) >= self%major_radius) return
        call evaluate_manufactured_torus(position(1), position(2), position(3), &
            self%major_radius, self%toroidal_flux, self%q_axis, &
            self%shear_shape, sample%bmod, sample%sqrtg, &
            sample%grad_log_b(1), sample%grad_log_b(2), sample%grad_log_b(3), &
            sample%hcov(1), sample%hcov(2), sample%hcov(3), &
            sample%hcon(1), sample%hcon(2), sample%hcon(3), &
            sample%curl_h(1), sample%curl_h(2), sample%curl_h(3), &
            sample%psi, sample%grad_psi(1), sample%grad_psi(2), &
            sample%grad_psi(3), q_unused, shear_unused)
        if (sample%bmod <= 0.0_dp .or. sample%sqrtg <= 0.0_dp) return
        status = GC_MODEL_SUCCESS
    end subroutine evaluate_field

    subroutine radial_profiles(self, radius, q_safety, magnetic_shear)
        class(manufactured_torus_field_t), intent(in) :: self
        real(dp), intent(in) :: radius
        real(dp), intent(out) :: q_safety, magnetic_shear
        real(dp) :: unused(18)

        call evaluate_manufactured_torus(radius, 0.0_dp, 0.0_dp, &
            self%major_radius, self%toroidal_flux, self%q_axis, &
            self%shear_shape, unused(1), unused(2), unused(3), unused(4), &
            unused(5), unused(6), unused(7), unused(8), unused(9), unused(10), &
            unused(11), unused(12), unused(13), unused(14), unused(15), &
            unused(16), unused(17), unused(18), q_safety, magnetic_shear)
    end subroutine radial_profiles

end module manufactured_gc_field
