module neort_gc_cylindrical_dynamics
    !! Hamiltonian guiding-center equations in physical cylindrical space.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: GC_CYL_INVALID_INPUT, &
        GC_CYL_SINGULAR_BSTAR, GC_CYL_SUCCESS, gc_cylindrical_field_sample_t, &
        gc_cylindrical_section_rate, gc_cylindrical_section_t, &
        gc_cylindrical_state_t

    implicit none
    private

    public :: gc_cylindrical_bstar_quantities
    public :: gc_cylindrical_rhs
    public :: gc_cylindrical_section_flux_density

contains

    pure subroutine gc_cylindrical_bstar_quantities(field, state, charge, c_light, &
            b_star, b_parallel_star, cylindrical_measure, status)
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        type(gc_cylindrical_state_t), intent(in) :: state
        real(dp), intent(in) :: charge, c_light
        real(dp), intent(out) :: b_star(3), b_parallel_star, cylindrical_measure
        integer, intent(out) :: status

        real(dp) :: alpha

        b_star = 0.0_dp
        b_parallel_star = 0.0_dp
        cylindrical_measure = 0.0_dp
        status = GC_CYL_INVALID_INPUT
        if (charge == 0.0_dp .or. c_light <= 0.0_dp) return
        if (state%R <= 0.0_dp .or. field%bmod <= 0.0_dp) return
        alpha = c_light*state%p_parallel/charge
        b_star = field%b + alpha*field%curl_bhat
        b_parallel_star = dot_product(field%bhat, b_star)
        if (abs(b_parallel_star) <= 100.0_dp*epsilon(b_parallel_star)) then
            status = GC_CYL_SINGULAR_BSTAR
            return
        end if
        cylindrical_measure = state%R*b_parallel_star
        status = GC_CYL_SUCCESS
    end subroutine gc_cylindrical_bstar_quantities

    pure subroutine gc_cylindrical_rhs(field, potential_gradient, mass, charge, &
            c_light, state, derivative, status)
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(in) :: potential_gradient(3), mass, charge, c_light
        type(gc_cylindrical_state_t), intent(in) :: state
        real(dp), intent(out) :: derivative(5)
        integer, intent(out) :: status

        real(dp) :: b_star(3), force(3), velocity(3), b_parallel_star
        real(dp) :: cylindrical_measure
        integer :: star_status

        derivative = 0.0_dp
        status = GC_CYL_INVALID_INPUT
        if (mass <= 0.0_dp) return
        call gc_cylindrical_bstar_quantities(field, state, charge, c_light, b_star, &
            b_parallel_star, cylindrical_measure, star_status)
        if (star_status /= GC_CYL_SUCCESS) then
            status = star_status
            return
        end if
        force = state%mu*field%grad_b + charge*potential_gradient
        velocity = (state%p_parallel/mass*b_star &
            +c_light/charge*cylindrical_cross(field%bhat, force)) &
            /b_parallel_star
        derivative(1) = velocity(1)
        derivative(2) = velocity(3)
        derivative(3) = velocity(2)/state%R
        derivative(4) = -dot_product(b_star, force)/b_parallel_star
        derivative(5) = 0.0_dp
        status = GC_CYL_SUCCESS
    end subroutine gc_cylindrical_rhs

    pure subroutine gc_cylindrical_section_flux_density(field, state, &
            potential_gradient, mass, charge, c_light, section, flux_density, status)
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        type(gc_cylindrical_state_t), intent(in) :: state
        real(dp), intent(in) :: potential_gradient(3)
        real(dp), intent(in) :: mass, charge, c_light
        type(gc_cylindrical_section_t), intent(in) :: section
        real(dp), intent(out) :: flux_density
        integer, intent(out) :: status

        real(dp) :: derivative(5), b_star(3), b_parallel_star, cylindrical_measure
        integer :: rhs_status, star_status

        flux_density = 0.0_dp
        call gc_cylindrical_rhs(field, potential_gradient, mass, charge, &
            c_light, state, derivative, rhs_status)
        if (rhs_status /= GC_CYL_SUCCESS) then
            status = rhs_status
            return
        end if
        call gc_cylindrical_bstar_quantities(field, state, charge, c_light, b_star, &
            b_parallel_star, cylindrical_measure, star_status)
        if (star_status /= GC_CYL_SUCCESS) then
            status = star_status
            return
        end if
        flux_density = cylindrical_measure*abs(gc_cylindrical_section_rate(section, &
            derivative))
        status = GC_CYL_SUCCESS
    end subroutine gc_cylindrical_section_flux_density

    pure function cylindrical_cross(left, right) result(value)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: value(3)

        value(1) = left(2)*right(3) - left(3)*right(2)
        value(2) = left(3)*right(1) - left(1)*right(3)
        value(3) = left(1)*right(2) - left(2)*right(1)
    end function cylindrical_cross

end module neort_gc_cylindrical_dynamics
