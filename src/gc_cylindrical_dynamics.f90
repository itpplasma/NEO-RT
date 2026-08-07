module neort_gc_cylindrical_dynamics
    !! Hamiltonian guiding-center equations in physical cylindrical space.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_cylindrical_bstar_symbolic, only: &
        evaluate_neort_cylindrical_bstar
    use neort_cylindrical_crossing_symbolic, only: &
        evaluate_neort_cylindrical_crossing_density
    use neort_cylindrical_littlejohn_symbolic, only: &
        evaluate_neort_cylindrical_littlejohn
    use neort_gc_cylindrical_model, only: GC_CYL_INVALID_INPUT, &
        GC_CYL_MEASURE_UNAVAILABLE, GC_CYL_SINGULAR_BSTAR, GC_CYL_SUCCESS, &
        gc_cylindrical_field_sample_t, &
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

        b_star = 0.0_dp
        b_parallel_star = 0.0_dp
        cylindrical_measure = 0.0_dp
        status = GC_CYL_INVALID_INPUT
        if (.not. all(ieee_is_finite([charge, c_light, state%R, &
                state%p_parallel, field%bmod, field%b, field%bhat, &
                field%curl_bhat]))) return
        if (charge == 0.0_dp .or. c_light <= 0.0_dp) return
        if (state%R <= 0.0_dp .or. field%bmod <= 0.0_dp) return
        call evaluate_neort_cylindrical_bstar(charge, c_light, state%R, &
            state%p_parallel, field%b(1), field%b(2), field%b(3), &
            field%bhat(1), field%bhat(2), field%bhat(3), &
            field%curl_bhat(1), field%curl_bhat(2), field%curl_bhat(3), &
            b_star(1), b_star(2), b_star(3), b_parallel_star, &
            cylindrical_measure)
        if (.not. all(ieee_is_finite([b_star, b_parallel_star, &
                cylindrical_measure]))) return
        if (abs(b_parallel_star) <= 100.0_dp*epsilon(b_parallel_star)) then
            status = GC_CYL_SINGULAR_BSTAR
            return
        end if
        status = GC_CYL_SUCCESS
    end subroutine gc_cylindrical_bstar_quantities

    pure subroutine gc_cylindrical_rhs(field, potential, potential_gradient, &
            mass, charge, c_light, state, derivative, status, &
            hamiltonian_value)
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(in) :: potential, potential_gradient(3)
        real(dp), intent(in) :: mass, charge, c_light
        type(gc_cylindrical_state_t), intent(in) :: state
        real(dp), intent(out) :: derivative(5)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: hamiltonian_value

        real(dp) :: b_star(3), force(3), cross_force(3), velocity(3)
        real(dp) :: b_parallel_star, cylindrical_measure
        real(dp) :: d_hamiltonian_dt, hamiltonian, canonical_p_phi

        derivative = 0.0_dp
        if (present(hamiltonian_value)) hamiltonian_value = 0.0_dp
        status = GC_CYL_INVALID_INPUT
        if (.not. all(ieee_is_finite([potential, potential_gradient, mass, charge, &
                c_light, state%R, state%Z, state%phi, state%p_parallel, &
                state%mu, field%bmod, field%psi, field%b, field%bhat, &
                field%curl_bhat, field%grad_b]))) return
        if (mass <= 0.0_dp .or. charge == 0.0_dp .or. c_light <= 0.0_dp) return
        if (state%mu < 0.0_dp) return
        if (state%R <= 0.0_dp .or. field%bmod <= 0.0_dp) return
        call evaluate_neort_cylindrical_littlejohn(mass, charge, c_light, &
            state%mu, field%bmod, potential, state%R, state%p_parallel, &
            field%psi, field%b(1), field%b(2), field%b(3), field%bhat(1), &
            field%bhat(2), field%bhat(3), field%curl_bhat(1), &
            field%curl_bhat(2), field%curl_bhat(3), field%grad_b(1), &
            field%grad_b(2), field%grad_b(3), potential_gradient(1), &
            potential_gradient(2), potential_gradient(3), b_star(1), &
            b_star(2), b_star(3), b_parallel_star, force(1), force(2), &
            force(3), cross_force(1), cross_force(2), cross_force(3), &
            velocity(1), velocity(2), velocity(3), derivative(1), &
            derivative(2), derivative(3), derivative(4), derivative(5), &
            cylindrical_measure, d_hamiltonian_dt, hamiltonian, &
            canonical_p_phi)
        if (.not. all(ieee_is_finite([derivative, b_star, force, cross_force, &
                velocity, b_parallel_star, cylindrical_measure, &
                d_hamiltonian_dt, hamiltonian, canonical_p_phi]))) then
            derivative = 0.0_dp
            return
        end if
        if (abs(b_parallel_star) <= 100.0_dp*epsilon(b_parallel_star)) then
            derivative = 0.0_dp
            status = GC_CYL_SINGULAR_BSTAR
            return
        end if
        if (present(hamiltonian_value)) hamiltonian_value = hamiltonian
        status = GC_CYL_SUCCESS
    end subroutine gc_cylindrical_rhs

    pure subroutine gc_cylindrical_section_flux_density(field, state, &
            potential, potential_gradient, mass, charge, c_light, section, &
            flux_density, status)
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        type(gc_cylindrical_state_t), intent(in) :: state
        real(dp), intent(in) :: potential, potential_gradient(3)
        real(dp), intent(in) :: mass, charge, c_light
        type(gc_cylindrical_section_t), intent(in) :: section
        real(dp), intent(out) :: flux_density
        integer, intent(out) :: status

        real(dp) :: derivative(5), b_star(3), b_parallel_star
        real(dp) :: signed_measure, generated_signed_measure
        real(dp) :: positive_crossing_density, signed_density
        real(dp) :: absolute_measure, absolute_cdot, cdot
        integer :: rhs_status, star_status

        flux_density = 0.0_dp
        call gc_cylindrical_rhs(field, potential, potential_gradient, mass, &
            charge, c_light, state, derivative, rhs_status)
        if (rhs_status /= GC_CYL_SUCCESS) then
            status = rhs_status
            return
        end if
        call gc_cylindrical_bstar_quantities(field, state, charge, c_light, b_star, &
            b_parallel_star, signed_measure, star_status)
        if (star_status /= GC_CYL_SUCCESS) then
            status = star_status
            return
        end if
        ! Eq. 14 uses a positive phase-space crossing density.  Keep the
        ! signed R*B_parallel* measure available above for orientation
        ! diagnostics; never use its sign to alter force/torque quantities.
        cdot = gc_cylindrical_section_rate(section, derivative)
        call evaluate_neort_cylindrical_crossing_density(signed_measure, cdot, &
            generated_signed_measure, signed_density, absolute_measure, &
            absolute_cdot, positive_crossing_density)
        if (.not. all(ieee_is_finite([generated_signed_measure, &
                signed_density, absolute_measure, absolute_cdot, &
                positive_crossing_density]))) then
            status = GC_CYL_MEASURE_UNAVAILABLE
            return
        end if
        flux_density = positive_crossing_density
        status = GC_CYL_SUCCESS
    end subroutine gc_cylindrical_section_flux_density

end module neort_gc_cylindrical_dynamics
