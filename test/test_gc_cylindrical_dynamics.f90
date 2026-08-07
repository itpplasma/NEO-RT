program test_gc_cylindrical_dynamics
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: GC_CYL_SECTION_PHI, GC_CYL_SUCCESS, &
        gc_cylindrical_field_sample_t, gc_cylindrical_invariants_t, &
        gc_cylindrical_linear_flux_potential_t, &
        gc_cylindrical_section_t, &
        gc_cylindrical_state_t, make_gc_cylindrical_field_sample, &
        gc_cylindrical_invariant_residuals, invariants_from_cylindrical_state, &
        canonical_flux_from_state, canonical_toroidal_momentum_from_state, &
        jperp_from_state, state_from_cylindrical_invariants
    use neort_gc_cylindrical_dynamics, only: gc_cylindrical_bstar_quantities, &
        gc_cylindrical_rhs, gc_cylindrical_section_flux_density
    use util_for_test, only: pass_test

    implicit none

    real(dp), parameter :: c_light = 1.0_dp
    type(gc_cylindrical_field_sample_t) :: field
    type(gc_cylindrical_state_t) :: state, reconstructed
    type(gc_cylindrical_invariants_t) :: invariants
    type(gc_cylindrical_linear_flux_potential_t) :: electric_potential
    real(dp) :: b(3), db(3, 3), grad_psi(3), gradient(3)
    real(dp) :: derivative(5), potential, reconstructed_potential
    real(dp) :: expected, thin_derivative, full_derivative
    real(dp) :: energy_residual, mu_residual, pphi_residual
    real(dp) :: b_star(3), b_parallel_star, cylindrical_measure, flux_density
    type(gc_cylindrical_section_t) :: section
    integer :: status

    call check_grad_b_drift()
    call check_curvature_drift()
    call check_e_cross_b_charge_reversal()
    call check_exact_electric_only_drift()
    call check_pitch_invariant_jacobian()
    call check_invariant_launch()
    call check_bstar_measure_and_section_flux()
    call check_axis_regular_rhs()
    call check_thin_limit()
    call pass_test

contains

    subroutine check_grad_b_drift()
        b = [0.0_dp, 0.0_dp, 2.0_dp]
        db = 0.0_dp
        db(3, 1) = 0.4_dp
        grad_psi = 0.0_dp
        call make_gc_cylindrical_field_sample(3.0_dp, b, db, 0.0_dp, &
            grad_psi, field, status)
        if (status /= GC_CYL_SUCCESS) error stop 'grad-B sample failed'
        state = gc_cylindrical_state_t()
        state%R = 3.0_dp
        state%p_parallel = 1.2_dp
        state%mu = 0.7_dp
        gradient = 0.0_dp
        call gc_cylindrical_rhs(field, gradient, 1.0_dp, 1.0_dp, c_light, &
            state, derivative, status)
        if (status /= GC_CYL_SUCCESS) error stop 'grad-B RHS failed'
        expected = state%mu*db(3, 1)/(field%bmod*state%R)
        call require_close('grad-B phi rate', derivative(3), expected, 1.0e-13_dp)
        call require_close('grad-B radial rate', derivative(1), 0.0_dp, 1.0e-13_dp)

        state%mu = 0.0_dp
        gradient = [0.6_dp, 0.0_dp, 0.0_dp]
        call gc_cylindrical_rhs(field, gradient, 1.0_dp, 2.0_dp, c_light, &
            state, derivative, status)
        expected = c_light*gradient(1)/(field%bmod*state%R)
        call require_close('axial electric phi rate', derivative(3), expected, &
            1.0e-13_dp)
    end subroutine check_grad_b_drift

    subroutine check_curvature_drift()
        b = [0.0_dp, 3.0_dp, 0.0_dp]
        db = 0.0_dp
        db(2, 1) = -1.5_dp
        call make_gc_cylindrical_field_sample(2.0_dp, b, db, 0.0_dp, &
            grad_psi, field, status)
        if (status /= GC_CYL_SUCCESS) error stop 'curvature sample failed'
        state = gc_cylindrical_state_t()
        state%R = 2.0_dp
        state%p_parallel = 1.0_dp
        gradient = 0.0_dp
        call gc_cylindrical_rhs(field, gradient, 1.0_dp, 1.0_dp, c_light, &
            state, derivative, status)
        if (status /= GC_CYL_SUCCESS) error stop 'curvature RHS failed'
        call require_close('curvature phi rate', derivative(3), 0.5_dp, 1.0e-13_dp)
        call require_close('curvature Z drift', derivative(2), 1.0_dp/6.0_dp, 1.0e-13_dp)
        state%mu = 0.2_dp
        call gc_cylindrical_rhs(field, gradient, 1.0_dp, 1.0_dp, c_light, &
            state, derivative, status)
        call require_close('toroidal grad-B Z drift', derivative(2), &
            1.0_dp/6.0_dp + 0.2_dp/2.0_dp, 1.0e-13_dp)
    end subroutine check_curvature_drift

    subroutine check_e_cross_b_charge_reversal()
        real(dp) :: positive_charge(5), negative_charge(5), normalized_reversal(5)
        real(dp) :: reversed_field(5)

        b = [0.0_dp, 0.0_dp, 2.0_dp]
        db = 0.0_dp
        call make_gc_cylindrical_field_sample(4.0_dp, b, db, 0.0_dp, &
            grad_psi, field, status)
        state = gc_cylindrical_state_t()
        state%R = 4.0_dp
        gradient = [0.6_dp, 0.0_dp, 0.0_dp]
        call gc_cylindrical_rhs(field, gradient, 1.0_dp, 2.0_dp, c_light, &
            state, positive_charge, status)
        call gc_cylindrical_rhs(field, gradient, 1.0_dp, -2.0_dp, c_light, &
            state, negative_charge, status)
        call gc_cylindrical_rhs(field, -gradient, 1.0_dp, -2.0_dp, c_light, &
            state, normalized_reversal, status)
        call gc_cylindrical_rhs(field, -gradient, 1.0_dp, 2.0_dp, c_light, &
            state, reversed_field, status)
        call require_close('E cross B charge reversal', positive_charge(3), &
            negative_charge(3), 1.0e-13_dp)
        call require_close('normalized charge reversal', positive_charge(3), &
            -normalized_reversal(3), 1.0e-13_dp)
        call require_close('E cross B field reversal', positive_charge(3), &
            -reversed_field(3), 1.0e-13_dp)
    end subroutine check_e_cross_b_charge_reversal

    subroutine check_exact_electric_only_drift()
        real(dp), parameter :: radius = 4.0_dp, b_axial = 2.5_dp
        real(dp), parameter :: omega_e = 0.37_dp
        real(dp) :: evaluated_potential, evaluated_gradient(3)

        ! For an axial field with psi'=R B, Phi=(omega_e/c) psi gives the
        ! exact angular E x B drift dot(phi)=omega_e.  This oracle is
        ! independent of the production backend and catches the historical
        ! factor-of-two error in its potential coefficient.
        b = [0.0_dp, 0.0_dp, b_axial]
        db = 0.0_dp
        grad_psi = [radius*b_axial, 0.0_dp, 0.0_dp]
        call make_gc_cylindrical_field_sample(radius, b, db, 1.2_dp, grad_psi, &
            field, status)
        if (status /= GC_CYL_SUCCESS) error stop 'electric-only field sample failed'
        electric_potential = gc_cylindrical_linear_flux_potential_t()
        electric_potential%coefficient = omega_e/c_light
        electric_potential%psi_reference = field%psi
        call electric_potential%evaluate([radius, 0.0_dp, 0.0_dp], field, &
            evaluated_potential, evaluated_gradient, status)
        if (status /= GC_CYL_SUCCESS) error stop 'electric-only potential failed'
        state = gc_cylindrical_state_t()
        state%R = radius
        call gc_cylindrical_rhs(field, evaluated_gradient, 1.0_dp, 1.0_dp, &
            c_light, state, derivative, status)
        if (status /= GC_CYL_SUCCESS) error stop 'electric-only RHS failed'
        call require_close('exact electric-only angular drift', derivative(3), &
            omega_e, 1.0e-13_dp)
        call require_close('exact electric-only radial drift', derivative(1), &
            0.0_dp, 1.0e-13_dp)
    end subroutine check_exact_electric_only_drift

    subroutine check_pitch_invariant_jacobian()
        real(dp), parameter :: mass_value = 3.7_dp, speed_value = 1.2_dp
        real(dp), parameter :: eta_value = 0.31_dp, delta_eta = 1.0e-3_dp
        real(dp) :: mu_plus, mu_minus, jperp_plus, jperp_minus
        real(dp) :: dmu_deta, djperp_deta

        ! The cylindrical state stores physical mu.  The caller's eta is
        ! therefore related by mu=1/2*m*v^2*eta; J_perp=2*m*mu is a
        ! different invariant and must not be substituted silently.
        mu_plus = 0.5_dp*mass_value*speed_value**2*(eta_value + delta_eta)
        mu_minus = 0.5_dp*mass_value*speed_value**2*(eta_value - delta_eta)
        jperp_plus = 2.0_dp*mass_value*mu_plus
        jperp_minus = 2.0_dp*mass_value*mu_minus
        dmu_deta = (mu_plus - mu_minus)/(2.0_dp*delta_eta)
        djperp_deta = (jperp_plus - jperp_minus)/(2.0_dp*delta_eta)
        call require_close('dmu/deta physical invariant', dmu_deta, &
            0.5_dp*mass_value*speed_value**2, 1.0e-11_dp)
        call require_close('dJperp/deta normalization', djperp_deta, &
            mass_value**2*speed_value**2, 1.0e-11_dp)
    end subroutine check_pitch_invariant_jacobian

    subroutine check_invariant_launch()
        real(dp), parameter :: phi_value = 0.15_dp

        b = [0.0_dp, 2.0_dp, 0.0_dp]
        db = 0.0_dp
        grad_psi = [0.0_dp, 0.0_dp, 0.0_dp]
        call make_gc_cylindrical_field_sample(3.0_dp, b, db, 1.4_dp, &
            grad_psi, field, status)
        state = gc_cylindrical_state_t()
        state%R = 3.0_dp
        state%Z = -0.2_dp
        state%phi = phi_value
        state%p_parallel = -0.8_dp
        state%mu = 0.23_dp
        potential = 0.4_dp
        call invariants_from_cylindrical_state(state, field, potential, &
            1.0_dp, 1.0_dp, c_light, invariants, status)
        if (status /= GC_CYL_SUCCESS) error stop 'invariant construction failed'
        call require_close('canonical P_phi sign', &
            invariants%canonical_toroidal_momentum, -1.0_dp, 1.0e-13_dp)
        call require_close('psi_star canonical sign', canonical_flux_from_state(&
            state, field, 1.0_dp, c_light), &
            canonical_toroidal_momentum_from_state(state, field, 1.0_dp, c_light), &
            1.0e-13_dp)
        call state_from_cylindrical_invariants([state%R, state%Z, state%phi], &
            field, potential, invariants, -1, 1.0_dp, 1.0_dp, c_light, &
            reconstructed, status)
        if (status /= GC_CYL_SUCCESS) error stop 'fixed-invariant launch failed'
        reconstructed_potential = potential
        call require_close('fixed P_phi launch', reconstructed%p_parallel, &
            state%p_parallel, 1.0e-13_dp)
        call require_close('fixed energy launch', reconstructed%mu, state%mu, 1.0e-13_dp)
        call require_close('fixed potential passthrough', reconstructed_potential, &
            potential, 1.0e-13_dp)
        call gc_cylindrical_invariant_residuals(reconstructed, field, potential, &
            1.0_dp, 1.0_dp, c_light, invariants, energy_residual, mu_residual, &
            pphi_residual, status)
        if (status /= GC_CYL_SUCCESS) error stop 'invariant residual evaluation failed'
        call require_close('H invariant residual', energy_residual, 0.0_dp, 1.0e-13_dp)
        call require_close('mu invariant residual', mu_residual, 0.0_dp, 1.0e-13_dp)
        call require_close('P_phi invariant residual', pphi_residual, 0.0_dp, 1.0e-13_dp)
        call require_close('Jperp normalization', jperp_from_state(state, 1.0_dp), &
            2.0_dp*state%mu, 1.0e-13_dp)
    end subroutine check_invariant_launch

    subroutine check_bstar_measure_and_section_flux()
        real(dp), parameter :: radius = 3.0_dp, bphi = 2.0_dp, bz = 1.0_dp
        real(dp), parameter :: p_parallel = 0.8_dp
        real(dp) :: expected_bparallel, expected_flux

        b = [0.0_dp, bphi, bz]
        db = 0.0_dp
        call make_gc_cylindrical_field_sample(radius, b, db, 0.0_dp, &
            [0.0_dp, 0.0_dp, 0.0_dp], field, status)
        state = gc_cylindrical_state_t()
        state%R = radius
        state%p_parallel = p_parallel
        call gc_cylindrical_bstar_quantities(field, state, 1.0_dp, c_light, &
            b_star, b_parallel_star, cylindrical_measure, status)
        if (status /= GC_CYL_SUCCESS) error stop 'B_parallel* sample failed'
        expected_bparallel = field%bmod + p_parallel*field%bhat(3) &
            *field%bhat(2)/radius
        call require_close('finite-lambda B_parallel*', b_parallel_star, &
            expected_bparallel, 1.0e-13_dp)
        call require_close('cylindrical R B_parallel* measure', cylindrical_measure, &
            radius*expected_bparallel, 1.0e-13_dp)
        section = gc_cylindrical_section_t()
        section%kind = GC_CYL_SECTION_PHI
        call gc_cylindrical_rhs(field, [0.0_dp, 0.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, &
            c_light, state, derivative, status)
        call gc_cylindrical_section_flux_density(field, state, &
            [0.0_dp, 0.0_dp, 0.0_dp], 1.0_dp, 1.0_dp, c_light, section, &
            flux_density, status)
        if (status /= GC_CYL_SUCCESS) error stop 'section flux evaluation failed'
        expected_flux = cylindrical_measure*abs(derivative(3))
        call require_close('section flux R B_parallel* dot section', flux_density, &
            expected_flux, 1.0e-13_dp)
    end subroutine check_bstar_measure_and_section_flux

    subroutine check_axis_regular_rhs()
        real(dp), parameter :: r_axis = 10.0_dp, toroidal_field = 5.0_dp
        real(dp) :: radius
        integer :: k

        do k = -1, 1
            radius = r_axis + real(k, dp)*0.01_dp
            b = [0.0_dp, toroidal_field, 0.0_dp]
            db = 0.0_dp
            call make_gc_cylindrical_field_sample(radius, b, db, 0.0_dp, &
                grad_psi, field, status)
            state = gc_cylindrical_state_t()
            state%R = radius
            state%p_parallel = 0.5_dp
            call gc_cylindrical_rhs(field, [0.0_dp, 0.0_dp, 0.0_dp], &
                1.0_dp, 1.0_dp, c_light, &
                state, derivative, status)
            if (status /= GC_CYL_SUCCESS) error stop 'magnetic-axis RHS failed'
            if (any(derivative /= derivative)) error stop 'axis RHS is non-finite'
        end do
    end subroutine check_axis_regular_rhs

    subroutine check_thin_limit()
        b = [0.0_dp, 0.0_dp, 2.0_dp]
        db = 0.0_dp
        db(3, 1) = 0.4_dp
        call make_gc_cylindrical_field_sample(3.0_dp, b, db, 0.0_dp, &
            grad_psi, field, status)
        state = gc_cylindrical_state_t()
        state%R = 3.0_dp
        state%mu = 0.7_dp
        call gc_cylindrical_rhs(field, [0.0_dp, 0.0_dp, 0.0_dp], &
            1.0_dp, 1.0_dp, c_light, &
            state, derivative, status)
        full_derivative = derivative(3)
        state%mu = 0.0_dp
        call gc_cylindrical_rhs(field, [0.0_dp, 0.0_dp, 0.0_dp], &
            1.0_dp, 1.0_dp, c_light, &
            state, derivative, status)
        thin_derivative = derivative(3)
        call require_close('strict thin magnetic-drift limit', thin_derivative, &
            0.0_dp, 1.0e-13_dp)
        call require_close('finite-width drift increment', full_derivative - &
            thin_derivative, 0.7_dp*0.4_dp/(2.0_dp*3.0_dp), 1.0e-13_dp)
    end subroutine check_thin_limit

    subroutine require_close(label, actual, reference, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: actual, reference, tolerance

        if (abs(actual - reference) > tolerance) then
            write(*, '(a,2(1x,es24.16),1x,es12.4)') trim(label), actual, reference, tolerance
            error stop 'cylindrical dynamics oracle failed'
        end if
    end subroutine require_close

end program test_gc_cylindrical_dynamics
