module test_gc_nonlocal_transport_fixture
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_CALLBACK_FAILURE, GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES, &
        GC_NONLOCAL_INVALID_INPUT, GC_NONLOCAL_SAMPLE_VALID, &
        GC_NONLOCAL_SAMPLE_WALL, GC_NONLOCAL_SUCCESS, &
        gc_nonlocal_component_t, gc_nonlocal_orbit_sample_t, &
        gc_nonlocal_transport_provider_t, gc_nonlocal_transport_quadrature_t, &
        gc_nonlocal_transport_reference_t
    implicit none
    private

    type, public, extends(gc_nonlocal_transport_provider_t) :: &
            manufactured_provider_t
        integer :: reference_calls = 0
        integer :: quadrature_calls = 0
        logical :: zero_width_mode = .false.
        logical :: unresolved = .false.
        logical :: fail_outer = .false.
        logical :: invalid_reference = .false.
        real(dp) :: width = 1.0_dp
        real(dp) :: tau0 = 3.5_dp
        real(dp) :: outer_factor = -0.25_dp
        real(dp) :: local_psi = 0.75_dp
    contains
        procedure :: get_section_reference => manufactured_get_reference
        procedure :: get_quadrature => manufactured_get_quadrature
        procedure :: get_components => manufactured_get_components
        procedure :: evaluate_orbit => manufactured_evaluate_orbit
        procedure :: evaluate_profiles => manufactured_evaluate_profiles
        procedure :: evaluate_outer_measure_factor => &
            manufactured_evaluate_outer_factor
    end type manufactured_provider_t

contains

    subroutine manufactured_get_reference(provider, reference, status)
        class(manufactured_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(out) :: reference
        integer, intent(out) :: status

        reference = gc_nonlocal_transport_reference_t()
        provider%reference_calls = provider%reference_calls + 1
        reference%section_id = 17
        reference%section_coordinate = 'Rc'
        reference%section_units = 'manufactured'
        reference%section_position = [1.7_dp, -0.2_dp, 0.0_dp]
        reference%section_flux = 0.25_dp
        reference%p_zeta_orientation = 1
        reference%frequency_semantics = GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES
        reference%hamiltonian_includes_phi = .true.
        reference%fixed = .true.
        if (provider%invalid_reference) then
            reference%frequency_semantics = 0
        end if
        status = GC_NONLOCAL_SUCCESS
    end subroutine manufactured_get_reference

    subroutine manufactured_get_quadrature(provider, h0_order, jk_order, &
            quadrature, status)
        class(manufactured_provider_t), intent(inout) :: provider
        integer, intent(in) :: h0_order, jk_order
        type(gc_nonlocal_transport_quadrature_t), intent(out) :: quadrature
        integer, intent(out) :: status

        integer :: i, j, node, n_nodes

        provider%quadrature_calls = provider%quadrature_calls + 1
        if (h0_order < 2 .or. jk_order < 2) then
            status = GC_NONLOCAL_INVALID_INPUT
            return
        end if
        quadrature = gc_nonlocal_transport_quadrature_t()
        n_nodes = h0_order*jk_order
        allocate(quadrature%h0(n_nodes), quadrature%j_k(n_nodes), &
            quadrature%weight(n_nodes), quadrature%j_k_upper_bound(n_nodes))
        node = 0
        do i = 1, h0_order
            do j = 1, jk_order
                node = node + 1
                quadrature%h0(node) = real(i, dp)/real(h0_order+1, dp)
                quadrature%j_k(node) = 2.0_dp*real(j, dp)/ &
                    real(jk_order+1, dp)
                quadrature%weight(node) = 0.5_dp
                quadrature%j_k_upper_bound(node) = 2.0_dp
            end do
        end do
        quadrature%h0_order = h0_order
        quadrature%jk_order = jk_order
        quadrature%n_nodes = n_nodes
        quadrature%paired_domain = .true.
        quadrature%domain_certified = .true.
        quadrature%converged = .true.
        quadrature%h0_min = 0.0_dp
        quadrature%h0_scale = 1.0_dp
        status = GC_NONLOCAL_SUCCESS
    end subroutine manufactured_get_quadrature

    subroutine manufactured_get_components(provider, reference, h0, jperp, &
            components, status)
        class(manufactured_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp
        type(gc_nonlocal_component_t), allocatable, intent(out) :: components(:)
        integer, intent(out) :: status

        if (reference%section_id /= 17) then
            status = GC_NONLOCAL_INVALID_INPUT
            return
        end if
        if (.not. all(ieee_is_finite([h0, jperp]))) then
            status = GC_NONLOCAL_INVALID_INPUT
            return
        end if
        allocate(components(2))
        if (provider%zero_width_mode) then
            if (provider%width <= 0.0_dp) then
                status = GC_NONLOCAL_INVALID_INPUT
                return
            end if
            components(1) = gc_nonlocal_component_t(101, 1, -provider%width, &
                provider%width)
            components(2) = gc_nonlocal_component_t(202, -1, -provider%width, &
                provider%width)
        else
            components(1) = gc_nonlocal_component_t(101, 1, -1.0_dp, 1.0_dp)
            components(2) = gc_nonlocal_component_t(202, -1, -1.0_dp, 1.0_dp)
        end if
        status = GC_NONLOCAL_SUCCESS
    end subroutine manufactured_get_components

    subroutine manufactured_evaluate_orbit(provider, reference, h0, jperp, x, &
            sigma, component_id, sample, status)
        class(manufactured_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
        integer, intent(out) :: status

        real(dp) :: tau_b, amplitude

        sample = gc_nonlocal_orbit_sample_t()
        status = GC_NONLOCAL_INVALID_INPUT
        if (reference%section_id /= 17) return
        if (.not. all(ieee_is_finite([h0, jperp, x]))) return
        if (abs(sigma) /= 1) return
        if (component_id /= 101 .and. component_id /= 202) return

        if (provider%unresolved) then
            if (component_id == 202) then
                if (x > 0.5_dp) then
                    sample%status = GC_NONLOCAL_SAMPLE_WALL
                    sample%component_id = component_id
                    sample%sigma = sigma
                    status = GC_NONLOCAL_SUCCESS
                    return
                end if
            end if
        end if

        tau_b = provider%tau0
        amplitude = 1.0_dp
        if (provider%zero_width_mode) then
            tau_b = provider%tau0*(1.0_dp + provider%width)
            amplitude = 1.0_dp + provider%width**2
        end if
        sample%status = GC_NONLOCAL_SAMPLE_VALID
        sample%component_id = component_id
        sample%sigma = sigma
        sample%psi_star = provider%local_psi + 2.0_dp*x
        sample%dpsi_star_dx = 2.0_dp
        sample%tau_b = tau_b
        sample%omega_b = 2.0_dp*acos(-1.0_dp)/tau_b
        sample%omega_phi = (x - 2.0_dp*acos(-1.0_dp)*2.0_dp/(-3.0_dp)) &
            /tau_b
        sample%domega_b_dx = 0.0_dp
        sample%domega_phi_dx = 1.0_dp/tau_b
        sample%h_m = cmplx(amplitude, 0.0_dp, kind=dp)
        sample%derivatives_available = .true.
        sample%nforce = 0
        status = GC_NONLOCAL_SUCCESS
    end subroutine manufactured_evaluate_orbit

    subroutine manufactured_evaluate_profiles(provider, reference, h0, jperp, &
            x, sigma, component_id, sample, force_count, force_values, status)
        class(manufactured_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id, force_count
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        real(dp), intent(out) :: force_values(force_count)
        integer, intent(out) :: status

        force_values = 0.0_dp
        status = GC_NONLOCAL_INVALID_INPUT
        if (reference%section_id /= 17) return
        if (force_count /= 2) return
        if (abs(sigma) /= 1) return
        if (component_id <= 0) return
        if (.not. all(ieee_is_finite([h0, jperp, x]))) return
        if (sample%status /= GC_NONLOCAL_SAMPLE_VALID) return
        force_values(1) = h0**2 + 2.0_dp*jperp + &
            sample%psi_star - provider%local_psi
        force_values(2) = 3.0_dp*h0 - jperp**2
        status = GC_NONLOCAL_SUCCESS
    end subroutine manufactured_evaluate_profiles

    subroutine manufactured_evaluate_outer_factor(provider, reference, h0, &
            jperp, x, sigma, component_id, sample, outer_factor, status)
        class(manufactured_provider_t), intent(inout) :: provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        real(dp), intent(out) :: outer_factor
        integer, intent(out) :: status

        outer_factor = 0.0_dp
        status = GC_NONLOCAL_INVALID_INPUT
        if (reference%section_id /= 17) return
        if (.not. all(ieee_is_finite([h0, jperp, x]))) return
        if (abs(sigma) /= 1) return
        if (component_id <= 0) return
        if (sample%status /= GC_NONLOCAL_SAMPLE_VALID) return
        if (provider%fail_outer) then
            status = GC_NONLOCAL_CALLBACK_FAILURE
            return
        end if
        outer_factor = provider%outer_factor
        status = GC_NONLOCAL_SUCCESS
    end subroutine manufactured_evaluate_outer_factor

end module test_gc_nonlocal_transport_fixture

program test_gc_nonlocal_transport_integral
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_nonlocal_transport_integral, only: &
        integrate_gc_nonlocal_transport
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_CALLBACK_FAILURE, GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES, &
        GC_NONLOCAL_INVALID_INPUT, GC_NONLOCAL_PARTIAL, GC_NONLOCAL_SUCCESS, &
        gc_nonlocal_resonance_options_t, gc_nonlocal_transport_options_t, &
        gc_nonlocal_transport_result_t
    use test_gc_nonlocal_transport_fixture, only: manufactured_provider_t
    implicit none

    type(manufactured_provider_t) :: provider
    type(gc_nonlocal_transport_options_t) :: options
    type(gc_nonlocal_transport_result_t) :: result
    real(dp) :: local_reference(2), phase_weight, frequency_weight
    integer :: status

    options = gc_nonlocal_transport_options_t()
    options%max_h0_nodes = 4
    options%max_jperp_nodes = 4
    options%max_total_nodes = 16
    options%h0_order = 2
    options%jk_order = 2
    options%require_converged = .false.
    options%resonance_options = gc_nonlocal_resonance_options_t()
    options%resonance_options%scan_intervals = 16
    options%resonance_options%max_root_iterations = 128
    options%resonance_options%max_roots = 32
    options%resonance_options%force_count = 2
    options%resonance_options%residual_tolerance = 1.0e-12_dp
    options%resonance_options%x_tolerance = 1.0e-12_dp
    options%resonance_options%derivative_tolerance = 1.0e-12_dp

    call integrate_gc_nonlocal_transport(provider, 2, -3, options, result, status)
    call require(status == GC_NONLOCAL_SUCCESS, &
        'manufactured tensor-product transport did not certify')
    call require(result%certified, 'certified transport result has false flag')
    call require(result%weighted_nodes == 4 .and. result%certified_nodes == 4, &
        'outer quadrature node accounting is wrong')
    call require(result%ncomponents == 8 .and. result%nroots == 8, &
        'disconnected class/root accounting is wrong')
    call require(result%reference%section_id == 17, &
        'fixed section reference was not retained')
    call require(result%reference%frequency_semantics == &
        GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES, &
        'frequency convention was not retained')
    call require(result%reference%hamiltonian_includes_phi, &
        'Hamiltonian Phi inclusion was not retained')

    call supplied_local_reference(0.0_dp, 0.0_dp, 2, local_reference, status)
    call require(status == GC_NONLOCAL_SUCCESS, 'local reference callback failed')
    call require_close(result%contribution(1), local_reference(1), &
        'frequency-form polynomial H,J force 1', 1.0e-8_dp)
    call require_close(result%contribution(2), local_reference(2), &
        'frequency-form polynomial H,J force 2', 1.0e-8_dp)

    call compare_frequency_and_phase_oracles(phase_weight, frequency_weight)
    call require_close(frequency_weight, phase_weight, &
        'phase/frequency root representations', 1.0e-12_dp)

    provider%zero_width_mode = .true.
    provider%width = 1.0e-2_dp
    call integrate_gc_nonlocal_transport(provider, 2, -3, options, result, status)
    call require(status == GC_NONLOCAL_SUCCESS, 'finite-width limit run failed')
    phase_weight = maxval(abs(result%contribution(1:2) - local_reference))

    provider%width = 1.0e-6_dp
    call integrate_gc_nonlocal_transport(provider, 2, -3, options, result, status)
    call require(status == GC_NONLOCAL_SUCCESS, 'zero-width limit run failed')
    frequency_weight = maxval(abs(result%contribution(1:2) - local_reference))
    call require(frequency_weight < phase_weight, &
        'zero-width result did not approach supplied local reference')
    call require(frequency_weight < 1.0e-3_dp, &
        'zero-width result remains outside local-reference tolerance')

    provider%zero_width_mode = .false.
    provider%unresolved = .true.
    call integrate_gc_nonlocal_transport(provider, 2, -3, options, result, status)
    call require(status == GC_NONLOCAL_PARTIAL, &
        'unresolved weighted class was not propagated as partial')
    call require(.not. result%certified, 'partial result was certified')
    call require(all(abs(result%contribution(1:2)) < 1.0e-14_dp), &
        'partial diagnostic values leaked into accepted contribution')
    call require(any(abs(result%diagnostic_contribution(1:2)) > 0.0_dp), &
        'partial diagnostic contribution was not retained')

    provider%unresolved = .false.
    provider%fail_outer = .true.
    call integrate_gc_nonlocal_transport(provider, 2, -3, options, result, status)
    call require(status == GC_NONLOCAL_CALLBACK_FAILURE, &
        'outer-factor callback status was not propagated')
    call require(.not. result%certified, 'callback failure was certified')
    call require(all(abs(result%contribution(1:2)) < 1.0e-14_dp), &
        'callback-failure diagnostic leaked into accepted contribution')

    provider%fail_outer = .false.
    provider%invalid_reference = .true.
    call integrate_gc_nonlocal_transport(provider, 2, -3, options, result, status)
    call require(status == GC_NONLOCAL_INVALID_INPUT, &
        'invalid fixed reference was accepted')

    provider%invalid_reference = .false.
    call integrate_gc_nonlocal_transport(provider, 2, 0, options, result, status)
    call require(status == GC_NONLOCAL_INVALID_INPUT, &
        'zero toroidal harmonic was accepted for Eq. 17')

    write (*, '(a)') 'test_gc_nonlocal_transport_integral OK'

contains

    subroutine supplied_local_reference(h0, jperp, force_count, values, status)
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: force_count
        real(dp), intent(out) :: values(force_count)
        integer, intent(out) :: status

        values = 0.0_dp
        status = GC_NONLOCAL_INVALID_INPUT
        if (.not. all(ieee_is_finite([h0, jperp]))) return
        if (force_count /= 2) return
        ! Two opposite-sigma branches, one |d psi_star/dx|=2 factor,
        ! signed W_outer=-1/4, |m3|=3, tau_b=3.5, and analytic H,J
        ! integrals 14/3 and 1/3.  This is only a zero-width seam oracle;
        ! it is not a proof of the full NEO local formula.
        values(1) = 2.0_dp*(-0.25_dp)*2.0_dp*3.0_dp*3.5_dp**2 &
            *(14.0_dp/3.0_dp)
        values(2) = 2.0_dp*(-0.25_dp)*2.0_dp*3.0_dp*3.5_dp**2 &
            *(1.0_dp/3.0_dp)
        status = GC_NONLOCAL_SUCCESS
    end subroutine supplied_local_reference

    subroutine compare_frequency_and_phase_oracles(phase_value, frequency_value)
        real(dp), intent(out) :: phase_value, frequency_value

        real(dp) :: tau_b, harmonic_n, dg_dx, dr_dx, base

        tau_b = 3.5_dp
        harmonic_n = -3.0_dp
        dg_dx = 1.0_dp
        dr_dx = harmonic_n/tau_b*dg_dx
        base = (-0.25_dp)*2.0_dp*1.0_dp
        frequency_value = harmonic_n**2*tau_b*base/abs(dr_dx)
        phase_value = abs(harmonic_n)*tau_b**2*base/abs(dg_dx)
    end subroutine compare_frequency_and_phase_oracles

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            write (*, '(a)') trim(message)
            error stop 1
        end if
    end subroutine require

    subroutine require_close(value, reference, label, tolerance)
        real(dp), intent(in) :: value, reference, tolerance
        character(*), intent(in) :: label

        if (abs(value - reference) > tolerance) then
            write (*, '(a,2(1x,es24.16))') trim(label), value, reference
            error stop 'manufactured oracle mismatch'
        end if
    end subroutine require_close

end program test_gc_nonlocal_transport_integral
