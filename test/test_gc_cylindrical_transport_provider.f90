module gc_cylindrical_transport_provider_test_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_class_adapter, only: &
        GC_CYL_CLASS_SUCCESS, gc_cylindrical_class_adapter_t, &
        gc_cylindrical_class_interval_t, gc_cylindrical_class_launch_t, &
        initialize_gc_cylindrical_class_adapter
    use neort_gc_cylindrical_model, only: &
        GC_CYL_SUCCESS, gc_cylindrical_field_sample_t, &
        gc_cylindrical_field_t, gc_cylindrical_potential_t, &
        make_gc_cylindrical_field_sample
    use neort_gc_cylindrical_nonlocal_provider, only: &
        GC_CYL_NONLOCAL_ORBIT_VALID, GC_CYL_NONLOCAL_ORBIT_WALL, &
        GC_CYL_NONLOCAL_SUCCESS, GC_CYL_NONLOCAL_WALL_CLEAR, &
        GC_CYL_NONLOCAL_WALL_HIT, gc_cylindrical_nonlocal_context_t, &
        gc_cylindrical_nonlocal_orbit_t, &
        initialize_gc_cylindrical_nonlocal_provider
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_SAMPLE_VALID, gc_nonlocal_component_t
    use neort_gc_nonlocal_transport_types, only: &
        gc_nonlocal_orbit_sample_t, gc_nonlocal_transport_reference_t
    implicit none
    private

    real(dp), parameter, public :: MASS = 1.0_dp
    real(dp), parameter, public :: CHARGE = 1.0_dp
    real(dp), parameter, public :: C_LIGHT = 2.0_dp
    real(dp), parameter, public :: H0_REFERENCE = 4.0_dp
    real(dp), parameter, public :: JPERP_REFERENCE = 0.5_dp
    real(dp), parameter, public :: PI = 3.14159265358979323846264338328_dp

    type, public :: transport_test_state_t
        real(dp) :: width = 1.0_dp
        logical :: omit_splitter = .false.
        logical :: wall_mode = .false.
        logical :: bad_section = .false.
        logical :: fail_outer = .false.
        integer :: factory_calls = 0
    end type transport_test_state_t

    type, extends(gc_cylindrical_field_t), public :: transport_test_field_t
    contains
        procedure :: evaluate => transport_field_evaluate
    end type transport_test_field_t

    type, extends(gc_cylindrical_potential_t), public :: transport_test_potential_t
    contains
        procedure :: evaluate => transport_potential_evaluate
    end type transport_test_potential_t

    public :: transport_node_factory
    public :: transport_outer_factor

    type(transport_test_field_t), target, save :: test_field
    type(transport_test_potential_t), target, save :: test_potential

contains

    subroutine transport_field_evaluate(self, position, sample, status)
        class(transport_test_field_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(out) :: sample
        integer, intent(out) :: status

        real(dp) :: b(3), db_dq(3, 3), grad_psi(3)

        associate (unused_self => self)
        end associate
        b = [0.0_dp, 1.0_dp + 0.1_dp*(position(1) - 2.0_dp), 0.0_dp]
        db_dq = 0.0_dp
        db_dq(2, 1) = 0.1_dp
        grad_psi = [1.0_dp, 0.0_dp, 0.0_dp]
        call make_gc_cylindrical_field_sample(position(1), b, db_dq, &
            position(1), grad_psi, sample, status)
    end subroutine transport_field_evaluate

    subroutine transport_potential_evaluate(self, position, field, potential, &
            gradient, status)
        class(transport_test_potential_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(out) :: potential, gradient(3)
        integer, intent(out) :: status

        associate (unused_self => self, unused_position => position, &
                unused_field => field)
        end associate
        potential = 0.0_dp
        gradient = 0.0_dp
        status = GC_CYL_SUCCESS
    end subroutine transport_potential_evaluate

    subroutine transport_cut_map(rc, user_data, position, dposition_drc, status)
        real(dp), intent(in) :: rc
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: position(3), dposition_drc(3)
        integer, intent(out) :: status

        associate (unused_user_data => user_data)
        end associate
        position = [2.0_dp + rc, 0.0_dp, 0.0_dp]
        dposition_drc = [1.0_dp, 0.0_dp, 0.0_dp]
        status = GC_CYL_SUCCESS
    end subroutine transport_cut_map

    subroutine transport_splitter(h0, jperp, sigma, candidate, user_data, &
            split_classes, certified, status)
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: sigma
        type(gc_cylindrical_class_interval_t), intent(in) :: candidate
        class(*), pointer, intent(inout) :: user_data
        type(gc_cylindrical_class_interval_t), allocatable, intent(out) :: &
            split_classes(:)
        logical, intent(out) :: certified
        integer, intent(out) :: status

        type(transport_test_state_t), pointer :: state

        certified = .false.
        status = GC_CYL_CLASS_SUCCESS
        if (.not. associated(user_data)) then
            status = 1
            return
        end if
        select type (user_data)
            type is (transport_test_state_t)
            state => user_data
        class default
            status = 1
            return
        end select
        associate (unused_h0 => h0, unused_jperp => jperp, unused_sigma => sigma)
        end associate
        if (state%omit_splitter) return
        allocate(split_classes(1))
        split_classes(1) = candidate
        split_classes(1)%orbit_return_certified = .true.
        split_classes(1)%topology_certified = .true.
        certified = .true.
    end subroutine transport_splitter

    subroutine transport_node_factory(h0, jperp, user_data, adapter, context, &
            status)
        real(dp), intent(in) :: h0, jperp
        class(*), pointer, intent(inout) :: user_data
        type(gc_cylindrical_class_adapter_t), intent(out) :: adapter
        type(gc_cylindrical_nonlocal_context_t), intent(out) :: context
        integer, intent(out) :: status

        type(transport_test_state_t), pointer :: state
        real(dp) :: section_reference(3)

        status = 1
        if (.not. associated(user_data)) return
        select type (user_data)
            type is (transport_test_state_t)
            state => user_data
        class default
            return
        end select
        state%factory_calls = state%factory_calls + 1
        call initialize_gc_cylindrical_class_adapter(test_field, test_potential, &
            h0, jperp, MASS, CHARGE, C_LIGHT, -0.5_dp*state%width, &
            0.5_dp*state%width, transport_cut_map, adapter, status, &
            splitter=transport_splitter, user_data=state)
        if (state%omit_splitter) then
            call initialize_gc_cylindrical_class_adapter(test_field, &
                test_potential, h0, jperp, MASS, CHARGE, C_LIGHT, &
                -0.5_dp*state%width, 0.5_dp*state%width, transport_cut_map, &
                adapter, status, user_data=state)
        end if
        if (status /= GC_CYL_CLASS_SUCCESS) return
        section_reference = [2.0_dp, 0.0_dp, 0.0_dp]
        if (state%bad_section) section_reference(1) = 2.25_dp
        call initialize_gc_cylindrical_nonlocal_provider(h0, jperp, context, &
            status, particle_charge=CHARGE, c_light=C_LIGHT, &
            component_provider=transport_components, &
            orbit_provider=transport_orbit, harmonic_provider=transport_harmonic, &
            force_provider=transport_force, section_coordinate='R_c', &
            section_reference=section_reference, section_reference_id='cut-R', &
            user_data=state)
    end subroutine transport_node_factory

    subroutine transport_components(h0, jperp, user_data, components, status)
        real(dp), intent(in) :: h0, jperp
        class(*), pointer, intent(inout) :: user_data
        type(gc_nonlocal_component_t), allocatable, intent(out) :: components(:)
        integer, intent(out) :: status

        type(transport_test_state_t), pointer :: state

        status = 1
        if (.not. associated(user_data)) return
        select type (user_data)
            type is (transport_test_state_t)
            state => user_data
        class default
            return
        end select
        if (h0 <= 0.0_dp .or. jperp < 0.0_dp) return
        allocate(components(2))
        components(1) = gc_nonlocal_component_t(component_id=1, sigma=1, &
            x_min=-0.5_dp*state%width, x_max=0.5_dp*state%width)
        components(2) = gc_nonlocal_component_t(component_id=1, sigma=-1, &
            x_min=-0.5_dp*state%width, x_max=0.5_dp*state%width)
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine transport_components

    subroutine transport_orbit(h0, jperp, x, sigma, component_id, user_data, &
            orbit, status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        class(*), pointer, intent(inout) :: user_data
        type(gc_cylindrical_nonlocal_orbit_t), intent(out) :: orbit
        integer, intent(out) :: status

        type(transport_test_state_t), pointer :: state
        real(dp) :: vparallel, dvparallel_dx, psi_star, dpsi_star_dx
        real(dp) :: p_phi, dp_phi_dx, vparallel_squared
        real(dp) :: tau_b, delta_phi_b

        orbit = gc_cylindrical_nonlocal_orbit_t()
        status = 1
        if (.not. associated(user_data)) return
        select type (user_data)
            type is (transport_test_state_t)
            state => user_data
        class default
            return
        end select
        if (component_id /= 1 .or. abs(sigma) /= 1) return
        vparallel_squared = 2.0_dp*h0 - jperp*(1.0_dp + 0.1_dp*x)
        if (vparallel_squared <= 0.0_dp) return
        vparallel = real(sigma, dp)*sqrt(vparallel_squared)
        dvparallel_dx = real(sigma, dp)*(-0.1_dp*jperp) &
            /(2.0_dp*sqrt(vparallel_squared))
        psi_star = (2.0_dp + x)*(1.0_dp + 2.0_dp*vparallel)
        dpsi_star_dx = 1.0_dp + 2.0_dp*vparallel &
            + 2.0_dp*(2.0_dp + x)*dvparallel_dx
        p_phi = 0.5_dp*psi_star
        dp_phi_dx = 0.5_dp*dpsi_star_dx
        tau_b = 2.0_dp + x
        delta_phi_b = PI + x
        orbit%component_id = component_id
        orbit%sigma = sigma
        orbit%winding = 1
        orbit%section_return_crossings = 2
        orbit%intersection_orientations = [-1, 1]
        orbit%intersection_times = [0.25_dp*tau_b, 0.75_dp*tau_b]
        orbit%intersection_rates = [1.0_dp, -1.0_dp]
        orbit%intersection_multiplicity_certified = .true.
        orbit%winding_available = .true.
        orbit%section_return_available = .true.
        orbit%complete_cycle_return = .true.
        orbit%p_phi_mapping_certified = .true.
        orbit%mapping_orientation = 1
        orbit%electric_potential_included = .true.
        orbit%delta_phi_b = delta_phi_b
        orbit%ddelta_phi_b_dx = 1.0_dp
        orbit%phase_advance_available = .true.
        orbit%period_derivative_available = .true.
        orbit%wall_checked = .true.
        orbit%wall_status = GC_CYL_NONLOCAL_WALL_CLEAR
        orbit%section%coordinate = 'R_c'
        orbit%section%reference = [2.0_dp, 0.0_dp, 0.0_dp]
        orbit%section%reference_id = 'cut-R'
        orbit%section%locked = .true.
        orbit%p_phi = p_phi
        orbit%dp_phi_dx = dp_phi_dx
        orbit%tau_b = tau_b
        orbit%dtau_b_dx = 1.0_dp
        orbit%omega_b = 2.0_dp*PI/tau_b
        orbit%omega_phi = delta_phi_b/tau_b
        orbit%domega_b_dx = -2.0_dp*PI/tau_b**2
        orbit%domega_phi_dx = (tau_b - delta_phi_b)/tau_b**2
        orbit%derivatives_available = .true.
        orbit%status = GC_CYL_NONLOCAL_ORBIT_VALID
        if (state%wall_mode) then
            orbit%status = GC_CYL_NONLOCAL_ORBIT_WALL
            orbit%wall_status = GC_CYL_NONLOCAL_WALL_HIT
        end if
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine transport_orbit

    subroutine transport_harmonic(h0, jperp, x, sigma, component_id, harmonic_m, &
            harmonic_n, orbit, user_data, h_m, status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id, harmonic_m, harmonic_n
        type(gc_cylindrical_nonlocal_orbit_t), intent(in) :: orbit
        class(*), pointer, intent(inout) :: user_data
        complex(dp), intent(out) :: h_m
        integer, intent(out) :: status

        h_m = cmplx(0.0_dp, 0.0_dp, kind=dp)
        status = 1
        if (.not. associated(user_data)) return
        if (component_id /= orbit%component_id .or. sigma /= orbit%sigma) return
        h_m = cmplx(1.0_dp + 0.02_dp*h0 + 0.1_dp*x, &
            0.25_dp*real(sigma, dp), kind=dp)
        associate (unused_jperp => jperp, unused_harmonic_m => harmonic_m, &
                unused_harmonic_n => harmonic_n)
        end associate
        associate (unused_user_data => user_data)
        end associate
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine transport_harmonic

    subroutine transport_force(h0, jperp, x, sigma, component_id, orbit, &
            user_data, force, status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_cylindrical_nonlocal_orbit_t), intent(in) :: orbit
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: force(:)
        integer, intent(out) :: status

        force = 0.0_dp
        status = 1
        if (.not. associated(user_data)) return
        if (component_id /= orbit%component_id .or. sigma /= orbit%sigma) return
        if (size(force) /= 2) return
        force(1) = h0 + 2.0_dp*x
        force(2) = real(sigma, dp)*(jperp + x)
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine transport_force

    subroutine transport_outer_factor(reference, h0, jperp, x, sigma, &
            component_id, launch, sample, user_data, outer_factor, status)
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_cylindrical_class_launch_t), intent(in) :: launch
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: outer_factor
        integer, intent(out) :: status

        type(transport_test_state_t), pointer :: state

        outer_factor = 0.0_dp
        status = 1
        if (.not. associated(user_data)) return
        select type (user_data)
            type is (transport_test_state_t)
            state => user_data
        class default
            return
        end select
        if (state%fail_outer) return
        associate (unused_x => x)
        end associate
        if (reference%section_id /= 17) return
        if (sample%status /= GC_NONLOCAL_SAMPLE_VALID) return
        if (launch%status /= GC_CYL_CLASS_SUCCESS) return
        if (component_id /= 1 .or. abs(sigma) /= 1) return
        outer_factor = -0.5_dp + 0.1_dp*h0 - 0.2_dp*jperp
        status = GC_CYL_NONLOCAL_SUCCESS
    end subroutine transport_outer_factor

end module gc_cylindrical_transport_provider_test_support

program test_gc_cylindrical_transport_provider
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_transport_provider, only: &
        GC_CYL_TRANSPORT_CLASS_UNCERTIFIED, &
        GC_CYL_TRANSPORT_NODE_MISMATCH, &
        GC_CYL_TRANSPORT_REFERENCE_MISMATCH, GC_CYL_TRANSPORT_SUCCESS, &
        GC_CYL_TRANSPORT_UNAVAILABLE, &
        gc_cylindrical_transport_provider_t, &
        initialize_gc_cylindrical_transport_provider
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_COMPONENT_IDENTITY, GC_NONLOCAL_PARTIAL, &
        GC_NONLOCAL_SAMPLE_VALID, GC_NONLOCAL_SAMPLE_WALL, &
        GC_NONLOCAL_SUCCESS, gc_nonlocal_component_t, &
        gc_nonlocal_orbit_sample_t, &
        gc_nonlocal_resonance_options_t
    use neort_gc_nonlocal_transport_integral, only: &
        integrate_gc_nonlocal_transport
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES, &
        gc_nonlocal_transport_options_t, gc_nonlocal_transport_quadrature_t, &
        gc_nonlocal_transport_reference_t, &
        gc_nonlocal_transport_result_t
    use gc_cylindrical_transport_provider_test_support, only: &
        H0_REFERENCE, JPERP_REFERENCE, PI, &
        transport_node_factory, transport_outer_factor, &
        transport_test_state_t
    implicit none

    type(gc_cylindrical_transport_provider_t) :: provider
    type(gc_nonlocal_transport_reference_t) :: reference
    type(gc_nonlocal_transport_options_t) :: options
    type(gc_nonlocal_transport_result_t) :: result
    type(gc_nonlocal_orbit_sample_t) :: sample
    type(transport_test_state_t), target :: state
    type(gc_nonlocal_transport_quadrature_t) :: quadrature
    type(gc_nonlocal_component_t), allocatable :: components(:)
    real(dp) :: h0_nodes(2), h0_weights(2), jperp_nodes(2), jperp_weights(2)
    real(dp) :: forces(2), outer_factor, expected, expected_second, expected_thin
    integer :: status, i, j, sigma

    reference = gc_nonlocal_transport_reference_t()
    reference%section_id = 17
    reference%section_coordinate = 'R_c'
    reference%section_units = 'cm'
    reference%section_position = [2.0_dp, 0.0_dp, 0.0_dp]
    reference%section_flux = 0.25_dp
    reference%p_zeta_orientation = 1
    reference%frequency_semantics = GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES
    reference%hamiltonian_includes_phi = .true.
    reference%fixed = .true.
    h0_nodes = [4.0_dp, 4.5_dp]
    h0_weights = [0.5_dp, 0.5_dp]
    jperp_nodes = [0.4_dp, 0.6_dp]
    jperp_weights = [1.0_dp, 1.0_dp]
    quadrature = gc_nonlocal_transport_quadrature_t()
    allocate(quadrature%h0(4), quadrature%j_k(4), quadrature%weight(4), &
        quadrature%j_k_upper_bound(4))
    quadrature%h0 = [4.0_dp, 4.0_dp, 4.5_dp, 4.5_dp]
    quadrature%j_k = [0.4_dp, 0.6_dp, 0.4_dp, 0.6_dp]
    quadrature%weight = 0.5_dp
    quadrature%j_k_upper_bound = 1.0_dp
    quadrature%h0_order = 2
    quadrature%jk_order = 2
    quadrature%n_nodes = 4
    quadrature%paired_domain = .true.
    quadrature%domain_certified = .true.
    quadrature%converged = .true.
    quadrature%h0_min = 0.0_dp
    quadrature%h0_scale = 1.0_dp
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

    call initialize_gc_cylindrical_transport_provider(provider, quadrature, &
        reference, 1, -2, 2, transport_node_factory, transport_outer_factor, &
        status=status, user_data=state, &
        section_reference_id='cut-R')
    call require(status == GC_CYL_TRANSPORT_SUCCESS, &
        'concrete provider initialization failed')
    call check_first_call_at_zero(reference, quadrature)
    call check_mutated_sample(reference, quadrature)

    call provider%get_components(reference, H0_REFERENCE, JPERP_REFERENCE, &
        components, status)
    call require(status == GC_CYL_TRANSPORT_SUCCESS, &
        'fresh certified component list failed')
    call require(size(components) == 2, &
        'opposite-sigma disconnected classes were not retained')
    call require(components(1)%component_id == 1 .and. &
        components(2)%component_id == 1 .and. &
        components(1)%sigma == 1 .and. components(2)%sigma == -1, &
        'component identity was not preserved')

    call provider%evaluate_orbit(reference, H0_REFERENCE, JPERP_REFERENCE, &
        0.0_dp, 1, 1, sample, status)
    call require(status == GC_CYL_TRANSPORT_SUCCESS, &
        'valid cylindrical orbit was rejected')
    call require(sample%status == GC_NONLOCAL_SAMPLE_VALID, &
        'valid orbit did not produce a valid sample')
    call require(sample%nforce == 0, &
        'orbit callback leaked force slots into the transport contract')
    call require_close(sample%omega_phi, PI/2.0_dp, &
        'full-cycle electric frequency', 1.0e-12_dp)
    call provider%evaluate_profiles(reference, H0_REFERENCE, JPERP_REFERENCE, &
        0.0_dp, 1, 1, sample, 2, forces, status)
    call require(status == GC_CYL_TRANSPORT_SUCCESS, &
        'signed native profile callback failed')
    call require_close(forces(1), H0_REFERENCE, 'force one', 1.0e-12_dp)
    call require_close(forces(2), JPERP_REFERENCE, 'force two sign', 1.0e-12_dp)
    call provider%evaluate_outer_measure_factor(reference, H0_REFERENCE, &
        JPERP_REFERENCE, 0.0_dp, 1, 1, sample, outer_factor, status)
    call require(status == GC_CYL_TRANSPORT_SUCCESS, &
        'signed outer factor callback failed')
    call require_close(outer_factor, -0.5_dp + 0.1_dp*H0_REFERENCE &
        - 0.2_dp*JPERP_REFERENCE, 'bare outer factor', 1.0e-12_dp)

    call provider%evaluate_orbit(reference, H0_REFERENCE, JPERP_REFERENCE, &
        0.0_dp, -1, 1, sample, status)
    call require(status == GC_CYL_TRANSPORT_SUCCESS, &
        'negative-sigma orbit was rejected')
    call provider%evaluate_profiles(reference, H0_REFERENCE, JPERP_REFERENCE, &
        0.0_dp, -1, 1, sample, 2, forces, status)
    call require_close(forces(2), -JPERP_REFERENCE, &
        'native negative-sigma force sign', 1.0e-12_dp)
    call provider%evaluate_orbit(reference, H0_REFERENCE, JPERP_REFERENCE, &
        0.0_dp, 1, 99, sample, status)
    call require(status == GC_NONLOCAL_COMPONENT_IDENTITY, &
        'stale component identity was accepted')

    call integrate_gc_nonlocal_transport(provider, 1, -2, options, result, &
        status)
    call require(status == GC_NONLOCAL_SUCCESS, &
        'full cylindrical transport integration did not certify')
    call require(result%certified, 'certified result has false flag')
    call require(result%weighted_nodes == 4 .and. result%certified_nodes == 4, &
        'H0,Jperp node accounting is wrong')
    call require(state%factory_calls >= 4, &
        'a fresh fixed-invariant node was not built for every quadrature node')
    expected = 0.0_dp
    expected_second = 0.0_dp
    do j = 1, size(jperp_nodes)
        do i = 1, size(h0_nodes)
            do sigma = -1, 1, 2
                expected = expected + h0_weights(i)*jperp_weights(j) &
                    *analytic_root_contribution(h0_nodes(i), jperp_nodes(j), &
                    sigma, 1)
                expected_second = expected_second &
                    + h0_weights(i)*jperp_weights(j) &
                    *analytic_root_contribution(h0_nodes(i), jperp_nodes(j), &
                    sigma, 2)
            end do
        end do
    end do
    call require_close(result%contribution(1), expected, &
        'independent transport oracle force one', 2.0e-9_dp)
    call require_close(result%contribution(2), expected_second, &
        'opposite-sigma signed force oracle', 2.0e-9_dp)

    state%width = 1.0e-6_dp
    state%factory_calls = 0
    call integrate_gc_nonlocal_transport(provider, 1, -2, options, result, &
        status)
    call require(status == GC_NONLOCAL_SUCCESS, 'thin-orbit limit failed')
    expected_thin = 0.0_dp
    do j = 1, size(jperp_nodes)
        do i = 1, size(h0_nodes)
            do sigma = -1, 1, 2
                expected_thin = expected_thin + h0_weights(i) &
                    *jperp_weights(j)*analytic_root_contribution( &
                    h0_nodes(i), jperp_nodes(j), sigma, 1)
            end do
        end do
    end do
    call require_close(result%contribution(1), expected_thin, &
        'thin-orbit limiting oracle', 2.0e-8_dp)

    state%width = 1.0_dp
    state%wall_mode = .true.
    call provider%evaluate_orbit(reference, H0_REFERENCE, JPERP_REFERENCE, &
        0.0_dp, 1, 1, sample, status)
    call require(status == GC_CYL_TRANSPORT_SUCCESS .and. &
        sample%status == GC_NONLOCAL_SAMPLE_WALL, &
        'wall status was not preserved')
    call integrate_gc_nonlocal_transport(provider, 1, -2, options, result, &
        status)
    call require(status == GC_NONLOCAL_PARTIAL .and. .not. result%certified, &
        'wall-containing transport was certified')
    state%wall_mode = .false.

    state%bad_section = .true.
    call provider%get_components(reference, H0_REFERENCE, JPERP_REFERENCE, &
        components, status)
    call require(status == GC_CYL_TRANSPORT_REFERENCE_MISMATCH, &
        'section-reference mismatch was accepted')
    state%bad_section = .false.

    state%omit_splitter = .true.
    call provider%get_components(reference, H0_REFERENCE, JPERP_REFERENCE, &
        components, status)
    call require(status == GC_CYL_TRANSPORT_CLASS_UNCERTIFIED, &
        'missing homoclinic certification was accepted')
    state%omit_splitter = .false.

    state%fail_outer = .true.
    call integrate_gc_nonlocal_transport(provider, 1, -2, options, result, &
        status)
    call require(status /= GC_NONLOCAL_SUCCESS .and. .not. result%certified, &
        'outer callback failure was certified')
    state%fail_outer = .false.

    call check_missing_callbacks(reference, h0_nodes, h0_weights, &
        jperp_nodes, jperp_weights, state)
    write (*, '(a)') 'test_gc_cylindrical_transport_provider OK'

contains

    subroutine check_first_call_at_zero(local_reference, local_quadrature)
        type(gc_nonlocal_transport_reference_t), intent(in) :: local_reference
        type(gc_nonlocal_transport_quadrature_t), intent(in) :: local_quadrature
        type(gc_cylindrical_transport_provider_t) :: local_provider
        type(gc_nonlocal_orbit_sample_t) :: local_sample
        type(transport_test_state_t), target :: local_state
        integer :: local_status

        local_state = transport_test_state_t()
        call initialize_gc_cylindrical_transport_provider(local_provider, &
            local_quadrature, local_reference, 1, -2, 2, transport_node_factory, &
            transport_outer_factor, status=local_status, user_data=local_state, &
            section_reference_id='cut-R')
        call require(local_status == GC_CYL_TRANSPORT_SUCCESS, &
            'zero-node regression provider initialization failed')
        call local_provider%evaluate_orbit(local_reference, 0.0_dp, 0.0_dp, &
            0.0_dp, 1, 1, local_sample, local_status)
        call require(local_state%factory_calls == 1, &
            'first zero node reused an uninitialized cache node')
        call require(local_status /= GC_CYL_TRANSPORT_SUCCESS, &
            'invalid first zero node was accepted')
    end subroutine check_first_call_at_zero

    subroutine check_mutated_sample(local_reference, local_quadrature)
        type(gc_nonlocal_transport_reference_t), intent(in) :: local_reference
        type(gc_nonlocal_transport_quadrature_t), intent(in) :: local_quadrature
        type(gc_cylindrical_transport_provider_t) :: local_provider
        type(gc_nonlocal_orbit_sample_t) :: local_sample
        type(transport_test_state_t), target :: local_state
        integer :: local_status

        local_state = transport_test_state_t()
        call initialize_gc_cylindrical_transport_provider(local_provider, &
            local_quadrature, local_reference, 1, -2, 2, transport_node_factory, &
            transport_outer_factor, status=local_status, user_data=local_state, &
            section_reference_id='cut-R')
        call require(local_status == GC_CYL_TRANSPORT_SUCCESS, &
            'mutation regression provider initialization failed')

        call local_provider%evaluate_orbit(local_reference, H0_REFERENCE, &
            JPERP_REFERENCE, 0.0_dp, 1, 1, local_sample, local_status)
        call require(local_status == GC_CYL_TRANSPORT_SUCCESS, &
            'class-kind mutation setup failed')
        local_sample%class_kind = local_sample%class_kind + 1
        call require_outer_mismatch(local_provider, local_reference, &
            local_sample, 'class_kind mutation was accepted')

        call local_provider%evaluate_orbit(local_reference, H0_REFERENCE, &
            JPERP_REFERENCE, 0.0_dp, 1, 1, local_sample, local_status)
        call require(local_status == GC_CYL_TRANSPORT_SUCCESS, &
            'nforce mutation setup failed')
        local_sample%nforce = 1
        call require_outer_mismatch(local_provider, local_reference, &
            local_sample, 'nforce mutation was accepted')

        call local_provider%evaluate_orbit(local_reference, H0_REFERENCE, &
            JPERP_REFERENCE, 0.0_dp, 1, 1, local_sample, local_status)
        call require(local_status == GC_CYL_TRANSPORT_SUCCESS, &
            'force mutation setup failed')
        local_sample%thermodynamic_force(1) = 1.0_dp
        call require_outer_mismatch(local_provider, local_reference, &
            local_sample, 'force mutation was accepted')

        call local_provider%evaluate_orbit(local_reference, H0_REFERENCE, &
            JPERP_REFERENCE, 0.0_dp, 1, 1, local_sample, local_status)
        call require(local_status == GC_CYL_TRANSPORT_SUCCESS, &
            'derivative mutation setup failed')
        local_sample%derivatives_available = .false.
        call require_outer_mismatch(local_provider, local_reference, &
            local_sample, 'derivative-availability mutation was accepted')

        call local_provider%evaluate_orbit(local_reference, H0_REFERENCE, &
            JPERP_REFERENCE, 0.0_dp, 1, 1, local_sample, local_status)
        call require(local_status == GC_CYL_TRANSPORT_SUCCESS, &
            'class-behavior mutation setup failed')
        local_sample%class_behavior_certified = .true.
        call require_outer_mismatch(local_provider, local_reference, &
            local_sample, 'class-behavior mutation was accepted')
    end subroutine check_mutated_sample

    subroutine require_outer_mismatch(local_provider, local_reference, sample, &
            message)
        type(gc_cylindrical_transport_provider_t), intent(inout) :: local_provider
        type(gc_nonlocal_transport_reference_t), intent(in) :: local_reference
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        character(*), intent(in) :: message
        real(dp) :: local_outer_factor
        integer :: local_status

        call local_provider%evaluate_outer_measure_factor(local_reference, &
            H0_REFERENCE, JPERP_REFERENCE, 0.0_dp, 1, 1, sample, &
            local_outer_factor, local_status)
        call require(local_status == GC_CYL_TRANSPORT_NODE_MISMATCH, message)
    end subroutine require_outer_mismatch

    real(dp) function analytic_root_contribution(h0, jperp, sigma, slot)
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: sigma, slot
        real(dp) :: vparallel, dpsi, hreal, himag, hsq, outer, force

        vparallel = real(sigma, dp)*sqrt(2.0_dp*h0 - jperp)
        dpsi = 1.0_dp + 2.0_dp*vparallel &
            - 0.1_dp*real(sigma, dp)*jperp*2.0_dp &
            /sqrt(2.0_dp*h0 - jperp)
        hreal = 1.0_dp + 0.02_dp*h0
        himag = 0.25_dp*real(sigma, dp)
        hsq = hreal**2 + himag**2
        outer = -0.5_dp + 0.1_dp*h0 - 0.2_dp*jperp
        force = h0
        if (slot == 2) force = real(sigma, dp)*jperp
        analytic_root_contribution = 4.0_dp*outer*abs(dpsi)*hsq*2.0_dp*force
    end function analytic_root_contribution

    subroutine check_missing_callbacks(local_reference, local_h0_nodes, &
            local_h0_weights, local_jperp_nodes, local_jperp_weights, local_state)
        type(gc_nonlocal_transport_reference_t), intent(in) :: local_reference
        real(dp), intent(in) :: local_h0_nodes(:), local_h0_weights(:)
        real(dp), intent(in) :: local_jperp_nodes(:), local_jperp_weights(:)
        type(transport_test_state_t), target, intent(inout) :: local_state
        type(gc_cylindrical_transport_provider_t) :: missing_provider
        type(gc_nonlocal_transport_quadrature_t) :: local_quadrature
        integer :: local_status

        local_quadrature = gc_nonlocal_transport_quadrature_t()
        allocate(local_quadrature%h0(4), local_quadrature%j_k(4), &
            local_quadrature%weight(4), local_quadrature%j_k_upper_bound(4))
        local_quadrature%h0 = [local_h0_nodes(1), local_h0_nodes(1), &
            local_h0_nodes(size(local_h0_nodes)), local_h0_nodes(size(local_h0_nodes))]
        local_quadrature%j_k = [local_jperp_nodes(1), &
            local_jperp_nodes(size(local_jperp_nodes)), local_jperp_nodes(1), &
            local_jperp_nodes(size(local_jperp_nodes))]
        local_quadrature%weight = 0.25_dp
        local_quadrature%j_k_upper_bound = maxval(local_jperp_nodes) + 1.0_dp
        local_quadrature%h0_order = 2
        local_quadrature%jk_order = 2
        local_quadrature%n_nodes = 4
        local_quadrature%paired_domain = .true.
        local_quadrature%domain_certified = .true.
        local_quadrature%converged = .true.
        local_quadrature%h0_min = 0.0_dp
        local_quadrature%h0_scale = 1.0_dp
        call initialize_gc_cylindrical_transport_provider(missing_provider, &
            local_quadrature, local_reference, 1, -2, 2, &
            status=local_status, user_data=local_state, &
            section_reference_id='cut-R')
        call require(local_status == GC_CYL_TRANSPORT_UNAVAILABLE, &
            'missing node and outer callbacks were not rejected')
    end subroutine check_missing_callbacks

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) then
            write (*, '(a)') trim(message)
            error stop 1
        end if
    end subroutine require

    subroutine require_close(value, expected_value, label, tolerance)
        real(dp), intent(in) :: value, expected_value, tolerance
        character(*), intent(in) :: label

        if (abs(value - expected_value) > tolerance) then
            write (*, '(a,2(1x,es24.16))') trim(label), value, expected_value
            error stop 1
        end if
    end subroutine require_close

end program test_gc_cylindrical_transport_provider
