module gc_cylindrical_class_adapter_test_support
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: &
        GC_CYL_EQUILIBRIUM_DOMAIN, GC_CYL_SUCCESS, &
        gc_cylindrical_field_sample_t, gc_cylindrical_field_t, &
        gc_cylindrical_potential_t, make_gc_cylindrical_field_sample
    use neort_gc_cylindrical_class_adapter, only: &
        GC_CYL_CLASS_INVALID_INPUT, GC_CYL_CLASS_SUCCESS, &
        gc_cylindrical_class_interval_t, &
        gc_cylindrical_class_cut_map_i, gc_cylindrical_class_splitter_i
    use neort_gc_cylindrical_topology, only: &
        gc_cylindrical_allowed_region_set_t

    implicit none
    private

    real(dp), parameter, public :: H0_REFERENCE = 10.0_dp
    real(dp), parameter, public :: JPERP_REFERENCE = 0.75_dp
    real(dp), parameter, public :: MASS = 3.0_dp
    real(dp), parameter, public :: CHARGE = 2.0_dp
    real(dp), parameter, public :: C_LIGHT = 10.0_dp

    type, extends(gc_cylindrical_field_t), public :: manufactured_field_t
        logical :: invalid_hole = .false.
    contains
        procedure :: evaluate => evaluate_manufactured_field
    end type manufactured_field_t

    type, extends(gc_cylindrical_potential_t), public :: manufactured_potential_t
        real(dp) :: profile_h0 = H0_REFERENCE
        real(dp) :: profile_jperp = JPERP_REFERENCE
        real(dp) :: profile_mass = MASS
        real(dp) :: profile_charge = CHARGE
        real(dp) :: profile_c = C_LIGHT
    contains
        procedure :: evaluate => evaluate_manufactured_potential
    end type manufactured_potential_t

    public :: manufactured_cut_map
    public :: identity_splitter
    public :: manufactured_topology_provider
    public :: target_vparallel_squared
    public :: target_dvparallel_squared

contains

    subroutine evaluate_manufactured_field(self, position, sample, status)
        class(manufactured_field_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(out) :: sample
        integer, intent(out) :: status

        real(dp) :: b(3), db_dq(3, 3), psi, grad_psi(3)
        integer :: local_status

        sample = gc_cylindrical_field_sample_t()
        status = GC_CYL_EQUILIBRIUM_DOMAIN
        if (position(1) <= 0.0_dp) return
        if (self%invalid_hole) then
            if (position(1) > 2.25_dp) then
                if (position(1) < 2.75_dp) return
            end if
        end if
        b = [0.0_dp, 0.4_dp, 0.0_dp]
        db_dq = 0.0_dp
        psi = 0.5_dp*position(1)**2
        grad_psi = [position(1), 0.0_dp, 0.0_dp]
        call make_gc_cylindrical_field_sample(position(1), b, db_dq, psi, &
            grad_psi, sample, local_status)
        status = local_status
    end subroutine evaluate_manufactured_field

    subroutine evaluate_manufactured_potential(self, position, field, potential, &
            gradient, status)
        class(manufactured_potential_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(out) :: potential, gradient(3)
        integer, intent(out) :: status

        real(dp) :: mu_phys, kinetic, derivative

        potential = 0.0_dp
        gradient = 0.0_dp
        status = GC_CYL_CLASS_INVALID_INPUT
        if (position(1) <= 0.0_dp) return
        mu_phys = self%profile_jperp*abs(self%profile_charge) &
            /(self%profile_mass*self%profile_c)
        kinetic = 0.5_dp*self%profile_mass &
            *target_vparallel_squared(position(1))
        derivative = 0.5_dp*self%profile_mass &
            *target_dvparallel_squared(position(1))
        potential = (self%profile_h0 - mu_phys*field%bmod - kinetic) &
            /self%profile_charge
        gradient(1) = -derivative/self%profile_charge
        status = GC_CYL_SUCCESS
    end subroutine evaluate_manufactured_potential

    subroutine manufactured_cut_map(rc, user_data, position, dposition_drc, &
            status)
        real(dp), intent(in) :: rc
        class(*), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: position(3)
        real(dp), intent(out) :: dposition_drc(3)
        integer, intent(out) :: status

        position = [rc, 0.0_dp, 0.0_dp]
        dposition_drc = [1.0_dp, 0.0_dp, 0.0_dp]
        status = GC_CYL_SUCCESS
    end subroutine manufactured_cut_map

    subroutine identity_splitter(h0, jperp, sigma, candidate, user_data, &
            split_classes, certified, status)
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: sigma
        type(gc_cylindrical_class_interval_t), intent(in) :: candidate
        class(*), pointer, intent(inout) :: user_data
        type(gc_cylindrical_class_interval_t), allocatable, intent(out) :: &
            split_classes(:)
        logical, intent(out) :: certified
        integer, intent(out) :: status

        integer :: base_id

        allocate(split_classes(1))
        split_classes(1) = candidate
        base_id = 1000 + candidate%component_id
        if (sigma < 0) base_id = base_id + 100
        split_classes(1)%component_id = base_id
        certified = .true.
        status = GC_CYL_CLASS_SUCCESS
    end subroutine identity_splitter

    subroutine manufactured_topology_provider(h0, jperp, sigma, rc_min, rc_max, &
            user_data, regions, status)
        real(dp), intent(in) :: h0, jperp, rc_min, rc_max
        integer, intent(in) :: sigma
        class(*), pointer, intent(inout) :: user_data
        type(gc_cylindrical_allowed_region_set_t), intent(out) :: regions
        integer, intent(out) :: status

        real(dp), parameter :: exact_roots(3) = [0.5_dp, 1.0_dp, 2.0_dp]
        integer :: i

        associate (unused_h0 => h0, unused_jperp => jperp, &
                unused_user_data => user_data)
        end associate
        regions = gc_cylindrical_allowed_region_set_t()
        status = GC_CYL_CLASS_INVALID_INPUT
        if (abs(sigma) /= 1) return
        if (rc_min /= 0.5_dp .or. rc_max /= 4.5_dp) return
        regions%nroots = size(exact_roots)
        allocate(regions%roots(regions%nroots), &
            regions%root_canonical(regions%nroots))
        regions%roots = exact_roots
        do i = 1, regions%nroots
            regions%root_canonical(i) = 0.5_dp*exact_roots(i)**2
        end do
        regions%ncomponents = 2
        allocate(regions%components(regions%ncomponents))
        call set_component(regions%components(1), 1, sigma, 0.5_dp, 1.0_dp, &
            .true., .true.)
        call set_component(regions%components(2), 2, sigma, 2.0_dp, 4.5_dp, &
            .true., .false.)
        regions%total_canonical_measure = &
            sum(regions%components%canonical_measure)
        regions%topology_certified = .true.
        regions%certificate_method = 'analytic-quartic-test-oracle'
        status = GC_CYL_SUCCESS

    contains

        subroutine set_component(component, component_id, branch_sigma, &
                x_begin, x_end, lower_root, upper_root)
            use neort_gc_cylindrical_model, only: &
                gc_cylindrical_allowed_component_t
            type(gc_cylindrical_allowed_component_t), intent(out) :: component
            integer, intent(in) :: component_id, branch_sigma
            real(dp), intent(in) :: x_begin, x_end
            logical, intent(in) :: lower_root, upper_root

            component%component_id = component_id
            component%sigma = branch_sigma
            component%x_begin = x_begin
            component%x_end = x_end
            component%canonical_begin = manufactured_psi_star(x_begin, &
                branch_sigma)
            component%canonical_end = manufactured_psi_star(x_end, branch_sigma)
            component%canonical_measure = abs(component%canonical_end &
                -component%canonical_begin)
            component%lower_root = lower_root
            component%upper_root = upper_root
        end subroutine set_component

        pure real(dp) function manufactured_psi_star(radius, branch_sigma)
            real(dp), intent(in) :: radius
            integer, intent(in) :: branch_sigma

            manufactured_psi_star = 0.5_dp*radius**2 &
                +MASS*C_LIGHT/CHARGE*real(branch_sigma, dp) &
                *sqrt(max(0.0_dp, target_vparallel_squared(radius)))*radius
        end function manufactured_psi_star

    end subroutine manufactured_topology_provider

    pure real(dp) function target_vparallel_squared(radius)
        real(dp), intent(in) :: radius

        target_vparallel_squared = 0.1_dp*(radius - 1.0_dp) &
            *(radius - 2.0_dp)*(radius - 0.5_dp)**2
    end function target_vparallel_squared

    pure real(dp) function target_dvparallel_squared(radius)
        real(dp), intent(in) :: radius

        target_dvparallel_squared = 0.1_dp*((radius - 2.0_dp) &
            *(radius - 0.5_dp)**2 + (radius - 1.0_dp) &
            *(radius - 0.5_dp)**2 + 2.0_dp*(radius - 1.0_dp) &
            *(radius - 2.0_dp)*(radius - 0.5_dp))
    end function target_dvparallel_squared

end module gc_cylindrical_class_adapter_test_support

program test_gc_cylindrical_class_adapter
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use gc_cylindrical_class_adapter_test_support, only: &
        C_LIGHT, CHARGE, H0_REFERENCE, JPERP_REFERENCE, MASS, &
        manufactured_cut_map, manufactured_field_t, manufactured_potential_t, &
        identity_splitter, manufactured_topology_provider, &
        target_dvparallel_squared, &
        target_vparallel_squared
    use neort_gc_cylindrical_class_adapter, only: &
        GC_CYL_CLASS_JK_UNITS, &
        GC_CYL_CLASS_INTERIOR_INVALID, GC_CYL_CLASS_SPLITTER_UNAVAILABLE, &
        GC_CYL_CLASS_NOT_ENUMERATED, GC_CYL_CLASS_SPLITTER_FAILURE, &
        GC_CYL_CLASS_SUCCESS, gc_cylindrical_class_adapter_t, &
        gc_cylindrical_class_interval_t, gc_cylindrical_class_launch_t, &
        gc_cylindrical_class_options_t, gc_cylindrical_class_point_t, &
        gc_cylindrical_class_result_t, &
        enumerate_gc_cylindrical_classes, &
        evaluate_gc_cylindrical_class_point, &
        initialize_gc_cylindrical_class_adapter, &
        launch_gc_cylindrical_class

    implicit none

    real(dp), parameter :: RC_MIN = 0.5_dp
    real(dp), parameter :: RC_MAX = 4.5_dp

    type(gc_cylindrical_class_adapter_t) :: adapter
    type(gc_cylindrical_class_adapter_t) :: split_adapter
    type(gc_cylindrical_class_adapter_t) :: thin_adapter
    type(gc_cylindrical_class_adapter_t) :: invalid_adapter
    type(gc_cylindrical_class_adapter_t) :: uncertified_adapter
    type(gc_cylindrical_class_result_t) :: result
    type(gc_cylindrical_class_result_t) :: split_result
    type(gc_cylindrical_class_result_t) :: invalid_result
    type(gc_cylindrical_class_result_t) :: uncertified_result
    type(gc_cylindrical_class_point_t) :: point_plus
    type(gc_cylindrical_class_point_t) :: point_minus
    type(gc_cylindrical_class_point_t) :: thin_point
    type(gc_cylindrical_class_launch_t) :: launch_plus
    type(gc_cylindrical_class_launch_t) :: launch_minus
    type(gc_cylindrical_class_launch_t) :: launch_tangent
    type(gc_cylindrical_class_options_t) :: options
    type(manufactured_field_t), target :: field
    type(manufactured_potential_t), target :: potential
    integer :: status, point_status, plus_id, minus_id, tangent_id
    real(dp) :: omega_c, expected_v2, expected_dv2, expected_v
    real(dp) :: expected_psi, expected_psi_star, expected_dpsi_star
    real(dp) :: thin_speed, thin_h0, thin_shift
    integer :: i, n_plus, n_minus

    options = gc_cylindrical_class_options_t()
    options%scan_points = 65
    options%validation_points = 257

    call initialize_gc_cylindrical_class_adapter(field, potential, &
        H0_REFERENCE, JPERP_REFERENCE, MASS, CHARGE, C_LIGHT, RC_MIN, RC_MAX, &
        manufactured_cut_map, adapter, status, options=options, &
        topology_provider=manufactured_topology_provider)
    call require(status == GC_CYL_CLASS_SUCCESS, &
        'adapter initialization without splitter failed')
    call enumerate_gc_cylindrical_classes(adapter, result, status)
    call require(status == GC_CYL_CLASS_SPLITTER_UNAVAILABLE, &
        'missing orbit-return splitter was not reported unavailable')
    call require(.not. result%class_complete, &
        'allowed intervals were claimed to be complete classes')
    call require(result%nallowed_intervals == 4, &
        'two disconnected allowed regions were not retained for both sigma')
    n_plus = count_sigma(result%allowed_intervals, 1)
    n_minus = count_sigma(result%allowed_intervals, -1)
    call require(n_plus == 2 .and. n_minus == 2, &
        'each parallel branch did not retain both allowed regions')
    call require(has_tangent_lower_endpoint(result%allowed_intervals, 1), &
        'the endpoint tangency was not classified')
    call launch_gc_cylindrical_class(adapter, 0.75_dp, 1, 1, launch_plus, &
        status)
    call require(status == GC_CYL_CLASS_SPLITTER_UNAVAILABLE, &
        'uncertified allowed interval was launchable')

    call evaluate_gc_cylindrical_class_point(adapter, 0.75_dp, 1, point_plus, &
        point_status)
    call require(point_status == GC_CYL_CLASS_SUCCESS, &
        'positive-sigma manufactured point failed')
    call evaluate_gc_cylindrical_class_point(adapter, 0.75_dp, -1, point_minus, &
        point_status)
    call require(point_status == GC_CYL_CLASS_SUCCESS, &
        'negative-sigma manufactured point failed')
    call require(point_plus%allowed .and. point_minus%allowed, &
        'both sigma values were not allowed in the first region')
    call require_close(point_plus%vparallel, -point_minus%vparallel, &
        'sigma did not reverse v_parallel', 1.0e-13_dp)
    call require_close(point_plus%potential_gradient(1), &
        point_plus%potential_gradient(1), 'radial Phi callback is nonfinite', &
        1.0e-13_dp)

    omega_c = abs(CHARGE)*0.4_dp/(MASS*C_LIGHT)
    expected_v2 = target_vparallel_squared(0.75_dp)
    expected_dv2 = target_dvparallel_squared(0.75_dp)
    expected_v = sqrt(expected_v2)
    expected_psi = 0.5_dp*0.75_dp**2
    expected_psi_star = expected_psi + MASS*C_LIGHT/CHARGE &
        *expected_v*0.75_dp
    expected_dpsi_star = 0.75_dp + MASS*C_LIGHT/CHARGE &
        *(expected_dv2/(2.0_dp*expected_v)*0.75_dp + expected_v)
    call require_close(point_plus%omega_c, omega_c, &
        'cyclotron frequency normalization', 1.0e-13_dp)
    call require_close(point_plus%vparallel_squared, expected_v2, &
        'Eq. 12 allowed v_parallel squared', 1.0e-12_dp)
    call require_close(point_plus%dvparallel_squared_drc, expected_dv2, &
        'radial-Phi derivative of v_parallel squared', 1.0e-11_dp)
    call require_close(point_plus%psi_star, expected_psi_star, &
        'non-unit c/q canonical normalization', 1.0e-12_dp)
    call require_close(point_plus%dpsi_star_drc, expected_dpsi_star, &
        'derivative of psi_star along Rc', 1.0e-10_dp)
    call require_close(point_plus%psi_star, point_plus%field%psi &
        +MASS*C_LIGHT/CHARGE*point_plus%vparallel &
        *point_plus%position(1)*point_plus%field%b(2) &
        /point_plus%field%bmod, 'cylindrical covariant B_phi factor', 1.0e-12_dp)

    call initialize_gc_cylindrical_class_adapter(field, potential, &
        H0_REFERENCE, JPERP_REFERENCE, MASS, CHARGE, C_LIGHT, RC_MIN, RC_MAX, &
        manufactured_cut_map, split_adapter, status, options=options, &
        splitter=identity_splitter, &
        topology_provider=manufactured_topology_provider)
    call require(status == GC_CYL_CLASS_SUCCESS, &
        'adapter initialization with splitter failed')
    call enumerate_gc_cylindrical_classes(split_adapter, split_result, status)
    call require(status == GC_CYL_CLASS_SUCCESS, &
        'certified identity splitter was not accepted')
    call require(split_result%class_complete, &
        'orbit-return splitter did not certify class completeness')
    call require(split_result%nclasses == 4, &
        'splitter changed the disconnected interval count unexpectedly')
    do i = 1, split_result%nclasses
        call require(split_result%classes(i)%orbit_return_certified, &
            'uncertified class escaped the splitter seam')
    end do

    plus_id = class_id_at(split_result%classes, 1, 0.75_dp)
    minus_id = class_id_at(split_result%classes, -1, 0.75_dp)
    tangent_id = class_id_at(split_result%classes, 1, RC_MIN)
    call require(plus_id > 0 .and. minus_id > 0 .and. tangent_id > 0, &
        'class IDs were not retained for both sigma and tangency')
    call launch_gc_cylindrical_class(split_adapter, 0.75_dp, 1, plus_id, &
        launch_plus, status)
    call require(status == GC_CYL_CLASS_SUCCESS, &
        'positive-sigma fixed-invariant launch failed')
    call launch_gc_cylindrical_class(split_adapter, 0.75_dp, -1, minus_id, &
        launch_minus, status)
    call require(status == GC_CYL_CLASS_SUCCESS, &
        'negative-sigma fixed-invariant launch failed')
    call require(launch_plus%class_certified .and. launch_minus%class_certified, &
        'launch states did not carry certified class identity')
    call require_close(launch_plus%energy_residual, 0.0_dp, &
        'H0 was not preserved by launch', 1.0e-11_dp)
    call require_close(launch_plus%jperp_residual, 0.0_dp, &
        'J_perp was not preserved by launch', 1.0e-11_dp)
    call require_close(launch_plus%canonical_residual, 0.0_dp, &
        'canonical p_phi was not preserved by launch', 1.0e-11_dp)
    call require(trim(launch_plus%psi_star_units) == 'psi_star=(c/q)*p_phi', &
        'launch units did not record psi_star normalization')
    call require(trim(launch_plus%jperp_units) == trim(GC_CYL_CLASS_JK_UNITS), &
        'launch units did not record J_K normalization')
    call require_close(launch_plus%state%p_parallel, MASS*expected_v, &
        'launch p_parallel is not m*v_parallel', 1.0e-12_dp)
    call require_close(launch_plus%state%mu, &
        JPERP_REFERENCE*abs(CHARGE)/(MASS*C_LIGHT), &
        'launch magnetic moment/action conversion', 1.0e-12_dp)
    call launch_gc_cylindrical_class(split_adapter, RC_MIN, 1, tangent_id, &
        launch_tangent, status)
    call require(status == GC_CYL_CLASS_SUCCESS, &
        'tangency launch was rejected')
    call require(launch_tangent%endpoint_tangent, &
        'endpoint tangency was not retained in launch metadata')

    thin_speed = 1.0e-5_dp
    thin_h0 = JPERP_REFERENCE*abs(CHARGE)/(MASS*C_LIGHT) &
        *point_plus%field%bmod &
        +CHARGE*point_plus%potential + 0.5_dp*MASS*thin_speed**2
    call initialize_gc_cylindrical_class_adapter(field, potential, thin_h0, &
        JPERP_REFERENCE, MASS, CHARGE, C_LIGHT, RC_MIN, RC_MAX, &
        manufactured_cut_map, thin_adapter, status, options=options)
    call require(status == GC_CYL_CLASS_SUCCESS, &
        'thin-width adapter initialization failed')
    call evaluate_gc_cylindrical_class_point(thin_adapter, 0.75_dp, 1, &
        thin_point, point_status)
    call require(point_status == GC_CYL_CLASS_SUCCESS, &
        'thin-width point evaluation failed')
    thin_shift = thin_point%psi_star - thin_point%field%psi
    call require(abs(thin_shift) < 1.0e-3_dp, &
        'psi_star did not approach psi in the thin-width limit')
    call require_close(thin_point%vparallel, thin_speed, &
        'thin-width speed oracle', 1.0e-10_dp)

    field%invalid_hole = .true.
    call initialize_gc_cylindrical_class_adapter(field, potential, &
        H0_REFERENCE, JPERP_REFERENCE, MASS, CHARGE, C_LIGHT, RC_MIN, RC_MAX, &
        manufactured_cut_map, invalid_adapter, status, options=options, &
        topology_provider=manufactured_topology_provider)
    call require(status == GC_CYL_CLASS_SUCCESS, &
        'invalid-hole adapter initialization failed')
    call enumerate_gc_cylindrical_classes(invalid_adapter, invalid_result, status)
    call require(status == GC_CYL_CLASS_INTERIOR_INVALID, &
        'interior invalid hole was not failed closed')
    call require(.not. invalid_result%class_complete, &
        'invalid-hole result claimed complete classes')

    field%invalid_hole = .false.
    call initialize_gc_cylindrical_class_adapter(field, potential, &
        H0_REFERENCE, JPERP_REFERENCE, MASS, CHARGE, C_LIGHT, RC_MIN, RC_MAX, &
        manufactured_cut_map, uncertified_adapter, status, options=options, &
        splitter=uncertified_splitter, &
        topology_provider=manufactured_topology_provider)
    call require(status == GC_CYL_CLASS_SUCCESS, &
        'uncertified-splitter adapter initialization failed')
    call enumerate_gc_cylindrical_classes(uncertified_adapter, &
        uncertified_result, status)
    call require(status == GC_CYL_CLASS_SPLITTER_FAILURE, &
        'uncertified homoclinic topology was accepted')
    call require(.not. uncertified_result%class_complete, &
        'uncertified splitter claimed complete classes')
    call launch_gc_cylindrical_class(uncertified_adapter, 0.75_dp, 1, 1, &
        launch_plus, status)
    call require(status == GC_CYL_CLASS_NOT_ENUMERATED, &
        'failed splitter left stale classes launchable')

    write (*, '(A)') 'test_gc_cylindrical_class_adapter OK'

contains

    subroutine uncertified_splitter(h0, jperp, sigma, candidate, user_data, &
            split_classes, certified, status)
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: sigma
        type(gc_cylindrical_class_interval_t), intent(in) :: candidate
        class(*), pointer, intent(inout) :: user_data
        type(gc_cylindrical_class_interval_t), allocatable, intent(out) :: &
            split_classes(:)
        logical, intent(out) :: certified
        integer, intent(out) :: status

        allocate(split_classes(1))
        split_classes(1) = candidate
        certified = .false.
        status = GC_CYL_CLASS_SUCCESS
    end subroutine uncertified_splitter

    integer function count_sigma(intervals, sigma)
        type(gc_cylindrical_class_interval_t), intent(in) :: intervals(:)
        integer, intent(in) :: sigma
        integer :: i

        count_sigma = 0
        do i = 1, size(intervals)
            if (intervals(i)%sigma == sigma) count_sigma = count_sigma + 1
        end do
    end function count_sigma

    logical function has_tangent_lower_endpoint(intervals, sigma)
        type(gc_cylindrical_class_interval_t), intent(in) :: intervals(:)
        integer, intent(in) :: sigma
        integer :: i

        has_tangent_lower_endpoint = .false.
        do i = 1, size(intervals)
            if (intervals(i)%sigma /= sigma) cycle
            if (intervals(i)%lower_tangent) then
                has_tangent_lower_endpoint = .true.
                return
            end if
        end do
    end function has_tangent_lower_endpoint

    integer function class_id_at(classes, sigma, rc)
        type(gc_cylindrical_class_interval_t), intent(in) :: classes(:)
        integer, intent(in) :: sigma
        real(dp), intent(in) :: rc
        integer :: i

        class_id_at = -1
        do i = 1, size(classes)
            if (classes(i)%sigma /= sigma) cycle
            if (rc < classes(i)%rc_min) cycle
            if (rc > classes(i)%rc_max) cycle
            class_id_at = classes(i)%component_id
            return
        end do
    end function class_id_at

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: message

        if (.not. condition) error stop trim(message)
    end subroutine require

    subroutine require_close(value, expected, message, tolerance)
        real(dp), intent(in) :: value, expected, tolerance
        character(len=*), intent(in) :: message
        real(dp) :: scale

        scale = max(1.0_dp, max(abs(value), abs(expected)))
        call require(abs(value - expected) <= tolerance*scale, message)
    end subroutine require_close

end program test_gc_cylindrical_class_adapter
