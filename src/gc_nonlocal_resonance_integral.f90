module neort_gc_nonlocal_resonance_integral
    !! Nonlocal Eq. 14--17 resonance kernel at fixed (H0,J_perp).
    !!
    !! For every disconnected class c and sigma this module evaluates
    !!
    !!   I_c = integral dx |d psi_star/dx| |H_m|^2 tau_b^2 A_*(psi_star)
    !!             delta(m omega_b + n omega_phi).
    !!
    !! The delta distribution is evaluated at simple roots.  H0 and J_perp
    !! are explicit arguments, so an outer quadrature can apply the phase
    !! space measure from Eq. 14 without hiding a local eta callback here.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_CALLBACK_FAILURE, GC_NONLOCAL_CAPACITY, &
        GC_NONLOCAL_COMPONENT_IDENTITY, GC_NONLOCAL_DERIVATIVE_MISSING, &
        GC_NONLOCAL_FORCE_CONTRACT, GC_NONLOCAL_INVALID_INPUT, &
        GC_NONLOCAL_MAX_COMPONENTS, GC_NONLOCAL_MAX_FORCE_VALUES, &
        GC_NONLOCAL_MAX_ROOTS, GC_NONLOCAL_MAX_SCAN_INTERVALS, &
        GC_NONLOCAL_MAX_ROOT_ITERATIONS, &
        GC_NONLOCAL_NONFINITE, GC_NONLOCAL_PARTIAL, &
        GC_NONLOCAL_ROOT_NOT_CONVERGED, GC_NONLOCAL_SAMPLE_INVALID, &
        GC_NONLOCAL_SAMPLE_UNRESOLVED, GC_NONLOCAL_SAMPLE_VALID, &
        GC_NONLOCAL_SAMPLE_WALL, GC_NONLOCAL_SINGULAR_RESONANCE, &
        GC_NONLOCAL_SUCCESS, gc_nonlocal_component_t, &
        gc_nonlocal_orbit_sample_t, gc_nonlocal_resonance_options_t, &
        gc_nonlocal_resonance_result_t

    implicit none
    private

    abstract interface
        subroutine gc_nonlocal_orbit_evaluator_i(h0, jperp, x, sigma, &
                component_id, sample, status)
            import :: dp, gc_nonlocal_orbit_sample_t
            real(dp), intent(in) :: h0, jperp, x
            integer, intent(in) :: sigma, component_id
            type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
            integer, intent(out) :: status
        end subroutine gc_nonlocal_orbit_evaluator_i
    end interface

    public :: gc_nonlocal_orbit_evaluator_i
    public :: integrate_gc_nonlocal_resonance
    public :: evaluate_gc_nonlocal_root_contribution

contains

    subroutine integrate_gc_nonlocal_resonance(evaluate, h0, jperp, harmonic_m, &
            harmonic_n, components, options, result, status)
        procedure(gc_nonlocal_orbit_evaluator_i) :: evaluate
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: harmonic_m, harmonic_n
        type(gc_nonlocal_component_t), intent(in) :: components(:)
        type(gc_nonlocal_resonance_options_t), intent(in) :: options
        type(gc_nonlocal_resonance_result_t), intent(out) :: result
        integer, intent(out) :: status

        integer :: component_index, component_status
        logical :: partial

        result = gc_nonlocal_resonance_result_t()
        result%h0 = h0
        result%jperp = jperp
        result%harmonic_m = harmonic_m
        result%harmonic_n = harmonic_n
        status = validate_request(h0, jperp, harmonic_m, harmonic_n, components, &
            options)
        if (status /= GC_NONLOCAL_SUCCESS) then
            result%status = status
            return
        end if
        call initialize_result(result, components, options)
        partial = .false.
        do component_index = 1, size(components)
            call integrate_component(evaluate, h0, jperp, harmonic_m, harmonic_n, &
                components(component_index), options, component_index, result, &
                component_status)
            if (.not. is_recoverable_status(component_status)) then
                status = component_status
                result%status = status
                return
            end if
            if (component_status == GC_NONLOCAL_PARTIAL) partial = .true.
        end do
        status = GC_NONLOCAL_SUCCESS
        if (partial) status = GC_NONLOCAL_PARTIAL
        result%status = status
        result%certified = status == GC_NONLOCAL_SUCCESS
    end subroutine integrate_gc_nonlocal_resonance

    pure integer function validate_request(h0, jperp, harmonic_m, harmonic_n, &
            components, options) result(status)
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: harmonic_m, harmonic_n
        type(gc_nonlocal_component_t), intent(in) :: components(:)
        type(gc_nonlocal_resonance_options_t), intent(in) :: options

        status = GC_NONLOCAL_INVALID_INPUT
        if (.not. ieee_is_finite(h0)) return
        if (.not. ieee_is_finite(jperp)) return
        if (jperp < 0.0_dp) return
        if (harmonic_m == 0 .and. harmonic_n == 0) return
        if (size(components) < 1) return
        if (size(components) > GC_NONLOCAL_MAX_COMPONENTS) return
        if (.not. valid_options(options)) return
        status = validate_components(components)
    end function validate_request

    pure logical function valid_options(options)
        type(gc_nonlocal_resonance_options_t), intent(in) :: options

        valid_options = .false.
        if (options%scan_intervals < 1) return
        if (options%scan_intervals > GC_NONLOCAL_MAX_SCAN_INTERVALS) return
        if (options%max_root_iterations < 1) return
        if (options%max_root_iterations > GC_NONLOCAL_MAX_ROOT_ITERATIONS) return
        if (options%max_roots < 1) return
        if (options%max_roots > GC_NONLOCAL_MAX_ROOTS) return
        if (options%force_count < 1) return
        if (options%force_count > GC_NONLOCAL_MAX_FORCE_VALUES) return
        if (.not. ieee_is_finite(options%residual_tolerance)) return
        if (.not. ieee_is_finite(options%x_tolerance)) return
        if (.not. ieee_is_finite(options%derivative_tolerance)) return
        if (options%residual_tolerance <= 0.0_dp) return
        if (options%x_tolerance <= 0.0_dp) return
        if (options%derivative_tolerance <= 0.0_dp) return
        valid_options = .true.
    end function valid_options

    pure integer function validate_components(components) result(status)
        type(gc_nonlocal_component_t), intent(in) :: components(:)

        integer :: i, j

        status = GC_NONLOCAL_INVALID_INPUT
        do i = 1, size(components)
            if (components(i)%component_id <= 0) return
            if (abs(components(i)%sigma) /= 1) return
            if (.not. ieee_is_finite(components(i)%x_min)) return
            if (.not. ieee_is_finite(components(i)%x_max)) return
            if (components(i)%x_max <= components(i)%x_min) return
            do j = 1, i - 1
                if (same_component_identity(components(i), components(j))) then
                    status = GC_NONLOCAL_COMPONENT_IDENTITY
                    return
                end if
                if (components(i)%sigma == components(j)%sigma) then
                    if (intervals_overlap(components(i), components(j))) return
                end if
            end do
        end do
        status = GC_NONLOCAL_SUCCESS
    end function validate_components

    pure logical function intervals_overlap(first, second)
        type(gc_nonlocal_component_t), intent(in) :: first, second
        real(dp) :: overlap_left, overlap_right

        overlap_left = max(first%x_min, second%x_min)
        overlap_right = min(first%x_max, second%x_max)
        intervals_overlap = overlap_right > overlap_left
    end function intervals_overlap

    pure logical function same_component_identity(first, second)
        type(gc_nonlocal_component_t), intent(in) :: first, second

        same_component_identity = first%sigma == second%sigma
        if (same_component_identity) then
            same_component_identity = first%component_id == second%component_id
        end if
    end function same_component_identity

    subroutine initialize_result(result, components, options)
        type(gc_nonlocal_resonance_result_t), intent(inout) :: result
        type(gc_nonlocal_component_t), intent(in) :: components(:)
        type(gc_nonlocal_resonance_options_t), intent(in) :: options

        integer :: i

        result%nforce = options%force_count
        result%ncomponents = size(components)
        allocate(result%component_id(result%ncomponents), result%sigma(result%ncomponents))
        allocate(result%component_contribution(result%nforce, result%ncomponents))
        allocate(result%root_x(options%max_roots), &
            result%root_residual_derivative(options%max_roots), &
            result%root_component_id(options%max_roots), &
            result%root_sigma(options%max_roots), &
            result%root_force_contribution(result%nforce, options%max_roots))
        result%component_contribution = 0.0_dp
        result%root_x = 0.0_dp
        result%root_residual_derivative = 0.0_dp
        result%root_component_id = 0
        result%root_sigma = 0
        result%root_force_contribution = 0.0_dp
        do i = 1, result%ncomponents
            result%component_id(i) = components(i)%component_id
            result%sigma(i) = components(i)%sigma
        end do
    end subroutine initialize_result

    subroutine integrate_component(evaluate, h0, jperp, harmonic_m, harmonic_n, &
            component, options, component_index, result, status)
        procedure(gc_nonlocal_orbit_evaluator_i) :: evaluate
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: harmonic_m, harmonic_n, component_index
        type(gc_nonlocal_component_t), intent(in) :: component
        type(gc_nonlocal_resonance_options_t), intent(in) :: options
        type(gc_nonlocal_resonance_result_t), intent(inout) :: result
        integer, intent(out) :: status

        type(gc_nonlocal_orbit_sample_t), allocatable :: samples(:)
        real(dp), allocatable :: residual(:), residual_derivative(:)
        logical, allocatable :: valid(:)
        integer :: nscan, i, local_status
        logical :: partial

        nscan = options%scan_intervals
        allocate(samples(0:nscan), residual(0:nscan), residual_derivative(0:nscan), &
            valid(0:nscan))
        partial = .false.
        do i = 0, nscan
            call evaluate_component_sample(evaluate, h0, jperp, harmonic_m, &
                harmonic_n, component, i, nscan, options, samples(i), residual(i), &
                residual_derivative(i), valid(i), result, local_status)
            if (.not. is_recoverable_status(local_status)) then
                status = local_status
                return
            end if
            if (local_status == GC_NONLOCAL_PARTIAL) partial = .true.
        end do
        call process_valid_segments(evaluate, h0, jperp, harmonic_m, harmonic_n, &
            component, component_index, options, samples, residual, &
            residual_derivative, valid, result, local_status)
        if (.not. is_recoverable_status(local_status)) then
            status = local_status
            return
        end if
        if (local_status == GC_NONLOCAL_PARTIAL) partial = .true.
        status = GC_NONLOCAL_SUCCESS
        if (partial) status = GC_NONLOCAL_PARTIAL
    end subroutine integrate_component

    subroutine evaluate_component_sample(evaluate, h0, jperp, harmonic_m, &
            harmonic_n, component, index, nscan, options, sample, residual, &
            residual_derivative, valid, result, status)
        procedure(gc_nonlocal_orbit_evaluator_i) :: evaluate
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: harmonic_m, harmonic_n, index, nscan
        type(gc_nonlocal_component_t), intent(in) :: component
        type(gc_nonlocal_resonance_options_t), intent(in) :: options
        type(gc_nonlocal_orbit_sample_t), intent(out) :: sample
        real(dp), intent(out) :: residual, residual_derivative
        logical, intent(out) :: valid
        type(gc_nonlocal_resonance_result_t), intent(inout) :: result
        integer, intent(out) :: status

        real(dp) :: x
        integer :: callback_status, sample_status

        x = component%x_min + real(index, dp)*(component%x_max - component%x_min) &
            / real(nscan, dp)
        call evaluate(h0, jperp, x, component%sigma, component%component_id, &
            sample, callback_status)
        if (callback_status /= GC_NONLOCAL_SUCCESS) then
            status = GC_NONLOCAL_CALLBACK_FAILURE
            valid = .false.
            return
        end if
        call classify_sample(sample, component, options%force_count, valid, &
            sample_status)
        if (sample_status == GC_NONLOCAL_PARTIAL) then
            result%unresolved_samples = result%unresolved_samples + 1
            if (sample%status == GC_NONLOCAL_SAMPLE_WALL) then
                result%wall_samples = result%wall_samples + 1
            end if
            status = GC_NONLOCAL_PARTIAL
            return
        end if
        if (sample_status /= GC_NONLOCAL_SUCCESS) then
            status = sample_status
            return
        end if
        residual = real(harmonic_m, dp)*sample%omega_b &
            +real(harmonic_n, dp)*sample%omega_phi
        residual_derivative = real(harmonic_m, dp)*sample%domega_b_dx &
            +real(harmonic_n, dp)*sample%domega_phi_dx
        if (.not. all(ieee_is_finite([residual, residual_derivative]))) then
            status = GC_NONLOCAL_NONFINITE
            valid = .false.
            return
        end if
        result%valid_samples = result%valid_samples + 1
        status = GC_NONLOCAL_SUCCESS
    end subroutine evaluate_component_sample

    subroutine classify_sample(sample, component, force_count, valid, status)
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        type(gc_nonlocal_component_t), intent(in) :: component
        integer, intent(in) :: force_count
        logical, intent(out) :: valid
        integer, intent(out) :: status

        valid = .false.
        status = GC_NONLOCAL_INVALID_INPUT
        if (sample%status /= GC_NONLOCAL_SAMPLE_VALID) then
            if (sample%status == GC_NONLOCAL_SAMPLE_WALL) status = GC_NONLOCAL_PARTIAL
            if (sample%status == GC_NONLOCAL_SAMPLE_UNRESOLVED) then
                status = GC_NONLOCAL_PARTIAL
            end if
            if (sample%status == GC_NONLOCAL_SAMPLE_INVALID) status = GC_NONLOCAL_PARTIAL
            return
        end if
        if (sample%component_id /= component%component_id) then
            status = GC_NONLOCAL_COMPONENT_IDENTITY
            return
        end if
        if (sample%sigma /= component%sigma) then
            status = GC_NONLOCAL_COMPONENT_IDENTITY
            return
        end if
        if (sample%nforce /= force_count) then
            status = GC_NONLOCAL_FORCE_CONTRACT
            return
        end if
        if (.not. sample%derivatives_available) then
            status = GC_NONLOCAL_DERIVATIVE_MISSING
            return
        end if
        if (.not. finite_sample(sample, force_count)) then
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        if (sample%tau_b <= 0.0_dp) return
        valid = .true.
        status = GC_NONLOCAL_SUCCESS
    end subroutine classify_sample

    pure logical function finite_sample(sample, force_count)
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        integer, intent(in) :: force_count
        real(dp) :: values(7)

        values = [sample%psi_star, sample%dpsi_star_dx, sample%tau_b, &
            sample%omega_b, sample%omega_phi, sample%domega_b_dx, &
            sample%domega_phi_dx]
        finite_sample = all(ieee_is_finite(values))
        if (.not. finite_sample) return
        if (.not. ieee_is_finite(real(sample%h_m, dp))) then
            finite_sample = .false.
            return
        end if
        if (.not. ieee_is_finite(aimag(sample%h_m))) then
            finite_sample = .false.
            return
        end if
        finite_sample = all(ieee_is_finite(sample%thermodynamic_force(1:force_count)))
    end function finite_sample

    subroutine process_valid_segments(evaluate, h0, jperp, harmonic_m, harmonic_n, &
            component, component_index, options, samples, residual, &
            residual_derivative, valid, result, status)
        procedure(gc_nonlocal_orbit_evaluator_i) :: evaluate
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: harmonic_m, harmonic_n, component_index
        type(gc_nonlocal_component_t), intent(in) :: component
        type(gc_nonlocal_resonance_options_t), intent(in) :: options
        type(gc_nonlocal_orbit_sample_t), intent(in) :: samples(0:)
        real(dp), intent(in) :: residual(0:), residual_derivative(0:)
        logical, intent(in) :: valid(0:)
        type(gc_nonlocal_resonance_result_t), intent(inout) :: result
        integer, intent(out) :: status

        integer :: i, segment_start, local_status
        logical :: partial

        status = GC_NONLOCAL_SUCCESS
        partial = .false.
        segment_start = -1
        do i = 0, ubound(valid, 1)
            if (valid(i)) then
                if (segment_start < 0) segment_start = i
            else
                if (segment_start >= 0) then
                    call process_segment(evaluate, h0, jperp, harmonic_m, harmonic_n, &
                        component, component_index, options, samples, residual, &
                        residual_derivative, segment_start, i - 1, result, &
                        local_status)
                    if (.not. is_recoverable_status(local_status)) then
                        status = local_status
                        return
                    end if
                    if (local_status == GC_NONLOCAL_PARTIAL) partial = .true.
                end if
                segment_start = -1
            end if
        end do
        if (segment_start >= 0) then
            call process_segment(evaluate, h0, jperp, harmonic_m, harmonic_n, &
                component, component_index, options, samples, residual, &
                residual_derivative, segment_start, ubound(valid, 1), result, &
                local_status)
            if (.not. is_recoverable_status(local_status)) then
                status = local_status
                return
            end if
            if (local_status == GC_NONLOCAL_PARTIAL) partial = .true.
        end if
        if (partial) status = GC_NONLOCAL_PARTIAL
    end subroutine process_valid_segments

    subroutine process_segment(evaluate, h0, jperp, harmonic_m, harmonic_n, &
            component, component_index, options, samples, residual, &
            residual_derivative, first, last, result, status)
        procedure(gc_nonlocal_orbit_evaluator_i) :: evaluate
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: harmonic_m, harmonic_n, component_index, first, last
        type(gc_nonlocal_component_t), intent(in) :: component
        type(gc_nonlocal_resonance_options_t), intent(in) :: options
        type(gc_nonlocal_orbit_sample_t), intent(in) :: samples(0:)
        real(dp), intent(in) :: residual(0:), residual_derivative(0:)
        type(gc_nonlocal_resonance_result_t), intent(inout) :: result
        integer, intent(out) :: status

        integer :: i, local_status
        logical :: partial
        type(gc_nonlocal_orbit_sample_t) :: root_sample
        real(dp) :: root_x

        status = GC_NONLOCAL_SUCCESS
        partial = .false.
        do i = first, last - 1
            if (is_exact_zero(residual(i))) then
                call record_root(samples(i), residual_derivative(i), component, &
                    component_index, options, result, &
                    grid_coordinate(component, i, options), local_status)
                if (local_status /= GC_NONLOCAL_SUCCESS) then
                    status = local_status
                    return
                end if
            end if
            if (residual(i)*residual(i + 1) < 0.0_dp) then
                call locate_root(evaluate, h0, jperp, harmonic_m, harmonic_n, &
                    component, options, grid_coordinate(component, i, options), &
                    grid_coordinate(component, i + 1, options), residual(i), &
                    residual(i + 1), root_x, &
                    root_sample, local_status)
                if (local_status /= GC_NONLOCAL_SUCCESS) then
                    status = local_status
                    return
                end if
                call record_root(root_sample, residual_derivative_at_root(root_sample, &
                    harmonic_m, harmonic_n), component, component_index, options, &
                    result, root_x, local_status)
                if (local_status /= GC_NONLOCAL_SUCCESS) then
                    status = local_status
                    return
                end if
            end if
            if (residual_derivative(i)*residual_derivative(i + 1) < 0.0_dp) then
                partial = .true.
            end if
        end do
        if (last >= first) then
            if (is_exact_zero(residual(last))) then
                call record_root(samples(last), residual_derivative(last), component, &
                    component_index, options, result, &
                    grid_coordinate(component, last, options), local_status)
                if (local_status /= GC_NONLOCAL_SUCCESS) then
                    status = local_status
                    return
                end if
            end if
        end if
        if (partial) status = GC_NONLOCAL_PARTIAL
    end subroutine process_segment

    subroutine locate_root(evaluate, h0, jperp, harmonic_m, harmonic_n, component, &
            options, left_coordinate, right_coordinate, left_residual, right_residual, &
            root_x, root_sample, status)
        procedure(gc_nonlocal_orbit_evaluator_i) :: evaluate
        real(dp), intent(in) :: h0, jperp, left_coordinate, right_coordinate, &
            left_residual, right_residual
        integer, intent(in) :: harmonic_m, harmonic_n
        type(gc_nonlocal_component_t), intent(in) :: component
        type(gc_nonlocal_resonance_options_t), intent(in) :: options
        real(dp), intent(out) :: root_x
        type(gc_nonlocal_orbit_sample_t), intent(out) :: root_sample
        integer, intent(out) :: status

        real(dp) :: left_x, right_x, left_value, right_value, middle_x, middle_value
        type(gc_nonlocal_orbit_sample_t) :: middle_sample
        integer :: iteration, callback_status, sample_status
        logical :: valid

        left_x = left_coordinate
        right_x = right_coordinate
        left_value = left_residual
        right_value = right_residual
        root_x = 0.5_dp*(left_x + right_x)
        root_sample = gc_nonlocal_orbit_sample_t()
        status = GC_NONLOCAL_ROOT_NOT_CONVERGED
        do iteration = 1, options%max_root_iterations
            middle_x = 0.5_dp*(left_x + right_x)
            call evaluate(h0, jperp, middle_x, component%sigma, &
                component%component_id, middle_sample, callback_status)
            if (callback_status /= GC_NONLOCAL_SUCCESS) then
                status = GC_NONLOCAL_CALLBACK_FAILURE
                return
            end if
            call classify_sample(middle_sample, component, options%force_count, &
                valid, sample_status)
            if (sample_status /= GC_NONLOCAL_SUCCESS) then
                status = sample_status
                return
            end if
            middle_value = real(harmonic_m, dp)*middle_sample%omega_b &
                +real(harmonic_n, dp)*middle_sample%omega_phi
            if (.not. ieee_is_finite(middle_value)) then
                status = GC_NONLOCAL_NONFINITE
                return
            end if
            root_x = middle_x
            root_sample = middle_sample
            if (abs(middle_value) <= options%residual_tolerance) then
                status = GC_NONLOCAL_SUCCESS
                return
            end if
            if (right_x - left_x <= options%x_tolerance) then
                status = GC_NONLOCAL_SUCCESS
                return
            end if
            if (left_value*middle_value < 0.0_dp) then
                right_x = middle_x
                right_value = middle_value
            else
                left_x = middle_x
                left_value = middle_value
            end if
        end do
    end subroutine locate_root

    pure real(dp) function grid_coordinate(component, index, options)
        type(gc_nonlocal_component_t), intent(in) :: component
        integer, intent(in) :: index
        type(gc_nonlocal_resonance_options_t), intent(in) :: options

        grid_coordinate = component%x_min + real(index, dp) &
            *(component%x_max - component%x_min)/real(options%scan_intervals, dp)
    end function grid_coordinate

    pure real(dp) function residual_derivative_at_root(sample, harmonic_m, harmonic_n)
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        integer, intent(in) :: harmonic_m, harmonic_n

        residual_derivative_at_root = real(harmonic_m, dp)*sample%domega_b_dx &
            +real(harmonic_n, dp)*sample%domega_phi_dx
    end function residual_derivative_at_root

    pure logical function is_exact_zero(value)
        real(dp), intent(in) :: value

        is_exact_zero = abs(value) <= tiny(value)
    end function is_exact_zero

    subroutine record_root(sample, residual_derivative, component, component_index, &
            options, result, root_x, status)
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        real(dp), intent(in) :: residual_derivative
        type(gc_nonlocal_component_t), intent(in) :: component
        integer, intent(in) :: component_index
        real(dp), intent(in) :: root_x
        type(gc_nonlocal_resonance_options_t), intent(in) :: options
        type(gc_nonlocal_resonance_result_t), intent(inout) :: result
        integer, intent(out) :: status

        real(dp) :: root_contribution(GC_NONLOCAL_MAX_FORCE_VALUES)
        integer :: root_number

        status = GC_NONLOCAL_SUCCESS
        if (root_already_recorded(result, component%component_id, component%sigma, &
            root_x, &
            options%x_tolerance)) return
        call evaluate_gc_nonlocal_root_contribution(sample, residual_derivative, &
            options%force_count, options%derivative_tolerance, root_contribution, &
            status)
        if (status /= GC_NONLOCAL_SUCCESS) return
        if (result%nroots >= size(result%root_x)) then
            status = GC_NONLOCAL_CAPACITY
            return
        end if
        root_number = result%nroots + 1
        result%nroots = root_number
        result%root_x(root_number) = root_x
        result%root_residual_derivative(root_number) = residual_derivative
        result%root_component_id(root_number) = component%component_id
        result%root_sigma(root_number) = component%sigma
        result%root_force_contribution(:, root_number) = &
            root_contribution(1:options%force_count)
        result%contribution(1:options%force_count) = &
            result%contribution(1:options%force_count) + &
            root_contribution(1:options%force_count)
        result%component_contribution(1:options%force_count, component_index) = &
            result%component_contribution(1:options%force_count, component_index) + &
            root_contribution(1:options%force_count)
    end subroutine record_root

    pure logical function root_already_recorded(result, component_id, sigma, x, &
            tolerance)
        type(gc_nonlocal_resonance_result_t), intent(in) :: result
        integer, intent(in) :: component_id, sigma
        real(dp), intent(in) :: x, tolerance
        integer :: i

        root_already_recorded = .false.
        do i = 1, result%nroots
            if (result%root_component_id(i) /= component_id) cycle
            if (result%root_sigma(i) /= sigma) cycle
            if (abs(result%root_x(i) - x) <= tolerance) then
                root_already_recorded = .true.
                return
            end if
        end do
    end function root_already_recorded

    pure subroutine evaluate_gc_nonlocal_root_contribution(sample, &
            residual_derivative, force_count, derivative_tolerance, contribution, &
            status)
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        real(dp), intent(in) :: residual_derivative, derivative_tolerance
        integer, intent(in) :: force_count
        real(dp), intent(out) :: contribution(:)
        integer, intent(out) :: status

        real(dp) :: hamiltonian_square, weight

        contribution = 0.0_dp
        status = GC_NONLOCAL_INVALID_INPUT
        if (force_count < 1 .or. force_count > GC_NONLOCAL_MAX_FORCE_VALUES) return
        if (size(contribution) < force_count) return
        if (sample%status /= GC_NONLOCAL_SAMPLE_VALID) return
        if (sample%nforce /= force_count) then
            status = GC_NONLOCAL_FORCE_CONTRACT
            return
        end if
        if (.not. sample%derivatives_available) then
            status = GC_NONLOCAL_DERIVATIVE_MISSING
            return
        end if
        if (.not. finite_sample(sample, force_count)) then
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        if (.not. ieee_is_finite(residual_derivative)) then
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        if (abs(residual_derivative) <= derivative_tolerance) then
            status = GC_NONLOCAL_SINGULAR_RESONANCE
            return
        end if
        hamiltonian_square = real(sample%h_m*conjg(sample%h_m), dp)
        if (.not. ieee_is_finite(hamiltonian_square)) then
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        weight = abs(sample%dpsi_star_dx)*hamiltonian_square*sample%tau_b**2 &
            /abs(residual_derivative)
        if (.not. ieee_is_finite(weight)) then
            status = GC_NONLOCAL_NONFINITE
            return
        end if
        contribution(1:force_count) = weight*sample%thermodynamic_force(1:force_count)
        if (.not. all(ieee_is_finite(contribution(1:force_count)))) then
            status = GC_NONLOCAL_NONFINITE
            contribution = 0.0_dp
            return
        end if
        status = GC_NONLOCAL_SUCCESS
    end subroutine evaluate_gc_nonlocal_root_contribution

    pure logical function is_recoverable_status(status)
        integer, intent(in) :: status

        is_recoverable_status = status == GC_NONLOCAL_SUCCESS
        if (status == GC_NONLOCAL_PARTIAL) is_recoverable_status = .true.
    end function is_recoverable_status

end module neort_gc_nonlocal_resonance_integral
