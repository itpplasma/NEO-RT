module neort_gc_frequency_splines
    !! Adapter from the invariant-preserving real-space frequency provider to
    !! NEO-RT's velocity/pitch frequency interface.  The orbit calculation is
    !! performed once at the reference speed; exact thin-limit scalings then
    !! give omega_b~v, Omega_B~v^2, and Omega_E~v^0.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use spline, only: spline_coeff, spline_val_0
    use neort_gc_frequency_provider, only: GC_FREQUENCY_SUCCESS, &
        gc_frequency_context_t, gc_frequency_result_t, &
        initialize_gc_frequency_context, evaluate_gc_frequency, &
        evaluate_gc_phase_average, evaluate_gc_full_orbit_frequency, &
        evaluate_gc_full_orbit_phase_average, gc_full_orbit_frequency_result_t
    use neort_gc_orbit_integrator, only: GC_ORBIT_TRAPPED, GC_ORBIT_PASSING, &
        gc_orbit_average_t, gc_orbit_perturbation_i
    use neort_thin_orbit_limit, only: THIN_LIMIT_RETURN_ERROR

    implicit none
    private

    integer, parameter, public :: GC_SPLINE_SUCCESS = 0
    integer, parameter, public :: GC_SPLINE_INVALID_INPUT = 1
    integer, parameter, public :: GC_SPLINE_PROVIDER_ERROR = 2
    integer, parameter, public :: GC_SPLINE_NOT_INITIALIZED = 3
    integer, parameter :: MAX_KNOTS = 100

    type :: gc_region_spline_t
        real(dp) :: omega_b_coeff(MAX_KNOTS - 1, 5) = 0.0_dp
        real(dp) :: omega_magnetic_coeff(MAX_KNOTS - 1, 5) = 0.0_dp
        real(dp) :: omega_electric_coeff(MAX_KNOTS - 1, 5) = 0.0_dp
        real(dp) :: eta_min = 0.0_dp
        real(dp) :: eta_max = 0.0_dp
        real(dp) :: eta_tp = 0.0_dp
        real(dp) :: period_log_slope = 0.0_dp
        real(dp) :: period_log_intercept = 0.0_dp
        real(dp) :: magnetic_ratio_log_slope = 0.0_dp
        real(dp) :: magnetic_ratio_log_intercept = 0.0_dp
        real(dp) :: electric_ratio_log_slope = 0.0_dp
        real(dp) :: electric_ratio_log_intercept = 0.0_dp
        integer :: knot_count = 0
        integer :: orbit_class = 0
        logical :: initialized = .false.
    end type gc_region_spline_t

    type, public :: gc_spline_diagnostics_t
        integer :: orbit_evaluations = 0
        integer :: maximum_refinements = 0
        real(dp) :: elapsed_seconds = 0.0_dp
        real(dp) :: maximum_magnetic_relative_error = 0.0_dp
        real(dp) :: maximum_electric_relative_error = 0.0_dp
        real(dp) :: minimum_observed_order = huge(1.0_dp)
        real(dp) :: minimum_lambda = huge(1.0_dp)
        integer :: failed_provider_status = 0
        integer :: failed_magnetic_status = 0
        integer :: failed_total_status = 0
        real(dp) :: failed_eta = 0.0_dp
        real(dp) :: failed_magnetic_error = 0.0_dp
        real(dp) :: failed_electric_error = 0.0_dp
        real(dp) :: failed_magnetic_order = 0.0_dp
        real(dp) :: failed_total_order = 0.0_dp
        integer :: boundary_extrapolations = 0
    end type gc_spline_diagnostics_t

    type(gc_frequency_context_t) :: surface_context
    type(gc_region_spline_t) :: trapped_spline, passing_spline
    type(gc_spline_diagnostics_t) :: diagnostics
    real(dp) :: reference_speed = 0.0_dp
    logical :: surface_initialized = .false.
    !$omp threadprivate (surface_context, trapped_spline, passing_spline)
    !$omp threadprivate (diagnostics, reference_speed, surface_initialized)

    public :: initialize_gc_spline_surface, fit_gc_frequency_region
    public :: evaluate_gc_spline, gc_spline_q, get_gc_spline_diagnostics
    public :: evaluate_gc_phase_average_surface
    public :: evaluate_gc_full_orbit_frequency_surface
    public :: evaluate_gc_full_orbit_phase_average_surface

contains

    subroutine initialize_gc_spline_surface(surface, reference_theta, &
            field_scale, omega_e, particle_mass, particle_charge, velocity, &
            status, selected_frequency_model, wall_path, wall_units)
        real(dp), intent(in) :: surface, reference_theta, field_scale, omega_e
        real(dp), intent(in) :: particle_mass, particle_charge, velocity
        integer, intent(out) :: status
        integer, intent(in), optional :: selected_frequency_model
        character(len=*), intent(in), optional :: wall_path, wall_units

        integer :: provider_status

        surface_initialized = .false.
        reference_speed = 0.0_dp
        trapped_spline = gc_region_spline_t()
        passing_spline = gc_region_spline_t()
        diagnostics = gc_spline_diagnostics_t()
        call initialize_gc_frequency_context(surface, reference_theta, &
            field_scale, omega_e, particle_mass, particle_charge, velocity, &
            surface_context, provider_status, selected_frequency_model, &
            wall_path, wall_units)
        if (provider_status /= GC_FREQUENCY_SUCCESS) then
            status = GC_SPLINE_PROVIDER_ERROR
            return
        end if
        reference_speed = velocity
        surface_initialized = .true.
        status = GC_SPLINE_SUCCESS
    end subroutine initialize_gc_spline_surface

    subroutine fit_gc_frequency_region(eta, period_estimate, eta_tp, &
            orbit_class, status, failed_knot)
        real(dp), intent(in) :: eta(:), period_estimate(:), eta_tp
        integer, intent(in) :: orbit_class
        integer, intent(out) :: status, failed_knot

        type(gc_region_spline_t) :: region
        type(gc_frequency_result_t) :: result
        real(dp) :: omega_b(MAX_KNOTS), omega_magnetic(MAX_KNOTS)
        real(dp) :: omega_electric(MAX_KNOTS)
        real(dp) :: start_time, end_time, magnetic_scale, electric_scale
        integer :: k, n, n_fit, provider_status
        logical :: boundary_extrapolation

        status = GC_SPLINE_INVALID_INPUT
        failed_knot = 0
        if (.not. surface_initialized) then
            status = GC_SPLINE_NOT_INITIALIZED
            return
        end if
        n = size(eta)
        if (n < 3 .or. n > MAX_KNOTS .or. size(period_estimate) /= n) return
        if (any(eta(2:n) <= eta(1:n - 1)) .or. any(period_estimate <= 0.0_dp)) return
        if (orbit_class /= GC_ORBIT_TRAPPED .and. &
            orbit_class /= GC_ORBIT_PASSING) return
        if (orbit_class == GC_ORBIT_TRAPPED .and. any(eta <= eta_tp)) return
        if (orbit_class == GC_ORBIT_PASSING .and. any(eta >= eta_tp)) return

        call cpu_time(start_time)
        n_fit = n
        boundary_extrapolation = .false.
        do k = 1, n
            call evaluate_gc_frequency(surface_context, eta(k), 1, orbit_class, &
                period_estimate(k), result, provider_status)
            if (provider_status /= GC_FREQUENCY_SUCCESS) then
                if (orbit_class == GC_ORBIT_PASSING .and. k >= n - 1 &
                    .and. result%magnetic_limit_status == THIN_LIMIT_RETURN_ERROR &
                    .and. result%total_limit_status == THIN_LIMIT_RETURN_ERROR) then
                    ! The last passing knot is deliberately close to the
                    ! trapped-passing separatrix.  A real-space return is
                    ! undefined there, while the logarithmic limit below is
                    ! the correct continuation of the last two converged
                    ! knots.  Interior failures remain fatal.
                    n_fit = k - 1
                    boundary_extrapolation = .true.
                    exit
                end if
                diagnostics%failed_provider_status = provider_status
                diagnostics%failed_magnetic_status = result%magnetic_limit_status
                diagnostics%failed_total_status = result%total_limit_status
                diagnostics%failed_eta = eta(k)
                diagnostics%failed_magnetic_error = result%magnetic_error
                diagnostics%failed_electric_error = result%electric_error
                diagnostics%failed_magnetic_order = result%magnetic_order
                diagnostics%failed_total_order = result%total_order
                failed_knot = k
                status = GC_SPLINE_PROVIDER_ERROR
                return
            end if
            omega_b(k) = abs(result%omega_b)
            omega_magnetic(k) = result%omega_magnetic
            omega_electric(k) = result%omega_electric
            diagnostics%orbit_evaluations = diagnostics%orbit_evaluations + 1
            diagnostics%maximum_refinements = max( &
                diagnostics%maximum_refinements, result%maximum_refinements)
            magnetic_scale = max(1.0_dp, abs(result%omega_magnetic))
            electric_scale = max(1.0_dp, abs(result%omega_electric))
            diagnostics%maximum_magnetic_relative_error = max( &
                diagnostics%maximum_magnetic_relative_error, &
                result%magnetic_error/magnetic_scale)
            diagnostics%maximum_electric_relative_error = max( &
                diagnostics%maximum_electric_relative_error, &
                result%electric_error/electric_scale)
            diagnostics%minimum_observed_order = min( &
                diagnostics%minimum_observed_order, result%magnetic_order, &
                result%total_order)
            diagnostics%minimum_lambda = min(diagnostics%minimum_lambda, &
                minval(result%lambda_used))
        end do
        call cpu_time(end_time)
        diagnostics%elapsed_seconds = diagnostics%elapsed_seconds &
            + max(0.0_dp, end_time - start_time)
        if (n_fit < 3) then
            status = GC_SPLINE_PROVIDER_ERROR
            return
        end if
        if (boundary_extrapolation) then
            diagnostics%boundary_extrapolations = &
                diagnostics%boundary_extrapolations + 1
        end if

        region = gc_region_spline_t()
        region%knot_count = n_fit
        region%orbit_class = orbit_class
        region%eta_min = eta(1)
        region%eta_max = eta(n_fit)
        region%eta_tp = eta_tp
        region%omega_b_coeff(1:n_fit - 1, :) = spline_coeff( &
            eta(1:n_fit), omega_b(1:n_fit))
        region%omega_magnetic_coeff(1:n_fit - 1, :) = &
            spline_coeff(eta(1:n_fit), omega_magnetic(1:n_fit))
        region%omega_electric_coeff(1:n_fit - 1, :) = &
            spline_coeff(eta(1:n_fit), omega_electric(1:n_fit))
        call fit_log_extrapolation(region, eta(1:n_fit), omega_b, &
            omega_magnetic, omega_electric)
        region%initialized = .true.
        if (orbit_class == GC_ORBIT_TRAPPED) then
            trapped_spline = region
        else
            passing_spline = region
        end if
        status = GC_SPLINE_SUCCESS
    end subroutine fit_gc_frequency_region

    subroutine evaluate_gc_spline(velocity, eta, parallel_direction, &
            orbit_class, omega_b, domega_b_dv, domega_b_deta, &
            omega_magnetic, domega_magnetic_dv, domega_magnetic_deta, &
            omega_electric, domega_electric_deta, status)
        real(dp), intent(in) :: velocity, eta
        integer, intent(in) :: parallel_direction, orbit_class
        real(dp), intent(out) :: omega_b, domega_b_dv, domega_b_deta
        real(dp), intent(out) :: omega_magnetic, domega_magnetic_dv
        real(dp), intent(out) :: domega_magnetic_deta
        real(dp), intent(out) :: omega_electric, domega_electric_deta
        integer, intent(out) :: status

        real(dp) :: omega_b_ref, domega_b_ref_deta
        real(dp) :: omega_magnetic_ref, domega_magnetic_ref_deta
        real(dp) :: omega_electric_ref, domega_electric_ref_deta, speed_ratio

        omega_b = 0.0_dp
        domega_b_dv = 0.0_dp
        domega_b_deta = 0.0_dp
        omega_magnetic = 0.0_dp
        domega_magnetic_dv = 0.0_dp
        domega_magnetic_deta = 0.0_dp
        omega_electric = 0.0_dp
        domega_electric_deta = 0.0_dp
        status = GC_SPLINE_INVALID_INPUT
        if (.not. surface_initialized .or. velocity <= 0.0_dp &
            .or. eta <= 0.0_dp .or. abs(parallel_direction) /= 1) return

        if (orbit_class == GC_ORBIT_TRAPPED) then
            call evaluate_region(trapped_spline, eta, omega_b_ref, &
                domega_b_ref_deta, omega_magnetic_ref, &
                domega_magnetic_ref_deta, omega_electric_ref, &
                domega_electric_ref_deta, status)
        else if (orbit_class == GC_ORBIT_PASSING) then
            call evaluate_region(passing_spline, eta, omega_b_ref, &
                domega_b_ref_deta, omega_magnetic_ref, &
                domega_magnetic_ref_deta, omega_electric_ref, &
                domega_electric_ref_deta, status)
        else
            return
        end if
        if (status /= GC_SPLINE_SUCCESS) return

        speed_ratio = velocity/reference_speed
        omega_b = real(parallel_direction, dp)*omega_b_ref*speed_ratio
        domega_b_dv = real(parallel_direction, dp)*omega_b_ref/reference_speed
        domega_b_deta = real(parallel_direction, dp) &
            *domega_b_ref_deta*speed_ratio
        omega_magnetic = omega_magnetic_ref*speed_ratio**2
        domega_magnetic_dv = 2.0_dp*omega_magnetic_ref*speed_ratio &
            /reference_speed
        domega_magnetic_deta = domega_magnetic_ref_deta*speed_ratio**2
        omega_electric = omega_electric_ref
        domega_electric_deta = domega_electric_ref_deta
    end subroutine evaluate_gc_spline

    function gc_spline_q() result(q_value)
        real(dp) :: q_value

        q_value = 0.0_dp
        if (surface_initialized) q_value = surface_context%q_fieldline
    end function gc_spline_q

    subroutine evaluate_gc_full_orbit_frequency_surface(velocity, eta, &
            parallel_direction, orbit_class, period_estimate, result, status)
        real(dp), intent(in) :: velocity, eta, period_estimate
        integer, intent(in) :: parallel_direction, orbit_class
        type(gc_full_orbit_frequency_result_t), intent(out) :: result
        integer, intent(out) :: status

        result = gc_full_orbit_frequency_result_t()
        status = GC_SPLINE_NOT_INITIALIZED
        if (.not. surface_initialized) return
        call evaluate_gc_full_orbit_frequency(surface_context, eta, &
            parallel_direction, orbit_class, period_estimate, result, status, &
            velocity)
    end subroutine evaluate_gc_full_orbit_frequency_surface

    subroutine evaluate_gc_full_orbit_phase_average_surface(velocity, eta, &
            parallel_direction, orbit_class, period_estimate, omega_b, omega_phi, mth, mph, &
            perturbation, result, status)
        real(dp), intent(in) :: velocity, eta, period_estimate, omega_b, omega_phi
        integer, intent(in) :: parallel_direction, orbit_class, mth, mph
        procedure(gc_orbit_perturbation_i) :: perturbation
        type(gc_orbit_average_t), intent(out) :: result
        integer, intent(out) :: status

        result = gc_orbit_average_t()
        status = GC_SPLINE_NOT_INITIALIZED
        if (.not. surface_initialized) return
        call evaluate_gc_full_orbit_phase_average(surface_context, eta, &
            parallel_direction, orbit_class, period_estimate, omega_b, omega_phi, mth, mph, &
            perturbation, result, status, velocity)
    end subroutine evaluate_gc_full_orbit_phase_average_surface

    subroutine get_gc_spline_diagnostics(value)
        type(gc_spline_diagnostics_t), intent(out) :: value

        value = diagnostics
    end subroutine get_gc_spline_diagnostics

    subroutine evaluate_gc_phase_average_surface(velocity, eta, &
            parallel_direction, orbit_class, period_estimate, omega_b, omega_phi, &
            q_fieldline, mth, mph, &
            perturbation, result, status)
        real(dp), intent(in) :: velocity, eta, period_estimate, omega_b
        real(dp), intent(in) :: omega_phi, q_fieldline
        integer, intent(in) :: parallel_direction, orbit_class, mth, mph
        procedure(gc_orbit_perturbation_i) :: perturbation
        type(gc_orbit_average_t), intent(out) :: result
        integer, intent(out) :: status

        result = gc_orbit_average_t()
        status = GC_SPLINE_NOT_INITIALIZED
        if (.not. surface_initialized) return
        call evaluate_gc_phase_average(surface_context, velocity, eta, &
            parallel_direction, orbit_class, period_estimate, omega_b, omega_phi, &
            q_fieldline, mth, mph, perturbation, result, status)
    end subroutine evaluate_gc_phase_average_surface

    subroutine fit_log_extrapolation(region, eta, omega_b, omega_magnetic, &
            omega_electric)
        type(gc_region_spline_t), intent(inout) :: region
        real(dp), intent(in) :: eta(:), omega_b(:), omega_magnetic(:)
        real(dp), intent(in) :: omega_electric(:)

        real(dp) :: log_distance(2), period_length(2)
        real(dp) :: magnetic_ratio(2), electric_ratio(2)
        integer :: index(2)

        if (region%orbit_class == GC_ORBIT_TRAPPED) then
            index = [1, 2]
        else
            index = [size(eta) - 1, size(eta)]
        end if
        log_distance = log(abs(eta(index) - region%eta_tp))
        period_length = 2.0_dp*acos(-1.0_dp)*reference_speed/omega_b(index)
        magnetic_ratio = omega_magnetic(index)/omega_b(index)
        electric_ratio = omega_electric(index)/omega_b(index)
        call fit_line(log_distance, period_length, region%period_log_slope, &
            region%period_log_intercept)
        call fit_line(log_distance, magnetic_ratio, &
            region%magnetic_ratio_log_slope, &
            region%magnetic_ratio_log_intercept)
        call fit_line(log_distance, electric_ratio, &
            region%electric_ratio_log_slope, &
            region%electric_ratio_log_intercept)
    end subroutine fit_log_extrapolation

    pure subroutine fit_line(x, y, slope, intercept)
        real(dp), intent(in) :: x(2), y(2)
        real(dp), intent(out) :: slope, intercept

        slope = (y(2) - y(1))/(x(2) - x(1))
        intercept = y(1) - slope*x(1)
    end subroutine fit_line

    subroutine evaluate_region(region, eta, omega_b, domega_b_deta, &
            omega_magnetic, domega_magnetic_deta, omega_electric, &
            domega_electric_deta, status)
        type(gc_region_spline_t), intent(in) :: region
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: omega_b, domega_b_deta
        real(dp), intent(out) :: omega_magnetic, domega_magnetic_deta
        real(dp), intent(out) :: omega_electric, domega_electric_deta
        integer, intent(out) :: status

        real(dp) :: value(3), log_distance, period_length
        real(dp) :: ratio, dratio_deta, distance
        logical :: extrapolate_to_boundary

        omega_b = 0.0_dp
        domega_b_deta = 0.0_dp
        omega_magnetic = 0.0_dp
        domega_magnetic_deta = 0.0_dp
        omega_electric = 0.0_dp
        domega_electric_deta = 0.0_dp
        status = GC_SPLINE_NOT_INITIALIZED
        if (.not. region%initialized) return
        if (region%orbit_class == GC_ORBIT_TRAPPED) then
            if (eta <= region%eta_tp) then
                status = GC_SPLINE_INVALID_INPUT
                return
            end if
            extrapolate_to_boundary = eta < region%eta_min
        else
            if (eta >= region%eta_tp) then
                status = GC_SPLINE_INVALID_INPUT
                return
            end if
            extrapolate_to_boundary = eta > region%eta_max
        end if

        if (.not. extrapolate_to_boundary) then
            value = spline_val_0(region%omega_b_coeff( &
                1:region%knot_count - 1, :), eta)
            omega_b = value(1)
            domega_b_deta = value(2)
            value = spline_val_0(region%omega_magnetic_coeff( &
                1:region%knot_count - 1, :), eta)
            omega_magnetic = value(1)
            domega_magnetic_deta = value(2)
            value = spline_val_0(region%omega_electric_coeff( &
                1:region%knot_count - 1, :), eta)
            omega_electric = value(1)
            domega_electric_deta = value(2)
            status = GC_SPLINE_SUCCESS
            return
        end if

        distance = eta - region%eta_tp
        log_distance = log(abs(distance))
        period_length = region%period_log_slope*log_distance &
            + region%period_log_intercept
        if (period_length <= 0.0_dp) then
            status = GC_SPLINE_INVALID_INPUT
            return
        end if
        omega_b = 2.0_dp*acos(-1.0_dp)*reference_speed/period_length
        domega_b_deta = -omega_b**2/(2.0_dp*acos(-1.0_dp)*reference_speed) &
            *region%period_log_slope/distance

        ratio = region%magnetic_ratio_log_slope*log_distance &
            + region%magnetic_ratio_log_intercept
        dratio_deta = region%magnetic_ratio_log_slope/distance
        omega_magnetic = ratio*omega_b
        domega_magnetic_deta = dratio_deta*omega_b &
            + ratio*domega_b_deta
        ratio = region%electric_ratio_log_slope*log_distance &
            + region%electric_ratio_log_intercept
        dratio_deta = region%electric_ratio_log_slope/distance
        omega_electric = ratio*omega_b
        domega_electric_deta = dratio_deta*omega_b &
            + ratio*domega_b_deta
        status = GC_SPLINE_SUCCESS
    end subroutine evaluate_region

end module neort_gc_frequency_splines
