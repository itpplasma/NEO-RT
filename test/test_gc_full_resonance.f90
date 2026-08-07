program test_gc_full_resonance
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_full_resonance, only: GC_RESONANCE_SUCCESS, GC_RESONANCE_INVALID_INPUT, &
        GC_RESONANCE_PARTIAL, &
        GC_RESONANCE_BOUNDARY_INVALID, GC_RESONANCE_SAMPLE_VALID, &
        GC_RESONANCE_SAMPLE_UNCONFINED, GC_RESONANCE_SAMPLE_WALL, &
        GC_RESONANCE_SAMPLE_RADIAL_DOMAIN, &
        GC_RESONANCE_SAMPLE_INVALID, gc_resonance_diagnostics_t, &
        find_gc_resonances

    implicit none

    real(dp) :: roots(8), derivatives(8)
    integer :: nroots, status
    type(gc_resonance_diagnostics_t) :: diagnostics

    call find_gc_resonances(invalid_open_endpoints, 0.0_dp, 1.0_dp, 8, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 1) &
        error stop "invalid open endpoints were treated as a partial scan"
    if (abs(roots(1) - 0.35_dp) > 1.0e-9_dp) &
        error stop "root in the shrunken open interval was not recovered"

    call find_gc_resonances(positive_local_minimum, 0.0_dp, 1.0_dp, 3, &
        1.0e-12_dp, 1.0e-10_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 0) &
        error stop "positive local minimum was reported as a failed scan/root"

    call find_gc_resonances(two_roots, 0.0_dp, 1.0_dp, 40, 1.0e-12_dp, &
        1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 2) &
        error stop "manufactured roots not recovered"
    if (maxval(abs(roots(1:2) - [0.25_dp, 0.75_dp])) > 1.0e-10_dp) &
        error stop "manufactured root positions mismatch"
    if (abs(derivatives(1) + 0.5_dp) > 1.0e-6_dp .or. &
        abs(derivatives(2) - 0.5_dp) > 1.0e-6_dp) &
        error stop "manufactured root derivatives mismatch"

    call find_gc_resonances(tangent_root, 0.0_dp, 1.0_dp, 3, &
        1.0e-12_dp, 1.0e-10_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 1) &
        error stop "tangent root was not recovered"
    if (abs(roots(1) - 0.37_dp) > 1.0e-9_dp) &
        error stop "tangent root position mismatch"
    if (abs(derivatives(1)) > 1.0e-7_dp) &
        error stop "tangent root derivative mismatch"

    call find_gc_resonances(close_sign_roots, 0.0_dp, 1.0_dp, 4, &
        1.0e-12_dp, 1.0e-10_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 2) &
        error stop "closely spaced roots were not both recovered"
    if (maxval(abs(roots(1:2) - [0.401_dp, 0.409_dp])) > 1.0e-9_dp) &
        error stop "closely spaced root positions mismatch"

    call find_gc_resonances(root_hidden_by_loss, 0.0_dp, 1.0_dp, 20, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_PARTIAL .or. nroots /= 0) &
        error stop "no-return interval was not exposed"

    call find_gc_resonances(one_bad_one_good, 0.0_dp, 1.0_dp, 4, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status)
    if (status /= GC_RESONANCE_PARTIAL .or. nroots /= 2) &
        error stop "valid root was discarded after a failed bracket"
    if (maxval(abs(roots(1:2) - [0.2_dp, 0.8_dp])) > 1.0e-10_dp) &
        error stop "retained root position mismatch"

    call find_gc_resonances(wall_region, 0.0_dp, 1.0_dp, 20, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status, &
        diagnostics, classify_statuses, canonical_measure)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 2) &
        error stop "physical loss components were treated as a partial scan"
    if (maxval(abs(roots(1:2) - [0.2_dp, 0.8_dp])) > 1.0e-10_dp) &
        error stop "root in a confined component was lost across physical loss"
    if (diagnostics%lost_orbits <= 0 .or. diagnostics%wall_orbits <= 0 .or. &
        diagnostics%radial_domain_orbits /= 0) &
        error stop "certified wall-loss status counters were not recorded"
    if (diagnostics%confined_coverage_fraction <= 0.0_dp .or. &
        diagnostics%physical_coverage_fraction <= &
        diagnostics%confined_coverage_fraction) &
        error stop "physical coverage fractions were not recorded"
    if (.not. diagnostics%canonical_measure_certified .or. &
        abs(diagnostics%canonical_scan_measure - 1.5_dp) > 1.0e-10_dp) &
        error stop "canonical Jacobian measure was not applied"

    call find_gc_resonances(numerical_gap, 0.0_dp, 1.0_dp, 20, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status, &
        diagnostics, classify_statuses, canonical_measure)
    if (status /= GC_RESONANCE_PARTIAL .or. nroots /= 2) &
        error stop "numerical interior gap was not fail-closed"
    if (diagnostics%numerical_samples <= 0 .or. diagnostics%lost_orbits /= 0) &
        error stop "numerical gap was misclassified as physical loss"
    if (diagnostics%canonical_unresolved_measure <= 0.0_dp) &
        error stop "numerical gap had no unresolved measure"

    call find_gc_resonances(two_roots, 0.0_dp, 1.0_dp, 20, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status, &
        diagnostics, classify_statuses, measure_failure)
    if (status /= GC_RESONANCE_PARTIAL .or. nroots /= 2) &
        error stop "unknown canonical measure was not fail-closed"
    if (diagnostics%measure_failures <= 0 .or. &
        diagnostics%unknown_measure_cells <= 0 .or. &
        diagnostics%unknown_measure_coordinate_span <= 0.0_dp .or. &
        diagnostics%canonical_measure_certified .or. &
        diagnostics%physical_coverage_fraction /= 0.0_dp) &
        error stop "unknown measure was treated as zero coverage"

    call find_gc_resonances(unconfined_region, 0.0_dp, 1.0_dp, 20, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status, &
        diagnostics, classify_statuses, canonical_measure)
    if (status /= GC_RESONANCE_PARTIAL .or. diagnostics%unconfined_samples <= 0) &
        error stop "unconfined timeout was not fail-closed"

    call find_gc_resonances(radial_region, 0.0_dp, 1.0_dp, 20, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status, &
        diagnostics, classify_statuses, canonical_measure)
    if (status /= GC_RESONANCE_PARTIAL .or. diagnostics%radial_domain_orbits <= 0) &
        error stop "computational radial-domain gap was not fail-closed"
    if (diagnostics%lost_orbits /= 0 .or. diagnostics%physical_coverage_fraction >= 1.0_dp) &
        error stop "computational radial-domain gap was certified as physical loss"

    call find_gc_resonances(component_split, 0.0_dp, 1.0_dp, 20, &
        1.0e-12_dp, 1.0e-12_dp, roots, derivatives, nroots, status, &
        diagnostics, classify_statuses, canonical_measure, component_key)
    if (status /= GC_RESONANCE_SUCCESS .or. nroots /= 2 .or. &
        diagnostics%component_count /= 2 .or. &
        .not. diagnostics%component_identity_certified) &
        error stop "component-key segmentation did not preserve both roots"

    write (*, '(A)') "test_gc_full_resonance OK"

contains

    subroutine invalid_open_endpoints(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = eta - 0.35_dp
        local_status = merge(GC_RESONANCE_BOUNDARY_INVALID, 0, &
            eta == 0.0_dp .or. eta == 1.0_dp)
    end subroutine invalid_open_endpoints

    subroutine positive_local_minimum(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.37_dp)**2 + 1.0e-4_dp
        local_status = 0
    end subroutine positive_local_minimum

    subroutine two_roots(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.25_dp)*(eta - 0.75_dp)
        local_status = 0
    end subroutine two_roots

    subroutine tangent_root(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.37_dp)**2
        local_status = 0
    end subroutine tangent_root

    subroutine close_sign_roots(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.401_dp)*(eta - 0.409_dp)
        local_status = 0
    end subroutine close_sign_roots

    subroutine root_hidden_by_loss(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = eta - 0.5_dp
        local_status = merge(1, 0, eta >= 0.45_dp .and. eta <= 0.55_dp)
    end subroutine root_hidden_by_loss

    subroutine one_bad_one_good(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.2_dp)*(eta - 0.8_dp)
        local_status = merge(1, 0, eta > 0.1_dp .and. eta < 0.15_dp)
    end subroutine one_bad_one_good

    subroutine wall_region(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.2_dp)*(eta - 0.8_dp)
        local_status = GC_RESONANCE_SUCCESS
        if (eta >= 0.4_dp .and. eta <= 0.6_dp) local_status = 12
    end subroutine wall_region

    subroutine numerical_gap(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.2_dp)*(eta - 0.8_dp)
        local_status = merge(13, GC_RESONANCE_SUCCESS, &
            eta >= 0.4_dp .and. eta <= 0.6_dp)
    end subroutine numerical_gap

    subroutine unconfined_region(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.2_dp)*(eta - 0.8_dp)
        local_status = merge(14, GC_RESONANCE_SUCCESS, &
            eta >= 0.4_dp .and. eta <= 0.6_dp)
    end subroutine unconfined_region

    subroutine radial_region(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.2_dp)*(eta - 0.8_dp)
        local_status = merge(11, GC_RESONANCE_SUCCESS, &
            eta >= 0.4_dp .and. eta <= 0.6_dp)
    end subroutine radial_region

    subroutine component_split(eta, residual, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: residual
        integer, intent(out) :: local_status

        residual = (eta - 0.2_dp)*(eta - 0.8_dp)
        local_status = GC_RESONANCE_SUCCESS
    end subroutine component_split

    subroutine canonical_measure(eta, density, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: density
        integer, intent(out) :: local_status

        density = 1.0_dp + eta
        local_status = GC_RESONANCE_SUCCESS
    end subroutine canonical_measure

    subroutine measure_failure(eta, density, local_status)
        real(dp), intent(in) :: eta
        real(dp), intent(out) :: density
        integer, intent(out) :: local_status

        density = 1.0_dp
        local_status = merge(GC_RESONANCE_INVALID_INPUT, GC_RESONANCE_SUCCESS, &
            eta >= 0.4_dp .and. eta <= 0.6_dp)
    end subroutine measure_failure

    integer function component_key(eta, local_status)
        real(dp), intent(in) :: eta
        integer, intent(in) :: local_status

        associate (unused_status => local_status)
        end associate
        component_key = merge(1, 2, eta < 0.5_dp)
    end function component_key

    integer function classify_statuses(local_status)
        integer, intent(in) :: local_status

        select case (local_status)
        case (GC_RESONANCE_SUCCESS)
            classify_statuses = GC_RESONANCE_SAMPLE_VALID
        case (11)
            classify_statuses = GC_RESONANCE_SAMPLE_RADIAL_DOMAIN
        case (12)
            classify_statuses = GC_RESONANCE_SAMPLE_WALL
        case (13)
            classify_statuses = GC_RESONANCE_SAMPLE_INVALID
        case (14)
            classify_statuses = GC_RESONANCE_SAMPLE_UNCONFINED
        case default
            classify_statuses = GC_RESONANCE_SAMPLE_INVALID
        end select
    end function classify_statuses

end program test_gc_full_resonance
