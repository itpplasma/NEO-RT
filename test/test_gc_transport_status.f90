program test_gc_transport_status
    use neort_transport, only: gc_transport_failure_t, gc_transport_failure_code, &
        GC_TRANSPORT_SUCCESS, GC_TRANSPORT_FULL_ORBIT_FAILURE
    use neort_gc_orbit_integrator, only: GC_ORBIT_FIELD_ERROR, GC_ORBIT_STATE_ERROR, &
        GC_ORBIT_START_ROOT_ERROR, GC_ORBIT_INTEGRATOR_ERROR, GC_ORBIT_NO_RETURN, &
        GC_ORBIT_PERTURBATION_ERROR, GC_ORBIT_UNCONFINED, GC_ORBIT_RADIAL_DOMAIN, &
        GC_ORBIT_WALL_LOSS, gc_orbit_status_is_physical_loss, &
        gc_orbit_status_is_unconfined

    implicit none

    type(gc_transport_failure_t) :: failures

    failures = gc_transport_failure_t()
    if (gc_transport_failure_code(failures) /= GC_TRANSPORT_SUCCESS) &
        error stop 'empty full-orbit status is not successful'

    failures = gc_transport_failure_t()
    failures%lost_orbits = 1
    failures%wall_orbits = 1
    failures%confined_coverage_fraction = 0.75
    failures%physical_coverage_fraction = 1.0
    if (gc_transport_failure_code(failures) /= GC_TRANSPORT_SUCCESS) &
        error stop 'physical loss counters were treated as numerical failure'

    if (gc_orbit_status_is_physical_loss(GC_ORBIT_NO_RETURN) .or. &
        .not. gc_orbit_status_is_physical_loss(GC_ORBIT_WALL_LOSS) .or. &
        gc_orbit_status_is_physical_loss(GC_ORBIT_RADIAL_DOMAIN)) &
        error stop 'orbit physical-loss classification is incomplete'
    if (gc_orbit_status_is_unconfined(GC_ORBIT_NO_RETURN) .or. &
        .not. gc_orbit_status_is_unconfined(GC_ORBIT_UNCONFINED)) &
        error stop 'timeout and explicit unconfined statuses were aliased'
    if (gc_orbit_status_is_physical_loss(GC_ORBIT_FIELD_ERROR) .or. &
        gc_orbit_status_is_physical_loss(GC_ORBIT_STATE_ERROR) .or. &
        gc_orbit_status_is_physical_loss(GC_ORBIT_START_ROOT_ERROR) .or. &
        gc_orbit_status_is_physical_loss(GC_ORBIT_INTEGRATOR_ERROR) .or. &
        gc_orbit_status_is_physical_loss(GC_ORBIT_PERTURBATION_ERROR)) &
        error stop 'numerical orbit status was classified as physical loss'

    failures = gc_transport_failure_t()
    failures%unconfined_samples = 1
    if (gc_transport_failure_code(failures) /= GC_TRANSPORT_FULL_ORBIT_FAILURE) &
        error stop 'unconfined timeout was not fail-closed'

    failures = gc_transport_failure_t()
    failures%radial_domain_orbits = 1
    if (gc_transport_failure_code(failures) /= GC_TRANSPORT_FULL_ORBIT_FAILURE) &
        error stop 'unresolved radial-domain sample was not fail-closed'

    failures = gc_transport_failure_t()
    failures%scanned_samples = 1
    failures%canonical_measure_certified = .false.
    if (gc_transport_failure_code(failures) /= GC_TRANSPORT_FULL_ORBIT_FAILURE) &
        error stop 'uncertified canonical measure was accepted'

    failures = gc_transport_failure_t()
    failures%scanned_samples = 1
    failures%canonical_measure_certified = .true.
    failures%component_identity_certified = .false.
    if (gc_transport_failure_code(failures) /= GC_TRANSPORT_FULL_ORBIT_FAILURE) &
        error stop 'uncertified disconnected topology was accepted'

    failures = gc_transport_failure_t()
    failures%resonance_partial = 1
    if (gc_transport_failure_code(failures) /= GC_TRANSPORT_FULL_ORBIT_FAILURE) &
        error stop 'partial resonance scan was not machine-verifiable failure'

    failures%frequency_failures = 2
    failures%phase_failures = 3
    failures%orbit_failures = 4
    if (gc_transport_failure_code(failures) /= GC_TRANSPORT_FULL_ORBIT_FAILURE) &
        error stop 'full-orbit failures were not aggregated'

    write (*, '(A)') 'test_gc_transport_status OK'
end program test_gc_transport_status
