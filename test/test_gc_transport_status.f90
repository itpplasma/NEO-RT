program test_gc_transport_status
    use neort_transport, only: gc_transport_failure_t, gc_transport_failure_code, &
        GC_TRANSPORT_SUCCESS, GC_TRANSPORT_FULL_ORBIT_FAILURE

    implicit none

    type(gc_transport_failure_t) :: failures

    failures = gc_transport_failure_t()
    if (gc_transport_failure_code(failures) /= GC_TRANSPORT_SUCCESS) &
        error stop 'empty full-orbit status is not successful'

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
