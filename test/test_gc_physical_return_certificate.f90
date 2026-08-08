program test_gc_physical_return_certificate
    !! The fixture below is deliberately a hand-built numerical observation.
    !! The certificate provider is a separate test oracle; the observation's
    !! event count is never used to manufacture its theorem result.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: GC_CYL_SUCCESS
    use neort_gc_cylindrical_physical_return, only: &
        GC_CYL_PHYSICAL_EVENT_RETURN, &
        GC_CYL_PHYSICAL_RETURN_CERTIFICATE_UNAVAILABLE, &
        GC_CYL_PHYSICAL_RETURN_CERTIFICATE_INVALID, &
        GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_CERTIFIED, &
        GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_UNKNOWN, &
        attach_gc_cylindrical_physical_return_certificate, &
        certify_gc_cylindrical_physical_return, &
        gc_cylindrical_physical_return_certificate_t, &
        gc_cylindrical_physical_return_t
    implicit none

    type(gc_cylindrical_physical_return_t) :: evidence
    type(gc_cylindrical_physical_return_certificate_t) :: certificate
    type(gc_cylindrical_physical_return_t) :: provider_evidence
    integer :: status
    character(len=256) :: message

    call test_provider_statuses()

    evidence = numerical_two_event_evidence()
    certificate = valid_certificate()
    call attach_gc_cylindrical_physical_return_certificate(evidence, certificate, &
        status, message)
    call require(status == GC_CYL_SUCCESS, 'valid independent certificate rejected')
    call require(evidence%intersection_multiplicity_certified, &
        'valid certificate did not mark multiplicity certified')
    call require(evidence%multiplicity_status == &
        GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_CERTIFIED, &
        'certified multiplicity status is not explicit')

    evidence = numerical_two_event_evidence()
    certificate = gc_cylindrical_physical_return_certificate_t()
    call attach_gc_cylindrical_physical_return_certificate(evidence, certificate, &
        status, message)
    call require(status == GC_CYL_PHYSICAL_RETURN_CERTIFICATE_INVALID, &
        'empty certificate was accepted')
    call require(.not. evidence%intersection_multiplicity_certified .and. &
        evidence%multiplicity_status == GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_UNKNOWN, &
        'invalid certificate did not fail closed')

    evidence = numerical_two_event_evidence()
    evidence%intersection_rates(1) = 0.0_dp
    certificate = valid_certificate()
    call attach_gc_cylindrical_physical_return_certificate(evidence, certificate, &
        status, message)
    call require(status == GC_CYL_PHYSICAL_RETURN_CERTIFICATE_INVALID, &
        'non-transverse evidence was accepted')

    write (*, '(a)') 'PASS: physical return certificate contract'

contains

    subroutine test_provider_statuses()
        provider_evidence = numerical_two_event_evidence()
        call certify_gc_cylindrical_physical_return(provider_evidence, &
            unavailable_provider, status, message)
        call require(status == GC_CYL_PHYSICAL_RETURN_CERTIFICATE_UNAVAILABLE, &
            'unavailable provider did not remain unavailable')
        call require(provider_evidence%multiplicity_status == &
            GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_UNKNOWN, &
            'unavailable provider did not fail closed')

        provider_evidence = numerical_two_event_evidence()
        call certify_gc_cylindrical_physical_return(provider_evidence, &
            invalid_provider, status, message)
        call require(status == GC_CYL_PHYSICAL_RETURN_CERTIFICATE_INVALID, &
            'malformed provider certificate was accepted')
        call require(provider_evidence%multiplicity_status == &
            GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_UNKNOWN, &
            'malformed provider did not fail closed')

        provider_evidence = numerical_two_event_evidence()
        call certify_gc_cylindrical_physical_return(provider_evidence, &
            valid_provider, status, message)
        call require(status == GC_CYL_SUCCESS .and. &
            provider_evidence%multiplicity_status == &
            GC_CYL_PHYSICAL_RETURN_MULTIPLICITY_CERTIFIED, &
            'valid provider certificate was not attached')
    end subroutine test_provider_statuses

    subroutine unavailable_provider(evidence, certificate, provider_status)
        type(gc_cylindrical_physical_return_t), intent(in) :: evidence
        type(gc_cylindrical_physical_return_certificate_t), intent(out) :: certificate
        integer, intent(out) :: provider_status

        associate (unused_evidence => evidence)
        end associate
        certificate = gc_cylindrical_physical_return_certificate_t()
        provider_status = GC_CYL_PHYSICAL_RETURN_CERTIFICATE_UNAVAILABLE
    end subroutine unavailable_provider

    subroutine invalid_provider(evidence, certificate, provider_status)
        type(gc_cylindrical_physical_return_t), intent(in) :: evidence
        type(gc_cylindrical_physical_return_certificate_t), intent(out) :: certificate
        integer, intent(out) :: provider_status

        associate (unused_evidence => evidence)
        end associate
        certificate = gc_cylindrical_physical_return_certificate_t()
        provider_status = GC_CYL_SUCCESS
    end subroutine invalid_provider

    subroutine valid_provider(evidence, certificate, provider_status)
        type(gc_cylindrical_physical_return_t), intent(in) :: evidence
        type(gc_cylindrical_physical_return_certificate_t), intent(out) :: certificate
        integer, intent(out) :: provider_status

        associate (unused_evidence => evidence)
        end associate
        certificate = valid_certificate()
        provider_status = GC_CYL_SUCCESS
    end subroutine valid_provider

    function numerical_two_event_evidence() result(value)
        type(gc_cylindrical_physical_return_t) :: value

        value = gc_cylindrical_physical_return_t()
        value%status = GC_CYL_SUCCESS
        value%event_kind = GC_CYL_PHYSICAL_EVENT_RETURN
        value%return_orientation = 1
        value%period = 2.0_dp
        value%intersection_count = 2
        value%intersection_orientations = [-1, 1]
        value%intersection_times = [1.0_dp, 2.0_dp]
        value%intersection_rates = [-1.0_dp, 1.0_dp]
        value%physical_return_found = .true.
        value%numerical_failure = .false.
    end function numerical_two_event_evidence

    function valid_certificate() result(value)
        type(gc_cylindrical_physical_return_certificate_t) :: value

        value = gc_cylindrical_physical_return_certificate_t()
        value%certificate_id = 901
        value%crossing_count = 2
        value%exactly_two_proved = .true.
    end function valid_certificate

    subroutine require(condition, text)
        logical, intent(in) :: condition
        character(len=*), intent(in) :: text

        if (.not. condition) error stop text
    end subroutine require

end program test_gc_physical_return_certificate
