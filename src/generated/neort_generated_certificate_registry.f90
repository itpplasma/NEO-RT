module neort_generated_certificate_registry
    implicit none
    character(*), parameter :: fortsym_revision = &
        'fortsym@58a0e06c95ecc943dfdcb044b7ca6a9964c1c55d'
    character(*), parameter :: regenerate_command = &
        'cd tools/gc_symbolics && fo exec gen_full_fow_physics ../../src/generated'
    integer, parameter :: certificate_count = 8
    character(len=32), parameter :: certificate_id(certificate_count) = &
        [character(len=32) :: 'geometry', 'littlejohn', 'eq13_cdot', 'boundary_limits', &
        'root_enclosures', 'interpolation', 'profile_endpoints', &
        'refinement' ]
    character(len=64), parameter :: certificate_fingerprint(certificate_count) = &
        [character(len=64) :: 'neort-cert-v1:geometry:19:fortsym-58a0e06', &
        'neort-cert-v1:littlejohn:22:fortsym-58a0e06', &
        'neort-cert-v1:eq13_cdot:3:fortsym-58a0e06', &
        'neort-cert-v1:boundary_limits:13:fortsym-58a0e06', &
        'neort-cert-v1:root_enclosures:3:fortsym-58a0e06', &
        'neort-cert-v1:interpolation:7:fortsym-58a0e06', &
        'neort-cert-v1:profile_endpoints:8:fortsym-58a0e06', &
        'neort-cert-v1:refinement:14:fortsym-58a0e06' ]
    ! Fingerprints are provenance/arity manifests, not algebraic proofs.
    ! Root multiplicity and crossing counts require interval/theorem gates.
contains
    pure integer function certificate_index(id) result(index)
        character(*), intent(in) :: id
        integer :: k
        index = 0
        do k = 1, certificate_count
            if (id == certificate_id(k)) index = k
        end do
    end function certificate_index
    pure logical function certificate_matches(id, fingerprint)
        character(*), intent(in) :: id, fingerprint
        integer :: k
        k = certificate_index(id)
        certificate_matches = .false.
        if (k > 0) then
            certificate_matches = fingerprint == &
                certificate_fingerprint(k)
        end if
    end function certificate_matches
end module neort_generated_certificate_registry
