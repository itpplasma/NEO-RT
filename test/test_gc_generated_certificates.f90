program test_gc_generated_certificates
    use neort_generated_certificate_registry, only: certificate_matches, &
        fortsym_revision, regenerate_command
    use util_for_test, only: pass_test
    implicit none

    if (.not. certificate_matches('geometry', &
            'neort-cert-v1:geometry:19:fortsym-a2b9bb3')) error stop &
        'geometry certificate registry mismatch'
    if (.not. certificate_matches('littlejohn', &
            'neort-cert-v1:littlejohn:22:fortsym-a2b9bb3')) error stop &
        'Littlejohn certificate registry mismatch'
    if (.not. certificate_matches('eq13_cdot', &
            'neort-cert-v1:eq13_cdot:3:fortsym-a2b9bb3')) error stop &
        'Cdot certificate registry mismatch'
    if (.not. certificate_matches('harmonic_integrand', &
            'neort-cert-v1:harmonic_integrand:8:fortsym-a2b9bb3')) error stop &
        'harmonic-integrand certificate registry mismatch'
    if (.not. certificate_matches('simple_root_force', &
            'neort-cert-v1:simple_root_force:3:fortsym-a2b9bb3')) error stop &
        'simple-root-force certificate registry mismatch'
    if (index(fortsym_revision, 'a2b9bb353816b03445282fc71d0584b34a787549') == 0) &
        error stop 'Fortsym revision provenance missing'
    if (index(regenerate_command, 'gen_full_fow_physics') == 0) &
        error stop 'regeneration provenance missing'

    write (*, '(A)') 'test_gc_generated_certificates OK'
    call pass_test
end program test_gc_generated_certificates
