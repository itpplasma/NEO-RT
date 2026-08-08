program test_gc_generated_certificates
    use neort_generated_certificate_registry, only: certificate_matches, &
        fortsym_revision, regenerate_command
    use util_for_test, only: pass_test
    implicit none

    if (.not. certificate_matches('geometry', &
            'neort-cert-v1:geometry:19:fortsym-5457884')) error stop &
        'geometry certificate registry mismatch'
    if (.not. certificate_matches('littlejohn', &
            'neort-cert-v1:littlejohn:22:fortsym-5457884')) error stop &
        'Littlejohn certificate registry mismatch'
    if (.not. certificate_matches('eq13_cdot', &
            'neort-cert-v1:eq13_cdot:3:fortsym-5457884')) error stop &
        'Cdot certificate registry mismatch'
    if (.not. certificate_matches('harmonic_integrand', &
            'neort-cert-v1:harmonic_integrand:8:fortsym-5457884')) error stop &
        'harmonic-integrand certificate registry mismatch'
    if (.not. certificate_matches('simple_root_force', &
            'neort-cert-v1:simple_root_force:3:fortsym-5457884')) error stop &
        'simple-root-force certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cell_jet', &
            'neort-cert-v1:eqdsk_cell_jet:10:fortsym-5457884')) error stop &
        'EQDSK cell-jet certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cell_fourth_jet', &
            'neort-cert-v1:eqdsk_cell_fourth_jet:5:fortsym-5457884')) &
        error stop 'EQDSK fourth-jet certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_profile_jet', &
            'neort-cert-v1:eqdsk_profile_jet:4:fortsym-5457884')) error stop &
        'EQDSK profile-jet certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cut_jet', &
            'neort-cert-v1:eqdsk_cut_jet:7:fortsym-5457884')) error stop &
        'EQDSK cut-jet certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cut_numerator_jet', &
            'neort-cert-v1:eqdsk_cut_numerator_jet:3:fortsym-5457884')) &
        error stop 'EQDSK cut-numerator certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cut_r_chart', &
            'neort-cert-v1:eqdsk_cut_r_chart:2:fortsym-5457884')) error stop &
        'EQDSK R-chart certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cut_z_chart', &
            'neort-cert-v1:eqdsk_cut_z_chart:2:fortsym-5457884')) error stop &
        'EQDSK Z-chart certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cut_numerator_hessian', &
            'neort-cert-v1:eqdsk_cut_numerator_hessian:3:fortsym-5457884')) &
        error stop 'EQDSK numerator-Hessian certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cut_r_flux_curvature', &
            'neort-cert-v1:eqdsk_cut_r_flux_curvature:4:fortsym-5457884')) &
        error stop 'EQDSK flux-curvature certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cut_axis_curvature', &
            'neort-cert-v1:eqdsk_cut_axis_curvature:2:fortsym-5457884')) &
        error stop 'EQDSK axis-curvature certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cut_axis_limit', &
            'neort-cert-v1:eqdsk_cut_axis_limit:3:fortsym-5457884')) &
        error stop 'EQDSK axis-limit certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_rho_tor_map', &
            'neort-cert-v1:eqdsk_rho_tor_map:2:fortsym-5457884')) &
        error stop 'EQDSK rho-tor-map certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_flux_profile_rho_chain', &
            'neort-cert-v1:eqdsk_flux_profile_rho_chain:1:fortsym-5457884')) &
        error stop 'EQDSK flux-profile rho-chain certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cut_flux_coordinate', &
            'neort-cert-v1:eqdsk_cut_flux_coordinate:3:fortsym-5457884')) &
        error stop 'EQDSK flux-coordinate certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_cut_axis_rho_limit', &
            'neort-cert-v1:eqdsk_cut_axis_rho_limit:3:fortsym-5457884')) &
        error stop 'EQDSK axis-rho-limit certificate registry mismatch'
    if (.not. certificate_matches('eqdsk_flux_profile_segment', &
            'neort-cert-v1:eqdsk_flux_profile_segment:3:fortsym-5457884')) &
        error stop 'EQDSK flux-profile certificate registry mismatch'
    if (index(fortsym_revision, '545788453a204d58705f735b519c3863c2f734c8') == 0) &
        error stop 'Fortsym revision provenance missing'
    if (index(regenerate_command, 'gen_full_fow_physics') == 0) &
        error stop 'regeneration provenance missing'

    write (*, '(A)') 'test_gc_generated_certificates OK'
    call pass_test
end program test_gc_generated_certificates
