module neort_generated_certificate_registry
    implicit none
    character(*), parameter :: fortsym_revision = &
        'fortsym@545788453a204d58705f735b519c3863c2f734c8'
    character(*), parameter :: regenerate_command = &
        'cd tools/gc_symbolics && fo exec gen_full_fow_physics ../../src/generated'
    integer, parameter :: certificate_count = 20
    character(len=32), parameter :: certificate_id(certificate_count) = &
        [character(len=32) :: 'geometry', 'littlejohn', 'eq13_cdot', 'boundary_limits', &
        'root_enclosures', 'interpolation', 'profile_endpoints', &
        'refinement', 'harmonic_integrand', 'simple_root_force', &
        'eqdsk_cell_jet', 'eqdsk_profile_jet', 'eqdsk_cut_jet', &
        'eqdsk_cut_numerator_jet', 'eqdsk_cut_r_chart', &
        'eqdsk_cut_z_chart', 'eqdsk_cut_r_flux_chart', &
        'eqdsk_cut_mean_value', 'eqdsk_cut_axis_curvature', &
        'eqdsk_cut_axis_limit' ]
    character(len=64), parameter :: certificate_fingerprint(certificate_count) = &
        [character(len=64) :: 'neort-cert-v1:geometry:19:fortsym-5457884', &
        'neort-cert-v1:littlejohn:22:fortsym-5457884', &
        'neort-cert-v1:eq13_cdot:3:fortsym-5457884', &
        'neort-cert-v1:boundary_limits:13:fortsym-5457884', &
        'neort-cert-v1:root_enclosures:3:fortsym-5457884', &
        'neort-cert-v1:interpolation:9:fortsym-5457884', &
        'neort-cert-v1:profile_endpoints:8:fortsym-5457884', &
        'neort-cert-v1:refinement:14:fortsym-5457884', &
        'neort-cert-v1:harmonic_integrand:8:fortsym-5457884', &
        'neort-cert-v1:simple_root_force:3:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cell_jet:10:fortsym-5457884', &
        'neort-cert-v1:eqdsk_profile_jet:4:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_jet:7:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_numerator_jet:3:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_r_chart:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_z_chart:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_r_flux_chart:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_mean_value:1:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_axis_curvature:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_axis_limit:3:fortsym-5457884' ]
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
