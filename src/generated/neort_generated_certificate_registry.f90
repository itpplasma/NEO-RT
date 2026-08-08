module neort_generated_certificate_registry
    implicit none
    character(*), parameter :: fortsym_revision = &
        'fortsym@545788453a204d58705f735b519c3863c2f734c8'
    character(*), parameter :: regenerate_command = &
        'cd tools/gc_symbolics && fo exec gen_full_fow_physics ../../src/generated'
    integer, parameter :: certificate_count = 55
    character(len=32), parameter :: certificate_id(certificate_count) = &
        [character(len=32) :: 'geometry', 'littlejohn', 'eq13_cdot', 'boundary_limits', &
        'root_enclosures', 'interpolation', 'profile_endpoints', &
        'refinement', 'harmonic_integrand', 'simple_root_force', &
        'eqdsk_cell_jet', 'eqdsk_cell_fourth_jet', &
        'eqdsk_profile_jet', 'eqdsk_cut_jet', &
        'eqdsk_cut_numerator_jet', 'eqdsk_cut_r_chart', &
        'eqdsk_cut_z_chart', 'eqdsk_cut_r_flux_chart', &
        'eqdsk_cut_numerator_hessian', &
        'eqdsk_cut_r_flux_curvature', 'eqdsk_cut_mean_value', &
        'eqdsk_cut_axis_curvature', &
        'eqdsk_cut_axis_limit', &
        'eqdsk_rho_tor_map', &
        'eqdsk_s_tor_to_rho', &
        'eqdsk_flux_profile_rho_chain', &
        'eqdsk_cut_flux_coordinate', &
        'eqdsk_flux_coordinate_interval', &
        'eqdsk_cut_axis_rho_limit', &
        'eqdsk_allowed_rho_chain', &
        'eqdsk_allowed_axis_rho_chain', &
        'eqdsk_cut_endpoint_system', &
        'eqdsk_cut_endpoint_newton', &
        'eqdsk_cut_endpoint_krawczyk', &
        'eqdsk_axis_stationarity_system', &
        'eqdsk_axis_stationarity_newton', &
        'eqdsk_axis_stationarity_krawczyk', &
        'eqdsk_flux_profile_segment', &
        'eqdsk_scaled_flux_normalization', &
        'eqdsk_physical_flux_norm', &
        'gauss_interval_map', &
        'profile_potential_map', &
        'eqdsk_allowed_energy', &
        'eqdsk_canonical_cut', &
        'eqdsk_turning_chart', &
        'eqdsk_physical_flux_map', &
        'full_fow_canonical_numerator', &
        'full_fow_canonical_turning', &
        'full_fow_canonical_symmetry', &
        'full_fow_cycle_frequency', 'full_fow_resonance_scalar', &
        'full_fow_cycle_average', 'full_fow_mode_mapping', &
        'full_fow_torque_assembly', 'full_fow_cut_topology' ]
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
        'neort-cert-v1:eqdsk_cell_fourth_jet:5:fortsym-5457884', &
        'neort-cert-v1:eqdsk_profile_jet:4:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_jet:7:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_numerator_jet:3:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_r_chart:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_z_chart:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_r_flux_chart:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_numerator_hessian:3:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_r_flux_curvature:4:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_mean_value:1:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_axis_curvature:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_axis_limit:3:fortsym-5457884', &
        'neort-cert-v1:eqdsk_rho_tor_map:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_s_tor_to_rho:1:fortsym-5457884', &
        'neort-cert-v1:eqdsk_flux_profile_rho_chain:1:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_flux_coordinate:3:fortsym-5457884', &
        'neort-cert-v1:eqdsk_flux_coordinate_interval:3:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_axis_rho_limit:3:fortsym-5457884', &
        'neort-cert-v1:eqdsk_allowed_rho_chain:7:fortsym-5457884', &
        'neort-cert-v1:eqdsk_allowed_axis_rho_chain:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_endpoint_system:6:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_endpoint_newton:7:fortsym-5457884', &
        'neort-cert-v1:eqdsk_cut_endpoint_krawczyk:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_axis_stationarity_system:6:fortsym-5457884', &
        'neort-cert-v1:eqdsk_axis_stationarity_newton:7:fortsym-5457884', &
        'neort-cert-v1:eqdsk_axis_stationarity_krawczyk:2:fortsym-5457884', &
        'neort-cert-v1:eqdsk_flux_profile_segment:3:fortsym-5457884', &
        'neort-cert-v1:eqdsk_scaled_flux_normalization:1:fortsym-5457884', &
        'neort-cert-v1:eqdsk_physical_flux_norm:1:fortsym-5457884', &
        'neort-cert-v1:gauss_interval_map:2:fortsym-5457884', &
        'neort-cert-v1:profile_potential_map:3:fortsym-5457884', &
        'neort-cert-v1:eqdsk_allowed_energy:17:fortsym-5457884', &
        'neort-cert-v1:eqdsk_canonical_cut:4:fortsym-5457884', &
        'neort-cert-v1:eqdsk_turning_chart:5:fortsym-5457884', &
        'neort-cert-v1:eqdsk_physical_flux_map:2:fortsym-5457884', &
        'neort-cert-v1:full_fow_canonical_numerator:3:fortsym-5457884', &
        'neort-cert-v1:full_fow_canonical_turning:2:fortsym-5457884', &
        'neort-cert-v1:full_fow_canonical_symmetry:21:fortsym-5457884', &
        'neort-cert-v1:full_fow_cycle_frequency:6:fortsym-5457884', &
        'neort-cert-v1:full_fow_resonance_scalar:2:fortsym-5457884', &
        'neort-cert-v1:full_fow_cycle_average:4:fortsym-5457884', &
        'neort-cert-v1:full_fow_mode_mapping:4:fortsym-5457884', &
        'neort-cert-v1:full_fow_torque_assembly:1:fortsym-5457884', &
        'neort-cert-v1:full_fow_cut_topology:6:fortsym-5457884' ]
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
