program test_gc_eqdsk_flux_profile_map
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_eqdsk_flux_profile_map, only: &
        EQDSK_FLUX_MAP_ENDPOINT_MISMATCH, EQDSK_FLUX_MAP_NOT_MONOTONE, &
        EQDSK_FLUX_MAP_OUT_OF_RANGE, EQDSK_FLUX_MAP_SUCCESS, &
        EQDSK_FLUX_MAP_CERTIFICATE_ID, &
        eqdsk_flux_profile_map_t, initialize_eqdsk_flux_profile_map, &
        map_eqdsk_psihat_interval_to_s_tor, map_eqdsk_psihat_to_s_tor, &
        map_eqdsk_rho_tor_to_psihat, map_eqdsk_s_tor_to_psihat, &
        map_eqdsk_scaled_psi_to_s_tor, &
        validate_eqdsk_flux_profile_map
    use neort_gc_eqdsk_cut_endpoint_certificate, only: &
        EQDSK_ENDPOINT_CERTIFICATE_ID
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_t
    implicit none

    type(eqdsk_flux_profile_map_t) :: map
    real(dp), parameter :: s_nodes(3) = [0.0_dp, 0.25_dp, 1.0_dp]
    real(dp), parameter :: scaled_psi(3) = [0.0_dp, 0.8_dp, 8.0_dp]
    real(dp) :: s_tor, psihat, derivative
    type(gc_outward_interval_t) :: s_tor_interval
    integer :: status, segments_covered

    ! Independent nonuniform oracle: field_scale*psi_sep=8, so the profile
    ! nodes are psi_hat=[0,0.1,1].  The two exact slopes are 0.4 and 1.2.
    call initialize_eqdsk_flux_profile_map(s_nodes, scaled_psi, 2.0_dp, &
        4.0_dp, map, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'valid nonuniform flux map was rejected')
    call validate_eqdsk_flux_profile_map(map, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'fresh flux-map certificate failed validation')
    call require(EQDSK_FLUX_MAP_CERTIFICATE_ID /= &
        EQDSK_ENDPOINT_CERTIFICATE_ID, &
        'flux-map and endpoint certificate IDs collide')
    call require(map%certificate_id == EQDSK_FLUX_MAP_CERTIFICATE_ID, &
        'flux-map certificate used the wrong named ID')
    call require_close(map%minimum_dpsihat_dstor, 0.4_dp, &
        'certified minimum profile slope')

    call map_eqdsk_s_tor_to_psihat(map, 0.125_dp, psihat, derivative, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, 'forward map failed')
    call require_close(psihat, 0.05_dp, 'forward nonuniform map')
    call require_close(derivative, 0.4_dp, 'forward nonuniform derivative')

    call map_eqdsk_psihat_to_s_tor(map, 0.7_dp, s_tor, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, 'inverse map failed')
    call require_close(s_tor, 0.75_dp, 'inverse nonuniform map')
    call map_eqdsk_psihat_interval_to_s_tor(map, &
        gc_outward_interval(0.05_dp, 0.7_dp), s_tor_interval, &
        segments_covered, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'interval inverse map failed')
    call require(segments_covered == 2, &
        'interval inverse did not retain both profile segments')
    call require_contains(s_tor_interval, 0.125_dp, &
        'interval inverse lower endpoint')
    call require_contains(s_tor_interval, 0.75_dp, &
        'interval inverse upper endpoint')
    call map_eqdsk_psihat_interval_to_s_tor(map, &
        gc_outward_interval(0.1_dp, 0.1_dp), s_tor_interval, &
        segments_covered, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'profile-knot interval inverse failed')
    call require(segments_covered == 2, &
        'profile knot did not retain both adjacent segments')
    call require_contains(s_tor_interval, 0.25_dp, &
        'profile-knot interval inverse')
    call map_eqdsk_psihat_interval_to_s_tor(map, &
        gc_outward_interval(-epsilon(1.0_dp), 0.1_dp), s_tor_interval, &
        segments_covered, status)
    call require(status == EQDSK_FLUX_MAP_OUT_OF_RANGE, &
        'out-of-range interval was silently clipped')
    call map_eqdsk_scaled_psi_to_s_tor(map, 5.6_dp, s_tor, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'scaled physical-flux inverse failed')
    call require_close(s_tor, 0.75_dp, 'scaled physical-flux inverse')

    call map_eqdsk_rho_tor_to_psihat(map, 0.25_dp, s_tor, psihat, &
        derivative, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, 'rho_tor map failed')
    call require_close(s_tor, 0.0625_dp, 'rho_tor squared oracle')
    call require_close(psihat, 0.025_dp, 'rho_tor to psi_hat oracle')
    call require_close(derivative, 0.2_dp, 'rho_tor derivative oracle')

    map%scaled_psi(2) = map%scaled_psi(2)+0.1_dp
    call validate_eqdsk_flux_profile_map(map, status)
    call require(status /= EQDSK_FLUX_MAP_SUCCESS, &
        'mutated profile retained its original certificate')
    map%scaled_psi(2) = scaled_psi(2)
    call validate_eqdsk_flux_profile_map(map, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'restored profile did not recover its certificate')

    call initialize_eqdsk_flux_profile_map(s_nodes, &
        [-epsilon(1.0_dp)*8.0_dp, 0.8_dp, &
        (1.0_dp-epsilon(1.0_dp))*8.0_dp], 2.0_dp, 4.0_dp, map, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'roundoff-sized endpoint residuals lost endpoint certification')
    call require(map%endpoints_certified, &
        'accepted endpoint residuals were not explicitly certified')

    call initialize_eqdsk_flux_profile_map(s_nodes, &
        [0.0_dp, 1.0_dp, 0.5_dp], 1.0_dp, 0.5_dp, map, status)
    call require(status == EQDSK_FLUX_MAP_NOT_MONOTONE, &
        'decreasing profile segment retained a monotonicity certificate')
    call initialize_eqdsk_flux_profile_map(s_nodes, scaled_psi, 2.0_dp, &
        5.0_dp, map, status)
    call require(status == EQDSK_FLUX_MAP_ENDPOINT_MISMATCH, &
        'incorrect separatrix normalization retained endpoint certification')

    write (*, '(a)') 'test_gc_eqdsk_flux_profile_map OK'

contains

    subroutine require_close(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label

        if (abs(actual-expected) > 5.0e-13_dp*max(1.0_dp, abs(expected))) then
            write (*, '(a,2(1x,es24.16))') trim(label), actual, expected
            error stop 'flux-profile map oracle failed'
        end if
    end subroutine require_close

    subroutine require_contains(interval, expected, label)
        type(gc_outward_interval_t), intent(in) :: interval
        real(dp), intent(in) :: expected
        character(*), intent(in) :: label

        call require(interval%lo <= expected .and. interval%hi >= expected, &
            trim(label)//' excluded the independent oracle')
    end subroutine require_contains

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_flux_profile_map
