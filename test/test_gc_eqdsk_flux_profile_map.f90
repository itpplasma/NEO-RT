program test_gc_eqdsk_flux_profile_map
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_eqdsk_flux_profile_map, only: &
        EQDSK_FLUX_MAP_ENDPOINT_MISMATCH, EQDSK_FLUX_MAP_NOT_MONOTONE, &
        EQDSK_FLUX_MAP_SUCCESS, eqdsk_flux_profile_map_t, &
        initialize_eqdsk_flux_profile_map, map_eqdsk_psihat_to_s_tor, &
        map_eqdsk_rho_tor_to_psihat, map_eqdsk_s_tor_to_psihat, &
        validate_eqdsk_flux_profile_map
    implicit none

    type(eqdsk_flux_profile_map_t) :: map
    real(dp), parameter :: s_nodes(3) = [0.0_dp, 0.25_dp, 1.0_dp]
    real(dp), parameter :: scaled_psi(3) = [0.0_dp, 0.8_dp, 8.0_dp]
    real(dp) :: s_tor, psihat, derivative
    integer :: status

    ! Independent nonuniform oracle: field_scale*psi_sep=8, so the profile
    ! nodes are psi_hat=[0,0.1,1].  The two exact slopes are 0.4 and 1.2.
    call initialize_eqdsk_flux_profile_map(s_nodes, scaled_psi, 2.0_dp, &
        4.0_dp, map, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'valid nonuniform flux map was rejected')
    call validate_eqdsk_flux_profile_map(map, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, &
        'fresh flux-map certificate failed validation')
    call require_close(map%minimum_dpsihat_dstor, 0.4_dp, &
        'certified minimum profile slope')

    call map_eqdsk_s_tor_to_psihat(map, 0.125_dp, psihat, derivative, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, 'forward map failed')
    call require_close(psihat, 0.05_dp, 'forward nonuniform map')
    call require_close(derivative, 0.4_dp, 'forward nonuniform derivative')

    call map_eqdsk_psihat_to_s_tor(map, 0.7_dp, s_tor, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, 'inverse map failed')
    call require_close(s_tor, 0.75_dp, 'inverse nonuniform map')

    call map_eqdsk_rho_tor_to_psihat(map, 0.25_dp, s_tor, psihat, &
        derivative, status)
    call require(status == EQDSK_FLUX_MAP_SUCCESS, 'rho_tor map failed')
    call require_close(s_tor, 0.0625_dp, 'rho_tor squared oracle')
    call require_close(psihat, 0.025_dp, 'rho_tor to psi_hat oracle')
    call require_close(derivative, 0.2_dp, 'rho_tor derivative oracle')

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

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_flux_profile_map
