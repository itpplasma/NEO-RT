program test_signed_toroidal_bound
    use resonance_mode_bounds_mod, only: resonant_delphi_bound, &
        canonical_flux_outside_lcfs
    implicit none

    integer, parameter :: m_modes(7) = [-3, -2, -1, 0, 1, 2, 3]
    integer, parameter :: n_positive(7) = 3
    integer, parameter :: n_negative(7) = -3
    double precision :: positive_bound, negative_bound

    positive_bound = resonant_delphi_bound(m_modes, n_positive)
    negative_bound = resonant_delphi_bound(m_modes, n_negative)

    if (positive_bound <= 0.d0) error stop "positive n produced nonpositive bound"
    if (negative_bound <= 0.d0) error stop "negative n produced nonpositive bound"
    if (positive_bound /= negative_bound) &
        error stop "search bound depends on toroidal-mode sign"

    if (canonical_flux_outside_lcfs(0.5d0, 0.d0, 1.d0)) &
        error stop "inside point rejected for increasing flux"
    if (.not. canonical_flux_outside_lcfs(1.1d0, 0.d0, 1.d0)) &
        error stop "outside point accepted for increasing flux"
    if (canonical_flux_outside_lcfs(-0.5d0, 0.d0, -1.d0)) &
        error stop "inside point rejected for decreasing flux"
    if (.not. canonical_flux_outside_lcfs(-1.1d0, 0.d0, -1.d0)) &
        error stop "outside point accepted for decreasing flux"
end program test_signed_toroidal_bound
