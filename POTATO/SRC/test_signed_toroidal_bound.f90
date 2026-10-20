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

    ! |m| <= 3 and |n| = 3 give 2*pi*(1 + 1/3), regardless of sign.
    if (abs(positive_bound - 8.d0*acos(-1.d0)/3.d0) > 1.d-13) &
        error stop "incorrect search extent"
    if (positive_bound <= 0.d0) error stop "positive n produced nonpositive bound"
    if (negative_bound <= 0.d0) error stop "negative n produced nonpositive bound"
    if (abs(positive_bound-negative_bound) > 1d-14) &
        error stop "search bound depends on toroidal-mode sign"

    if (canonical_flux_outside_lcfs(0.5d0, 0.d0, 1.d0)) &
        error stop "inside point rejected for increasing flux"
    if (.not. canonical_flux_outside_lcfs(1.1d0, 0.d0, 1.d0)) &
        error stop "outside point accepted for increasing flux"
    if (canonical_flux_outside_lcfs(-0.5d0, 0.d0, -1.d0)) &
        error stop "inside point rejected for decreasing flux"
    if (.not. canonical_flux_outside_lcfs(-1.1d0, 0.d0, -1.d0)) &
        error stop "outside point accepted for decreasing flux"
    ! A native nonzero axis distinguishes this from the old psi_star/psi_edge
    ! comparison: the interior point exceeds the edge in a decreasing gauge.
    if (canonical_flux_outside_lcfs(5.5d0, 6.d0, 5.d0)) &
        error stop "native decreasing-flux inside point rejected"
    if (.not. canonical_flux_outside_lcfs(4.9d0, 6.d0, 5.d0)) &
        error stop "native decreasing-flux outside point accepted"
    if (canonical_flux_outside_lcfs(-5.5d0, -6.d0, -5.d0)) &
        error stop "native increasing-flux inside point rejected"
    if (.not. canonical_flux_outside_lcfs(-4.9d0, -6.d0, -5.d0)) &
        error stop "native increasing-flux outside point accepted"
    if (canonical_flux_outside_lcfs(5.d0, 6.d0, 5.d0)) &
        error stop "LCFS point rejected for decreasing flux"
    if (canonical_flux_outside_lcfs(-5.d0, -6.d0, -5.d0)) &
        error stop "LCFS point rejected for increasing flux"
end program test_signed_toroidal_bound
