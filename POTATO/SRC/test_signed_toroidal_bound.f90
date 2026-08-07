program test_signed_toroidal_bound
    use resonance_mode_bounds_mod, only: resonant_delphi_bound, &
        canonical_flux_outside_lcfs, set_resonance_harmonic_guard, &
        resonance_harmonic_values, resonance_no_root_for_any, &
        harmonic_guard_success,harmonic_guard_invalid_set
    implicit none

    integer, parameter :: m_modes(3) = [-2, 0, 3]
    integer, parameter :: n_positive(3) = 3
    integer, parameter :: n_negative(3) = -3
    integer, parameter :: n_zero(3) = [3, 0, -3]
    double precision :: positive_bound, negative_bound,guard_bound
    double precision :: target,extent,resonance_g,no_root_margin
    logical :: no_root
    integer :: ierr

    positive_bound = resonant_delphi_bound(m_modes, n_positive)
    negative_bound = resonant_delphi_bound(m_modes, n_negative)

    if (positive_bound <= 0.d0) error stop "positive n produced nonpositive bound"
    if (negative_bound <= 0.d0) error stop "negative n produced nonpositive bound"
    if (positive_bound /= negative_bound) &
        error stop "search bound depends on toroidal-mode sign"

    call set_resonance_harmonic_guard(m_modes,n_negative,guard_bound,ierr)
    if(ierr /= harmonic_guard_success) error stop "harmonic guard rejected valid set"
    if(guard_bound /= negative_bound) error stop "finite-set envelope changed"
    call set_resonance_harmonic_guard(m_modes,n_zero,guard_bound,ierr)
    if(ierr /= harmonic_guard_invalid_set) error stop "zero toroidal mode was accepted"
    call set_resonance_harmonic_guard(m_modes,n_negative,guard_bound,ierr)
    if(ierr /= harmonic_guard_success) error stop "valid harmonic guard was not restored"

    ! The signed target for (m,n)=(-2,-3) is -4*pi/3.  The extent is its
    ! magnitude only, and equality is retained as a possible root.
    call resonance_harmonic_values(-2,-3,0.d0,target,extent,resonance_g, &
                                   no_root_margin,ierr)
    if(ierr /= harmonic_guard_success) error stop "harmonic kernel failed"
    if(abs(target + 4.d0*acos(-1.d0)/3.d0) > 1.d-13) &
        error stop "signed harmonic target was erased"
    if(abs(extent + target) > 1.d-13) error stop "harmonic extent is not |target|"
    if(resonance_g /= -target) error stop "signed g was not preserved"

    call resonance_no_root_for_any(10.d0*acos(-1.d0),no_root,ierr)
    if(ierr /= harmonic_guard_success .or. .not.no_root) &
        error stop "strict per-harmonic no-root guard failed"
    call resonance_no_root_for_any(-4.d0*acos(-1.d0)/3.d0,no_root,ierr)
    if(ierr /= harmonic_guard_success .or. no_root) &
        error stop "root equality was incorrectly clipped"

    if (canonical_flux_outside_lcfs(0.5d0, 0.d0, 1.d0)) &
        error stop "inside point rejected for increasing flux"
    if (.not. canonical_flux_outside_lcfs(1.1d0, 0.d0, 1.d0)) &
        error stop "outside point accepted for increasing flux"
    if (canonical_flux_outside_lcfs(-0.5d0, 0.d0, -1.d0)) &
        error stop "inside point rejected for decreasing flux"
    if (.not. canonical_flux_outside_lcfs(-1.1d0, 0.d0, -1.d0)) &
        error stop "outside point accepted for decreasing flux"
end program test_signed_toroidal_bound
