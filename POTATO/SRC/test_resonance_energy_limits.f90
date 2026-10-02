program test_resonance_energy_limits
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan
    use resonance_energy_limits_mod, only: resonance_energy_limits
    implicit none
    real(dp) :: low, high, nan
    integer :: first
    logical :: valid

    ! A gauge offset changes total energies while leaving their interval intact.
    call resonance_energy_limits(-10.0_dp, -8.0_dp, 2.0_dp, 6.0_dp, &
        0.0_dp, low, high, first, valid)
    if (.not. valid) error stop 'Legacy energy domain rejected'
    if (low /= -10.0_dp .or. high /= -2.0_dp .or. first /= 2) &
        error stop 'Legacy energy grid changed'
    ! Every positive-cutoff midpoint must exceed the maximum physical potential.
    call resonance_energy_limits(-10.0_dp, -8.0_dp, 2.0_dp, 6.0_dp, &
        0.5_dp, low, high, first, valid)
    if (.not. valid) error stop 'Positive physical kinetic cutoff rejected'
    if (low /= -7.0_dp .or. high /= -2.0_dp .or. first /= 1) &
        error stop 'Physical energy cutoff mismatch'
    call resonance_energy_limits(0.0_dp, 1.0_dp, 2.0_dp, 6.0_dp, &
        3.0_dp, low, high, first, valid)
    if (valid) error stop 'Empty energy interval accepted'
    call resonance_energy_limits(0.0_dp, 1.0_dp, 2.0_dp, 6.0_dp, &
        4.0_dp, low, high, first, valid)
    if (valid) error stop 'Negative energy interval accepted'
    call resonance_energy_limits(1.0_dp, 0.0_dp, 2.0_dp, 6.0_dp, &
        0.0_dp, low, high, first, valid)
    if (valid) error stop 'Reversed potential extrema accepted'
    call resonance_energy_limits(-huge(1.0_dp), 0.0_dp, 2.0_dp, &
        huge(1.0_dp)/2.0_dp, 0.0_dp, low, high, first, valid)
    if (valid) error stop 'Overflowing finite-endpoint energy interval accepted'
    nan = ieee_value(0.0_dp, ieee_quiet_nan)
    call resonance_energy_limits(0.0_dp, nan, 2.0_dp, 6.0_dp, &
        0.0_dp, low, high, first, valid)
    if (valid) error stop 'Nonfinite energy interval accepted'
    print *, 'Physical energy-limit oracle PASS'
end program test_resonance_energy_limits
