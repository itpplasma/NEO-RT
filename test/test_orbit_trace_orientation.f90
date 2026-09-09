program test_orbit_trace_orientation
    use iso_fortran_env, only: dp => real64
    use diag_orbit_trace, only: native_torque_mode_factors, physical_orientation, &
        toroidal_velocity_from_components
    implicit none

    logical :: valid
    real(dp) :: covector, drive_mode, product

    if (physical_orientation(2.0_dp, 3.0_dp) /= 1) error stop "co-positive orientation"
    if (physical_orientation(-2.0_dp, 3.0_dp) /= -1) error stop "co-negative orientation"
    if (physical_orientation(2.0_dp, -3.0_dp) /= -1) error stop "reversed-chart orientation"
    if (physical_orientation(-2.0_dp, -3.0_dp) /= 1) error stop "double-reversed orientation"
    if (physical_orientation(0.0_dp, 3.0_dp) /= 0) error stop "zero state orientation"
    if (physical_orientation(2.0_dp, 0.0_dp) /= 0) error stop "zero chart-factor orientation"

    if (toroidal_velocity_from_components(2.0_dp, 3.0_dp, 5.0_dp, 6.0_dp) &
        /= 17.0_dp) error stop "toroidal velocity decomposition"
    if (toroidal_velocity_from_components(-2.0_dp, 3.0_dp, -5.0_dp, 6.0_dp) &
        /= -5.0_dp) error stop "signed toroidal velocity decomposition"

    call native_torque_mode_factors(3, -2.0_dp, -4.0_dp, -1.0_dp, &
        covector, drive_mode, product, valid)
    if (.not. valid) error stop "valid native torque factors rejected"
    if (covector /= -3.0_dp .or. drive_mode /= 3.0_dp .or. product /= -9.0_dp) &
        error stop "native torque factor signs"
    call native_torque_mode_factors(-3, -2.0_dp, -4.0_dp, 1.0_dp, &
        covector, drive_mode, product, valid)
    if (.not. valid) error stop "valid reversed-mode torque factors rejected"
    if (covector /= -3.0_dp .or. drive_mode /= -3.0_dp .or. product /= 9.0_dp) &
        error stop "reversed-mode torque factors"
    call native_torque_mode_factors(3, 0.0_dp, -4.0_dp, -1.0_dp, &
        covector, drive_mode, product, valid)
    if (valid) error stop "zero flux derivative admitted as a torque orientation"

    print *, "orbit trace orientation map: PASS"
end program test_orbit_trace_orientation
