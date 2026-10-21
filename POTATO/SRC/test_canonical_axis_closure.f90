program test_canonical_axis_closure
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_value, ieee_quiet_nan, &
        ieee_positive_inf
    use potato_input_mod, only: canonical_axis_constant_extension, &
        zero_canonical_axis_drive, bounce_period_admitted
    implicit none
    real(dp), parameter :: delta = 1.0e-6_dp, axis = 0.7_dp
    real(dp) :: edge, psi, derivative, value_left, value_right
    integer :: flux_sign, charge_sign

    if (bounce_period_admitted(ieee_value(0.0_dp,ieee_quiet_nan))) &
        error stop 'NaN bounce period admitted'
    if (bounce_period_admitted(ieee_value(0.0_dp,ieee_positive_inf))) &
        error stop 'infinite bounce period admitted'
    if (bounce_period_admitted(0.0_dp)) error stop 'zero bounce period admitted'
    if (bounce_period_admitted(-1.0_dp)) error stop 'negative period admitted'
    if (.not. bounce_period_admitted(1.0_dp)) error stop 'valid period rejected'
    if (canonical_axis_constant_extension) error stop 'closure must default off'
    do flux_sign = -1, 1, 2
        edge = axis + real(flux_sign,dp)*3.5_dp
        psi = axis - 0.1_dp*(edge-axis)
        if (zero_canonical_axis_drive(psi,axis,edge,.true.)) &
            error stop 'strict default changed'
    end do
    canonical_axis_constant_extension = .true.
    do flux_sign = -1, 1, 2
        edge = axis + real(flux_sign,dp)*3.5_dp
        do charge_sign = -1, 1, 2
            psi = axis - 0.1_dp*(edge-axis)
            if (.not. zero_canonical_axis_drive(psi,axis,edge,.true.)) &
                error stop 'constant-axis drive was not suppressed'
            if (zero_canonical_axis_drive(psi,axis,edge,.false.)) &
                error stop 'closure applied before orbit admission'
            derivative = (log_f0(psi+delta,axis,edge,charge_sign) &
                        - log_f0(psi-delta,axis,edge,charge_sign))/(2*delta)
            if (derivative /= 0.0_dp) error stop 'extended F0 derivative not zero'
            psi = axis + 0.2_dp*(edge-axis)
            if (zero_canonical_axis_drive(psi,axis,edge,.true.)) &
                error stop 'in-domain canonical drive changed'
            derivative = (log_f0(psi+delta,axis,edge,charge_sign) &
                        - log_f0(psi-delta,axis,edge,charge_sign))/(2*delta)
            if (abs(derivative) < 1.0e-3_dp) &
                error stop 'manufactured in-domain drive is degenerate'
            value_left = log_f0(axis-delta,axis,edge,charge_sign)
            value_right = log_f0(axis+delta,axis,edge,charge_sign)
            if (abs(value_right-value_left) > 2*delta) &
                error stop 'extended F0 is discontinuous at the axis'
        end do
        if (zero_canonical_axis_drive(axis,axis,edge,.true.)) &
            error stop 'axis endpoint must retain strict evaluation'
    end do
contains
    pure real(dp) function log_f0(psi_star, psi_axis, psi_edge, charge) result(f)
        real(dp), intent(in) :: psi_star, psi_axis, psi_edge
        integer, intent(in) :: charge
        real(dp) :: s, density, temperature, potential

        ! Independent manufactured F0 at fixed H=1.7, with all three profiles
        ! varying in-domain and continuous constant extensions below the axis.
        s = max(0.0_dp,(psi_star-psi_axis)/(psi_edge-psi_axis))
        density = exp(0.4_dp*s)
        temperature = 2.0_dp+0.7_dp*s
        potential = 0.3_dp*s+0.1_dp*s*s
        f = log(density)-1.5_dp*log(temperature) &
            -(1.7_dp-real(charge,dp)*potential)/temperature
    end function log_f0
end program test_canonical_axis_closure
