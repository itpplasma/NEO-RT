program test_gc_circular_flux_continuation
    !! Independent geometric and quadrature oracle for the circular fixture.
    !!
    !! In particular, this test does not reproduce the Fortsym atanh primitive.
    !! It integrates dpsi/dr and the field-line pitch directly on the circular
    !! contour.  The pitch check fails for the old missing-sqrt normalization.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_circular_flux_continuation_limit_symbolic, only: &
        evaluate_neort_circular_flux_continuation_limit
    use neort_circular_flux_continuation_symbolic, only: &
        evaluate_neort_circular_flux_continuation_raw => &
        evaluate_neort_circular_flux_continuation
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: major_radius = 1.60_dp
    real(dp), parameter :: edge_radius = 0.50_dp
    real(dp), parameter :: psi_edge = 0.4_dp
    real(dp), parameter :: toroidal_flux = 3.2_dp
    real(dp), parameter :: q_axis = 1.5_dp
    real(dp), parameter :: delta_q = 2.5_dp
    real(dp), parameter :: mu0 = 4.0e-7_dp*acos(-1.0_dp)
    real(dp), parameter :: h = 1.0e-5_dp
    integer, parameter :: quadrature_steps = 20000
    real(dp) :: psi, dpsi, q, psi_left, psi_right, dpsi_left
    real(dp) :: dpsi_right, psi_limit, dpsi_limit, q_limit
    real(dp) :: psi_limit_left, psi_limit_right
    real(dp) :: radius_value, oracle_value, psi_axis_oracle
    real(dp) :: q_pitch, toroidal_flux_value, toroidal_flux_edge
    real(dp) :: rho_tor, dvolume_ds, current_oracle
    real(dp) :: psi_tor_emitted, s_tor_emitted, rho_tor_emitted
    real(dp) :: dvolume_ds_emitted, dvolume_ds_cgs_emitted, current_emitted
    integer :: k

    call evaluate_neort_circular_flux_continuation(edge_radius, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, major_radius, psi, dpsi, q)
    call require_close('LCFS value', psi, psi_edge, 2.0e-13_dp)
    call evaluate_neort_circular_flux_continuation(edge_radius-h, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, major_radius, psi_left, &
        dpsi_left, q)
    call evaluate_neort_circular_flux_continuation(edge_radius+h, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, major_radius, psi_right, &
        dpsi_right, q)
    call require_close('left first derivative', (psi-psi_left)/h, dpsi, &
        2.0e-5_dp)
    call require_close('right first derivative', (psi_right-psi)/h, dpsi, &
        2.0e-5_dp)
    call require_close('first derivative continuity', dpsi_left, dpsi_right, &
        4.0e-5_dp)

    psi_axis_oracle = numerical_flux(0.0_dp)
    call evaluate_neort_circular_flux_continuation(0.0_dp, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, major_radius, psi, dpsi, q)
    call require_close('axis value from quadrature', psi, psi_axis_oracle, &
        2.0e-9_dp)

    do k = 1, 8
        radius_value = 0.05_dp + 1.20_dp*real(k, dp)/8.0_dp
        call evaluate_neort_circular_flux_continuation( &
            radius_value, edge_radius, psi_edge, toroidal_flux, q_axis, &
            delta_q, major_radius, psi, dpsi, q)
        if (dpsi <= 0.0_dp) error stop 'exterior flux is not monotone'
        oracle_value = numerical_flux(radius_value)
        call require_close('numerical flux oracle', psi, oracle_value, &
            2.0e-9_dp)
        call require_close('quadrature field-line pitch', &
            field_line_pitch(radius_value, dpsi), target_q(radius_value), &
            2.0e-9_dp)
        call require_close('generated q equals pitch', q, &
            field_line_pitch(radius_value, dpsi), 2.0e-9_dp)
    end do

    call evaluate_neort_circular_flux_continuation(edge_radius, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, major_radius, psi, dpsi, q)
    call require_close('defining LCFS slope', dpsi, &
        toroidal_flux*edge_radius/(target_q(edge_radius)* &
        sqrt(major_radius**2-edge_radius**2)), 2.0e-13_dp)

    call evaluate_neort_circular_flux_continuation_limit(edge_radius, &
        edge_radius, psi_edge, toroidal_flux, q_axis, major_radius, &
        psi_limit, dpsi_limit, &
        q_limit)
    call require_close('zero-shear LCFS value', psi_limit, psi_edge, &
        2.0e-13_dp)
    call evaluate_neort_circular_flux_continuation_limit(edge_radius-h, &
        edge_radius, psi_edge, toroidal_flux, q_axis, major_radius, &
        psi_limit_left, &
        dpsi_left, q_limit)
    call evaluate_neort_circular_flux_continuation_limit(edge_radius+h, &
        edge_radius, psi_edge, toroidal_flux, q_axis, major_radius, &
        psi_limit_right, &
        dpsi_right, q_limit)
    call require_close('zero-shear first derivative', &
        (psi_limit_right-psi_limit_left)/(2.0_dp*h), dpsi_limit, 2.0e-10_dp)
    call require_close('zero-shear q', q_limit, q_axis, 2.0e-13_dp)
    call require_close('zero-shear finite-aspect derivative', dpsi_limit, &
        toroidal_flux*edge_radius/(q_axis* &
        sqrt(major_radius**2-edge_radius**2)), 2.0e-13_dp)

    call evaluate_neort_circular_flux_continuation(edge_radius, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, major_radius, psi, dpsi, q)
    toroidal_flux_edge = direct_toroidal_flux(edge_radius)
    toroidal_flux_value = integrated_toroidal_flux(edge_radius)
    call require_close('toroidal-flux differential', q*dpsi, &
        toroidal_flux*edge_radius/sqrt(major_radius**2-edge_radius**2), &
        2.0e-12_dp)
    call require_close('toroidal-flux edge normalization', &
        toroidal_flux_value/toroidal_flux_edge, 1.0_dp, 4.0e-9_dp)

    do k = 1, 7
        radius_value = 0.05_dp + 0.40_dp*real(k, dp)/7.0_dp
        call evaluate_neort_circular_flux_continuation( &
            radius_value, edge_radius, psi_edge, toroidal_flux, q_axis, &
            delta_q, major_radius, psi, dpsi, q)
        call evaluate_full(radius_value, psi, dpsi, q, psi_tor_emitted, &
            s_tor_emitted, rho_tor_emitted, dvolume_ds_emitted, &
            dvolume_ds_cgs_emitted, current_emitted)
        toroidal_flux_value = integrated_toroidal_flux(radius_value)
        call require_close('emitted toroidal flux', psi_tor_emitted, &
            direct_toroidal_flux(radius_value), 4.0e-9_dp)
        call require_close('quadrature toroidal flux', toroidal_flux_value, &
            direct_toroidal_flux(radius_value), 4.0e-9_dp)
        rho_tor = sqrt(direct_toroidal_flux(radius_value)/toroidal_flux_edge)
        call require_close('emitted s_tor map', s_tor_emitted, &
            direct_toroidal_flux(radius_value)/toroidal_flux_edge, 4.0e-9_dp)
        call require_close('emitted rho_tor map', rho_tor_emitted, rho_tor, &
            4.0e-9_dp)
        dvolume_ds = direct_volume_jacobian(radius_value)
        call require_close('emitted dV/ds_tor', dvolume_ds_emitted, &
            dvolume_ds, 2.0e-7_dp)
        call require_close('emitted dV/ds_tor CGS', dvolume_ds_cgs_emitted, &
            1.0e6_dp*dvolume_ds, 2.0e-7_dp)
        call require_close('finite-difference dV/ds_tor', dvolume_ds, &
            finite_difference_volume_jacobian(radius_value), 2.0e-7_dp)
    end do

    call evaluate_neort_circular_flux_continuation(edge_radius, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, major_radius, psi, dpsi, q)
    call evaluate_full(edge_radius, psi, dpsi, q, psi_tor_emitted, &
        s_tor_emitted, rho_tor_emitted, dvolume_ds_emitted, &
        dvolume_ds_cgs_emitted, current_emitted)
    current_oracle = ampere_contour_current(dpsi)
    call require_close('emitted Ampere contour current', current_emitted, &
        current_oracle, 2.0e-10_dp)

    write (*, '(A)') 'test_gc_circular_flux_continuation OK'
    call pass_test

contains

    subroutine evaluate_neort_circular_flux_continuation(radius_value, &
            edge_radius_value, psi_edge_value, toroidal_flux_value, &
            q_axis_value, delta_q_value, major_radius_value, psi_value, &
            dpsi_value, q_value)
        real(dp), intent(in) :: radius_value, edge_radius_value, psi_edge_value
        real(dp), intent(in) :: toroidal_flux_value, q_axis_value, delta_q_value
        real(dp), intent(in) :: major_radius_value
        real(dp), intent(out) :: psi_value, dpsi_value, q_value
        real(dp) :: psi_tor_dummy, s_tor_dummy, rho_tor_dummy
        real(dp) :: dvolume_dummy, dvolume_cgs_dummy, current_dummy

        call evaluate_neort_circular_flux_continuation_raw(radius_value, &
            edge_radius_value, psi_edge_value, toroidal_flux_value, &
            q_axis_value, delta_q_value, major_radius_value, psi_value, &
            dpsi_value, q_value, psi_tor_dummy, s_tor_dummy, rho_tor_dummy, &
            dvolume_dummy, dvolume_cgs_dummy, current_dummy)
    end subroutine evaluate_neort_circular_flux_continuation

    subroutine evaluate_full(radius_value, psi_value, dpsi_value, q_value, &
            psi_tor_value, s_tor_value, rho_tor_value, dvolume_value, &
            dvolume_cgs_value, current_value)
        real(dp), intent(in) :: radius_value
        real(dp), intent(out) :: psi_value, dpsi_value, q_value
        real(dp), intent(out) :: psi_tor_value, s_tor_value, rho_tor_value
        real(dp), intent(out) :: dvolume_value, dvolume_cgs_value, current_value

        call evaluate_neort_circular_flux_continuation_raw(radius_value, &
            edge_radius, psi_edge, toroidal_flux, q_axis, delta_q, &
            major_radius, psi_value, &
            dpsi_value, q_value, psi_tor_value, s_tor_value, rho_tor_value, &
            dvolume_value, dvolume_cgs_value, current_value)
    end subroutine evaluate_full

    subroutine require_close(label, got, expected, tolerance)
        character(*), intent(in) :: label
        real(dp), intent(in) :: got, expected, tolerance

        if (abs(got-expected) > tolerance*max(1.0_dp, abs(expected))) then
            write (*, '(A,1X,ES24.16,1X,ES24.16)') trim(label), got, expected
            error stop 'circular continuation oracle failed'
        end if
    end subroutine require_close

    real(dp) function numerical_flux(radius_value)
        real(dp), intent(in) :: radius_value
        real(dp) :: lower, step, left, right
        integer :: j

        lower = edge_radius
        step = (radius_value-lower)/real(quadrature_steps, dp)
        numerical_flux = psi_edge
        left = oracle_dpsi(lower)
        do j = 1, quadrature_steps
            right = oracle_dpsi(lower+real(j, dp)*step)
            numerical_flux = numerical_flux + 0.5_dp*step*(left+right)
            left = right
        end do
    end function numerical_flux

    pure real(dp) function oracle_dpsi(radius_value)
        real(dp), intent(in) :: radius_value

        oracle_dpsi = toroidal_flux*radius_value/(target_q(radius_value)* &
            sqrt(major_radius**2-radius_value**2))
    end function oracle_dpsi

    pure real(dp) function target_q(radius_value)
        real(dp), intent(in) :: radius_value

        target_q = q_axis+delta_q*(radius_value/edge_radius)**2
    end function target_q

    real(dp) function field_line_pitch(radius_value, dpsi_value)
        real(dp), intent(in) :: radius_value, dpsi_value
        real(dp) :: step, theta, total, major_radius_local
        integer :: j

        step = 2.0_dp*acos(-1.0_dp)/real(quadrature_steps, dp)
        total = 0.0_dp
        do j = 0, quadrature_steps-1
            theta = step*real(j, dp)
            major_radius_local = major_radius+radius_value*cos(theta)
            total = total + toroidal_flux*radius_value/ &
                (major_radius_local*dpsi_value)
        end do
        field_line_pitch = total/real(quadrature_steps, dp)
    end function field_line_pitch

    real(dp) function direct_toroidal_flux(radius_value)
        real(dp), intent(in) :: radius_value

        direct_toroidal_flux = toroidal_flux*(major_radius- &
            sqrt(major_radius**2-radius_value**2))
    end function direct_toroidal_flux

    real(dp) function integrated_toroidal_flux(radius_value)
        real(dp), intent(in) :: radius_value
        real(dp) :: step, left, right
        integer :: j

        step = radius_value/real(quadrature_steps, dp)
        integrated_toroidal_flux = 0.0_dp
        left = target_q(0.0_dp)*oracle_dpsi(0.0_dp)
        do j = 1, quadrature_steps
            right = target_q(real(j, dp)*step)* &
                oracle_dpsi(real(j, dp)*step)
            integrated_toroidal_flux = integrated_toroidal_flux + &
                0.5_dp*step*(left+right)
            left = right
        end do
    end function integrated_toroidal_flux

    real(dp) function direct_volume_jacobian(radius_value)
        real(dp), intent(in) :: radius_value

        direct_volume_jacobian = 4.0_dp*acos(-1.0_dp)**2*major_radius* &
            (major_radius-sqrt(major_radius**2-edge_radius**2))* &
            sqrt(major_radius**2-radius_value**2)
    end function direct_volume_jacobian

    real(dp) function finite_difference_volume_jacobian(radius_value)
        real(dp), intent(in) :: radius_value
        real(dp) :: dh

        dh = 1.0e-6_dp
        finite_difference_volume_jacobian = &
            (direct_volume(radius_value+dh)-direct_volume(radius_value-dh))/ &
            (direct_s_tor(radius_value+dh)-direct_s_tor(radius_value-dh))
    end function finite_difference_volume_jacobian

    real(dp) function direct_volume(radius_value)
        real(dp), intent(in) :: radius_value

        direct_volume = 2.0_dp*acos(-1.0_dp)**2*major_radius*radius_value**2
    end function direct_volume

    real(dp) function direct_s_tor(radius_value)
        real(dp), intent(in) :: radius_value

        direct_s_tor = direct_toroidal_flux(radius_value)/toroidal_flux_edge
    end function direct_s_tor

    real(dp) function ampere_contour_current(dpsi_edge_value)
        real(dp), intent(in) :: dpsi_edge_value
        real(dp) :: step, theta, total, local_R
        integer :: j

        step = 2.0_dp*acos(-1.0_dp)/real(quadrature_steps, dp)
        total = 0.0_dp
        do j = 0, quadrature_steps-1
            theta = step*real(j, dp)
            local_R = major_radius+edge_radius*cos(theta)
            total = total+dpsi_edge_value*edge_radius/local_R
        end do
        ampere_contour_current = total*step/mu0
    end function ampere_contour_current

end program test_gc_circular_flux_continuation
