program test_gc_circular_flux_continuation
    !! Independent numerical oracle for the Fortsym-emitted continuation.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_circular_flux_continuation_limit_symbolic, only: &
        evaluate_neort_circular_flux_continuation_limit
    use neort_circular_flux_continuation_symbolic, only: &
        evaluate_neort_circular_flux_continuation
    use util_for_test, only: pass_test
    implicit none

    real(dp), parameter :: edge_radius = 0.5_dp
    real(dp), parameter :: psi_edge = 0.4_dp
    real(dp), parameter :: toroidal_flux = 3.2_dp
    real(dp), parameter :: q_axis = 1.5_dp
    real(dp), parameter :: delta_q = 2.5_dp
    real(dp), parameter :: h = 1.0e-5_dp
    integer, parameter :: quadrature_steps = 20000
    real(dp) :: psi, dpsi, q, psi_left, psi_right, dpsi_left
    real(dp) :: dpsi_right, psi_limit, dpsi_limit, q_limit
    real(dp) :: psi_limit_left, psi_limit_right
    real(dp) :: radius_value, oracle_value
    integer :: k

    call evaluate_neort_circular_flux_continuation(edge_radius, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, psi, dpsi, q)
    call require_close('LCFS value', psi, psi_edge, 2.0e-13_dp)
    call evaluate_neort_circular_flux_continuation(edge_radius-h, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, psi_left, dpsi_left, q)
    call evaluate_neort_circular_flux_continuation(edge_radius+h, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, psi_right, dpsi_right, q)
    call require_close('left first derivative', (psi-psi_left)/h, dpsi, &
        2.0e-5_dp)
    call require_close('right first derivative', (psi_right-psi)/h, dpsi, &
        2.0e-5_dp)
    call require_close('first derivative continuity', dpsi_left, dpsi_right, &
        4.0e-5_dp)

    do k = 0, 256
        radius_value = edge_radius + 1.5_dp*edge_radius*real(k, dp)/256.0_dp
        call evaluate_neort_circular_flux_continuation( &
            radius_value, &
            edge_radius, psi_edge, toroidal_flux, q_axis, delta_q, psi, dpsi, q)
        if (dpsi <= 0.0_dp) error stop 'exterior flux is not strictly monotone'
        if (k > 0) then
            if (psi <= psi_left) error stop 'exterior flux values are not increasing'
        end if
        oracle_value = numerical_flux(radius_value)
        call require_close('numerical flux oracle', psi, oracle_value, &
            2.0e-9_dp)
        psi_left = psi
    end do
    call evaluate_neort_circular_flux_continuation(edge_radius, edge_radius, &
        psi_edge, toroidal_flux, q_axis, delta_q, psi, dpsi, q)
    call require_close('defining LCFS slope', dpsi, &
        toroidal_flux*edge_radius/(q_axis+delta_q), 2.0e-13_dp)

    call evaluate_neort_circular_flux_continuation_limit(edge_radius, &
        edge_radius, psi_edge, toroidal_flux, q_axis, psi_limit, dpsi_limit, &
        q_limit)
    call require_close('zero-shear LCFS value', psi_limit, psi_edge, &
        2.0e-13_dp)
    call evaluate_neort_circular_flux_continuation_limit(edge_radius-h, &
        edge_radius, psi_edge, toroidal_flux, q_axis, psi_limit_left, &
        dpsi_left, q_limit)
    call evaluate_neort_circular_flux_continuation_limit(edge_radius+h, &
        edge_radius, psi_edge, toroidal_flux, q_axis, psi_limit_right, &
        dpsi_right, q_limit)
    call require_close('zero-shear first derivative', &
        (psi_limit_right-psi_limit_left)/(2.0_dp*h), dpsi_limit, 2.0e-10_dp)
    call require_close('zero-shear q', q_limit, q_axis, 2.0e-13_dp)

    write (*, '(A)') 'test_gc_circular_flux_continuation OK'
    call pass_test

contains

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

        oracle_dpsi = toroidal_flux*radius_value/(q_axis + delta_q* &
            (radius_value/edge_radius)**2)
    end function oracle_dpsi

end program test_gc_circular_flux_continuation
