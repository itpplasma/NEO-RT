program test_gc_eqdsk_cut_axis_limit
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_eqdsk_cut_axis_curvature_interval_symbolic, only: &
        evaluate_neort_eqdsk_cut_axis_curvature_interval
    use neort_eqdsk_cut_axis_limit_symbolic, only: &
        evaluate_neort_eqdsk_cut_axis_limit
    use neort_eqdsk_cut_axis_rho_limit_symbolic, only: &
        evaluate_neort_eqdsk_cut_axis_rho_limit
    use neort_eqdsk_axis_stationarity_krawczyk_interval_symbolic, only: &
        evaluate_neort_eqdsk_axis_stationarity_krawczyk_interval
    use neort_eqdsk_axis_stationarity_newton_symbolic, only: &
        evaluate_neort_eqdsk_axis_stationarity_newton
    use neort_eqdsk_cut_flux_coordinate_symbolic, only: &
        evaluate_neort_eqdsk_cut_flux_coordinate
    use neort_eqdsk_rho_tor_map_symbolic, only: &
        evaluate_neort_eqdsk_rho_tor_map
    use neort_gc_outward_interval, only: gc_outward_interval, &
        gc_outward_interval_t
    implicit none

    real(dp), parameter :: slope = -0.2_dp
    real(dp), parameter :: psi_rr = 3.0_dp
    real(dp), parameter :: psi_rz = 0.4_dp
    real(dp), parameter :: psi_zz = 2.0_dp
    real(dp), parameter :: psi_sep = 5.0_dp
    real(dp), parameter :: delta_psihat = 0.03_dp
    real(dp), parameter :: tolerance = 2.0e-14_dp
    real(dp) :: curvature, absolute_derivative, delta_r
    real(dp) :: oracle_curvature, oracle_delta_r, oracle_derivative
    real(dp) :: s_tor, dstor_drho, dpsihat_drho, dR_drho, dZ_drho
    real(dp) :: oracle_dR_drho
    real(dp) :: axis_determinant, axis_newton_R, axis_newton_Z, axis_inverse(4)
    type(gc_outward_interval_t) :: interval_curvature, interval_determinant
    type(gc_outward_interval_t) :: axis_krawczyk_R, axis_krawczyk_Z

    ! Independent quadratic axis oracle.  Along Z=slope*R, evaluate
    ! psi_hat(R) directly and recover its second derivative from a centered
    ! second difference; no generated cut expression is reused here.
    oracle_curvature = (quadratic_psihat(1.0_dp) &
        - 2.0_dp*quadratic_psihat(0.0_dp) &
        + quadratic_psihat(-1.0_dp))
    oracle_delta_r = sqrt(2.0_dp*delta_psihat/oracle_curvature)
    oracle_derivative = oracle_curvature*oracle_delta_r

    call evaluate_neort_eqdsk_cut_axis_limit(slope, psi_rr, psi_rz, psi_zz, &
        psi_sep, delta_psihat, -1.0_dp, curvature, absolute_derivative, delta_r)
    call require_close(curvature, oracle_curvature, tolerance, &
        'axis curvature')
    call require_close(absolute_derivative, oracle_derivative, tolerance, &
        'axis derivative limit')
    call require_close(delta_r, -oracle_delta_r, tolerance, &
        'inboard axis inverse limit')

    call evaluate_neort_eqdsk_cut_axis_limit(slope, psi_rr, psi_rz, psi_zz, &
        psi_sep, delta_psihat, 1.0_dp, curvature, absolute_derivative, delta_r)
    call require_close(delta_r, oracle_delta_r, tolerance, &
        'outboard axis inverse limit')

    call evaluate_neort_eqdsk_rho_tor_map(0.4_dp, s_tor, dstor_drho)
    call require_close(s_tor, 0.16_dp, tolerance, 'rho_tor squared map')
    call require_close(dstor_drho, 0.8_dp, tolerance, &
        'rho_tor derivative map')
    call evaluate_neort_eqdsk_cut_flux_coordinate(dstor_drho, 0.7_dp, &
        0.2_dp, -0.3_dp, dpsihat_drho, dR_drho, dZ_drho)
    call require_close(dpsihat_drho, 0.56_dp, tolerance, &
        'rho_tor flux chain')
    call require_close(dR_drho, 2.8_dp, tolerance, &
        'inverse flux-chart chain')
    call require_close(dZ_drho, -0.84_dp, tolerance, &
        'inverse vertical-chart chain')

    oracle_dR_drho = sqrt(2.0_dp*0.7_dp/oracle_curvature)
    call evaluate_neort_eqdsk_cut_axis_rho_limit(slope, psi_rr, psi_rz, &
        psi_zz, psi_sep, 0.7_dp, -1.0_dp, curvature, dR_drho, dZ_drho)
    call require_close(dR_drho, -oracle_dR_drho, tolerance, &
        'inboard axis rho_tor derivative')
    call require_close(dZ_drho, -slope*oracle_dR_drho, tolerance, &
        'inboard axis vertical rho_tor derivative')
    call evaluate_neort_eqdsk_cut_axis_rho_limit(slope, psi_rr, psi_rz, &
        psi_zz, psi_sep, 0.7_dp, 1.0_dp, curvature, dR_drho, dZ_drho)
    call require_close(dR_drho, oracle_dR_drho, tolerance, &
        'outboard axis rho_tor derivative')

    call evaluate_neort_eqdsk_cut_axis_curvature_interval( &
        gc_outward_interval(-0.21_dp, -0.19_dp), &
        gc_outward_interval(2.99_dp, 3.01_dp), &
        gc_outward_interval(0.39_dp, 0.41_dp), &
        gc_outward_interval(1.99_dp, 2.01_dp), &
        gc_outward_interval(4.99_dp, 5.01_dp), interval_curvature, &
        interval_determinant)
    call require(interval_curvature%lo <= oracle_curvature .and. &
        interval_curvature%hi >= oracle_curvature, &
        'axis curvature interval excluded its quadratic oracle')
    call require(interval_curvature%lo > 0.0_dp, &
        'nondegenerate positive axis curvature was not certified')
    call require(interval_determinant%lo > 0.0_dp, &
        'positive-definite axis Hessian determinant was not certified')

    ! Independent quadratic stationarity oracle.  The Hessian is
    ! [[3,0.5],[0.5,2]] and the exact root is (1,2).  At (1.2,1.7),
    ! grad(psi)=(0.45,-0.5), so one generated Newton step must reach the root.
    call evaluate_neort_eqdsk_axis_stationarity_newton(1.2_dp, 1.7_dp, &
        0.45_dp, -0.5_dp, 3.0_dp, 0.5_dp, 2.0_dp, axis_determinant, &
        axis_newton_R, axis_newton_Z, axis_inverse(1), axis_inverse(2), &
        axis_inverse(3), axis_inverse(4))
    call require_close(axis_determinant, 5.75_dp, tolerance, &
        'axis Hessian determinant')
    call require_close(axis_newton_R, 1.0_dp, tolerance, 'axis Newton R')
    call require_close(axis_newton_Z, 2.0_dp, tolerance, 'axis Newton Z')
    call evaluate_neort_eqdsk_axis_stationarity_krawczyk_interval( &
        point(1.2_dp), point(1.7_dp), point(0.45_dp), point(-0.5_dp), &
        point(3.0_dp), point(0.5_dp), point(2.0_dp), point(3.0_dp), &
        point(0.5_dp), point(2.0_dp), &
        gc_outward_interval(-0.2_dp, 0.2_dp), &
        gc_outward_interval(-0.3_dp, 0.3_dp), axis_krawczyk_R, &
        axis_krawczyk_Z)
    call require(axis_krawczyk_R%lo <= 1.0_dp .and. &
        axis_krawczyk_R%hi >= 1.0_dp, &
        'axis Krawczyk R enclosure excluded the exact root')
    call require(axis_krawczyk_Z%lo <= 2.0_dp .and. &
        axis_krawczyk_Z%hi >= 2.0_dp, &
        'axis Krawczyk Z enclosure excluded the exact root')

    write (*, '(a)') 'test_gc_eqdsk_cut_axis_limit OK'

contains

    pure function point(value) result(interval)
        real(dp), intent(in) :: value
        type(gc_outward_interval_t) :: interval

        interval = gc_outward_interval(value, value)
    end function point

    pure real(dp) function quadratic_psihat(radius_offset)
        real(dp), intent(in) :: radius_offset
        real(dp) :: height_offset

        height_offset = slope*radius_offset
        quadratic_psihat = 0.5_dp*(psi_rr*radius_offset**2 &
            + 2.0_dp*psi_rz*radius_offset*height_offset &
            + psi_zz*height_offset**2)/psi_sep
    end function quadratic_psihat

    subroutine require_close(actual, expected, local_tolerance, label)
        real(dp), intent(in) :: actual, expected, local_tolerance
        character(*), intent(in) :: label

        call require(abs(actual-expected) <= local_tolerance &
            *max(1.0_dp, abs(expected)), trim(label)//' mismatch')
    end subroutine require_close

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_cut_axis_limit
