program test_eqdsk_axis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s, do_magfie, psi_pr, q, sign_theta
    use do_magfie_pert_mod, only: inp_swi_pert, read_boozer_pert_file, &
        init_magfie_pert_at_s, do_magfie_pert_amp, set_mph
    use neort_orbit, only: fieldline_label_component
    use neort_magfie, only: init_flux_surface_average
    use driftorbit, only: dVds
    use util, only: pi
    use util_for_test, only: pass_test

    implicit none

    integer, parameter :: ntheta = 4096
    character(len=1024) :: eqdsk_file, perturbation_file
    real(dp) :: x(3), bmod, sqrtg, bder(3), hcov(3), hcon(3), hcurl(3)
    real(dp) :: theta, dtheta, q_integral
    real(dp) :: dVds_inner, dVds_outer, dVds_expected
    real(dp) :: dVds_scan(151), dVds_curvature, dVds_max_curvature
    real(dp), parameter :: fixture_rmajor_cm = 160.0_dp
    real(dp), parameter :: fixture_aminor_cm = 50.0_dp
    real(dp), parameter :: s_tor_inner = (0.001_dp/fixture_aminor_cm)**2
    real(dp), parameter :: s_tor_outer = (0.002_dp/fixture_aminor_cm)**2
    complex(dp) :: bamp
    integer :: k

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    call get_environment_variable('RZ_PERT_FILE', perturbation_file)
    if (len_trim(eqdsk_file) == 0 .or. len_trim(perturbation_file) == 0) then
        error stop "GEQDSK test fixtures are not configured"
    end if

    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))

    ! The synthetic fixture has s_tor = r**2/a**2 and concentric circular
    ! surfaces. Its independent geometric oracle is therefore
    ! dV/ds_tor = 2*pi**2*R0*a**2, including the finite magnetic-axis limit.
    dVds_expected = 2.0_dp*pi**2*fixture_rmajor_cm*fixture_aminor_cm**2
    call set_s(s_tor_inner)
    call init_magfie_at_s()
    call init_flux_surface_average(s_tor_inner)
    dVds_inner = dVds
    call set_s(s_tor_outer)
    call init_magfie_at_s()
    call init_flux_surface_average(s_tor_outer)
    dVds_outer = dVds
    if (abs(dVds_inner/dVds_expected - 1.0_dp) > 1.0e-2_dp .or. &
        abs(dVds_outer/dVds_expected - 1.0_dp) > 1.0e-2_dp) then
        write(*,*) 'Circular-axis dV/ds_tor normalization failed:', &
            dVds_inner, dVds_outer, dVds_expected
        error stop "GEQDSK near-axis volume derivative failed"
    end if

    do k = 10, 160
        call set_s(real(k, dp)/250.0_dp)
        call init_magfie_at_s()
        call init_flux_surface_average(real(k, dp)/250.0_dp)
        dVds_scan(k - 9) = dVds
    end do
    dVds_max_curvature = 0.0_dp
    do k = 2, size(dVds_scan) - 1
        dVds_curvature = abs(dVds_scan(k + 1) - 2.0_dp*dVds_scan(k) &
            + dVds_scan(k - 1))/dVds_expected
        dVds_max_curvature = max(dVds_max_curvature, dVds_curvature)
    end do
    if (dVds_max_curvature > 8.0e-4_dp) then
        write(*,*) 'Circular-equilibrium dV/ds_tor radial metric rings:', &
            dVds_max_curvature
        error stop "GEQDSK radial metric rings between grid cells"
    end if

    call set_s(0.25_dp)
    call init_magfie_at_s()

    dtheta = 2.0_dp*pi/real(ntheta, dp)
    q_integral = 0.0_dp
    do k = 0, ntheta - 1
        theta = -pi + (real(k, dp) + 0.5_dp)*dtheta
        x = [0.25_dp, 0.0_dp, theta]
        call do_magfie(x, bmod, sqrtg, bder, hcov, hcon, hcurl)
        if (.not. all(ieee_is_finite([bmod, sqrtg, bder, hcov, hcon, hcurl]))) then
            error stop "GEQDSK field contains non-finite values"
        end if
        if (bmod <= 0.0_dp .or. abs(sqrtg) <= tiny(sqrtg) .or. &
            abs(hcon(3)) <= tiny(hcon(3))) then
            error stop "GEQDSK field is degenerate"
        end if
        q_integral = q_integral + hcon(2)/hcon(3)*dtheta/(2.0_dp*pi)
    end do
    if (q*q_integral <= 0.0_dp .or. abs(q_integral/q) < 0.5_dp .or. &
        abs(q_integral/q) > 2.0_dp) then
        write(*,*) 'GEQDSK q and integrated field-line pitch are inconsistent:', q, q_integral
        error stop "GEQDSK pitch oracle failed"
    end if
    if (sign_theta*sqrtg <= 0.0_dp) then
        write(*,*) 'GEQDSK handedness and coordinate Jacobian are inconsistent:', &
            sign_theta, sqrtg
        error stop "GEQDSK handedness oracle failed"
    end if
    if (psi_pr <= 0.0_dp) then
        write(*,*) 'GEQDSK toroidal flux disagrees with the analytic B_phi:', psi_pr
        error stop "GEQDSK toroidal-flux oracle failed"
    end if
    if (abs(fieldline_label_component(5.0_dp, 1.0_dp, 6.0_dp, 2.0_dp) - &
        2.0_dp) > 10.0_dp*epsilon(1.0_dp)) then
        error stop "field-line-label component failed its closed-form oracle"
    end if
    if (abs(fieldline_label_component(hcon(2), hcon(3), hcon(2), hcon(3))) > &
        10.0_dp*epsilon(1.0_dp)*abs(hcon(2))) then
        error stop "parallel motion changed the direct field-line label"
    end if

    inp_swi_pert = 11
    call set_mph(-3)
    call read_boozer_pert_file(trim(perturbation_file))
    call init_magfie_pert_at_s()
    x = [0.25_dp, 0.0_dp, 0.0_dp]
    call do_magfie_pert_amp(x, bamp)
    if (.not. ieee_is_finite(real(bamp)) .or. &
        .not. ieee_is_finite(aimag(bamp))) error stop "R-Z perturbation is non-finite"

    call pass_test
end program test_eqdsk_axis
