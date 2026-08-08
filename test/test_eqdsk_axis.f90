program test_eqdsk_axis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s, do_magfie, s, psi_pr, q, dqds, iota, &
        Bthcov, Bphcov, sign_theta
    use neort_gc_dynamics, only: gc_field_sample_t
    use neort_gc_eqdsk_adapter, only: eqdsk_gc_field_t
    use neort_gc_models, only: GC_MODEL_SUCCESS
    use neort_circular_flux_continuation_symbolic, only: &
        evaluate_neort_circular_flux_continuation
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
    real(dp) :: theta, dtheta, q_integral, q_adapter_integral
    real(dp) :: q_adapter, dqds_adapter, psi_adapter, dpsi_adapter
    real(dp) :: psi_edge_adapter, h_radial
    real(dp) :: legacy_state(7)
    real(dp) :: dVds_inner, dVds_outer
    real(dp) :: dVds_expected, dVds_max_relative_error
    real(dp), parameter :: s_inner = (0.001_dp/64.0_dp)**2
    real(dp), parameter :: s_outer = (0.002_dp/64.0_dp)**2
    ! The 65x65 generated grid gives 1.775e-3 against this pointwise oracle;
    ! the next nested 129x129 grid gives 7.559e-4.
    real(dp), parameter :: dVds_interpolation_tolerance = 1.0e-3_dp
    complex(dp) :: bamp
    type(eqdsk_gc_field_t) :: gc_field
    type(gc_field_sample_t) :: gc_sample, gc_sample_minus, gc_sample_plus
    type(gc_field_sample_t) :: gc_sample_theta_minus, gc_sample_theta_plus
    integer :: k, gc_status

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    call get_environment_variable('RZ_PERT_FILE', perturbation_file)
    if (len_trim(eqdsk_file) == 0 .or. len_trim(perturbation_file) == 0) then
        error stop "GEQDSK test fixtures are not configured"
    end if

    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))

    call set_s(s_inner)
    call init_magfie_at_s()
    call init_flux_surface_average(s_inner)
    dVds_inner = dVds
    call set_s(s_outer)
    call init_magfie_at_s()
    call init_flux_surface_average(s_outer)
    dVds_outer = dVds
    if (abs(dVds_outer/dVds_inner - 1.0_dp) > 1.0e-2_dp) then
        write(*,*) 'Circular-axis dV/ds must have a finite constant limit:', &
            dVds_inner, dVds_outer
        error stop "GEQDSK near-axis volume derivative failed"
    end if

    dVds_max_relative_error = 0.0_dp
    do k = 10, 160
        call set_s(real(k, dp)/250.0_dp)
        call init_magfie_at_s()
        call init_flux_surface_average(real(k, dp)/250.0_dp)
        call circular_dvolume_oracle(real(k, dp)/250.0_dp, dVds_expected)
        dVds_max_relative_error = max(dVds_max_relative_error, &
            abs(dVds/dVds_expected - 1.0_dp))
    end do
    if (dVds_max_relative_error > dVds_interpolation_tolerance) then
        write(*,*) 'Circular-equilibrium dV/ds_tor interpolation error: ', &
            dVds_max_relative_error
        error stop "GEQDSK radial metric differs from generated oracle"
    end if

    call set_s(0.25_dp)
    call init_magfie_at_s()
    dtheta = 2.0_dp*pi/real(ntheta, dp)

    ! The real-space orbit adapter is a read-only view.  Its poloidal flux is
    ! checked both against an independent finite difference and against the
    ! coordinate identity sqrt(g) B^theta=d(psi_pol)/ds.
    legacy_state = [s, q, dqds, iota, psi_pr, Bthcov, Bphcov]
    gc_field%field_scale = 1.0_dp
    call gc_field%radial_profiles(0.64_dp, q_adapter, dqds_adapter, &
        psi_adapter, dpsi_adapter, psi_edge_adapter, gc_status)
    if (gc_status /= GC_MODEL_SUCCESS) error stop 'GEQDSK GC profiles failed'
    call gc_field%evaluate([0.64_dp, 0.0_dp, 0.37_dp], gc_sample, gc_status)
    if (gc_status /= GC_MODEL_SUCCESS) error stop 'GEQDSK GC sample failed'
    h_radial = 1.0e-5_dp
    call gc_field%evaluate([0.64_dp - h_radial, 0.0_dp, 0.37_dp], &
        gc_sample_minus, gc_status)
    if (gc_status /= GC_MODEL_SUCCESS) error stop 'GEQDSK GC minus sample failed'
    call gc_field%evaluate([0.64_dp + h_radial, 0.0_dp, 0.37_dp], &
        gc_sample_plus, gc_status)
    if (gc_status /= GC_MODEL_SUCCESS) error stop 'GEQDSK GC plus sample failed'
    if (abs((gc_sample_plus%psi - gc_sample_minus%psi)/(2.0_dp*h_radial) &
            - gc_sample%grad_psi(1)) > 2.0e-6_dp*abs(gc_sample%grad_psi(1))) then
        error stop 'GEQDSK GC poloidal-flux derivative failed'
    end if
    call gc_field%evaluate([0.64_dp, 0.0_dp, 0.37_dp - h_radial], &
        gc_sample_theta_minus, gc_status)
    if (gc_status /= GC_MODEL_SUCCESS) error stop 'GEQDSK GC theta-minus sample failed'
    call gc_field%evaluate([0.64_dp, 0.0_dp, 0.37_dp + h_radial], &
        gc_sample_theta_plus, gc_status)
    if (gc_status /= GC_MODEL_SUCCESS) error stop 'GEQDSK GC theta-plus sample failed'
    if (abs((gc_sample_theta_plus%psi - gc_sample_theta_minus%psi) &
            /(2.0_dp*h_radial) - gc_sample%grad_psi(3)) &
            > 2.0e-6_dp*max(abs(gc_sample%grad_psi(1)), 1.0_dp)) then
        error stop 'GEQDSK GC poloidal flux has wrong theta derivative'
    end if
    if (abs(gc_sample%sqrtg*gc_sample%bmod*gc_sample%hcon(3) &
            - gc_sample%grad_psi(1)) > 2.0e-12_dp*abs(gc_sample%grad_psi(1))) then
        write(*,*) 'sqrt(g) B^theta, dpsi_pol/ds:', &
            gc_sample%sqrtg*gc_sample%bmod*gc_sample%hcon(3), &
            gc_sample%grad_psi(1)
        error stop 'GEQDSK GC flux/field orientation failed'
    end if
    if (abs(gc_sample%sqrtg*gc_sample%bmod*gc_sample%hcon(1) &
            + gc_sample%grad_psi(3)) &
            > 2.0e-12_dp*abs(gc_sample%grad_psi(1))) then
        error stop 'GEQDSK GC radial field/flux identity failed'
    end if

    q_adapter_integral = 0.0_dp
    do k = 0, ntheta - 1
        theta = -pi + (real(k, dp) + 0.5_dp)*dtheta
        call gc_field%evaluate([0.64_dp, 0.0_dp, theta], gc_sample, gc_status)
        if (gc_status /= GC_MODEL_SUCCESS) error stop 'GEQDSK GC q sample failed'
        q_adapter_integral = q_adapter_integral &
            + gc_sample%hcon(2)/gc_sample%hcon(3)*dtheta/(2.0_dp*pi)
    end do
    if (abs(q_adapter_integral/q_adapter - 1.0_dp) > 1.0e-2_dp) then
        error stop 'GEQDSK GC field-line pitch failed'
    end if
    if (any([s, q, dqds, iota, psi_pr, Bthcov, Bphcov] /= legacy_state)) then
        error stop 'GEQDSK GC sampling changed the legacy surface cache'
    end if

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

contains

    subroutine circular_dvolume_oracle(s_tor_target, dvolume_cgs)
        !! The Fortsym kernel owns both the normalized toroidal-flux map and
        !! dV/ds_tor.  Bisection is numerical orchestration only; the focused
        !! circular-continuation test independently validates both quantities
        !! against contour quadrature and a finite-difference volume oracle.
        real(dp), intent(in) :: s_tor_target
        real(dp), intent(out) :: dvolume_cgs
        real(dp), parameter :: edge_radius = 0.5_dp
        real(dp), parameter :: major_radius = 1.6_dp
        real(dp), parameter :: toroidal_flux = 3.2_dp
        real(dp), parameter :: q_axis = 1.5_dp
        real(dp), parameter :: delta_q = 2.5_dp
        integer, parameter :: bisection_steps = 80
        real(dp) :: radius_left, radius_right, radius_midpoint
        real(dp) :: psi_value, dpsi_value, q_value, psi_tor_value
        real(dp) :: s_tor_value, rho_tor_value, dvolume_value
        real(dp) :: dvolume_cgs_value, current_value
        integer :: j

        radius_left = 0.0_dp
        radius_right = edge_radius
        do j = 1, bisection_steps
            radius_midpoint = 0.5_dp*(radius_left + radius_right)
            call evaluate_neort_circular_flux_continuation(radius_midpoint, &
                edge_radius, 1.0_dp, toroidal_flux, q_axis, delta_q, &
                major_radius, psi_value, dpsi_value, q_value, psi_tor_value, &
                s_tor_value, rho_tor_value, dvolume_value, dvolume_cgs_value, &
                current_value)
            if (s_tor_value < s_tor_target) then
                radius_left = radius_midpoint
            else
                radius_right = radius_midpoint
            end if
        end do
        radius_midpoint = 0.5_dp*(radius_left + radius_right)
        call evaluate_neort_circular_flux_continuation(radius_midpoint, &
            edge_radius, 1.0_dp, toroidal_flux, q_axis, delta_q, major_radius, &
            psi_value, dpsi_value, q_value, psi_tor_value, s_tor_value, &
            rho_tor_value, dvolume_value, dvolume_cgs_value, current_value)
        dvolume_cgs = dvolume_cgs_value
    end subroutine circular_dvolume_oracle
end program test_eqdsk_axis
