program test_gc_cylindrical_topology
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_gc_cylindrical_model, only: GC_CYL_SUCCESS, &
        gc_cylindrical_allowed_value_i
    use neort_gc_cylindrical_topology, only: gc_cylindrical_allowed_region_set_t, &
        find_gc_cylindrical_allowed_regions, canonical_measure_from_curve
    use util_for_test, only: pass_test

    implicit none

    type(gc_cylindrical_allowed_region_set_t) :: regions
    type(gc_cylindrical_allowed_region_set_t) :: negative_sigma_regions
    real(dp), parameter :: lower = 0.0_dp, upper = 6.0_dp
    integer :: status

    call find_gc_cylindrical_allowed_regions(disconnected_allowed_value, lower, &
        upper, 601, 1, regions, status)
    if (status /= GC_CYL_SUCCESS) error stop 'allowed-region scan failed'
    if (regions%ncomponents /= 2) error stop 'disconnected components were merged'
    if (regions%components(1)%component_id /= 1 .or. &
        regions%components(2)%component_id /= 2) then
        error stop 'component ids are not data-oriented'
    end if
    if (any(regions%components(:)%sigma /= 1)) error stop 'sigma was not retained'
    if (regions%nroots /= 4) error stop 'all v_parallel roots were not retained'
    if (abs(regions%components(1)%x_begin - 1.0_dp) > 2.0e-5_dp .or. &
        abs(regions%components(1)%x_end - 2.0_dp) > 2.0e-5_dp) then
        error stop 'first allowed component bounds are wrong'
    end if
    if (abs(regions%components(2)%x_begin - 4.0_dp) > 2.0e-5_dp .or. &
        abs(regions%components(2)%x_end - 5.0_dp) > 2.0e-5_dp) then
        error stop 'second allowed component bounds are wrong'
    end if
    if (abs(regions%components(1)%canonical_measure - 1.0_dp) > 2.0e-3_dp) then
        error stop 'canonical component measure is not |d psi*/dx|'
    end if
    if (abs(regions%components(2)%canonical_measure - 1.0_dp) > 2.0e-3_dp) then
        error stop 'second canonical component measure is wrong'
    end if
    call find_gc_cylindrical_allowed_regions(disconnected_allowed_value, lower, &
        upper, 601, -1, negative_sigma_regions, status)
    if (status /= GC_CYL_SUCCESS) error stop 'negative-sigma scan failed'
    if (any(negative_sigma_regions%components(:)%sigma /= -1)) then
        error stop 'negative sigma was not retained'
    end if

    call check_curve_measure()
    call check_tangent_root()
    call pass_test

contains

    subroutine disconnected_allowed_value(x, value, derivative, psi_star, &
            dpsi_star_dx, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value, derivative, psi_star, dpsi_star_dx
        integer, intent(out) :: status

        value = -(x - 1.0_dp)*(x - 2.0_dp)*(x - 4.0_dp)*(x - 5.0_dp)
        derivative = -4.0_dp*x**3 + 36.0_dp*x**2 - 98.0_dp*x + 78.0_dp
        psi_star = x
        dpsi_star_dx = 1.0_dp
        status = GC_CYL_SUCCESS
    end subroutine disconnected_allowed_value

    subroutine check_curve_measure()
        real(dp) :: x(5), psi(5), measure

        x = [0.0_dp, 0.5_dp, 1.0_dp, 1.5_dp, 2.0_dp]
        psi = [0.0_dp, 1.0_dp, 0.0_dp, -1.0_dp, 0.0_dp]
        call canonical_measure_from_curve(x, psi, measure, status)
        if (status /= GC_CYL_SUCCESS) error stop 'curve measure failed'
        if (abs(measure - 4.0_dp) > 1.0e-13_dp) then
            error stop 'curve canonical measure omitted an extremum'
        end if
    end subroutine check_curve_measure

    subroutine check_tangent_root()
        call find_gc_cylindrical_allowed_regions(tangent_allowed_value, 0.0_dp, &
            5.0_dp, 501, -1, regions, status)
        if (status /= GC_CYL_SUCCESS) error stop 'tangent-root scan failed'
        if (regions%nroots /= 3) error stop 'even root was omitted'
        if (regions%ncomponents /= 3) error stop 'tangent root merged components'
        if (abs(regions%roots(1) - 1.003_dp) > 2.0e-5_dp) then
            error stop 'tangent root location is wrong'
        end if
        if (abs(regions%components(1)%x_end - 1.003_dp) > 2.0e-5_dp) then
            error stop 'tangent root did not delimit first component'
        end if
        if (abs(regions%components(2)%x_begin - 1.003_dp) > 2.0e-5_dp) then
            error stop 'tangent root did not delimit second component'
        end if
    end subroutine check_tangent_root

    subroutine tangent_allowed_value(x, value, derivative, psi_star, &
            dpsi_star_dx, status)
        real(dp), intent(in) :: x
        real(dp), intent(out) :: value, derivative, psi_star, dpsi_star_dx
        integer, intent(out) :: status
        real(dp) :: first_factor, second_factor, third_factor

        first_factor = x - 1.003_dp
        second_factor = x - 2.0_dp
        third_factor = x - 4.0_dp
        value = first_factor**2*second_factor*third_factor
        derivative = 2.0_dp*first_factor*second_factor*third_factor &
            +first_factor**2*(second_factor + third_factor)
        psi_star = x
        dpsi_star_dx = 1.0_dp
        status = GC_CYL_SUCCESS
    end subroutine tangent_allowed_value

end program test_gc_cylindrical_topology
