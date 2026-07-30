program test_boozer_chartmap_geometry
    ! Independent analytic oracle for the circular-torus chart
    !
    ! x = ((R0 + a*rho*cos(theta))*cos(phi),
    !      (R0 + a*rho*cos(theta))*sin(phi),
    !       a*rho*sin(theta)),  rho=sqrt(s).
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_chartmap_geometry, only: boozer_chartmap_geometry_t, &
        evaluate_boozer_metric, read_boozer_chartmap_geometry
    implicit none

    character(1024) :: path
    type(boozer_chartmap_geometry_t) :: geometry
    real(dp), allocatable :: metric(:, :, :), jacobian(:)
    real(dp), parameter :: s = 0.25_dp
    real(dp), parameter :: rho = sqrt(s)
    real(dp), parameter :: major_radius = 3.0_dp
    real(dp), parameter :: minor_scale = 0.5_dp
    real(dp), parameter :: tolerance = 2.0e-12_dp
    real(dp) :: theta, cylindrical_radius
    integer :: index

    if (command_argument_count() /= 1) error stop "chartmap fixture path required"
    call get_command_argument(1, path)
    call read_boozer_chartmap_geometry(trim(path), geometry)
    allocate (metric(3, 3, geometry%theta_count))
    allocate (jacobian(geometry%theta_count))
    call evaluate_boozer_metric(geometry, s, metric, jacobian)

    do index = 1, geometry%theta_count
        theta = geometry%theta(index)
        cylindrical_radius = major_radius + minor_scale*rho*cos(theta)
        call assert_close( &
            metric(1, 1, index), minor_scale**2/(4.0_dp*rho**2), "g_ss")
        call assert_close(metric(1, 2, index), 0.0_dp, "g_sphi")
        call assert_close(metric(1, 3, index), 0.0_dp, "g_stheta")
        call assert_close( &
            metric(2, 2, index), cylindrical_radius**2, "g_phiphi")
        call assert_close(metric(2, 3, index), 0.0_dp, "g_phitheta")
        call assert_close( &
            metric(3, 3, index), (minor_scale*rho)**2, "g_thetatheta")
        call assert_close( &
            jacobian(index), &
            minor_scale**2*cylindrical_radius/2.0_dp, "Jacobian")
    end do
    print *, "test_boozer_chartmap_geometry: all checks passed"

contains

    subroutine assert_close(actual, expected, label)
        real(dp), intent(in) :: actual, expected
        character(*), intent(in) :: label

        if (abs(actual - expected) > tolerance*max(1.0_dp, abs(expected))) then
            print *, trim(label), ": expected", expected, "got", actual
            error stop
        end if
    end subroutine assert_close

end program test_boozer_chartmap_geometry
