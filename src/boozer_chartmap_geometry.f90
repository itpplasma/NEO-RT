module boozer_chartmap_geometry
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use netcdf
    implicit none
    private

    type, public :: boozer_chartmap_geometry_t
        integer :: radial_count = 0
        integer :: theta_count = 0
        integer :: zeta_count = 0
        integer :: field_periods = 1
        real(dp), allocatable :: rho(:)
        real(dp), allocatable :: theta(:)
        real(dp), allocatable :: zeta(:)
        real(dp), allocatable :: position(:, :, :, :)
        real(dp), allocatable :: position_theta(:, :, :, :)
        real(dp), allocatable :: position_zeta(:, :, :, :)
    end type boozer_chartmap_geometry_t

    public :: evaluate_boozer_metric
    public :: read_boozer_chartmap_geometry

contains

    subroutine read_boozer_chartmap_geometry(path, geometry)
        character(*), intent(in) :: path
        type(boozer_chartmap_geometry_t), intent(out) :: geometry
        integer :: ncid, component

        call check_netcdf(nf90_open(trim(path), nf90_nowrite, ncid), &
            "open "//trim(path))
        call read_dimension(ncid, "rho", geometry%radial_count)
        call read_dimension(ncid, "theta", geometry%theta_count)
        call read_dimension(ncid, "zeta", geometry%zeta_count)
        if (geometry%radial_count < 3) error stop "chartmap radial grid too short"
        if (geometry%theta_count < 3) error stop "chartmap theta grid too short"
        if (geometry%zeta_count < 3) error stop "chartmap zeta grid too short"
        allocate (geometry%rho(geometry%radial_count))
        allocate (geometry%theta(geometry%theta_count))
        allocate (geometry%zeta(geometry%zeta_count))
        allocate (geometry%position( &
            3, geometry%radial_count, geometry%theta_count, geometry%zeta_count))
        allocate (geometry%position_theta, mold=geometry%position)
        allocate (geometry%position_zeta, mold=geometry%position)
        call read_real_vector(ncid, "rho", geometry%rho)
        call read_real_vector(ncid, "theta", geometry%theta)
        call read_real_vector(ncid, "zeta", geometry%zeta)
        call read_position(ncid, "x", geometry%position(1, :, :, :))
        call read_position(ncid, "y", geometry%position(2, :, :, :))
        call read_position(ncid, "z", geometry%position(3, :, :, :))
        call read_integer_scalar( &
            ncid, "num_field_periods", geometry%field_periods)
        call require_attribute(ncid, "rho_convention", "rho_tor")
        call require_attribute(ncid, "zeta_convention", "boozer")
        call require_variable_attribute(ncid, "x", "units", "cm")
        call require_variable_attribute(ncid, "y", "units", "cm")
        call require_variable_attribute(ncid, "z", "units", "cm")
        call check_netcdf(nf90_close(ncid), "close "//trim(path))

        geometry%position = 0.01_dp*geometry%position
        call validate_grid(geometry%rho, "rho", 0.0_dp)
        call validate_periodic_grid(geometry%theta, "theta", 2.0_dp*acos(-1.0_dp))
        call validate_periodic_grid(geometry%zeta, "zeta", &
            2.0_dp*acos(-1.0_dp)/real(geometry%field_periods, dp))
        do component = 1, 3
            call differentiate_angles(geometry, component)
        end do
    end subroutine read_boozer_chartmap_geometry

    subroutine evaluate_boozer_metric(geometry, s, metric, jacobian)
        type(boozer_chartmap_geometry_t), intent(in) :: geometry
        real(dp), intent(in) :: s
        real(dp), intent(out) :: metric(3, 3, geometry%theta_count)
        real(dp), intent(out) :: jacobian(geometry%theta_count)
        real(dp) :: rho, e_s(3), e_phi(3), e_theta(3)
        integer :: component, theta_index

        if (s <= 0.0_dp) error stop "chartmap metric is singular at s <= 0"
        rho = sqrt(s)
        if (rho < geometry%rho(1) .or. &
            rho > geometry%rho(geometry%radial_count)) then
            error stop "requested surface is outside chartmap geometry"
        end if
        do theta_index = 1, geometry%theta_count
            do component = 1, 3
                e_s(component) = interpolate_radial_derivative( &
                    geometry%rho, &
                    geometry%position(component, :, theta_index, 1), rho)/ &
                    (2.0_dp*rho)
                e_phi(component) = interpolate_radial_value( &
                    geometry%rho, &
                    geometry%position_zeta(component, :, theta_index, 1), rho)
                e_theta(component) = interpolate_radial_value( &
                    geometry%rho, &
                    geometry%position_theta(component, :, theta_index, 1), rho)
            end do
            metric(1, 1, theta_index) = dot_product(e_s, e_s)
            metric(1, 2, theta_index) = dot_product(e_s, e_phi)
            metric(1, 3, theta_index) = dot_product(e_s, e_theta)
            metric(2, 1, theta_index) = metric(1, 2, theta_index)
            metric(2, 2, theta_index) = dot_product(e_phi, e_phi)
            metric(2, 3, theta_index) = dot_product(e_phi, e_theta)
            metric(3, 1, theta_index) = metric(1, 3, theta_index)
            metric(3, 2, theta_index) = metric(2, 3, theta_index)
            metric(3, 3, theta_index) = dot_product(e_theta, e_theta)
            jacobian(theta_index) = dot_product(e_s, cross(e_phi, e_theta))
        end do
        if (any(abs(jacobian) <= tiny(1.0_dp))) then
            error stop "chartmap metric has a singular Jacobian"
        end if
    end subroutine evaluate_boozer_metric

    subroutine differentiate_angles(geometry, component)
        type(boozer_chartmap_geometry_t), intent(inout) :: geometry
        integer, intent(in) :: component
        real(dp) :: theta_values(geometry%theta_count)
        real(dp) :: zeta_values(geometry%zeta_count)
        integer :: radial_index, theta_index, zeta_index

        do radial_index = 1, geometry%radial_count
            do zeta_index = 1, geometry%zeta_count
                theta_values = geometry%position( &
                    component, radial_index, :, zeta_index)
                geometry%position_theta( &
                    component, radial_index, :, zeta_index) = &
                    periodic_spectral_derivative(theta_values, 2.0_dp*acos(-1.0_dp))
            end do
            do theta_index = 1, geometry%theta_count
                zeta_values = geometry%position( &
                    component, radial_index, theta_index, :)
                geometry%position_zeta( &
                    component, radial_index, theta_index, :) = &
                    periodic_spectral_derivative( &
                    zeta_values, &
                    2.0_dp*acos(-1.0_dp)/real(geometry%field_periods, dp))
            end do
        end do
    end subroutine differentiate_angles

    pure function periodic_spectral_derivative(values, period) result(derivative)
        real(dp), intent(in) :: values(:), period
        real(dp) :: derivative(size(values))
        real(dp) :: angle, wave_number
        integer :: destination, mode, source, count, maximum_mode

        count = size(values)
        maximum_mode = (count - 1)/2
        derivative = 0.0_dp
        do destination = 1, count
            do source = 1, count
                angle = 2.0_dp*acos(-1.0_dp)* &
                    real(destination - source, dp)/real(count, dp)
                do mode = 1, maximum_mode
                    wave_number = 2.0_dp*acos(-1.0_dp)* &
                        real(mode, dp)/period
                    derivative(destination) = derivative(destination) - &
                        2.0_dp*wave_number*sin(real(mode, dp)*angle)* &
                        values(source)/real(count, dp)
                end do
            end do
        end do
    end function periodic_spectral_derivative

    pure function interpolate_radial_value(x, values, target) result(value)
        real(dp), intent(in) :: x(:), values(:), target
        real(dp) :: value
        integer :: first

        first = interpolation_window(x, target)
        value = interpolant_value(x(first:first + 2), values(first:first + 2), target)
    end function interpolate_radial_value

    pure function interpolate_radial_derivative(x, values, target) result(value)
        real(dp), intent(in) :: x(:), values(:), target
        real(dp) :: value
        integer :: first

        first = interpolation_window(x, target)
        value = interpolant_derivative( &
            x(first:first + 2), values(first:first + 2), target)
    end function interpolate_radial_derivative

    pure function interpolation_window(x, target) result(first)
        real(dp), intent(in) :: x(:), target
        integer :: first, nearest

        nearest = minloc(abs(x - target), dim=1)
        first = max(1, min(nearest - 1, size(x) - 2))
    end function interpolation_window

    pure function interpolant_value(x, values, target) result(value)
        real(dp), intent(in) :: x(3), values(3), target
        real(dp) :: value, weight
        integer :: j, k

        value = 0.0_dp
        do j = 1, 3
            weight = 1.0_dp
            do k = 1, 3
                if (k /= j) weight = weight*(target - x(k))/(x(j) - x(k))
            end do
            value = value + weight*values(j)
        end do
    end function interpolant_value

    pure function interpolant_derivative(x, values, target) result(value)
        real(dp), intent(in) :: x(3), values(3), target
        real(dp) :: value, weight
        integer :: j, k, l

        value = 0.0_dp
        do j = 1, 3
            weight = 0.0_dp
            do k = 1, 3
                if (k == j) cycle
                l = 6 - j - k
                weight = weight + (target - x(l))/ &
                    ((x(j) - x(k))*(x(j) - x(l)))
            end do
            value = value + weight*values(j)
        end do
    end function interpolant_derivative

    pure function cross(left, right) result(product)
        real(dp), intent(in) :: left(3), right(3)
        real(dp) :: product(3)

        product = [ &
            left(2)*right(3) - left(3)*right(2), &
            left(3)*right(1) - left(1)*right(3), &
            left(1)*right(2) - left(2)*right(1)]
    end function cross

    subroutine read_dimension(ncid, name, length)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        integer, intent(out) :: length
        integer :: dimid

        call check_netcdf(nf90_inq_dimid(ncid, name, dimid), "dimension "//name)
        call check_netcdf( &
            nf90_inquire_dimension(ncid, dimid, len=length), "length "//name)
    end subroutine read_dimension

    subroutine read_real_vector(ncid, name, values)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        real(dp), intent(out) :: values(:)
        integer :: varid

        call check_netcdf(nf90_inq_varid(ncid, name, varid), "variable "//name)
        call check_netcdf(nf90_get_var(ncid, varid, values), "read "//name)
    end subroutine read_real_vector

    subroutine read_position(ncid, name, values)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        real(dp), intent(out) :: values(:, :, :)
        integer :: varid

        call check_netcdf(nf90_inq_varid(ncid, name, varid), "variable "//name)
        call check_netcdf(nf90_get_var(ncid, varid, values), "read "//name)
    end subroutine read_position

    subroutine read_integer_scalar(ncid, name, value)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        integer, intent(out) :: value
        integer :: varid

        call check_netcdf(nf90_inq_varid(ncid, name, varid), "variable "//name)
        call check_netcdf(nf90_get_var(ncid, varid, value), "read "//name)
    end subroutine read_integer_scalar

    subroutine require_attribute(ncid, name, expected)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name, expected
        character(128) :: actual

        actual = ""
        call check_netcdf( &
            nf90_get_att(ncid, nf90_global, name, actual), "attribute "//name)
        if (trim(actual) /= expected) error stop "unsupported chartmap attribute"
    end subroutine require_attribute

    subroutine require_variable_attribute(ncid, variable, name, expected)
        integer, intent(in) :: ncid
        character(*), intent(in) :: variable, name, expected
        character(128) :: actual
        integer :: varid

        actual = ""
        call check_netcdf( &
            nf90_inq_varid(ncid, variable, varid), "variable "//variable)
        call check_netcdf( &
            nf90_get_att(ncid, varid, name, actual), &
            "attribute "//variable//":"//name)
        if (trim(actual) /= expected) error stop "unsupported chartmap variable unit"
    end subroutine require_variable_attribute

    subroutine validate_grid(values, name, minimum)
        real(dp), intent(in) :: values(:), minimum
        character(*), intent(in) :: name

        if (values(1) < minimum) error stop name//" grid starts below its domain"
        if (any(values(2:) <= values(:size(values) - 1))) then
            error stop name//" grid is not increasing"
        end if
    end subroutine validate_grid

    subroutine validate_periodic_grid(values, name, period)
        real(dp), intent(in) :: values(:), period
        character(*), intent(in) :: name
        real(dp) :: spacing
        integer :: component

        call validate_grid(values, name, 0.0_dp)
        spacing = period/real(size(values), dp)
        if (abs(values(1)) > 100.0_dp*epsilon(1.0_dp)) then
            error stop name//" periodic grid does not start at zero"
        end if
        if (maxval(abs(values - spacing* &
            [(real(component, dp), component=0, size(values) - 1)])) > &
            1.0e-10_dp*period) then
            error stop name//" periodic grid is not uniform"
        end if
    end subroutine validate_periodic_grid

    subroutine check_netcdf(status, context)
        integer, intent(in) :: status
        character(*), intent(in) :: context

        if (status /= nf90_noerr) then
            print *, trim(context)//": "//trim(nf90_strerror(status))
            error stop "Boozer chartmap geometry NetCDF error"
        end if
    end subroutine check_netcdf

end module boozer_chartmap_geometry
