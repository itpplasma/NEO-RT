module boozer_vector_io
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use boozer_ampere, only: boozer_harmonic_jacobian_current
    use netcdf
    implicit none
    private

    type, public :: boozer_vector_harmonics_t
        integer :: radial_count = 0
        integer :: mode_count = 0
        integer :: toroidal_mode = 0
        real(dp), allocatable :: s(:)
        integer, allocatable :: modes(:)
        complex(dp), allocatable :: b_covariant(:, :, :)
    end type boozer_vector_harmonics_t

    public :: read_boozer_vector_harmonics
    public :: compute_boozer_weighted_current

contains

    subroutine read_boozer_vector_harmonics(path, field)
        character(*), intent(in) :: path
        type(boozer_vector_harmonics_t), intent(out) :: field
        integer :: ncid

        call check_netcdf(nf90_open(trim(path), nf90_nowrite, ncid), &
            "open "//trim(path))
        call read_dimensions(ncid, field)
        allocate (field%s(field%radial_count), field%modes(field%mode_count))
        allocate (field%b_covariant(3, field%radial_count, field%mode_count))
        call read_real_vector(ncid, "s", field%s)
        call read_integer_vector(ncid, "m", field%modes)
        call check_netcdf(nf90_get_att( &
            ncid, nf90_global, "toroidal_mode", field%toroidal_mode), &
            "attribute toroidal_mode")
        call require_attribute(ncid, "coordinate_order", "s,phi,theta")
        call require_attribute(ncid, "component_variance", "covariant")
        call require_attribute(ncid, "radial_coordinate", "s_tor")
        call require_attribute(ncid, "magnetic_component_units", "T m")
        call require_attribute(ncid, "fourier_convention", &
            "B_hat(s) exp(i*(n*phi+m*theta))")
        call read_complex_component(ncid, "B_s", field%b_covariant(1, :, :))
        call read_complex_component(ncid, "B_phi", field%b_covariant(2, :, :))
        call read_complex_component(ncid, "B_theta", field%b_covariant(3, :, :))
        call check_netcdf(nf90_close(ncid), "close "//trim(path))
        call validate_field(field)
    end subroutine read_boozer_vector_harmonics

    subroutine compute_boozer_weighted_current(field, mu0, weighted_current)
        type(boozer_vector_harmonics_t), intent(in) :: field
        real(dp), intent(in) :: mu0
        complex(dp), allocatable, intent(out) :: weighted_current(:, :, :)
        complex(dp), allocatable :: derivative(:, :)
        integer :: radial_index, mode_index

        if (mu0 == 0.0_dp) error stop "zero magnetic permeability"
        allocate (weighted_current(3, field%radial_count, field%mode_count))
        allocate (derivative(3, field%radial_count))
        do mode_index = 1, field%mode_count
            derivative = cmplx(0.0_dp, 0.0_dp, dp)
            call differentiate_radial( &
                field%s, field%b_covariant(2, :, mode_index), derivative(2, :))
            call differentiate_radial( &
                field%s, field%b_covariant(3, :, mode_index), derivative(3, :))
            do radial_index = 1, field%radial_count
                call boozer_harmonic_jacobian_current( &
                    mu0, field%modes(mode_index), field%toroidal_mode, &
                    field%b_covariant(:, radial_index, mode_index), &
                    derivative(:, radial_index), &
                    weighted_current(:, radial_index, mode_index))
            end do
        end do
    end subroutine compute_boozer_weighted_current

    subroutine differentiate_radial(x, values, derivative)
        real(dp), intent(in) :: x(:)
        complex(dp), intent(in) :: values(:)
        complex(dp), intent(out) :: derivative(:)
        integer :: center, first, last

        if (size(x) /= size(values) .or. size(x) /= size(derivative)) then
            error stop "radial derivative shape mismatch"
        end if
        if (size(x) < 3) error stop "radial derivative needs at least three points"
        do center = 1, size(x)
            first = max(1, min(center - 1, size(x) - 2))
            last = first + 2
            derivative(center) = derivative_of_interpolant( &
                x(first:last), values(first:last), x(center))
        end do
    end subroutine differentiate_radial

    pure function derivative_of_interpolant(x, values, target) result(derivative)
        real(dp), intent(in) :: x(3), target
        complex(dp), intent(in) :: values(3)
        complex(dp) :: derivative
        real(dp) :: weight
        integer :: j, k, l

        derivative = cmplx(0.0_dp, 0.0_dp, dp)
        do j = 1, 3
            weight = 0.0_dp
            do k = 1, 3
                if (k == j) cycle
                l = 6 - j - k
                weight = weight + (target - x(l))/ &
                    ((x(j) - x(k))*(x(j) - x(l)))
            end do
            derivative = derivative + weight*values(j)
        end do
    end function derivative_of_interpolant

    subroutine read_dimensions(ncid, field)
        integer, intent(in) :: ncid
        type(boozer_vector_harmonics_t), intent(inout) :: field
        integer :: dimid

        call check_netcdf(nf90_inq_dimid(ncid, "s", dimid), "dimension s")
        call check_netcdf(nf90_inquire_dimension( &
            ncid, dimid, len=field%radial_count), "length s")
        call check_netcdf(nf90_inq_dimid(ncid, "mode", dimid), "dimension mode")
        call check_netcdf(nf90_inquire_dimension( &
            ncid, dimid, len=field%mode_count), "length mode")
    end subroutine read_dimensions

    subroutine read_real_vector(ncid, name, values)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        real(dp), intent(out) :: values(:)
        integer :: varid

        call check_netcdf(nf90_inq_varid(ncid, name, varid), "variable "//name)
        call check_netcdf(nf90_get_var(ncid, varid, values), "read "//name)
    end subroutine read_real_vector

    subroutine read_integer_vector(ncid, name, values)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        integer, intent(out) :: values(:)
        integer :: varid

        call check_netcdf(nf90_inq_varid(ncid, name, varid), "variable "//name)
        call check_netcdf(nf90_get_var(ncid, varid, values), "read "//name)
    end subroutine read_integer_vector

    subroutine read_complex_component(ncid, prefix, values)
        integer, intent(in) :: ncid
        character(*), intent(in) :: prefix
        complex(dp), intent(out) :: values(:, :)
        real(dp), allocatable :: real_part(:, :), imaginary_part(:, :)

        allocate (real_part(size(values, 1), size(values, 2)))
        allocate (imaginary_part(size(values, 1), size(values, 2)))
        call read_real_matrix(ncid, prefix//"_real", real_part)
        call read_real_matrix(ncid, prefix//"_imag", imaginary_part)
        values = cmplx(real_part, imaginary_part, dp)
    end subroutine read_complex_component

    subroutine read_real_matrix(ncid, name, values)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name
        real(dp), intent(out) :: values(:, :)
        integer :: varid

        call check_netcdf(nf90_inq_varid(ncid, name, varid), "variable "//name)
        call check_netcdf(nf90_get_var(ncid, varid, values), "read "//name)
    end subroutine read_real_matrix

    subroutine require_attribute(ncid, name, expected)
        integer, intent(in) :: ncid
        character(*), intent(in) :: name, expected
        character(128) :: actual

        actual = ""
        call check_netcdf(nf90_get_att(ncid, nf90_global, name, actual), &
            "attribute "//name)
        if (trim(actual) /= expected) error stop "unsupported Boozer vector attribute"
    end subroutine require_attribute

    subroutine validate_field(field)
        type(boozer_vector_harmonics_t), intent(in) :: field
        integer :: index

        if (field%radial_count < 3) error stop "Boozer vector radial grid too short"
        if (field%mode_count < 1) error stop "Boozer vector has no modes"
        if (any(field%s(2:) <= field%s(:field%radial_count - 1))) then
            error stop "Boozer vector radial grid is not increasing"
        end if
        do index = 1, field%mode_count
            if (count(field%modes == field%modes(index)) /= 1) then
                error stop "Boozer vector modes are not unique"
            end if
        end do
    end subroutine validate_field

    subroutine check_netcdf(status, context)
        integer, intent(in) :: status
        character(*), intent(in) :: context

        if (status /= nf90_noerr) then
            print *, trim(context)//": "//trim(nf90_strerror(status))
            error stop "Boozer vector NetCDF error"
        end if
    end subroutine check_netcdf

end module boozer_vector_io
