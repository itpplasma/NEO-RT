module mars_response_io
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    type, public :: mars_vector_harmonics_t
        integer :: mode_count = 0
        integer :: radial_count = 0
        integer :: toroidal_mode = 0
        integer, allocatable :: modes(:)
        complex(dp), allocatable :: c1(:, :)
        complex(dp), allocatable :: c2(:, :)
        complex(dp), allocatable :: c3(:, :)
    end type mars_vector_harmonics_t

    public :: read_mars_vector

contains

    subroutine read_mars_vector(path, vector)
        character(*), intent(in) :: path
        type(mars_vector_harmonics_t), intent(out) :: vector
        integer :: status, unit

        open (newunit=unit, file=path, status="old", action="read", iostat=status)
        if (status /= 0) error stop "Cannot open MARS vector file"
        call read_header(unit, vector)
        call allocate_vector(vector)
        call read_modes(unit, vector%modes)
        call read_components(unit, vector)
        close (unit)
    end subroutine read_mars_vector

    subroutine read_header(unit, vector)
        integer, intent(in) :: unit
        type(mars_vector_harmonics_t), intent(inout) :: vector
        real(dp) :: toroidal_mode
        integer :: status

        read (unit, *, iostat=status) vector%mode_count, vector%radial_count, toroidal_mode
        if (status /= 0) error stop "Invalid MARS vector header"
        if (vector%mode_count <= 0) error stop "Invalid MARS mode count"
        if (vector%radial_count <= 0) error stop "Invalid MARS radial count"
        vector%toroidal_mode = nint(toroidal_mode)
    end subroutine read_header

    subroutine allocate_vector(vector)
        type(mars_vector_harmonics_t), intent(inout) :: vector

        allocate (vector%modes(vector%mode_count))
        allocate (vector%c1(vector%radial_count, vector%mode_count))
        allocate (vector%c2(vector%radial_count, vector%mode_count))
        allocate (vector%c3(vector%radial_count, vector%mode_count))
    end subroutine allocate_vector

    subroutine read_modes(unit, modes)
        integer, intent(in) :: unit
        integer, intent(out) :: modes(:)
        real(dp) :: row(6)
        integer :: index, status

        do index = 1, size(modes)
            read (unit, *, iostat=status) row
            if (status /= 0) error stop "Invalid MARS mode row"
            modes(index) = nint(row(1))
        end do
    end subroutine read_modes

    subroutine read_components(unit, vector)
        integer, intent(in) :: unit
        type(mars_vector_harmonics_t), intent(inout) :: vector
        real(dp) :: row(6)
        integer :: mode_index, radial_index, status

        do mode_index = 1, vector%mode_count
            do radial_index = 1, vector%radial_count
                read (unit, *, iostat=status) row
                if (status /= 0) error stop "Invalid MARS component row"
                vector%c1(radial_index, mode_index) = cmplx(row(1), row(2), dp)
                vector%c2(radial_index, mode_index) = cmplx(row(3), row(4), dp)
                vector%c3(radial_index, mode_index) = cmplx(row(5), row(6), dp)
            end do
        end do
    end subroutine read_components

end module mars_response_io
