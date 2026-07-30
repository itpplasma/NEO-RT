program neo_rt_boozer_current
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use boozer_vector_io, only: boozer_vector_harmonics_t, &
        compute_boozer_weighted_current, read_boozer_vector_harmonics
    implicit none

    character(1024) :: input_path, output_path
    type(boozer_vector_harmonics_t) :: field
    complex(dp), allocatable :: current(:, :, :)
    real(dp) :: mu0

    call read_arguments(input_path, output_path, mu0)
    call read_boozer_vector_harmonics(trim(input_path), field)
    call compute_boozer_weighted_current(field, mu0, current)
    call write_current(trim(output_path), field, current)
    write (*, '(A,I0,A,I0)') "computed weighted current on ", &
        field%radial_count, " surfaces for ", field%mode_count

contains

    subroutine read_arguments(input, output, permeability)
        character(*), intent(out) :: input, output
        real(dp), intent(out) :: permeability
        character(64) :: argument
        real(dp), parameter :: pi = acos(-1.0_dp)

        if (command_argument_count() < 2) then
            write (error_unit, '(A)') &
                "usage: neo_rt_boozer_current.x VECTOR.nc OUTPUT [MU0]"
            error stop
        end if
        call get_command_argument(1, input)
        call get_command_argument(2, output)
        permeability = 4.0e-7_dp*pi
        if (command_argument_count() >= 3) then
            call get_command_argument(3, argument)
            read (argument, *) permeability
        end if
    end subroutine read_arguments

    subroutine write_current(path, vector, values)
        character(*), intent(in) :: path
        type(boozer_vector_harmonics_t), intent(in) :: vector
        complex(dp), intent(in) :: values(:, :, :)
        integer :: mode_index, radial_index, unit

        open (newunit=unit, file=path, status="replace", action="write")
        write (unit, '(A)') &
            "# s m Re(JJ_s) Im(JJ_s) Re(JJ_phi) Im(JJ_phi) Re(JJ_theta) Im(JJ_theta)"
        do mode_index = 1, vector%mode_count
            do radial_index = 1, vector%radial_count
                write (unit, '(ES24.16,1X,I0,1X,6(ES24.16,1X))') &
                    vector%s(radial_index), vector%modes(mode_index), &
                    real(values(1, radial_index, mode_index), dp), &
                    aimag(values(1, radial_index, mode_index)), &
                    real(values(2, radial_index, mode_index), dp), &
                    aimag(values(2, radial_index, mode_index)), &
                    real(values(3, radial_index, mode_index), dp), &
                    aimag(values(3, radial_index, mode_index))
            end do
        end do
        close (unit)
    end subroutine write_current

end program neo_rt_boozer_current
