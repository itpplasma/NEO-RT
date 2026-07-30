program neo_rt_jxb_profile
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use response_jxb, only: integrate_mars_profile

    implicit none

    character(1024) :: edge_path, profile_path, scale_text
    real(dp), allocatable :: density(:), edges(:)
    real(dp) :: scale

    call get_command_argument(1, edge_path)
    call get_command_argument(2, profile_path)
    call get_command_argument(3, scale_text)
    if (len_trim(profile_path) == 0) call usage

    scale = 1.0_dp
    if (len_trim(scale_text) > 0) read (scale_text, *) scale
    call read_first_column(trim(edge_path), edges)
    call read_second_column(trim(profile_path), density)
    if (size(edges) /= size(density) + 1) then
        write (error_unit, '(A,I0,A,I0)') "Expected one more edge than profile rows; got ", &
            size(edges), " and ", size(density)
        error stop
    end if
    write (*, '(ES24.16)') integrate_mars_profile(edges, density, scale)

contains

    subroutine usage
        write (error_unit, '(A)') &
            "usage: neo_rt_jxb_profile.x PROFEQ.OUT TORQUEJXB.OUT [SI_SCALE]"
        error stop
    end subroutine usage

    subroutine read_first_column(path, values)
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: values(:)
        integer :: index, number_of_rows, status, unit

        number_of_rows = count_rows(path)
        allocate (values(number_of_rows))
        open (newunit=unit, file=path, status="old", action="read")
        do index = 1, number_of_rows
            read (unit, *, iostat=status) values(index)
            if (status /= 0) error stop "Invalid edge-profile row"
        end do
        close (unit)
    end subroutine read_first_column

    subroutine read_second_column(path, values)
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: values(:)
        real(dp) :: coordinate
        integer :: index, number_of_rows, status, unit

        number_of_rows = count_rows(path)
        allocate (values(number_of_rows))
        open (newunit=unit, file=path, status="old", action="read")
        do index = 1, number_of_rows
            read (unit, *, iostat=status) coordinate, values(index)
            if (status /= 0) error stop "Invalid torque-profile row"
        end do
        close (unit)
    end subroutine read_second_column

    function count_rows(path) result(number_of_rows)
        character(*), intent(in) :: path
        character(4096) :: line
        integer :: number_of_rows, status, unit

        number_of_rows = 0
        open (newunit=unit, file=path, status="old", action="read")
        do
            read (unit, '(A)', iostat=status) line
            if (status /= 0) exit
            if (len_trim(line) > 0) number_of_rows = number_of_rows + 1
        end do
        close (unit)
    end function count_rows

end program neo_rt_jxb_profile
