program neo_rt_jxb_from_mars
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use mars_response_io, only: mars_vector_harmonics_t, read_mars_vector
    use response_jxb, only: integrate_mars_profile, reconstruct_mars_profile

    implicit none

    character(1024) :: b_path, edge_path, j_path, output_path
    type(mars_vector_harmonics_t) :: magnetic, current
    real(dp), allocatable :: edges(:), half_mesh(:), processed(:), raw(:)
    real(dp) :: edge_cutoff
    integer :: first_torque_cell, smoothing_passes

    call read_arguments( &
        edge_path, b_path, j_path, output_path, first_torque_cell, &
        smoothing_passes, edge_cutoff)
    call read_first_column(trim(edge_path), edges)
    call read_mars_vector(trim(b_path), magnetic)
    call read_mars_vector(trim(j_path), current)
    call check_vectors(magnetic, current, size(edges))
    half_mesh = 0.5_dp*(edges(:size(edges) - 1) + edges(2:))
    call reconstruct_mars_profile( &
        current%c1, current%c2, magnetic%c1, magnetic%c2, half_mesh, &
        first_torque_cell, smoothing_passes, edge_cutoff, raw, processed)
    call write_profiles(trim(output_path), half_mesh, raw, processed)
    write (*, '(A,ES24.16)') "raw native total: ", &
        integrate_mars_profile(edges, raw, 1.0_dp)
    write (*, '(A,ES24.16)') "processed native total: ", &
        integrate_mars_profile(edges, processed, 1.0_dp)

contains

    subroutine read_arguments( &
            edges, magnetic, current, output, first_cell, passes, cutoff)
        character(*), intent(out) :: edges, magnetic, current, output
        integer, intent(out) :: first_cell, passes
        real(dp), intent(out) :: cutoff
        character(64) :: argument

        if (command_argument_count() < 4) call usage
        call get_command_argument(1, edges)
        call get_command_argument(2, magnetic)
        call get_command_argument(3, current)
        call get_command_argument(4, output)
        first_cell = 1
        passes = 0
        cutoff = 1.0_dp
        if (command_argument_count() >= 5) then
            call get_command_argument(5, argument)
            read (argument, *) first_cell
        end if
        if (command_argument_count() >= 6) then
            call get_command_argument(6, argument)
            read (argument, *) passes
        end if
        if (command_argument_count() >= 7) then
            call get_command_argument(7, argument)
            read (argument, *) cutoff
        end if
    end subroutine read_arguments

    subroutine usage
        write (error_unit, '(A)') "usage: neo_rt_jxb_from_mars.x PROFEQ BPLASMA "// &
            "JPLASMA OUTPUT [NTORQ] [SMOOTHING_PASSES] [CTEDGE]"
        error stop
    end subroutine usage

    subroutine check_vectors(magnetic, current, edge_count)
        type(mars_vector_harmonics_t), intent(in) :: magnetic, current
        integer, intent(in) :: edge_count

        if (magnetic%radial_count < edge_count) error stop "BPLASMA radial grid too short"
        if (current%radial_count < edge_count) error stop "JPLASMA radial grid too short"
        if (magnetic%mode_count /= current%mode_count) error stop "MARS mode-count mismatch"
        if (magnetic%toroidal_mode /= current%toroidal_mode) error stop "MARS n mismatch"
        if (any(magnetic%modes /= current%modes)) error stop "MARS m mismatch"
    end subroutine check_vectors

    subroutine read_first_column(path, values)
        character(*), intent(in) :: path
        real(dp), allocatable, intent(out) :: values(:)
        character(4096) :: line
        integer :: count, index, status, unit

        count = 0
        open (newunit=unit, file=path, status="old", action="read")
        do
            read (unit, '(A)', iostat=status) line
            if (status /= 0) exit
            if (len_trim(line) > 0) count = count + 1
        end do
        rewind (unit)
        allocate (values(count))
        do index = 1, count
            read (unit, *, iostat=status) values(index)
            if (status /= 0) error stop "Invalid MARS radial-grid row"
        end do
        close (unit)
    end subroutine read_first_column

    subroutine write_profiles(path, radius, raw, processed)
        character(*), intent(in) :: path
        real(dp), intent(in) :: radius(:), raw(:), processed(:)
        integer :: index, unit

        open (newunit=unit, file=path, status="replace", action="write")
        do index = 1, size(radius)
            write (unit, '(3(ES24.16,1X))') radius(index), raw(index), processed(index)
        end do
        close (unit)
    end subroutine write_profiles

end program neo_rt_jxb_from_mars
