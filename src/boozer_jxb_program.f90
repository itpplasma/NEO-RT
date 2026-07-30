program neo_rt_boozer_jxb
    use, intrinsic :: iso_fortran_env, only: dp => real64, error_unit
    use boozer_chartmap_geometry, only: boozer_chartmap_geometry_t, &
        read_boozer_chartmap_geometry
    use boozer_jxb_profile, only: compute_boozer_jxb_profile
    use boozer_vector_io, only: boozer_vector_harmonics_t, &
        read_boozer_vector_harmonics
    implicit none

    character(1024) :: axisymmetric_path, perturbation_path, output_path
    type(boozer_chartmap_geometry_t) :: geometry
    type(boozer_vector_harmonics_t) :: field
    real(dp), allocatable :: torque_density(:)
    real(dp) :: mu0

    call read_arguments( &
        axisymmetric_path, perturbation_path, output_path, mu0)
    call read_boozer_chartmap_geometry(trim(axisymmetric_path), geometry)
    call read_boozer_vector_harmonics(trim(perturbation_path), field)
    call compute_boozer_jxb_profile(geometry, field, mu0, torque_density)
    call write_profile(trim(output_path), field%s, torque_density)
    write (*, '(A,I0,A)') &
        "computed Boozer jxb torque density on ", field%radial_count, " surfaces"

contains

    subroutine read_arguments(axisymmetric, perturbation, output, permeability)
        character(*), intent(out) :: axisymmetric, perturbation, output
        real(dp), intent(out) :: permeability
        character(64) :: argument
        real(dp), parameter :: pi = acos(-1.0_dp)

        if (command_argument_count() < 3) then
            write (error_unit, '(A)') &
                "usage: neo_rt_boozer_jxb.x AXISYMMETRIC.nc PERTURBATION.nc "// &
                "OUTPUT [MU0]"
            error stop
        end if
        call get_command_argument(1, axisymmetric)
        call get_command_argument(2, perturbation)
        call get_command_argument(3, output)
        permeability = 4.0e-7_dp*pi
        if (command_argument_count() >= 4) then
            call get_command_argument(4, argument)
            read (argument, *) permeability
        end if
    end subroutine read_arguments

    subroutine write_profile(path, s, density)
        character(*), intent(in) :: path
        real(dp), intent(in) :: s(:), density(:)
        integer :: radial_index, unit

        if (size(s) /= size(density)) error stop "Boozer jxb profile size mismatch"
        open (newunit=unit, file=path, status="replace", action="write")
        write (unit, '(A)') "# s_tor rho_tor dT_phi_ds_N_m"
        do radial_index = 1, size(s)
            write (unit, '(3(ES24.16,1X))') &
                s(radial_index), sqrt(s(radial_index)), density(radial_index)
        end do
        close (unit)
    end subroutine write_profile

end program neo_rt_boozer_jxb
