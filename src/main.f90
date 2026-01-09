module neort_main

    implicit none

    character(1024) :: runname

contains

    subroutine main
        use do_magfie_mod, only: s
        use neort, only: check_magfie, write_magfie_data_to_files, write_transport_data_to_files
        use neort_datatypes, only: magfie_data_t
        use neort_lib
        use util, only: files_exist

        character(len=*), parameter :: boozer_file = "in_file"
        character(len=*), parameter :: boozer_pert_file = "in_file_pert"
        character(len=*), parameter :: plasma_file = "plasma.in"
        character(len=*), parameter :: profile_file = "profile.in"

        type(magfie_data_t) :: magfie_data
        type(transport_data_t) :: transport_data

        call get_command_argument(1, runname)

        call neort_init(trim(runname), boozer_file, boozer_pert_file)

        if (files_exist(plasma_file, profile_file)) then
            call neort_prepare_splines(plasma_file, profile_file)
            call neort_compute_at_s(s, transport_data)
        else
            call neort_compute_no_splines(transport_data)
        end if

        call check_magfie(magfie_data)  ! diagnostics

        call write_magfie_data_to_files(magfie_data, trim(runname))
        call write_transport_data_to_files(transport_data, trim(runname))
    end subroutine main

end module neort_main
