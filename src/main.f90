module neort_main

    implicit none

    character(1024) :: runname

contains

    subroutine main
        use do_magfie_mod, only: s
        use neort, only: check_magfie
        use neort_config, only: output_format
        use neort_datatypes, only: magfie_data_t
        use neort_lib, only: transport_data_t, mixture_result_t, neort_init, &
            neort_prepare_splines, neort_compute_at_s, neort_compute_no_splines, &
            neort_compute_mixture_at_s
        use neort_output, only: write_output, write_mixture_output
        use neort_species, only: mixture_enabled
        use util, only: files_exist

        character(len=*), parameter :: boozer_file = "in_file"
        character(len=*), parameter :: boozer_pert_file = "in_file_pert"
        character(len=*), parameter :: plasma_file = "plasma.in"
        character(len=*), parameter :: profile_file = "profile.in"

        type(magfie_data_t) :: magfie_data
        type(transport_data_t) :: transport_data
        type(mixture_result_t) :: mixture_result

        call get_command_argument(1, runname)

        call neort_init(trim(runname)//".in", boozer_file, boozer_pert_file)

        if (mixture_enabled) then
            call neort_compute_mixture_at_s(s, mixture_result)
            call write_mixture_output(mixture_result, trim(runname), output_format)
            return
        end if

        if (files_exist(plasma_file, profile_file)) then
            call neort_prepare_splines(plasma_file, profile_file)
            call neort_compute_at_s(s, transport_data)
        else
            call neort_compute_no_splines(transport_data)
        end if

        call check_magfie(magfie_data)  ! diagnostics

        call write_output(magfie_data, transport_data, trim(runname), output_format)
    end subroutine main

end module neort_main
