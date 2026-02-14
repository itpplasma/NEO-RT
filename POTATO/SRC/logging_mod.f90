module logging_mod
    use iso_fortran_env, only: output_unit
    implicit none

    integer, save :: log_unit = -1
    logical, save :: logging_enabled = .false.

contains

    subroutine init_logging(filename)
        character(len=*), intent(in) :: filename
        integer :: ierr

        open(newunit=log_unit, file=filename, status='replace', action='write', iostat=ierr)
        if (ierr == 0) then
            logging_enabled = .true.
            write(log_unit, '(A)') '=== POTATO Log Started ==='
            flush(log_unit)
        else
            log_unit = -1
            logging_enabled = .false.
        endif
    end subroutine init_logging

    subroutine close_logging()
        if (logging_enabled) then
            write(log_unit, '(A)') '=== POTATO Log Ended ==='
            close(log_unit)
            logging_enabled = .false.
            log_unit = -1
        endif
    end subroutine close_logging

    subroutine flush_log()
        if (logging_enabled) flush(log_unit)
    end subroutine flush_log

    function get_log_unit() result(unit_num)
        integer :: unit_num
        unit_num = log_unit
    end function get_log_unit

    function is_logging_enabled() result(enabled)
        logical :: enabled
        enabled = logging_enabled
    end function is_logging_enabled

    subroutine tee_message(msg)
        character(len=*), intent(in) :: msg

        write(output_unit, '(A)') msg
        if (logging_enabled) then
            write(log_unit, '(A)') msg
            flush(log_unit)
        endif
    end subroutine tee_message

end module logging_mod
