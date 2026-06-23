module logging_mod
    use iso_fortran_env, only: output_unit
    implicit none

    integer, save :: log_unit = -1
    logical, save :: logging_enabled = .false.

contains

    function timestamp() result(ts)
        character(len=23) :: ts
        character(len=8)  :: date_str
        character(len=10) :: time_str

        call date_and_time(date=date_str, time=time_str)
        write(ts, '(A4,"-",A2,"-",A2," ",A2,":",A2,":",A2,".",A3)') &
            date_str(1:4), date_str(5:6), date_str(7:8), &
            time_str(1:2), time_str(3:4), time_str(5:6), time_str(8:10)
    end function timestamp

    subroutine init_logging(filename)
        character(len=*), intent(in) :: filename
        integer :: ierr

        open(newunit=log_unit, file=filename, status='replace', action='write', iostat=ierr)
        if (ierr == 0) then
            logging_enabled = .true.
            write(log_unit, '(A,": ",A)') timestamp(), '=== POTATO Log Started ==='
            flush(log_unit)
        else
            log_unit = -1
            logging_enabled = .false.
        endif
    end subroutine init_logging

    subroutine close_logging()
        if (logging_enabled) then
            write(log_unit, '(A,": ",A)') timestamp(), '=== POTATO Log Ended ==='
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

    subroutine log_message(msg)
        character(len=*), intent(in) :: msg

        if (logging_enabled) then
            write(log_unit, '(A,": ",A)') timestamp(), msg
            flush(log_unit)
        endif
    end subroutine log_message

    subroutine tee_message(msg)
        character(len=*), intent(in) :: msg

        write(output_unit, '(A)') msg
        if (logging_enabled) then
            write(log_unit, '(A,": ",A)') timestamp(), msg
            flush(log_unit)
        endif
    end subroutine tee_message

end module logging_mod
