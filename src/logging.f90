module logger
  implicit none
  private
  integer, public, parameter :: LOG_TRACE=0, LOG_DEBUG=1, LOG_INFO=2, LOG_WARNING=3, LOG_ERROR=4
  integer :: current_level = LOG_INFO
  public :: set_log_level, get_log_level, trace, debug, info, warning, error
contains
  subroutine set_log_level(level)
    integer, intent(in) :: level
    current_level = max(LOG_TRACE, min(LOG_ERROR, level))
  end subroutine set_log_level

  function get_log_level() result(level)
    integer :: level
    level = current_level
  end function get_log_level

  subroutine trace(msg)
    character(*), intent(in) :: msg
    if (current_level <= LOG_TRACE) then
      write(6,'(A)') '[TRACE] ' // trim(msg)
    end if
  end subroutine trace

  subroutine debug(msg)
    character(*), intent(in) :: msg
    if (current_level <= LOG_DEBUG) then
      write(6,'(A)') '[DEBUG] ' // trim(msg)
    end if
  end subroutine debug

  subroutine info(msg)
    character(*), intent(in) :: msg
    if (current_level <= LOG_INFO) then
      write(6,'(A)') '[INFO ] ' // trim(msg)
    end if
  end subroutine info

  subroutine warning(msg)
    character(*), intent(in) :: msg
    if (current_level <= LOG_WARNING) then
      write(0,'(A)') '[WARN ] ' // trim(msg)
    end if
  end subroutine warning

  subroutine error(msg)
    character(*), intent(in) :: msg
    write(0,'(A)') '[ERROR] ' // trim(msg)
    error stop
  end subroutine error
end module logger

