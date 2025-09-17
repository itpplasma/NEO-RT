module logger
  implicit none
  private
  ! Expose levels 1..5 as requested: 1=ERROR, 2=WARNING, 3=INFO, 4=DEBUG, 5=TRACE
  integer, public, parameter :: LVL_ERROR=1, LVL_WARNING=2, LVL_INFO=3, LVL_DEBUG=4, &
  LVL_TRACE=5
  ! Backward-compatible aliases used in imports elsewhere
  integer, public, parameter :: LOG_ERROR=LVL_ERROR, LOG_WARNING=LVL_WARNING, &
  LOG_INFO=LVL_INFO, LOG_DEBUG=LVL_DEBUG, LOG_TRACE=LVL_TRACE
  integer :: current_level = LVL_INFO
  public :: set_log_level, get_log_level, trace, debug, info, warning, error
contains
  subroutine set_log_level(level)
    integer, intent(in) :: level
    current_level = max(LVL_ERROR, min(LVL_TRACE, level))
  end subroutine set_log_level

  function get_log_level() result(level)
    integer :: level
    level = current_level
  end function get_log_level

  subroutine trace(msg)
    character(*), intent(in) :: msg
    if (current_level >= LVL_TRACE) write(6,'(A)') '[TRACE] ' // trim(msg)
  end subroutine trace

  subroutine debug(msg)
    character(*), intent(in) :: msg
    if (current_level >= LVL_DEBUG) write(6,'(A)') '[DEBUG] ' // trim(msg)
  end subroutine debug

  subroutine info(msg)
    character(*), intent(in) :: msg
    if (current_level >= LVL_INFO) write(6,'(A)') '[INFO ] ' // trim(msg)
  end subroutine info

  subroutine warning(msg)
    character(*), intent(in) :: msg
    if (current_level >= LVL_WARNING) write(0,'(A)') '[WARN ] ' // trim(msg)
  end subroutine warning

  subroutine error(msg)
    character(*), intent(in) :: msg
    write(0,'(A)') '[ERROR] ' // trim(msg)
    error stop
  end subroutine error
end module logger
