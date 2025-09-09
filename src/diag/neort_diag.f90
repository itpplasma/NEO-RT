program neort_diag
  use iso_fortran_env, only: int32
  use diag_bounce_nonlin, only: run_bounce_nonlin_diag
  implicit none
  character(len=256) :: diag, runname
  call get_command_argument(1, diag)
  call get_command_argument(2, runname)
  if (len_trim(diag) == 0 .or. len_trim(runname) == 0) then
     print *, "Usage: neo_rt_diag.x <diagnostic> <runname>"
     print *, "Diagnostics: bounce_nonlin"
     stop 1
  end if

  select case (trim(adjustl(diag)))
  case ("bounce_nonlin")
     call run_bounce_nonlin_diag(trim(adjustl(runname)))
  case default
     print *, "Unknown diagnostic:", trim(diag)
     stop 2
  end select

end program neort_diag
