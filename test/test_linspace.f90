program test_linspace
  use, intrinsic :: iso_fortran_env, only: dp => real64
  implicit none
  real(dp) :: arr(5)
  interface
    subroutine linspace(a, b, cnt, out)
      import dp
      real(dp), intent(in) :: a, b
      integer, intent(in) :: cnt
      real(dp), intent(out) :: out(cnt)
    end subroutine linspace
  end interface
  call linspace(0.0_dp, 1.0_dp, 5, arr)
  if (abs(arr(3)-0.5_dp) > 1d-12) then
     error stop "linspace midpoint incorrect"
  endif
end program test_linspace
