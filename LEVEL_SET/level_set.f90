!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Implementation
!

subroutine level_set_step_2D(fun, h, x)
!
! Input:
!   fun ... subroutine that yields values and
!           gradients of the input function F:
!           fun(x, F, dF)
!   h ...   step length
!
! Input/Output:
!   x ...   position that is propagated along level line
!
  external :: fun
  real(8), intent(in) :: h
  real(8), intent(inout) :: x(2)

  real(8) :: xold(2), fvec(2)
  integer(4) :: info
  real(8), parameter :: tol = 1d-13

  xold = x
  call hybrd1(fcn, 2, x, fvec, tol, info)
  x = 2.0d0*x - xold  ! Update with midpoint: (xold + xnew)/2.0 = xmid
  !  xnew - xold = xmid*2.0 - xold*2.0 = h*xdot

  contains

  subroutine fcn(n, xmid, fvec2, iflag)
      integer(4), intent(in) :: n
      real(8), intent(in) :: xmid(n)
      real(8) :: fvec2(n)
      integer(4) :: iflag

      real(8) :: F, dF(2), jac

      if (iflag == 1) then
        call fun(xmid, F, dF)
        jac = abs(sqrt(sum(dF**2)))
        fvec2(1) = xmid(1) - xold(1) + 0.5d0*h*dF(2)/jac
        fvec2(2) = xmid(2) - xold(2) - 0.5d0*h*dF(1)/jac
      end if
    end subroutine fcn
end subroutine level_set_step_2D
