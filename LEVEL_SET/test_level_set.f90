!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! For testing
!

program test_level_set
  integer(4), parameter :: nstep = 1000
  integer(4) :: kstep
  real(8) :: h = 0.1d0
  real(8) :: xiso(2, nstep+1), xtemp(2), F, F0, dF(2)

  xiso(1,1) = 2.0d0
  xiso(2,1) = 0.0d0

  call fun_test(xiso(:,1), F, dF)
  F0 = F

  do kstep = 1, nstep
    xtemp = xiso(:, kstep)
    call level_set_step_2D(fun_test, h, xtemp)
    xiso(:, kstep+1) = xtemp
    call fun_test(xtemp, F, dF)
    print *, kstep, F - F0, xtemp
  end do

  contains

  subroutine fun_test(x, Fout, dFout)
    real(8), intent(in)  :: x(2)
    real(8), intent(out) :: Fout
    real(8), intent(out) :: dFout(2)

    Fout = 0.5d0*x(2)**2 - cos(x(1))
    dFout(1) = sin(x(1))
    dFout(2) = x(2)
  end subroutine

  subroutine fun_test2(x, Fout, dFout)
    real(8), intent(in)  :: x(2)
    real(8), intent(out) :: Fout
    real(8), intent(out) :: dFout(2)

    Fout = x(1)**2 + x(2)**2 - 1.0d0
    dFout(1) = 2.0d0*x(1)
    dFout(2) = 2.0d0*x(2)
  end subroutine

end program test_level_set