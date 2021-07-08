!
  subroutine integrator_stw(N,x0,xn,func,antideriv,ierr)
!
! Computes integral of function "func" given at the equidistant grid of N intervals between x0 and xn
! as function of upper limit (lower limit is x0). Uses cubic polynomial interpolation of function.
! Input: 
!        N              - number of intervals
!        x0             - lower limit
!        xn             - upper limt
!        func(0:N)      - values of function at grid nodes x0,x1,..,xn
! Output:
!        antideriv(0:N) - values of integral (antiderivative) at grid nodes
!        ierr           - error code (0 - normal work, 1 - number of intervals is less than 3
!
  implicit none
!
  double precision, dimension(0:3) :: wl=(/3.d0/8.d0, 19.d0/24.d0, -5.d0/24.d0, 1.d0/24.d0/),    &
                                      wm=(/-1.d0/24.d0, 13.d0/24.d0, 13.d0/24.d0, -1.d0/24.d0/), &
                                      wr=(/1.d0/24.d0, -5.d0/24.d0, 19.d0/24.d0, 3.d0/8.d0/)
!
  double precision :: x0,xn,integral,h
  integer :: N,i,ierr
  double precision, dimension(0:N), intent(in)  :: func
  double precision, dimension(0:N), intent(out) :: antideriv
!
  if(N.lt.3) then
    print *,'integrator_stw: N < 3'
    ierr=1
  else
    ierr=0
  endif
!
  h=(xn-x0)/N
!
  antideriv(0)=0.d0
  antideriv(1)=sum(wl*func(0:3))
!  
  do i=2,N-1
    antideriv(i)=antideriv(i-1)+sum(wm*func(i-2:i+1))
  enddo
!
  antideriv(N)=antideriv(N-1)+sum(wr*func(N-3:N))
!
  antideriv=antideriv*h
!
  end subroutine integrator_stw
!
