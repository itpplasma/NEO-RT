!-----------------------------------------------------------------------
function sigmav_ion(Te)
  implicit none
  double precision, parameter :: Tm=2.d0 !eV
  double precision :: sigmav_ion, Te, x,y
  if(Te .lt. Tm) then
     sigmav_ion = 0.d0
     return
  endif
  x = log(Te)
  y=-3.271396786375d+01+x*(1.353655609057d+01+x*(-5.739328757388d+00            &
    +x*(1.563154982022d+00+x*(-2.877056004391d-01+x*(3.482559773737d-02         &
    +x*(-2.631976175590d-03+x*(1.119543953861d-04+x*(-2.039149852002d-06))))))))
  sigmav_ion = exp(y)
 return
end function sigmav_ion
!-----------------------------------------------------------------------
!function sigmav_ion(Te)
!  implicit none
!  double precision :: sigmav_ion, Te, x, Tm=2.d0,  T0=20.d0 ! eV
!  if(Te .lt. Tm) then
!     sigmav_ion = 0.d0
!     return
!  endif
!  x = log10(Te)
!  if(Te .lt. T0) then
!    sigmav_ion = 1.d1**(-3.054d0*x - 1.572d1*exp(-x) + 1.603d0*exp(-x**2))
! else
!    sigmav_ion = 1.d1**(-5.151d-1*x -2.563d0/x - 5.231d0)
! endif
! return
!end function sigmav_ion
!-----------------------------------------------------------------------
function sigmav_cx(T)
  implicit none
  double precision :: sigmav_cx, T
  double precision, save ::  Tmin, Tmax, hm1, dt1, dt2
  logical, save :: firstcall=.true.
  integer, save :: nmin, nmax, ih
  integer :: i
  double precision, save, dimension(:), allocatable :: svcx

  if(firstcall) then
     firstcall = .false.
     open(1,file='cx.dat')
     read(1,*)nmin, nmax, ih
     Tmin = dble(nmin)
     Tmax = dble(nmax)
     hm1  = 1.d0/dble(ih)
     allocate(svcx(nmin:nmax))
     do i=nmin,nmax,ih
        read(1,*) svcx(i)
     enddo
     close(1)
  endif
  if(T .lt. Tmin) T = Tmin
  if(T .gt. Tmax) T = Tmax
  dt1 = (T - Tmin)*hm1
  i = int(dt1) + 1
  dt2 = dt1  + dble(1 - i)
  sigmav_cx = svcx(i) + (svcx(i+1) - svcx(i))*dt2 
  return
end function sigmav_cx
!--------------------------------------------------------------
