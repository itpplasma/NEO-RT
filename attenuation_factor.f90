!
  module thetadata_mod
    logical :: prop=.true.
    integer :: npoiarg,iunitarg=741
    double precision :: argmin,argstepinv
    double precision, dimension(:), allocatable :: argvals,thetavals
  end module thetadata_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine attenuation_factor(D,Theta)
!
  use thetadata_mod
  USE polylag_3,  ONLY : mp,indef, plag1d
!
  implicit none
!
  integer :: i
  double precision :: D,Theta,arglog,fun,der
  integer,          dimension(mp) :: indu
  double precision, dimension(mp) :: xp,fp
!
  if(prop) then
    prop=.false.
    open(iunitarg,file='thetafun_inp.dat')
    npoiarg=0
    do
      read(iunitarg,*,end=1) arglog
      npoiarg=npoiarg+1
    enddo
1   close(iunitarg)
    allocate(argvals(npoiarg),thetavals(npoiarg))
    open(iunitarg,file='thetafun_inp.dat')
    do i=1,npoiarg
      read(iunitarg,*) arglog,Theta
      argvals(i)=log(arglog)
      thetavals(i)=log(Theta)
    enddo
    close(iunitarg)
    argmin=argvals(1)
    argstepinv=1.d0/(argvals(2)-argvals(1))
  endif
!
  arglog=log(D)
!
  call indef(arglog,argmin,argstepinv,npoiarg,indu)
!
  xp=argvals(indu)
  fp=thetavals(indu)
!
  call plag1d(arglog,fp,argstepinv,xp,fun,der)
!
  Theta=exp(fun)
!
  end subroutine attenuation_factor
