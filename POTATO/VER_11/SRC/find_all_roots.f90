!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Modules:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module find_all_roots_mod
    logical :: customgrid=.false.
    integer :: nroots, nsearch_min=100, ncustom, niter=100
    double precision :: relerr_allroots=1.d-12
    double precision, dimension(:), allocatable :: xcustom,roots
  end module find_all_roots_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Routines:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_all_roots(fun,x1in,x2in,ierr)
!
! Finds all roots of function "fun" at the interval between x1 and x2.
! Subroutine fun(x,f,df) should provide the value of function f and
! its derivative df for a given input argument x.
! Output via module find_all_roots_mod:
! nroots          - number of roots
! roots(1:nroots) - array with roots
!
  use find_all_roots_mod, only : nsearch_min,niter,relerr_allroots,  &
                                 customgrid,ncustom,xcustom,         &
                                 nroots,roots
!
  implicit none
!
  integer :: i,iter,nsearch,ierr,ndummy,kx,kxc,k
  double precision :: x,x1in,x2in,f,df,hx,dx,errdist,xb,xe,dfb,dfe,xxtr,fxtr
  double precision :: x1,x2
  double precision, dimension(:), allocatable :: xarr,farr,dfarr,dummy1d
  external :: fun
!
  ierr=0
  x1=x1in
  x2=x2in
!
  errdist=relerr_allroots*abs(x2-x1)
!
! equidistant grid for primary search:
  nsearch=nsearch_min
  allocate(xarr(0:nsearch))
  hx=(x2-x1)/dble(nsearch)
!
  do i=0,nsearch
    xarr(i)=x1+hx*dble(i)
  enddo
!
! Merge the equidistant grid with externally provided custom (generally non-equidistand) grid:
!
  if(customgrid) then
    ndummy=nsearch+ncustom
    allocate(dummy1d(0:ndummy))
    if(xarr(0).lt.xcustom(1)) then
      dummy1d(0)=xarr(0)
      kx=1
      kxc=1
    else
      dummy1d(0)=xcustom(1)
      kx=0
      kxc=2
    endif
!
    k=0
!
    do while(kx.le.nsearch .and. kxc.le.ncustom)
      if(xarr(kx).lt.xcustom(kxc)) then
        k=k+1
        dummy1d(k)=xarr(kx)
        kx=kx+1
      else
        k=k+1
        dummy1d(k)=xcustom(kxc)
        kxc=kxc+1
      endif
    enddo
!
    do while(kx.le.nsearch)
      k=k+1
      dummy1d(k)=xarr(kx)
      kx=kx+1
    enddo
!
    do while(kxc.le.ncustom)
      k=k+1
      dummy1d(k)=xcustom(kxc)
      kxc=kxc+1
    enddo
!
    deallocate(xarr)
    nsearch=ndummy
    allocate(xarr(0:nsearch))
    xarr=dummy1d
!
! eliminate small intervals:
    hx=0.5d0*min(hx,minval(xcustom(2:ncustom)-xcustom(1:ncustom-1)))
    nsearch=0
!
    do k=1,ndummy
      if(xarr(k)-dummy1d(nsearch).gt.hx) then
        nsearch=nsearch+1
        dummy1d(nsearch)=xarr(k)
      endif
    enddo
!
    deallocate(xarr)
    allocate(xarr(0:nsearch))
    xarr=dummy1d(0:nsearch)
    deallocate(dummy1d)
  endif
!
  x1=minval(xarr)
  x2=maxval(xarr)
!
! End merge the equidistant grid with externally provided custom grid
!
  allocate(farr(0:nsearch),dfarr(0:nsearch))
!
  do i=0,nsearch
!
    call fun(xarr(i),farr(i),dfarr(i))
!
  enddo
!
  nroots=0
  if(allocated(roots)) deallocate(roots)
!
  do i=1,nsearch
    if(dfarr(i)*dfarr(i-1).gt.0.d0) then
! first derivative does not change sign in the interval, look for the root:
      if(farr(i)*farr(i-1).le.0.d0) then
        if(farr(i).eq.0.d0) cycle
        x=(farr(i)*xarr(i-1)-farr(i-1)*xarr(i))/(farr(i)-farr(i-1))
!
        call addroot
!
      endif
    else
! first derivative changes sign in the interval, find an extremum first:
      xb=xarr(i-1)
      xe=xarr(i)
      dfb=dfarr(i-1)
      dfe=dfarr(i)
!
      do iter=1,niter
        x=0.5d0*(xb+xe)
!
        call fun(x,f,df)
!
        if(dfb*df.gt.0.d0) then
          xb=x
        else
          xe=x
        endif
        if(abs(xe-xb).lt.errdist) exit
      enddo
!
! split interval in two by the extremum point and look for roots:
      xxtr=x
      fxtr=f
      if(farr(i-1)*fxtr.lt.0.d0) then
        x=(fxtr*xarr(i-1)-farr(i-1)*xxtr)/(fxtr-farr(i-1))
!
        call addroot
!
      endif
      if(fxtr*farr(i).le.0.d0) then
        if(farr(i).eq.0.d0) cycle
        x=(farr(i)*xxtr-fxtr*xarr(i))/(farr(i)-fxtr)
!
        call addroot
!
      endif
    endif
  enddo
!
!------------  
!
  contains
!
!------------  
!
  subroutine addroot
!
! Root adjustment by Newton method
!
  double precision :: xx
!
! extend storage arragy for a new root:
  nroots=nroots+1
  if(allocated(roots)) then 
    allocate(dummy1d(nroots-1))
    dummy1d=roots
    deallocate(roots)
    allocate(roots(nroots))
    roots(1:nroots-1)=dummy1d
    deallocate(dummy1d)
  else
    allocate(roots(nroots))
  endif
!
! Newton iterations:
  do iter=1,niter
!
    call fun(x,f,df)
!
    dx=f/df
    xx=x-dx
! keep intermediate solution within search bounds x1<x<x2:
    if(xx.le.x1) then
      x=0.5d0*(x+x1)
    elseif(xx.ge.x2) then
      x=0.5d0*(x+x2)
    else
      x=xx
    endif
    if(abs(dx).lt.errdist) exit
  enddo
!
  if(iter.gt.niter) then
    print *,'No convergence of Newton in find_all_roots: error = ',sngl(abs(dx)), &
            ' tolerance = ',sngl(errdist)
    ierr=1
  endif
!
  roots(nroots)=x
!
  end subroutine addroot
!
!------------  
!
  end subroutine find_all_roots
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
