!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
! Modules:
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module find_all_roots_mod
    integer, parameter :: root_success=0
    integer, parameter :: root_nonconverged=1
    integer, parameter :: root_invalid_domain=2
    integer, parameter :: root_invalid_interval=3
    integer, parameter :: root_no_intersection=4
    logical :: customgrid=.false.
    logical :: root_eval_valid=.true.
    integer :: root_eval_error=0
    integer :: nroots, nsearch_min=100, ncustom, niter=100
    double precision :: relerr_allroots=1.d-12
    double precision, dimension(:), allocatable :: xcustom,roots
    !$omp threadprivate(customgrid,root_eval_valid,root_eval_error,ncustom,niter, &
    !$omp&                relerr_allroots,xcustom,nroots,roots)
  end module find_all_roots_mod
!
! The two public entry points intentionally share one bounded implementation.
! A Newton extrapolation is not a root certificate when the derivative is flat
! or when the physical function is undefined on part of the requested interval.
!
  subroutine find_all_roots(fun,x1in,x2in,ierr)
  external :: fun
  integer, intent(out) :: ierr
  double precision, intent(in) :: x1in,x2in

  call find_all_roots_bracketed(fun,x1in,x2in,ierr)

  end subroutine find_all_roots
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine find_all_roots_bracketed(fun,x1in,x2in,ierr)

  use find_all_roots_mod, only : nsearch_min,niter,relerr_allroots,       &
                                 customgrid,ncustom,xcustom,              &
                                 nroots,roots,root_eval_valid,            &
                                 root_eval_error,root_success,            &
                                 root_nonconverged,root_invalid_domain,   &
                                 root_invalid_interval

  implicit none

  integer, intent(out) :: ierr
  double precision, intent(in) :: x1in,x2in
  integer :: i,iter,nsearch,ndummy,kx,kxc,k,kxc_first,kxc_last
  integer :: ncustom_inside
  double precision :: x,hx,errdist,xb,xe,dfb,xext,fext
  double precision, dimension(:), allocatable :: xarr,farr,dfarr,dummy1d
  logical, dimension(:), allocatable :: validarr
  logical :: valid,ok
  external :: fun

  ierr=root_success
  nroots=0
  if(allocated(roots)) deallocate(roots)
  if(x2in.le.x1in) then
    ierr=root_invalid_interval
    return
  endif

  errdist=relerr_allroots*abs(x2in-x1in)
  if(errdist.le.0.d0) errdist=epsilon(1.d0)*abs(x2in-x1in)

  nsearch=nsearch_min
  if(nsearch.lt.1) nsearch=1
  allocate(xarr(0:nsearch))
  hx=(x2in-x1in)/dble(nsearch)
  do i=0,nsearch
    xarr(i)=x1in+hx*dble(i)
  enddo
  xarr(0)=x1in
  xarr(nsearch)=x2in

! Merge the equidistant grid with the caller's sorted custom grid.  Values
! outside the requested interval are ignored, never used as extrapolation.
  if(customgrid) then
    kxc_first=1
    do while(kxc_first.le.ncustom)
      if(xcustom(kxc_first).ge.x1in) exit
      kxc_first=kxc_first+1
    enddo
    kxc_last=ncustom
    do while(kxc_last.ge.kxc_first)
      if(xcustom(kxc_last).le.x2in) exit
      kxc_last=kxc_last-1
    enddo
    ncustom_inside=max(0,kxc_last-kxc_first+1)
    ndummy=nsearch+ncustom_inside
    allocate(dummy1d(0:ndummy))
    if(ncustom_inside.eq.0) then
      dummy1d(0)=xarr(0)
      kx=1
      kxc=kxc_first
    elseif(xarr(0).lt.xcustom(kxc_first)) then
      dummy1d(0)=xarr(0)
      kx=1
      kxc=kxc_first
    else
      dummy1d(0)=xcustom(kxc_first)
      kx=0
      kxc=kxc_first+1
    endif
    k=0
    do while(kx.le.nsearch .and. kxc.le.kxc_last)
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
    do while(kx.le.nsearch)
      k=k+1
      dummy1d(k)=xarr(kx)
      kx=kx+1
    enddo
    do while(kxc.le.kxc_last)
      k=k+1
      dummy1d(k)=xcustom(kxc)
      kxc=kxc+1
    enddo
    deallocate(xarr)
    nsearch=ndummy
    allocate(xarr(0:nsearch))
    xarr=dummy1d
    if(ncustom_inside.gt.1) then
      hx=0.5d0*min(hx,minval(xcustom(kxc_first+1:kxc_last) &
                             -xcustom(kxc_first:kxc_last-1)))
    endif
    nsearch=0
    do k=1,ndummy
      if(xarr(k)-dummy1d(nsearch).gt.hx) then
        nsearch=nsearch+1
        dummy1d(nsearch)=xarr(k)
      endif
    enddo
    deallocate(xarr)
    allocate(xarr(0:nsearch))
    xarr=dummy1d(0:nsearch)
    deallocate(dummy1d)
  endif

  allocate(farr(0:nsearch),dfarr(0:nsearch),validarr(0:nsearch))
  do i=0,nsearch
    call evaluate_at(xarr(i),farr(i),dfarr(i),valid)
    validarr(i)=valid
    if(.not.valid) then
      call fail_search(root_invalid_domain)
      return
    endif
  enddo

  do i=1,nsearch
    if(farr(i-1).eq.0.d0) then
      x=xarr(i-1)
      call addroot_bracketed
      if(ierr.ne.root_success) then
        call fail_search(ierr)
        return
      endif
    endif

    if(farr(i-1)*farr(i).lt.0.d0) then
      call refine_bracket(xarr(i-1),xarr(i),farr(i-1),farr(i))
      if(ierr.ne.root_success) then
        call fail_search(ierr)
        return
      endif
    endif

! A derivative-zero interval is not evidence of an extremum.  A sign change
! (including an endpoint zero) is, and is refined by derivative bisection.
    if(dfarr(i).eq.0.d0 .and. dfarr(i-1).eq.0.d0) cycle
    if(dfarr(i)*dfarr(i-1).le.0.d0) then
      xb=xarr(i-1)
      xe=xarr(i)
      dfb=dfarr(i-1)
      xext=xb
      fext=farr(i-1)
      call refine_derivative(xb,xe,dfb,dfarr(i),xext,fext,ok)
      if(.not.ok) then
        call fail_search(ierr)
        return
      endif
      call refine_bracket(xarr(i-1),xext,farr(i-1),fext)
      if(ierr.ne.root_success) then
        call fail_search(ierr)
        return
      endif
      call refine_bracket(xext,xarr(i),fext,farr(i))
      if(ierr.ne.root_success) then
        call fail_search(ierr)
        return
      endif
    endif
  enddo

  if(farr(nsearch).eq.0.d0) then
    x=xarr(nsearch)
    call addroot_bracketed
  endif
  if(ierr.ne.root_success) then
    call fail_search(ierr)
    return
  endif

  deallocate(xarr,farr,dfarr,validarr)

  contains

  subroutine evaluate_at(xin,fout,dfout,is_valid)
  double precision, intent(in) :: xin
  double precision, intent(out) :: fout,dfout
  logical, intent(out) :: is_valid

  root_eval_valid=.true.
  root_eval_error=root_success
  call fun(xin,fout,dfout)
  is_valid=root_eval_valid
  if(.not.is_valid) then
    fout=0.d0
    dfout=0.d0
    if(root_eval_error.eq.root_success) root_eval_error=root_invalid_domain
    ierr=root_eval_error
  endif
  end subroutine evaluate_at

  subroutine addroot_bracketed
  double precision :: rootdist

  if(x.lt.x1in .or. x.gt.x2in) then
    ierr=root_invalid_interval
    return
  endif
  rootdist=max(errdist,epsilon(1.d0)*max(1.d0,abs(x)))
  if(nroots.gt.0) then
    if(minval(abs(roots(1:nroots)-x)).le.rootdist) return
  endif
  nroots=nroots+1
  if(allocated(dummy1d)) deallocate(dummy1d)
  if(nroots.gt.1) then
    allocate(dummy1d(nroots-1))
    dummy1d=roots
    deallocate(roots)
    allocate(roots(nroots))
    roots(1:nroots-1)=dummy1d
    deallocate(dummy1d)
  else
    allocate(roots(nroots))
  endif
  roots(nroots)=x
  end subroutine addroot_bracketed

  subroutine refine_derivative(xlo,xhi,dflo,dfhi,xout,fout,ok)
  double precision, intent(in) :: xlo,xhi,dflo,dfhi
  double precision, intent(out) :: xout,fout
  logical, intent(out) :: ok
  double precision :: a,b,da,db,mid,fm,dm

  ok=.false.
  a=xlo
  b=xhi
  da=dflo
  db=dfhi
  if(da.eq.0.d0) then
    xout=a
    call evaluate_at(xout,fout,dm,ok)
    return
  endif
  if(db.eq.0.d0) then
    xout=b
    call evaluate_at(xout,fout,dm,ok)
    return
  endif
  if(da*db.gt.0.d0) then
    ierr=root_invalid_interval
    return
  endif
  do iter=1,niter
    mid=0.5d0*(a+b)
    call evaluate_at(mid,fm,dm,ok)
    if(.not.ok) return
    if(dm.eq.0.d0 .or. abs(b-a).le.errdist) then
      xout=mid
      fout=fm
      return
    endif
    if(da*dm.le.0.d0) then
      b=mid
      db=dm
    else
      a=mid
      da=dm
    endif
  enddo
  ierr=root_nonconverged
  end subroutine refine_derivative

  subroutine refine_bracket(xlo,xhi,flo,fhi)
  double precision, intent(in) :: xlo,xhi,flo,fhi
  double precision :: a,b,fa,fb,mid,fm,dfm
  logical :: valid_mid

  if(flo.eq.0.d0) then
    x=xlo
    call addroot_bracketed
    return
  endif
  if(fhi.eq.0.d0) then
    x=xhi
    call addroot_bracketed
    return
  endif
  if(flo*fhi.gt.0.d0) return

  a=xlo
  b=xhi
  fa=flo
  fb=fhi
  do iter=1,niter
    mid=0.5d0*(a+b)
    call evaluate_at(mid,fm,dfm,valid_mid)
    if(.not.valid_mid) return
    if(fm.eq.0.d0 .or. abs(b-a).le.errdist) then
      x=mid
      call addroot_bracketed
      return
    endif
    if((fa.lt.0.d0 .and. fm.gt.0.d0) .or. &
       (fa.gt.0.d0 .and. fm.lt.0.d0)) then
      b=mid
      fb=fm
    else
      a=mid
      fa=fm
    endif
  enddo
  ierr=root_nonconverged
  end subroutine refine_bracket

  subroutine fail_search(status)
  integer, intent(in) :: status

  ierr=status
  nroots=0
  if(allocated(roots)) deallocate(roots)
  if(allocated(xarr)) deallocate(xarr)
  if(allocated(farr)) deallocate(farr)
  if(allocated(dfarr)) deallocate(dfarr)
  if(allocated(validarr)) deallocate(validarr)
  if(allocated(dummy1d)) deallocate(dummy1d)
  end subroutine fail_search

  end subroutine find_all_roots_bracketed
