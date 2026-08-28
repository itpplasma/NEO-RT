program test_find_all_roots_bracketed
  use find_all_roots_mod, only : customgrid,ncustom,niter,nroots,nsearch_min, &
                                 relerr_allroots,roots,xcustom
  implicit none

  integer :: ierr
  logical :: evaluated_outside_bounds

  customgrid=.false.
  niter=100
  nsearch_min=10
  relerr_allroots=1.d-12

  call find_all_roots_bracketed(two_roots,0.d0,1.d0,ierr)
  call require(ierr.eq.0,'two_roots ierr')
  call require(nroots.eq.2,'two_roots nroots')
  call require(abs(roots(1)-0.25d0).lt.1.d-10,'two_roots first root')
  call require(abs(roots(2)-0.75d0).lt.1.d-10,'two_roots second root')

  call find_all_roots_bracketed(tangent_root,0.d0,1.d0,ierr)
  call require(ierr.eq.0,'tangent_root ierr')
  call require(nroots.eq.1,'tangent_root nroots')
  call require(abs(roots(1)-0.5d0).lt.1.d-10,'tangent_root root')

  allocate(xcustom(4))
  xcustom=[-2.d0,0.25d0,0.75d0,3.d0]
  ncustom=size(xcustom)
  customgrid=.true.
  evaluated_outside_bounds=.false.
  call find_all_roots(bounded_two_roots,0.d0,1.d0,ierr)
  call require(ierr.eq.0,'bounded legacy ierr')
  call require(.not.evaluated_outside_bounds, &
               'legacy root finder evaluated outside requested interval')
  call require(nroots.eq.2,'bounded legacy nroots')
  evaluated_outside_bounds=.false.
  call find_all_roots_bracketed(bounded_two_roots,0.d0,1.d0,ierr)
  call require(ierr.eq.0,'bounded bracketed ierr')
  call require(.not.evaluated_outside_bounds, &
               'bracketed root finder evaluated outside requested interval')
  call require(nroots.eq.2,'bounded bracketed nroots')

contains

  subroutine require(ok,msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg

    if(ok) return
    print *,trim(msg)
    error stop 1
  end subroutine require

  subroutine two_roots(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    f=(x-0.25d0)*(x-0.75d0)
    df=2.d0*x-1.d0
  end subroutine two_roots

  subroutine tangent_root(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    f=(x-0.5d0)**2
    df=2.d0*(x-0.5d0)
  end subroutine tangent_root

  subroutine bounded_two_roots(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    if(x.lt.0.d0 .or. x.gt.1.d0) evaluated_outside_bounds=.true.
    call two_roots(x,f,df)
  end subroutine bounded_two_roots

end program test_find_all_roots_bracketed
