program test_find_all_roots_bracketed
  use find_all_roots_mod, only : customgrid,ncustom,niter,nroots,nsearch_min, &
                                 relerr_allroots,roots,xcustom
  implicit none

  integer :: ierr
  logical :: evaluated_outside_bounds
  double precision :: xunusable

  customgrid=.false.
  niter=100
  nsearch_min=1
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

  nsearch_min=10
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

  xcustom=[-2.d0,-1.d0]
  call check_upper_root('custom nodes below interval')
  xcustom=[2.d0,3.d0]
  call check_upper_root('custom nodes above interval')
  xcustom=[-2.d0,3.d0]
  call check_upper_root('no custom nodes inside interval')
  xcustom=[-2.d0,0.35d0,3.d0]
  call check_upper_root('one custom node inside interval')

  call check_bracket_end('uniform grid',.false.)
  call check_bracket_end('custom grid',.true.)

contains

  subroutine check_bracket_end(label,custom)
! find_valid_class_roots ends a bracket one ulp below a class node whose
! interpolation stencil is unusable. The function has no root below that node
! and jumps to the other sign on it, so any reported root is spurious.
    character(len=*), intent(in) :: label
    logical, intent(in) :: custom
    integer :: k
    double precision :: x1,x2

    nsearch_min=400
    customgrid=custom
    do k=1,200
      x1=-dble(k)/7.d0
      xunusable=dble(k)/13.d0
      x2=nearest(xunusable,-1.d0)
      if(custom) then
        xcustom=[x1-1.d0,x1,0.5d0*(x1+x2),xunusable,xunusable+1.d0]
        ncustom=size(xcustom)
      endif
      evaluated_outside_bounds=.false.
      call find_all_roots_bracketed(jump_at_unusable_node,x1,x2,ierr)
      call require(ierr.eq.0,label//' bracket end ierr')
      call require(.not.evaluated_outside_bounds, &
                   label//' evaluated on the node above the bracket end')
      call require(nroots.eq.0,label//' spurious root at the bracket end')
    enddo
  end subroutine check_bracket_end

  subroutine jump_at_unusable_node(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    if(x.ge.xunusable) then
      evaluated_outside_bounds=.true.
      f=1.d0
      df=0.d0
    else
      f=-2.d0-sin(x)
      df=-cos(x)
    endif
  end subroutine jump_at_unusable_node

  subroutine check_upper_root(label)
    character(len=*), intent(in) :: label
    integer :: entry

    ncustom=size(xcustom)
    nsearch_min=1
    do entry=1,2
      evaluated_outside_bounds=.false.
      if(entry.eq.1) then
        call find_all_roots(bounded_upper_root,0.d0,1.d0,ierr)
      else
        call find_all_roots_bracketed(bounded_upper_root,0.d0,1.d0,ierr)
      endif
      call require(ierr.eq.0,label//' ierr')
      call require(.not.evaluated_outside_bounds,label//' evaluated outside interval')
      call require(nroots.eq.1,label//' lost upper-interval root')
      call require(abs(roots(1)-0.95d0).lt.1.d-10,label//' root')
    enddo
  end subroutine check_upper_root

  subroutine bounded_upper_root(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    if(x.lt.0.d0 .or. x.gt.1.d0) evaluated_outside_bounds=.true.
    f=x-0.95d0
    df=1.d0
  end subroutine bounded_upper_root

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
