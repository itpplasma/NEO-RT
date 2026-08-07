program test_find_all_roots_bracketed
  use find_all_roots_mod, only : customgrid,ncustom,niter,nroots,nsearch_min, &
                                 relerr_allroots,roots,xcustom,              &
                                 root_eval_valid,root_eval_error,            &
                                 root_invalid_domain,root_left_endpoint_contracted, &
                                 root_right_endpoint_contracted,root_search_left,   &
                                 root_search_right,root_left_invalid_bracket,       &
                                 root_right_invalid_bracket
  implicit none

  integer :: ierr
  logical :: evaluated_outside_bounds

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

  allocate(xcustom(4))
  xcustom=[-2.d0,0.2d0,0.8d0,3.d0]
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

  customgrid=.false.
  nsearch_min=100
  evaluated_outside_bounds=.false.
  call find_all_roots(decimal_interval_roots,0.1d0,0.9d0,ierr)
  call require(ierr.eq.0,'decimal legacy ierr')
  call require(.not.evaluated_outside_bounds, &
               'legacy equidistant grid exceeded decimal interval')
  evaluated_outside_bounds=.false.
  call find_all_roots_bracketed(decimal_interval_roots,0.1d0,0.9d0,ierr)
  call require(ierr.eq.0,'decimal bracketed ierr')
  call require(.not.evaluated_outside_bounds, &
               'bracketed equidistant grid exceeded decimal interval')

  customgrid=.false.
  nsearch_min=1
  call find_all_roots_bracketed(endpoint_roots,0.d0,1.d0,ierr)
  call require(ierr.eq.0,'endpoint roots ierr')
  call require(nroots.eq.2,'endpoint roots nroots')
  call require(abs(roots(1)).lt.1.d-14,'endpoint left root')
  call require(abs(roots(2)-0.75d0).lt.1.d-10,'endpoint interior root')

  call find_all_roots(flat_derivative_root,0.d0,1.d0,ierr)
  call require(ierr.eq.0,'flat derivative ierr')
  call require(nroots.eq.1,'flat derivative nroots')
  call require(abs(roots(1)-0.5d0).lt.1.d-10,'flat derivative root')

  nsearch_min=8
  call find_all_roots_bracketed(left_endpoint_invalid,0.d0,1.d0,ierr)
  call require(ierr.eq.0,'left endpoint contraction ierr')
  call require(nroots.eq.2,'left endpoint contraction nroots')
  call require(root_left_endpoint_contracted,'left endpoint was not recorded')
  call require(.not.root_right_endpoint_contracted,'unexpected right contraction')
  call require(root_left_invalid_bracket.lt.root_search_left, &
               'left contraction was not bracketed')
  call require(root_search_left.gt.0.d0 .and. root_search_left.lt.0.2d0, &
               'left contraction escaped endpoint domain')
  call require(abs(roots(1)-0.3d0).lt.1.d-10,'left endpoint first root')
  call require(abs(roots(2)-0.7d0).lt.1.d-10,'left endpoint second root')

  call find_all_roots_bracketed(right_endpoint_invalid,0.d0,1.d0,ierr)
  call require(ierr.eq.0,'right endpoint contraction ierr')
  call require(nroots.eq.2,'right endpoint contraction nroots')
  call require(.not.root_left_endpoint_contracted,'unexpected left contraction')
  call require(root_right_endpoint_contracted,'right endpoint was not recorded')
  call require(root_search_right.lt.root_right_invalid_bracket, &
               'right contraction was not bracketed')
  call require(root_search_right.gt.0.8d0 .and. root_search_right.lt.1.d0, &
               'right contraction escaped endpoint domain')
  call require(abs(roots(1)-0.3d0).lt.1.d-10,'right endpoint first root')
  call require(abs(roots(2)-0.7d0).lt.1.d-10,'right endpoint second root')

  call find_all_roots_bracketed(both_endpoints_invalid,0.d0,1.d0,ierr)
  call require(ierr.eq.0,'both endpoint contractions ierr')
  call require(nroots.eq.2,'both endpoint contractions nroots')
  call require(root_left_endpoint_contracted,'both endpoints left record')
  call require(root_right_endpoint_contracted,'both endpoints right record')
  call require(abs(roots(1)-0.3d0).lt.1.d-10,'both endpoints first root')
  call require(abs(roots(2)-0.7d0).lt.1.d-10,'both endpoints second root')

  call find_all_roots_bracketed(interior_invalid_hole,0.d0,1.d0,ierr)
  call require(ierr.eq.root_invalid_domain,'interior hole status')
  call require(nroots.eq.0,'interior hole fabricated root')
  call require(.not.root_left_endpoint_contracted,'interior hole left trim')
  call require(.not.root_right_endpoint_contracted,'interior hole right trim')

  call find_all_roots_bracketed(no_valid_domain,0.d0,1.d0,ierr)
  call require(ierr.eq.root_invalid_domain,'empty valid domain status')
  call require(nroots.eq.0,'empty valid domain fabricated root')

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

  subroutine decimal_interval_roots(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    if(x.lt.0.1d0 .or. x.gt.0.9d0) evaluated_outside_bounds=.true.
    f=(x-0.3d0)*(x-0.7d0)
    df=2.d0*x-1.d0
  end subroutine decimal_interval_roots

  subroutine endpoint_roots(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    f=x*(x-0.75d0)
    df=2.d0*x-0.75d0
  end subroutine endpoint_roots

  subroutine flat_derivative_root(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    f=(x-0.5d0)**3
    df=3.d0*(x-0.5d0)**2
  end subroutine flat_derivative_root

  subroutine left_endpoint_invalid(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    if(x.lt.0.08d0) then
      root_eval_valid=.false.
      root_eval_error=root_invalid_domain
      f=0.d0
      df=0.d0
      return
    endif
    f=(x-0.3d0)*(x-0.7d0)
    df=2.d0*x-1.d0
  end subroutine left_endpoint_invalid

  subroutine right_endpoint_invalid(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    if(x.gt.0.92d0) then
      root_eval_valid=.false.
      root_eval_error=root_invalid_domain
      f=0.d0
      df=0.d0
      return
    endif
    f=(x-0.3d0)*(x-0.7d0)
    df=2.d0*x-1.d0
  end subroutine right_endpoint_invalid

  subroutine both_endpoints_invalid(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    if(x.lt.0.08d0 .or. x.gt.0.92d0) then
      root_eval_valid=.false.
      root_eval_error=root_invalid_domain
      f=0.d0
      df=0.d0
      return
    endif
    f=(x-0.3d0)*(x-0.7d0)
    df=2.d0*x-1.d0
  end subroutine both_endpoints_invalid

  subroutine interior_invalid_hole(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    if(x.gt.0.44d0 .and. x.lt.0.56d0) then
      root_eval_valid=.false.
      root_eval_error=root_invalid_domain
      f=0.d0
      df=0.d0
      return
    endif
    f=(x-0.25d0)*(x-0.75d0)
    df=2.d0*x-1.d0
  end subroutine interior_invalid_hole

  subroutine no_valid_domain(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    root_eval_valid=.false.
    root_eval_error=root_invalid_domain
    f=0.d0
    df=0.d0
  end subroutine no_valid_domain

end program test_find_all_roots_bracketed
