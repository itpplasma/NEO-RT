program test_endpoint_stationary_root
  use find_all_roots_mod, only : nroots,roots
  use pitch_boundary_mod, only : pitch_boundary_success,resolve_pitch_squared
  use potato_topology_mod, only : choose_two_sided_step, &
                                  root_has_two_sided_neighborhood,root_is_open_interval
  implicit none

  integer :: ierr,i,nopen
  logical :: resolved_step
  double precision :: pitch_squared_zero,pitch_squared_plus
  double precision :: root,h_fixed,h_adaptive,qminus,qzero,qplus,curvature

  call find_all_roots_bracketed(manufactured_stationary_roots,0.d0,1.d0,ierr)
  call require(ierr.eq.0,'manufactured stationary-root scan ierr')
  call require(nroots.eq.3,'manufactured stationary-root count')

  nopen=0
  do i=1,nroots
    if(root_is_open_interval(roots(i),0.d0,1.d0)) nopen=nopen+1
  enddo

  call require(.not.root_is_open_interval(roots(1),0.d0,1.d0), &
               'left stationary boundary classified as an open-interval extremum')
  call require(root_is_open_interval(roots(2),0.d0,1.d0), &
               'interior stationary root rejected')
  call require(.not.root_is_open_interval(roots(3),0.d0,1.d0), &
               'right stationary boundary classified as an open-interval extremum')
  call require(nopen.eq.1,'wrong number of open-interval stationary roots')

  call require(.not.root_has_two_sided_neighborhood(0.d0,0.d0,1.d0,1.d-6), &
               'left endpoint admitted by two-sided classifier')
  call require(root_has_two_sided_neighborhood(0.5d0,0.d0,1.d0,1.d-6), &
               'strictly interior root rejected by two-sided classifier')
  call require(.not.root_has_two_sided_neighborhood(1.d0,0.d0,1.d0,1.d-6), &
               'right endpoint admitted by two-sided classifier')
  call require(.not.root_has_two_sided_neighborhood(1.d-7,0.d0,1.d0,1.d-6), &
               'one-sided-neighborhood root admitted by classifier')

! A fixed probe can cross a reflecting boundary even when the stationary root
! itself is strictly interior.  The adaptive probe remains centered and keeps
! the quadratic curvature oracle well-defined on both sides.
  root=8.d-4
  h_fixed=1.d-3
  call require(root-h_fixed.le.0.d0,'fixed probe did not cross manufactured boundary')
  call choose_two_sided_step(root,0.d0,1.d0,h_fixed,0.5d0,h_adaptive,resolved_step)
  call require(resolved_step,'adaptive near-boundary step rejected')
  call require(root-h_adaptive.gt.0.d0 .and. root+h_adaptive.lt.1.d0, &
               'adaptive step is not strictly centered in the domain')
  qminus=(root-h_adaptive-root)**2
  qzero=0.d0
  qplus=(root+h_adaptive-root)**2
  curvature=(qplus-2.d0*qzero+qminus)/h_adaptive**2
  call require(abs(curvature-2.d0).lt.1.d-12, &
               'adaptive centered quadratic curvature is wrong')
  call require(curvature.gt.0.d0,'adaptive quadratic minimum was misclassified')

  h_adaptive=1.d0
  call choose_two_sided_step(1.d0,1.d0-1.d-15,2.d0,1.d-3,0.5d0,h_adaptive,resolved_step)
  call require(.not.resolved_step,'roundoff-scale interval was not failed closed')
  call require(h_adaptive.eq.0.d0,'failed-closed step was not cleared')

! At fixed local B/p^2=1, the normalized J_perp argument is the pitch ratio.
! The exact J_perp=0 limit is parallel motion (lambda^2=1); 0+ must approach
! it from inside the allowed pitch domain without becoming an endpoint root.
  call resolve_pitch_squared(0.d0,pitch_squared_zero,ierr)
  call require(ierr.eq.pitch_boundary_success,'Jperp=0 pitch-limit status')
  call require(pitch_squared_zero.eq.1.d0,'Jperp=0 pitch-limit value')
  call resolve_pitch_squared(1.d-12,pitch_squared_plus,ierr)
  call require(ierr.eq.pitch_boundary_success,'Jperp=0+ pitch-limit status')
  call require(pitch_squared_plus.lt.pitch_squared_zero .and. &
               pitch_squared_plus.gt.0.d0,'Jperp=0+ does not approach from inside')

  print *,'endpoint stationary-root classification: passed'

contains

  subroutine require(ok,msg)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: msg

    if(ok) return
    print *,trim(msg)
    error stop 1
  end subroutine require

  subroutine manufactured_stationary_roots(x,f,df)
    double precision, intent(in) :: x
    double precision, intent(out) :: f,df

    f=x*(x-0.5d0)*(x-1.d0)
    df=3.d0*x*x-3.d0*x+0.5d0
  end subroutine manufactured_stationary_roots

end program test_endpoint_stationary_root
