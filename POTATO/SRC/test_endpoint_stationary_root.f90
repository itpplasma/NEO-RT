program test_endpoint_stationary_root
  use find_all_roots_mod, only : nroots,roots
  use pitch_boundary_mod, only : pitch_boundary_success,resolve_pitch_squared
  use sample_matrix_out_mod, only : nlagr,n1,n2,npoi,itermax,xbeg,xend,eps, &
                                    icount,topology_context_h, &
                                    sample_matrix_topology_transition, &
                                    topology_stencil_is_compatible,topology_arr, &
                                    xarr,amat_arr
  use potato_topology_mod, only : choose_two_sided_step, &
                                  root_has_two_sided_neighborhood,root_is_open_interval
  implicit none

  integer :: ierr,i,nopen
  integer :: signatures(4)
  logical :: resolved_step
  double precision :: pitch_squared_zero,pitch_squared_plus
  double precision :: root,h_fixed,h_adaptive,qminus,qzero,qplus,curvature
  external :: sample_matrix_out,sample_matrix_out_partitioned

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

! A sampler may converge numerically on each smooth branch, but the two
! branches are not one interpolation problem.  This is the independent
! convergence oracle for the outer J_perp topology gate.
  signatures=(/7,7,9,9/)
  call require(topology_stencil_is_compatible(signatures,1,2), &
               'same-topology interpolation stencil was rejected')
  call require(.not.topology_stencil_is_compatible(signatures,2,3), &
               'cross-topology interpolation stencil was admitted')
  call require(topology_stencil_is_compatible(signatures,3,4), &
               'second same-topology interpolation stencil was rejected')
  call test_topology_gate
  call test_partitioned_piecewise
  call test_endpoint_contraction

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

  subroutine test_topology_gate
    integer :: local_ierr

    nlagr=3
    n1=1
    n2=1
    npoi=6
    itermax=4
    xbeg=0.d0
    xend=1.d0
    eps=1.d-6
    icount=0
    topology_context_h=42.d0
    call sample_matrix_out(manufactured_topology_matrix,local_ierr)
    call require(local_ierr.eq.sample_matrix_topology_transition, &
                 'cross-topology sampler did not fail closed')
  end subroutine test_topology_gate

  subroutine manufactured_topology_matrix
    use sample_matrix_out_mod, only : x,amat,topology_signature,topology_error

    topology_error=0
    if(x.lt.0.5d0) then
      topology_signature=7
      amat(1,1)=x*x
    else
      topology_signature=9
      amat(1,1)=x*x+1.d0
    endif
  end subroutine manufactured_topology_matrix

  subroutine test_partitioned_piecewise
    integer :: local_ierr,j
    double precision :: integral,exact

    nlagr=3
    n1=1
    n2=1
    npoi=9
    itermax=4
    xbeg=0.d0
    xend=1.d0
    eps=1.d-10
    icount=0
    topology_context_h=43.d0
    call sample_matrix_out_partitioned(manufactured_piecewise_matrix,local_ierr)
    call require(local_ierr.eq.0,'piecewise sampler did not converge')

    integral=0.d0
    do j=2,npoi
      if(topology_arr(j).eq.topology_arr(j-1)) then
        integral=integral+0.5d0*(amat_arr(1,1,j)+amat_arr(1,1,j-1)) &
            *(xarr(j)-xarr(j-1))
      endif
    enddo
    exact=1.375d0
    call require(abs(integral-exact).lt.1.d-5, &
                 'piecewise trapezoid missed analytic branch integral')
  end subroutine test_partitioned_piecewise

  subroutine manufactured_piecewise_matrix
    use sample_matrix_out_mod, only : x,amat,topology_signature,topology_error

    topology_error=0
    if(x.le.0.5d0) then
      topology_signature=7
      amat(1,1)=x
    else
      topology_signature=9
      amat(1,1)=2.d0*x+1.d0
    endif
  end subroutine manufactured_piecewise_matrix

  subroutine test_endpoint_contraction
    integer :: local_ierr

    nlagr=3
    n1=1
    n2=1
    npoi=9
    itermax=4
    xbeg=0.d0
    xend=1.d0
    eps=1.d-10
    icount=0
    topology_context_h=44.d0
    call sample_matrix_out_partitioned(manufactured_invalid_endpoint,local_ierr)
    call require(local_ierr.eq.0,'invalid global endpoint was not contracted')
    call require(abs(xarr(1)-0.1d0).lt.1.d-6, &
                 'endpoint contraction was not refined to the valid boundary')
  end subroutine test_endpoint_contraction

  subroutine manufactured_invalid_endpoint
    use sample_matrix_out_mod, only : x,amat,topology_signature,topology_error

    if(x.lt.0.1d0) then
      topology_error=2
      topology_signature=0
      amat(1,1)=0.d0
    else
      topology_error=0
      topology_signature=11
      amat(1,1)=x
    endif
  end subroutine manufactured_invalid_endpoint

end program test_endpoint_stationary_root
