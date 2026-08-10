program test_endpoint_stationary_root
  use find_all_roots_mod, only : nroots,roots
  use pitch_boundary_mod, only : pitch_boundary_success,resolve_pitch_squared
  use sample_matrix_out_mod, only : nlagr,n1,n2,npoi,itermax,xbeg,xend,eps, &
                                    icount,topology_context_h, &
                                    sample_matrix_topology_transition, &
                                    sample_matrix_contribution_unresolved, &
                                    topology_stencil_is_compatible,topology_arr, &
                                    xarr,amat_arr,topology_gap_measure, &
                                    topology_gap_bound,topology_partition_tol, &
                                    topology_transition_count, &
                                    topology_contribution_error_certified
  use sample_matrix_mod, only : matrix_eval_starter_failure, &
                                matrix_eval_orbit_failure, &
                                matrix_eval_wall_loss, matrix_eval_nonfinite, &
                                matrix_eval_cut_domain, &
                                matrix_failure_is_open_boundary
  use potato_topology_mod, only : choose_two_sided_step, &
                                  root_has_two_sided_neighborhood,root_is_open_interval, &
                                  solve_fixedpoint_quadratic
  use resonance_status_mod, only : resonance_status_success, &
                                   resonance_status_root_failure, &
                                   resonance_status_tangent_root, &
                                   resonance_status_nonfinite_weight, &
                                   resonance_status_wall_not_zero, &
                                   classify_resonance_root, &
                                   resonance_status_from_root_error
  use potato_symbolic_kernel_mod, only : potato_jperp_kernel, potato_hm_kernel, &
                                         potato_root_jacobian_kernel, &
                                         potato_frequency_reduction_kernel, &
                                         potato_resonance_torque_kernel, &
                                         potato_gap_contribution_kernel, &
                                         potato_gap_square_map, &
                                         potato_gap_sqrt_integral, &
                                         potato_gap_sqrt_coefficient, &
                                         potato_class_linear_extrapolation, &
                                         potato_class_quadratic_extrapolation, &
                                         potato_limiting_kernel, &
                                         potato_limiting_kernel_checked, &
                                         potato_limiting_invalid_reference, &
                                         potato_limiting_invalid_distance
  use topology_gap_quadrature_mod, only : estimate_topology_gap_from_samples
  implicit none

  integer :: ierr,i,nopen
  integer :: signatures(4)
  logical :: resolved_step
  double precision :: pitch_squared_zero,pitch_squared_plus
  double precision :: root,h_fixed,h_adaptive,qminus,qzero,qplus,curvature
  double precision :: quadratic_roots(2)
  integer :: nquadratic
  logical :: quadratic_ok
  external :: sample_matrix_out,sample_matrix_out_partitioned, &
              sample_matrix_out_partitioned_certified

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
  call test_uncredentialed_narrow_island
  call test_partitioned_piecewise
  call test_certified_narrow_island
  call test_endpoint_contraction
  call test_endpoint_failure_classification
  call test_fixedpoint_quadratic
  call test_resonance_status
  call test_conjugate_mode_representation
  call test_gap_square_map
  call test_gap_sample_bound
  call test_symbolic_contract

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

! The same adaptive probe must work in a chart whose R coordinate is negative.
  root=-3.3819244005704327d2
  call choose_two_sided_step(root,-3.3819244984583725d2,-3.3799667783783626d2, &
                             1.d-6*max(1.d0,abs(root)),0.5d0,h_adaptive,resolved_step)
  call require(resolved_step,'negative-R adaptive step rejected')
  call require(root-h_adaptive.gt.-3.3819244984583725d2 .and. &
               root+h_adaptive.lt.-3.3799667783783626d2, &
               'negative-R adaptive step is not strictly centered')

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
    double precision :: integral,exact,integral_coarse,gap_coarse

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
    topology_partition_tol=1.d-2
    call sample_matrix_out_partitioned_certified( &
        manufactured_piecewise_matrix,manufactured_piecewise_boundaries,local_ierr, &
        manufactured_integrand_envelope)
    call require(local_ierr.eq.0,'piecewise sampler did not converge')
    call require(topology_contribution_error_certified, &
                 'piecewise contribution bound was not certified')

    integral=0.d0
    do j=2,npoi
      if(topology_arr(j).eq.topology_arr(j-1)) then
        integral=integral+0.5d0*(amat_arr(1,1,j)+amat_arr(1,1,j-1)) &
            *(xarr(j)-xarr(j-1))
      endif
    enddo
    exact=1.375d0
    integral_coarse=integral
    gap_coarse=topology_gap_measure
    call require(topology_gap_measure.le.topology_gap_bound+1.d-12, &
                 'coarse piecewise gap exceeded its certificate bound')

    topology_partition_tol=1.d-6
    call sample_matrix_out_partitioned_certified( &
        manufactured_piecewise_matrix,manufactured_piecewise_boundaries,local_ierr, &
        manufactured_integrand_envelope)
    call require(local_ierr.eq.0,'fine piecewise sampler did not converge')
    integral=0.d0
    do j=2,npoi
      if(topology_arr(j).eq.topology_arr(j-1)) then
        integral=integral+0.5d0*(amat_arr(1,1,j)+amat_arr(1,1,j-1)) &
            *(xarr(j)-xarr(j-1))
      endif
    enddo
    call require(topology_gap_measure.lt.gap_coarse, &
                 'piecewise omitted gap did not converge under refinement')
    call require(abs(integral-exact).lt.abs(integral_coarse-exact), &
                 'piecewise integral did not converge with bracket refinement')
    call require(abs(integral-exact).lt.1.d-5, &
                 'fine piecewise trapezoid missed analytic branch integral')
    topology_partition_tol=1.d-8
  end subroutine test_partitioned_piecewise

  subroutine test_uncredentialed_narrow_island
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
    topology_context_h=42.5d0
    call test_old_midpoint_misses_island
    call sample_matrix_out_partitioned(manufactured_narrow_matrix,local_ierr)
    call require(local_ierr.eq.4, &
                 'uncertified narrow A-B-A topology was accepted')
  end subroutine test_uncredentialed_narrow_island

  subroutine test_old_midpoint_misses_island
! Reproduce the 70c7f65 sampling geometry independently of the production
! sampler: 31 equidistant seeds plus one midpoint for an equal-signature
! interval.  The island is deliberately between the .1 and .1333 seeds and
! away from their .1167 midpoint, so that accepting this scan would silently
! lose a finite-measure component.
    integer, parameter :: nseed=31
    integer :: local_i
    double precision :: left,right,midpoint
    logical :: seen

    seen=.false.
    do local_i=0,nseed-1
      left=dble(local_i)/dble(nseed-1)
      right=dble(local_i+1)/dble(nseed-1)
      if(left.gt.0.1231d0 .and. left.lt.0.1233d0) seen=.true.
      if(right.gt.0.1231d0 .and. right.lt.0.1233d0) seen=.true.
      midpoint=0.5d0*(left+right)
      if(midpoint.gt.0.1231d0 .and. midpoint.lt.0.1233d0) seen=.true.
    enddo
    call require(.not.seen,'70c7f65 midpoint scan happened to see the island')
  end subroutine test_old_midpoint_misses_island

  subroutine test_certified_narrow_island
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
    topology_context_h=42.75d0
    topology_partition_tol=1.d-6
    call sample_matrix_out_partitioned_certified( &
        manufactured_narrow_matrix,manufactured_narrow_boundaries,local_ierr, &
        manufactured_integrand_envelope)
    call require(local_ierr.eq.0,'certified narrow A-B-A topology failed')
    call require(topology_transition_count.eq.2, &
                 'certified narrow topology transition count is wrong')
    integral=0.d0
    do j=2,npoi
      if(topology_arr(j).eq.topology_arr(j-1)) then
        integral=integral+0.5d0*(amat_arr(1,1,j)+amat_arr(1,1,j-1)) &
            *(xarr(j)-xarr(j-1))
      endif
    enddo
    exact=0.5d0+(0.1233d0-0.1231d0)
    call require(abs(integral-exact).lt.1.d-4, &
                 'certified narrow A-B-A analytic integral is wrong')
    call require(topology_gap_measure.le.topology_gap_bound+1.d-12, &
                 'narrow transition gap exceeded its certificate bound')
    call require(topology_contribution_error_certified, &
                 'narrow-island contribution bound was not certified')
    topology_partition_tol=1.d-8
  end subroutine test_certified_narrow_island

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

  subroutine manufactured_piecewise_boundaries(nmax,nfound,boundaries,local_ierr)
    integer, intent(in) :: nmax
    integer, intent(out) :: nfound,local_ierr
    double precision, intent(out) :: boundaries(nmax)

    nfound=2
    local_ierr=0
    boundaries=0.d0
    ! This first certified root is topology-inert: both open-side signatures
    ! remain 7.  It must not create a contribution gap because the sampler
    ! does not split the matrix interval there.
    boundaries(1)=1.d-8
    boundaries(2)=0.5d0
  end subroutine manufactured_piecewise_boundaries

  subroutine manufactured_narrow_matrix
    use sample_matrix_out_mod, only : x,amat,topology_signature,topology_error

    topology_error=0
    if(x.gt.0.1231d0 .and. x.lt.0.1233d0) then
      topology_signature=13
      amat(1,1)=x+1.d0
    else
      topology_signature=7
      amat(1,1)=x
    endif
  end subroutine manufactured_narrow_matrix

  subroutine manufactured_narrow_boundaries(nmax,nfound,boundaries,local_ierr)
    integer, intent(in) :: nmax
    integer, intent(out) :: nfound,local_ierr
    double precision, intent(out) :: boundaries(nmax)

    nfound=2
    local_ierr=0
    boundaries=0.d0
    boundaries(1)=0.1231d0
    boundaries(2)=0.1233d0
  end subroutine manufactured_narrow_boundaries

  subroutine manufactured_integrand_envelope(gap_lo,gap_hi,envelope,local_ierr)
    double precision, intent(in) :: gap_lo,gap_hi
    double precision, intent(out) :: envelope
    integer, intent(out) :: local_ierr

! The manufactured matrix values are bounded in absolute value by four on
! the full [0,1] domain.  This is an independent envelope oracle, not a
! sample-derived maximum.
    envelope=4.d0
    local_ierr=0
    if(gap_hi.lt.gap_lo) local_ierr=1
  end subroutine manufactured_integrand_envelope

  subroutine manufactured_unavailable_envelope(gap_lo,gap_hi,envelope,local_ierr)
    double precision, intent(in) :: gap_lo,gap_hi
    double precision, intent(out) :: envelope
    integer, intent(out) :: local_ierr

    envelope=0.d0
    local_ierr=sample_matrix_contribution_unresolved
  end subroutine manufactured_unavailable_envelope

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
    call sample_matrix_out_partitioned_certified( &
        manufactured_invalid_endpoint,manufactured_no_boundaries,local_ierr, &
        manufactured_integrand_envelope)
    call require(local_ierr.eq.0,'invalid global endpoint was not contracted')
    call require(abs(xarr(1)-0.1d0).lt.1.d-6, &
                 'endpoint contraction was not refined to the valid boundary')
  end subroutine test_endpoint_contraction

  subroutine test_endpoint_failure_classification
    call require(matrix_failure_is_open_boundary(matrix_eval_orbit_failure), &
                 'orbit-return failure was not classified as an open endpoint')
    call require(matrix_failure_is_open_boundary(matrix_eval_wall_loss), &
                 'wall loss was not classified as an open endpoint')
    call require(matrix_failure_is_open_boundary(matrix_eval_cut_domain), &
                 'physical cut/kinetic domain was not classified as an open endpoint')
    call require(.not.matrix_failure_is_open_boundary(matrix_eval_starter_failure), &
                 'starter failure was misclassified as an open endpoint')
    call require(.not.matrix_failure_is_open_boundary(matrix_eval_nonfinite), &
                 'nonfinite matrix was misclassified as an open endpoint')
  end subroutine test_endpoint_failure_classification

  subroutine test_fixedpoint_quadratic
! Independent algebraic oracle for the exact fixed-point discriminant used by
! the outer certificate.  It covers the two-root, linear, and no-real-root
! branches without depending on a field initialization or on the sampler.
    call solve_fixedpoint_quadratic(1.d0,-1.d0,0.1875d0,quadratic_roots, &
                                    nquadratic,quadratic_ok)
    call require(quadratic_ok .and. nquadratic.eq.2, &
                 'fixed-point quadratic did not return two roots')
    call require(abs(quadratic_roots(1)-0.25d0).lt.1.d-14 .and. &
                 abs(quadratic_roots(2)-0.75d0).lt.1.d-14, &
                 'fixed-point quadratic roots are wrong')
    call solve_fixedpoint_quadratic(0.d0,2.d0,-1.d0,quadratic_roots, &
                                    nquadratic,quadratic_ok)
    call require(quadratic_ok .and. nquadratic.eq.1 .and. &
                 abs(quadratic_roots(1)-0.5d0).lt.1.d-14, &
                 'fixed-point linear limit is wrong')
    call solve_fixedpoint_quadratic(1.d0,0.d0,1.d0,quadratic_roots, &
                                    nquadratic,quadratic_ok)
    call require(quadratic_ok .and. nquadratic.eq.0, &
                 'negative fixed-point discriminant was accepted')
  end subroutine test_fixedpoint_quadratic

  subroutine test_resonance_status
    integer :: local_status,local_ierr
    double precision :: local_factor

! A tangent root is an unresolved delta-function weight, not a zero-weight
! resonance.  The wall-zero branch is accepted only with a finite root factor
! and an exactly zero certified perturbation amplitude.
    call classify_resonance_root(0.d0,1.d0,1.d0,.false.,local_factor,local_status)
    call require(local_status.eq.resonance_status_tangent_root, &
                 'tangent resonance root was silently accepted')
    call classify_resonance_root(1.d0,2.d0,0.d0,.true.,local_factor,local_status)
    call require(local_status.eq.resonance_status_success .and. &
                 abs(local_factor-2.d0).lt.1.d-14, &
                 'finite certified wall-zero root was not accepted')
    call classify_resonance_root(1.d0,2.d0,1.d-30,.true.,local_factor,local_status)
    call require(local_status.eq.resonance_status_wall_not_zero, &
                 'nonzero wall amplitude was confused with certified zero')
    call classify_resonance_root(1.d-10,huge(1.d0),1.d0,.false.,local_factor,local_status)
    call require(local_status.eq.resonance_status_nonfinite_weight, &
                 'nonfinite resonance factor was silently zeroed')
    call require(resonance_status_from_root_error(0).eq.resonance_status_success, &
                 'successful root callback was marked failed')
    call require(resonance_status_from_root_error(2).eq.resonance_status_root_failure, &
                 'failed root callback was confused with a searched-zero result')

! The sampler must expose a callback failure rather than accepting a partial
! matrix; this is the independent outer-callback oracle.
    nlagr=3
    n1=1
    n2=1
    npoi=9
    itermax=4
    xbeg=0.d0
    xend=1.d0
    eps=1.d-10
    icount=0
    topology_context_h=45.d0
    call sample_matrix_out_partitioned_certified( &
        manufactured_callback_failure,manufactured_no_boundaries,local_ierr, &
        manufactured_integrand_envelope)
    call require(local_ierr.ne.0,'callback failure was accepted as a valid grid')
    call sample_matrix_out_partitioned_certified( &
        manufactured_piecewise_matrix,manufactured_piecewise_boundaries,local_ierr, &
        manufactured_unavailable_envelope)
    call require(local_ierr.eq.sample_matrix_contribution_unresolved, &
                 'unavailable contribution envelope was accepted')
  end subroutine test_resonance_status

  subroutine test_conjugate_mode_representation
    double precision :: hm_real,hm_imag,hm_squared,coefficient
    double precision :: hm_real_conj,hm_imag_conj,hm_squared_conj
    double precision :: delta_weight,torque_weight
    double precision :: delta_weight_conj,torque_weight_conj
    double precision :: omega_b,omega_t,resonance,resonance_conjugate

! Only the paired Fourier representation changes (m_2,m_3,A) to
! (-m_2,-m_3,conjg(A)).  This is not a claim that conjugating A at fixed
! (m_2,m_3) is a physical symmetry.  The generated H_m and Eq.(17) kernels
! must preserve the paired representation; no hand-normalized factor is used
! as the oracle here.
    call potato_hm_kernel(5.d0,1.d0,2.d0,0.75d0,2.d0,3.d0,0.4d0, &
                          hm_real,hm_imag,hm_squared,coefficient)
    call potato_hm_kernel(5.d0,1.d0,2.d0,0.75d0,2.d0,-3.d0,-0.4d0, &
                          hm_real_conj,hm_imag_conj,hm_squared_conj, &
                          coefficient)
    call require(abs(hm_squared-hm_squared_conj).lt.1.d-12, &
                 'conjugate Fourier representation changed H_m squared')
    call potato_resonance_torque_kernel(-2.d0,4.d0,3.d0,hm_squared,1.d0,1.d0, &
                                        1.d0,1.d0,1.d0,delta_weight, &
                                        torque_weight)
    call potato_resonance_torque_kernel(2.d0,-4.d0,-3.d0,hm_squared_conj,1.d0, &
                                        1.d0,1.d0,1.d0,1.d0, &
                                        delta_weight_conj,torque_weight_conj)
    call require(abs(delta_weight-delta_weight_conj).lt.1.d-12 .and. &
                 abs(torque_weight-torque_weight_conj).lt.1.d-12, &
                 'conjugate Fourier representation changed torque weight')
    omega_b=1.7d0
    omega_t=-0.2d0
    resonance=2.d0*omega_b+3.d0*omega_t
    resonance_conjugate=(-2.d0)*omega_b+(-3.d0)*omega_t
    call require(abs(resonance+resonance_conjugate).lt.1.d-14, &
                 'conjugate Fourier representation changed resonance zero set')
  end subroutine test_conjugate_mode_representation

  subroutine test_gap_square_map
    double precision :: coordinate,jacobian

! The transformed coordinate resolves an integrable square-root endpoint:
! x = x_b + s*d*u^2, dx/du = 2*s*d*u.  These values are an independent
! behavioral oracle for the generated map, including both orientations.
    call potato_gap_square_map(1.5d0,-1.d0,0.4d0,0.5d0,coordinate,jacobian)
    call require(abs(coordinate-1.4d0).lt.1.d-14 .and. &
                 abs(jacobian-0.4d0).lt.1.d-14, &
                 'negative-direction quadratic gap map is wrong')
    call potato_gap_square_map(1.5d0,1.d0,0.4d0,0.5d0,coordinate,jacobian)
    call require(abs(coordinate-1.6d0).lt.1.d-14 .and. &
                 abs(jacobian-0.4d0).lt.1.d-14, &
                 'positive-direction quadratic gap map is wrong')
  end subroutine test_gap_square_map

  subroutine test_gap_sample_bound
    double precision :: average,limit_integral,limit_coefficient
    integer :: j

! The production provider must use the completed branch grid, not launch new
! resonance solves.  A manufactured inverse-square-root branch gives the same
! exact average as the limiting expression.
    if(allocated(xarr)) deallocate(xarr)
    if(allocated(amat_arr)) deallocate(amat_arr)
    if(allocated(topology_arr)) deallocate(topology_arr)
    n1=1
    n2=1
    npoi=4
    xbeg=0.d0
    xend=1.d0
    allocate(xarr(npoi),amat_arr(n1,n2,npoi),topology_arr(npoi))
    xarr=(/0.25d0,0.5d0,0.75d0,1.d0/)
    topology_arr=7
    do j=1,npoi
      amat_arr(1,1,j)=1.d0/sqrt(xarr(j))
    enddo
    call estimate_topology_gap_from_samples(0.d0,0.25d0,average,ierr)
    call require(ierr.eq.0,'sampled branch gap bound failed')
    call require(abs(average-4.d0).lt.1.d-12, &
                 'sampled branch inverse-square-root bound is wrong')
    call potato_gap_sqrt_integral(1.d0,0.25d0,limit_integral)
    call require(abs(limit_integral-1.d0).lt.1.d-14, &
                 'inverse-square-root limiting integral is wrong')
    call potato_gap_sqrt_coefficient(2.d0,0.25d0,limit_coefficient)
    call require(abs(limit_coefficient-1.d0).lt.1.d-14, &
                 'inverse-square-root limiting coefficient is wrong')
    deallocate(xarr,amat_arr,topology_arr)
  end subroutine test_gap_sample_bound

  subroutine test_symbolic_contract
    double precision :: candidate,positive_bound,derivative
    double precision :: hm_real,hm_imag,hm_squared,coefficient
    double precision :: hm_real_sign,hm_imag_sign,hm_squared_sign
    double precision :: hm_real_twice,hm_imag_twice,hm_squared_twice
    double precision :: hm_real_phase,hm_imag_phase,hm_squared_phase
    double precision :: root_weight,gap_measure,gap_bound
    double precision :: delta_weight_base,torque_weight_base
    double precision :: delta_weight_sign,torque_weight_sign
    double precision :: delta_weight_twice,torque_weight_twice
    double precision :: delta_weight_phase,torque_weight_phase
    double precision :: hdet,lambda_local,c_tau,regular_offset,regular_jacobian
    double precision :: xpoint_offset,xpoint_jacobian,regular_tau
    double precision :: separatrix_tau,xpoint_tau
    double precision :: separatrix_tau_ref2,xpoint_tau_ref2
    double precision :: alpha,rotated_real,rotated_imag
    double precision :: frequency,frequency_prime,frequency_prime_root
    double precision :: frequency_remainder,frequency_delta_n2
    double precision :: collapsed_n_tau2
    double precision :: extrapolated_value,extrapolated_derivative
    integer :: limit_ierr

! Independent action-domain oracle with both qPhi' and B' nonzero.  elefie
! supplies qPhi' with its physical positive sign; the generated chain rule
! supplies the minus sign from d(H-qPhi)/dR.
    call potato_jperp_kernel(5.d0,1.d0,2.d0,3.d0,0.4d0,candidate, &
                             positive_bound,derivative)
    call require(abs(candidate-2.d0).lt.1.d-14, &
                 'generated Jperp candidate is wrong')
    call require(abs(positive_bound-2.d0).lt.1.d-14, &
                 'generated positive Jperp is wrong')
    call require(abs(derivative+1.9d0).lt.1.d-14, &
                 'generated Jperp chain rule has the wrong Phi/B-prime sign')
    call potato_jperp_kernel(5.d0,6.d0,2.d0,3.d0,0.4d0,candidate, &
                             positive_bound,derivative)
    call require(candidate.lt.0.d0 .and. positive_bound.eq.0.d0, &
                 'qPhi sign was silently clipped before the positive part')

! A sign, scale, phase, and paired-conjugation oracle for Eq. (4) H_m.
    alpha=0.73d0
    call potato_hm_kernel(5.d0,1.d0,2.d0,0.75d0,1.25d0,-0.75d0,0.31d0, &
                          hm_real,hm_imag,hm_squared,coefficient)
    call potato_hm_kernel(5.d0,1.d0,2.d0,0.75d0,-1.25d0,0.75d0,0.31d0, &
                          hm_real_sign,hm_imag_sign,hm_squared_sign, &
                          coefficient)
    call require(abs(hm_squared-hm_squared_sign).lt.1.d-12, &
                 'A -> -A changed H_m squared')
    call potato_hm_kernel(5.d0,1.d0,2.d0,0.75d0,2.5d0,-1.5d0,0.31d0, &
                          hm_real_twice,hm_imag_twice,hm_squared_twice, &
                          coefficient)
    call require(abs(hm_squared_twice-4.d0*hm_squared).lt.1.d-12, &
                 'A -> 2A did not scale H_m squared by four')
    rotated_real=1.25d0*cos(alpha)-(-0.75d0)*sin(alpha)
    rotated_imag=1.25d0*sin(alpha)+(-0.75d0)*cos(alpha)
    call potato_hm_kernel(5.d0,1.d0,2.d0,0.75d0,rotated_real,rotated_imag, &
                          0.31d0,hm_real_phase,hm_imag_phase,hm_squared_phase, &
                          coefficient)
    call require(abs(hm_squared_phase-hm_squared).lt.1.d-12, &
                 'global complex phase changed H_m squared')

    call potato_resonance_torque_kernel(-2.d0,4.d0,3.d0,hm_squared,1.d0, &
                                        1.d0,1.d0,1.d0,1.d0, &
                                        delta_weight_base,torque_weight_base)
    call potato_resonance_torque_kernel(-2.d0,4.d0,3.d0,hm_squared_sign,1.d0, &
                                        1.d0,1.d0,1.d0,1.d0, &
                                        delta_weight_sign,torque_weight_sign)
    call potato_resonance_torque_kernel(-2.d0,4.d0,3.d0,hm_squared_twice,1.d0, &
                                        1.d0,1.d0,1.d0,1.d0, &
                                        delta_weight_twice,torque_weight_twice)
    call potato_resonance_torque_kernel(-2.d0,4.d0,3.d0,hm_squared_phase,1.d0, &
                                        1.d0,1.d0,1.d0,1.d0, &
                                        delta_weight_phase,torque_weight_phase)
    call require(abs(delta_weight_base-delta_weight_sign).lt.1.d-12 .and. &
                 abs(torque_weight_base-torque_weight_sign).lt.1.d-12, &
                 'A -> -A changed torque weight')
    call require(abs(torque_weight_twice-4.d0*torque_weight_base).lt.1.d-12, &
                 'A -> 2A did not scale torque by four')
    call require(abs(delta_weight_base-delta_weight_phase).lt.1.d-12 .and. &
                 abs(torque_weight_base-torque_weight_phase).lt.1.d-12, &
                 'global complex phase changed torque weight')

    call potato_root_jacobian_kernel(-2.d0,4.d0,root_weight)
    call require(abs(root_weight-0.5d0).lt.1.d-14, &
                 'generated resonance Jacobian is wrong')
    call potato_gap_contribution_kernel(3.d0,0.1d0,0.4d0,gap_measure, &
                                        gap_bound)
    call require(abs(gap_measure-0.3d0).lt.1.d-14 .and. &
                 abs(gap_bound-0.9d0).lt.1.d-14, &
                 'generated topology contribution bound is wrong')

    call potato_class_linear_extrapolation(7.d0,4.d0,0.5d0,0.25d0,2.d0,3.d0, &
                                           extrapolated_value,extrapolated_derivative)
    call require(abs(extrapolated_value-6.5d0).lt.1.d-14 .and. &
                 abs(extrapolated_derivative-6.d0).lt.1.d-14, &
                 'generated linear class endpoint continuation is wrong')
    call potato_class_quadratic_extrapolation(7.d0,4.d0,0.75d0,0.25d0,2.d0,3.d0, &
                                              1.d0,extrapolated_value, &
                                              extrapolated_derivative)
    call require(abs(extrapolated_value-5.d0).lt.1.d-14 .and. &
                 abs(extrapolated_derivative-18.d0).lt.1.d-14, &
                 'generated quadratic class endpoint continuation is wrong')

! q=Delta_phi+2*pi*m/n and F_freq=n*q/tau.  At q=0 the q*tau' term
! vanishes, so the original n^2 delta factor reduces before the transport tau
! is applied.  The generated kernel must contain only one |n| in the result.
    call potato_frequency_reduction_kernel(3.d0,0.d0,2.d0,4.d0, &
                                           0.7d0,frequency,frequency_prime, &
                                           frequency_prime_root, &
                                           frequency_remainder,frequency_delta_n2, &
                                           collapsed_n_tau2)
    call require(abs(frequency).lt.1.d-14, &
                 'frequency resonance did not vanish at q=0')
    call require(abs(frequency_prime-frequency_prime_root).lt.1.d-14 .and. &
                 abs(frequency_remainder).lt.1.d-14, &
                 'F_freq prime did not reduce to n*q-prime/tau at the root')
    call require(abs(frequency_delta_n2-6.d0).lt.1.d-14, &
                 'original n-squared delta factor was not reduced correctly')
    call require(abs(collapsed_n_tau2-24.d0).lt.1.d-14, &
                 'transport weight has an extra or missing n factor')
    call potato_frequency_reduction_kernel(-3.d0,0.d0,2.d0,4.d0, &
                                           0.7d0,frequency,frequency_prime, &
                                           frequency_prime_root, &
                                           frequency_remainder,frequency_delta_n2, &
                                           collapsed_n_tau2)
    call require(abs(collapsed_n_tau2-24.d0).lt.1.d-14, &
                 'signed n changed the paired absolute reduction')

! The limit kernel takes the local Hessian and cut maps, not a fitted C.
! det(H'')=-6, lambda=sqrt(6), and C_tau=1/lambda in the supplied time unit.
    call potato_limiting_kernel(2.d0,0.d0,-3.d0,7.d0,2.d0,4.d0,1.d-3,1.d0, &
                                hdet,lambda_local,c_tau,regular_offset, &
                                regular_jacobian,xpoint_offset,xpoint_jacobian, &
                                regular_tau,separatrix_tau,xpoint_tau)
    call require(abs(hdet+6.d0).lt.1.d-14 .and. lambda_local.gt.0.d0, &
                 'local Hamiltonian Hessian saddle rate is wrong')
    call require(abs(c_tau*lambda_local-1.d0).lt.1.d-12, &
                 'homoclinic coefficient is not inverse saddle rate')
    call require(abs(regular_tau-7.d0).lt.1.d-14 .and. &
                 abs(regular_offset-2.d-3).lt.1.d-14 .and. &
                 abs(regular_jacobian-2.d0).lt.1.d-14, &
                 'regular boundary limit/map is wrong')
    call require(abs(xpoint_offset-4.d-6).lt.1.d-14 .and. &
                 abs(xpoint_jacobian-8.d-3).lt.1.d-14, &
                 'X-point action map/Jacobian is wrong')
    call require(abs(xpoint_tau-2.d0*separatrix_tau).lt.1.d-12, &
                 'X-point logarithmic limit is not twice separatrix limit')

    call potato_limiting_kernel_checked(2.d0,0.d0,-3.d0,7.d0,2.d0,4.d0, &
                                        1.d-3,2.d0,hdet,lambda_local,c_tau, &
                                        regular_offset,regular_jacobian, &
                                        xpoint_offset,xpoint_jacobian, &
                                        regular_tau,separatrix_tau_ref2, &
                                        xpoint_tau_ref2,limit_ierr)
    call require(limit_ierr.eq.0 .and. &
                 abs((separatrix_tau_ref2-separatrix_tau)-c_tau*log(2.d0)) &
                 .lt.1.d-12, 'limiting log did not use delta_R_ref ratio')

    call potato_limiting_kernel_checked(2.d0,0.d0,-3.d0,7.d0,2.d0,4.d0, &
                                        1.d-3,0.d0,hdet,lambda_local,c_tau, &
                                        regular_offset,regular_jacobian, &
                                        xpoint_offset,xpoint_jacobian, &
                                        regular_tau,separatrix_tau,xpoint_tau, &
                                        limit_ierr)
    call require(limit_ierr.eq.potato_limiting_invalid_reference, &
                 'nonpositive limiting reference scale was accepted')
    call potato_limiting_kernel_checked(2.d0,0.d0,-3.d0,7.d0,2.d0,4.d0, &
                                        0.d0,1.d0,hdet,lambda_local,c_tau, &
                                        regular_offset,regular_jacobian, &
                                        xpoint_offset,xpoint_jacobian, &
                                        regular_tau,separatrix_tau,xpoint_tau, &
                                        limit_ierr)
    call require(limit_ierr.eq.potato_limiting_invalid_distance, &
                 'zero limiting distance was accepted')
  end subroutine test_symbolic_contract

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

  subroutine manufactured_callback_failure
    use sample_matrix_out_mod, only : x,amat,topology_signature,topology_error

    if(x.gt.0.4d0 .and. x.lt.0.6d0) then
      topology_signature=0
      topology_error=17
      amat(1,1)=0.d0
    else
      topology_signature=11
      topology_error=0
      amat(1,1)=x
    endif
  end subroutine manufactured_callback_failure

  subroutine manufactured_no_boundaries(nmax,nfound,boundaries,local_ierr)
    integer, intent(in) :: nmax
    integer, intent(out) :: nfound,local_ierr
    double precision, intent(out) :: boundaries(nmax)

    nfound=0
    local_ierr=0
    boundaries=0.d0
  end subroutine manufactured_no_boundaries

end program test_endpoint_stationary_root
