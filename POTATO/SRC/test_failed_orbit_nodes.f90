program test_failed_orbit_nodes
! A class may contain orbits that leave the field domain.  get_matrix reports
! those through ierr_get_matrix and leaves no usable value behind, so the grid
! refinement must not read them as structure and keep bisecting towards them.
!
! The stub below is linear (exactly reproduced by the Lagrange interpolation
! used for the split test) outside a failed stretch, so every split the
! refinement performs is driven by the failed nodes alone.
  use sample_matrix_mod, only : nlagr,n1,n2,npoi,itermax,x,amat,xarr,xbeg,xend, &
                                eps,valid_arr,ierr_interpolate,find_valid_class_roots
  use find_all_roots_mod, only : nroots,roots,customgrid,nsearch_min
  use interp_cache_mod, only : interp_cache_reset
  implicit none

! Refinement starts from nlagr+2 nodes.  Bisecting towards the failed stretch
! adds nodes on every one of the itermax passes; staying under this bound means
! the failed nodes were excluded rather than chased.
  integer, parameter :: npoi_max_expected = 20
  integer :: ierr,i,nodes_in_gap
  character(len=64) :: msg
  double precision :: vec(1),dvec(1)

  external :: stub_matrix,stub_quartic

  nlagr   = 7
  n1      = 1
  n2      = 1
  itermax = 20
  eps     = 1.d-3
  xbeg    = 0.d0
  xend    = 1.d0

  call sample_matrix(stub_matrix,ierr)

  call require(ierr.eq.0,'sample_matrix reported failure')

  write(msg,'(A,I0)') 'grid grew chasing failed nodes, npoi = ',npoi
  call require(npoi.le.npoi_max_expected,trim(msg))

! No node may be inserted inside the failed stretch: nothing there is worth
! resolving, and the values carry no physical weight.
  nodes_in_gap = 0
  do i = 1,npoi
    if(xarr(i).gt.0.4d0 .and. xarr(i).lt.0.6d0) nodes_in_gap = nodes_in_gap+1
  end do
  write(msg,'(A,I0)') 'nodes inserted inside failed stretch = ',nodes_in_gap
  call require(nodes_in_gap.le.2,trim(msg))

! The interior failed node survives endpoint trimming. It must be rejected
! by the production interpolator, including before any cache lookup.
  call require(count(.not.valid_arr).gt.0,'stub did not fail an interior node')
  call interp_cache_reset
  call interpolate_class_doublecount(0.5d0,vec,dvec)
  call require(ierr_interpolate.ne.0,'interpolator consumed a failed stencil')

! Quadratic interpolation leaves valid windows on both sides of this same grid.
! Keep the genuine roots at .1 and .9; exclude a third root inside the gap.
  nlagr=2
  customgrid=.false.
  nsearch_min=100
  call find_valid_class_roots(residual,xbeg,xend,ierr)
  call require(ierr.eq.0,'valid-window root search failed')
  call require(nroots.eq.2,'root search bridged the failed interior region')
  call require(maxval(abs(roots-[0.1d0,0.9d0])).lt.1.d-8, &
      'valid resonances were lost or moved')

! Exercise rebuilding the mask and retaining the final endpoint through a
! refinement pass. A quartic requires a completed refinement pass.
  nlagr=2
  itermax=1
  eps=1.d-12
  call sample_matrix(stub_quartic,ierr)
  call require(npoi.gt.nlagr+2,'endpoint test did not refine')
  do i=1,npoi
    call require(valid_arr(i), &
        'refinement lost a retained node validity bit')
  enddo

  print *,'test_failed_orbit_nodes: OK  (npoi = ',npoi,')'

contains

  subroutine residual(xval,f,df)
    double precision, intent(in) :: xval
    double precision, intent(out) :: f,df
    double precision :: v(1),dv(1),u

    call interpolate_class_doublecount(xval,v,dv)
    call require(ierr_interpolate.eq.0,'root search evaluated an invalid window')
    call require(abs(v(1)-(1.d0+xval)).lt.1.d-12, &
        'valid interpolation changed the linear oracle')
    u=xval-0.5d0
    f=u*(u*u-0.16d0)
    df=3.d0*u*u-0.16d0
  end subroutine residual

  subroutine require(ok,text)
    logical, intent(in) :: ok
    character(len=*), intent(in) :: text

    if(ok) return
    print *,text
    error stop 1
  end subroutine require

end program test_failed_orbit_nodes

subroutine stub_matrix
! Stands in for get_matrix_doublecount: linear where the orbit closes, and a
! reported failure with no usable value where it leaves the domain.
  use sample_matrix_mod, only : x,amat,ierr_get_matrix

  implicit none

  if(x.gt.0.4d0 .and. x.lt.0.6d0) then
    amat = 0.d0
    ierr_get_matrix = 1
    return
  endif

  ierr_get_matrix = 0
  amat(1,1) = 1.d0+x
end subroutine stub_matrix

subroutine stub_quartic
  use sample_matrix_mod, only : x,amat,ierr_get_matrix
  implicit none

  ierr_get_matrix=0
  amat(1,1)=x**4
end subroutine stub_quartic
