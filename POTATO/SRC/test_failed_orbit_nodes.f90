program test_failed_orbit_nodes
! A class may contain orbits that leave the field domain.  get_matrix reports
! those through ierr_get_matrix and leaves no usable value behind, so the grid
! refinement must not read them as structure and keep bisecting towards them.
!
! The stub below is linear (exactly reproduced by the Lagrange interpolation
! used for the split test) outside a failed stretch, so every split the
! refinement performs is driven by the failed nodes alone.
  use sample_matrix_mod, only : nlagr,n1,n2,npoi,itermax,x,amat,xarr,xbeg,xend, &
                                eps
  implicit none

! Refinement starts from nlagr+2 nodes.  Bisecting towards the failed stretch
! adds nodes on every one of the itermax passes; staying under this bound means
! the failed nodes were excluded rather than chased.
  integer, parameter :: npoi_max_expected = 20
  integer :: ierr,i,nodes_in_gap
  character(len=64) :: msg

  external :: stub_matrix

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

  print *,'test_failed_orbit_nodes: OK  (npoi = ',npoi,')'

contains

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
