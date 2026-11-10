program test_interp_cache
! The hash index of interp_cache_mod must return exactly the entry the former
! linear scan (first i with xcache(i).eq.x) returned, across table growth,
! resets, signed zeros, NaN and keys that differ only in low mantissa bits.
  use, intrinsic :: ieee_arithmetic, only : ieee_value, ieee_quiet_nan
  use sample_matrix_mod, only : n1
  use interp_cache_mod,  only : ncache, xcache, interp_cache_find, interp_cache_reset
  implicit none

  integer, parameter :: nkey=3*2000
  double precision :: keys(nkey), vec(3), dvec(3), nan
  integer :: k

  n1=3
  vec=0.d0
  dvec=0.d0
  do k=1,2000
    keys(k)=dble(k)
    keys(2000+k)=1.d0+dble(k)*epsilon(1.d0)
    keys(4000+k)=dble(k)*1.d-300
  enddo
  keys(1)=0.d0

  do k=1,nkey
    call interp_cache_store(keys(k),vec,dvec)
  enddo
  if(ncache.ne.nkey) error stop 'not every key was stored'
  do k=1,nkey
    call require(keys(k), k)
  enddo
  call require(-0.d0, 1)
  call require(0.5d0, 0)
  call require(-1.d0, 0)
  call require(nearest(keys(4500),2.d0), 0)

  nan=ieee_value(nan, ieee_quiet_nan)
  call interp_cache_store(nan,vec,dvec)
  call require(nan, 0)

  call interp_cache_reset
  call require(keys(5), 0)
  call interp_cache_store(keys(5),vec,dvec)
  call require(keys(5), 1)
  call require(keys(6), 0)

contains

  subroutine require(x, expected)
    double precision, intent(in) :: x
    integer, intent(in) :: expected
    integer :: i, linear
    linear=0
    do i=1,ncache
      if(xcache(i).eq.x) then
        linear=i
        exit
      endif
    enddo
    if(linear.ne.expected) error stop 'test reference disagrees with expectation'
    if(interp_cache_find(x).ne.linear) error stop 'hash lookup differs from linear scan'
  end subroutine require

end program test_interp_cache
