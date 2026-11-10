program test_int_powers
! int_powers must reproduce x**i (integer variable i, i.e. libgcc __powidf2)
! bit for bit, so the profile polynomials keep their exact values.  The
! exponent is volatile so the reference is the runtime __powidf2 call, as in
! the original polynomial loops, not a compile-time expansion of x**i.
  use phielec_of_psi_mod, only : npolyphi, int_powers
  implicit none

  double precision :: x, pw(0:npolyphi)
  integer :: i, k
  integer, volatile :: iv

  do k=0,20000
    x=-1.3d0+dble(k)*1.37d-4
    call int_powers(x,pw)
    do i=0,npolyphi
      iv=i
      if(transfer(pw(i),1_8).ne.transfer(x**iv,1_8)) error stop 'int_powers differs from x**i'
    enddo
  enddo
end program test_int_powers
