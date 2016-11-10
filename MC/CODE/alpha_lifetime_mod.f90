!
  module chamb_mod
    double precision :: rbig,rcham2
  end module chamb_mod
!
  module parmot_mod
    double precision :: rmu,ro0
  end module parmot_mod
!
  module collis_alp
    integer, parameter :: nsorts=3, ns=10000 !original: ns=10000
    integer :: iswmod
    logical :: swcoll=.false.
    double precision, dimension(nsorts)    :: efcolf,velrat,enrat
    double precision, dimension(nsorts,ns) :: efcolf_arr,velrat_arr,enrat_arr
  end module collis_alp
!
