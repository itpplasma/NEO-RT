!
  module parmot_mod
    logical :: gradpsiast=.false.             !<=NEW in version 4
    double precision :: rmu,ro0
    double precision :: dpsiast_dR,dpsiast_dZ !<=NEW in version 4
! gradpsiast and dpsiast_dR/dZ are per-orbit scratch: starter_doublecount sets
! gradpsiast=.true. around a velo call that writes dpsiast_dR/dZ, then resets it.
! The grid-build node-fill loops run starter+velo in parallel, so each thread
! needs its own copy.  rmu,ro0 are physical constants set once and read only,
! so they stay shared.
    !$omp threadprivate(gradpsiast,dpsiast_dR,dpsiast_dZ)
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
  module elefie_mod
    double precision :: rbig
    double precision, dimension(:,:), allocatable :: Mtprofile
    double precision, dimension(:,:), allocatable :: plasma
    double precision :: amb,am1,am2,Zb,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,v0
    double precision :: escale=1.0,bscale=1.0
  end module elefie_mod

  module constants
    double precision, parameter :: pi=3.14159265358979d0
    double precision,parameter  :: c=2.9979d10
    double precision,parameter  :: e_charge=4.8032d-10
    double precision,parameter  :: e_mass=9.1094d-28
    double precision,parameter  :: p_mass=1.6726d-24
    double precision,parameter  :: ev=1.6022d-12
  end module constants
