module common
  ! Contains common constants in CGS units,
  ! and some commonly used parameters and subroutines.

  implicit none
  save

  complex(8), parameter :: imun=(0d0,1d0)

  real(8), parameter, public   :: pi = 4*atan(1d0)

  real(8), parameter, public ::  &
    qe  = 4.803204d-10,       &  ! Elementary charge
    me  = 9.109382d-28,       &  ! Electron mass
    mu  = 1.660538d-24,       &  ! 1u in g
    c   = 2.997925d+10,       &  ! Speed of light in cm/s
    kb  = 1.381649d-16,       &  ! Boltzmann constant
    eV  = 1.602176d-12           ! 1 electron volt

  real(8), public :: qi, mi

  integer(8), parameter :: neqm=5  ! Number of equations of motion
end module common
