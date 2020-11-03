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

  real(8), public :: qi, mi      ! Ion charge and mass, must be initialized !

  integer(4), parameter :: neqm=5  ! Number of equations of motion

  contains

  function linspace(a, b, cnt)
    integer :: cnt
    real(8) :: linspace(cnt)
    integer :: i
    real(8) :: a, b, delta

    delta = (b-a)/(cnt-1)
    linspace = a + delta*(/(i,i=0,cnt-1)/)
  end function linspace

end module common
