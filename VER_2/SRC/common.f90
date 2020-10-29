module common
    implicit none
    save

    complex(8), parameter :: imun=(0d0,1d0)

    real(8), parameter, public   :: pi = 4*atan(1d0)

    real(8), parameter, public ::  &
         qe  = 4.803204d-10,       & ! elementary charge
         me  = 9.109382d-28,       & ! electron mass,
         mu  = 1.660538d-24,       & ! 1u
         c   = 2.997925d+10,       & ! speed of light
         kb  = 1.381649d-16,       & ! Boltzmann constant
         eV  = 1.602176d-12          ! 1 electron volt

    real(8), public :: qi, mi

    contains

    subroutine bounce_average(n, fcn)
        integer(4), intent(in) :: n
        external fcn
    end subroutine bounce_average
end module common
