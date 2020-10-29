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

  subroutine bounce_average(n, integrand, out)
    integer(4), intent(in) :: n
    external integrand
    real(8), intent(out) :: out(n)

    ! TODO
    out = 0.0d0

    contains

    subroutine timestep_ext(neq, t, y, dydt)
      use orbit, only: timestep

      integer(4), intent(in) :: neq
      real(8), intent(in) :: t
      real(8), intent(in) :: y(neq)
      real(8), intent(out) :: dydt(neq)

      call timestep(t, y(1:5), dydt(1:5))
      call integrand(y(1:5), dydt(5:(5+n)))
    end subroutine timestep_ext

  end subroutine bounce_average
end module common
