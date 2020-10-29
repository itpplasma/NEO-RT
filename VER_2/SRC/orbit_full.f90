module orbit

  contains

  subroutine timestep(tau, z, vz)
  ! See velo for details
    integer(4) :: n  ! Number of equations
    real(8) :: tau   ! Time
    real(8) :: z(5)  ! Phase-position
    real(8) :: vz(5) ! Phase-velocity

    call velo(tau, z, vz)
  end subroutine timestep

end module orbit
