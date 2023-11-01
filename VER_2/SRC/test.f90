!
!   use dvode_f90_m, only: dvode_f90, vode_opts, set_normal_opts

!   real(8) :: atol(neqm), rtol, tstart, tend
!   integer :: neq, itask, istate
!   type (vode_opts) :: options

!   integer(4) :: k
!
!   neq = neqm
!   rtol = 1d-9
!   atol = 1d-10
!   itask = 1
!   istate = 1
!   options = set_normal_opts(abserr_vector=atol, relerr=rtol, nevents=0)

!   tstart = 0.0d0
!   tend = 0.0d0

!   write(99, *) tend, z

!   do k = 0, 10
!     tstart = tend
!     tend = tstart + 10.0

!     call dvode_f90(timestep2, neq, z, tstart, tend, itask, istate, options)
!     !call odeint_allroutines(z, neqm, tstart, tend, atol, timestep)

!     write(99, *) tend, z

!   end do

!   contains

!   subroutine timestep2(neq, t, y, ydot)
!     ! Wrapper routine for timestep to work with VODE
!         integer, intent (in) :: neq
!         real(8), intent (in) :: t
!         real(8), intent (in) :: y(neq)
!         real(8), intent (out) :: ydot(neq)

!         call timestep(t, y, ydot)
!   end subroutine timestep2