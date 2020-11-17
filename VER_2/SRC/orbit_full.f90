submodule (orbit) orbit_full
  implicit none

  contains

  module procedure timestep
    call velo(tau, z, vz)
  end procedure timestep

  module procedure bounce_average
    real(8) :: dtau

    dtau = 2.0d0  ! TODO: time step is for strongly passing

    ret = 0.0d0

    call find_bounce(n, timestep_ext, dtau, z, taub, delphi, ret)

    contains

    subroutine timestep_ext(t, y, dydt)
      real(8), intent(in) :: t
      real(8), intent(in) :: y(neqm+n)
      real(8), intent(out) :: dydt(neqm+n)

      call timestep(t, y(1:neqm), dydt(1:neqm))
      call integrand(t, y(1:neqm), dydt(neqm+1:(neqm+n)))
    end subroutine timestep_ext

  end procedure bounce_average

  module procedure bounce_integral_box
    use dvode_f90_m
    use field_eq_mod, only: psif, psi_axis, psi_sep

    real(8) :: dtau, taub, delphi
    real(8) :: sprev, snext ! Previous and next flux radius box
    real(8) :: y(neqm+n)

    integer(4) :: nmax  ! Maximum loop iterations

    integer(4) :: k
    real(8) :: ti

    real(8) :: atol(neqm+n), rtol, tout, rstats(22)
    integer(4) :: neq, itask, istate, istats(31), numevents
    type (vode_opts) :: options
    real(8) :: bmod, phi_elec, s
    real(8) :: sold, told, yold(neqm+n)
    integer(4) :: sind, sind0 ! s index

    integer(4) :: jroots(2)

    ! TODO make this more efficient
    dtau = 2.0d0  ! TODO: time step is for strongly passing
    ret = 0.0d0
    y = 0.0d0
    y(1:neqm) = z
    call find_bounce(n, timestep_ext_forbounce, dtau, z, taub, delphi, ret)

    neq = neqm + n
    rtol = 1d-12
    atol = 1d-13
    itask = 1
    istate = 1
    numevents = 2
    options = set_normal_opts(abserr_vector=atol, relerr=rtol, nevents=numevents)

    nmax = 3*size(sbox)
    tau = 0d0
    ti = 0d0

    call get_bmod_and_Phi(z(1:3), bmod, phi_elec)

    s = abs((psif-psi_axis)/(psi_sep-psi_axis))
    sind = size(sbox)+1
    sprev = -1e5
    snext = 1e5
    do k = 1,size(sbox)
        if(sbox(k)>s) then
          sind = k
          snext = sbox(k)
          if(k>1) sprev = sbox(k-1)
          exit
        end if
    enddo
    sind0 = sind

    told = 0d0
    do k = 1,nmax
        yold = y
        sold = s
        tout = taub
        call dvode_f90( &
          timestep_ext, neq, y, ti, tout, itask, istate, options, g_fcn = sroots &
        )
        if (istate == 2) exit
        if (istate == 3) then
          tau(sind) = tau(sind) + ti-told
          told = ti
          call get_stats(rstats, istats, numevents, jroots)
          if (jroots(2).ne.0) then
              sind = sind + 1 ! moving outwards
              sprev = snext
              if (sind == size(sbox)+1) then
                snext = 1e5
              else
                snext = sbox(sind)
              end if
          end if
          if (jroots(1).ne.0) then
              sind = sind - 1 ! moving inwards
              snext = sprev
              if (sind == 1) then
                sprev = -1e5
              else
                sprev = sbox(sind-1)
              end if
          end if
        end if
    end do

    tau(sind) = tau(sind) + taub-told

    contains

    subroutine timestep_ext_forbounce(t, y2, dydt)
      real(8), intent(in) :: t
      real(8), intent(in) :: y2(neqm+n)
      real(8), intent(out) :: dydt(neqm+n)

      call timestep(t, y2(1:neqm), dydt(1:neqm))
      call integrand(y2(1:neqm), dydt(neqm+1:(neqm+n)))
    end subroutine timestep_ext_forbounce

    subroutine timestep_ext(neq2, t, y2, dydt)
      integer(4), intent(in) :: neq2
      real(8), intent(in) :: t
      real(8), intent(in) :: y2(neq2)
      real(8), intent(out) :: dydt(neq2)

      call timestep(t, y2(1:neqm), dydt(1:neqm))
      call integrand(y2(1:neqm), dydt(neqm+1:neq2))
    end subroutine timestep_ext

    subroutine sroots (neq2, t, y2, ng, gout)
    ! For finding roots between boxes

      integer, intent(in) :: neq2, ng
      real(8), intent(in) :: t, y2(neq)
      real(8), intent(out) :: gout(ng)

      call get_bmod_and_Phi(y2(1:3), bmod, phi_elec)

      s = abs((psif-psi_axis)/(psi_sep-psi_axis))

      gout(1) = s - sprev
      gout(2) = s - snext
    end subroutine sroots
  end procedure bounce_integral_box


  module procedure bounce_harmonic
    real(8) :: omb, omphi
    ! TODO
    call bounce_average(next, z, fn, taub, delphi, ret)

    omb = 2d0*pi/taub
    omphi = delphi/taub

    call bounce_average(next, z, integrand, taub, delphi, ret)

    contains

    subroutine integrand(t, y, res)
      ! Integrand for Fourier integral: real and imaginary part
        real(8), intent(in)  :: t           ! Orbit time parameter
        real(8), intent(in)  :: y(5)        ! Orbit phase-space variables
        real(8), intent(out) :: res(2)      ! Output

        real(8) :: bmod, sqrtg
        real(8), dimension(3) :: bder,hcovar,hctrvr,hcurl

        real(8) :: fnres(2)
        complex(8) :: fnval, expfac, resval

        call fn(t, y, fnres)
        fnval = cmplx(fnres(1), fnres(2))
        expfac = exp(-imun*(mb*omb*t + nph*(y(2) - omphi*t)))
        resval = fnval*expfac

        res(1) = real(resval)
        res(2) = aimag(resval)
      end subroutine integrand
  end procedure bounce_harmonic


end submodule orbit_full
