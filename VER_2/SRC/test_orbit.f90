program test_orbit
    use util
    use parmot_mod, only: rmu, ro0
    use orbit, only: timestep, neqm, bounce_average, time_in_box

    implicit none

    real(8) :: z(neqm), zdot(neqm)
    real(8) :: bmod_ref, bmod00, tempi1, am, Zb, v0, rlarm
    real(8) :: bmod, phi_elec
    real(8) :: mth, mph  ! Poloidal and toroidal harmonic numbers
    real(8) :: alpha(3)  ! Invariants
    real(8) :: HmReIm(2), taub, delphi

    integer(4), parameter :: nbox = 101
    real(8) :: sbox(nbox)
    real(8) :: taubox(nbox)
    real(8) :: HmReImbox(2, nbox)

    call main

    contains

    subroutine main

      real(8) :: z(neqm), zdot(neqm)
      real(8) :: bmod_ref, bmod00, tempi1, am, Zb, v0, rlarm
      real(8) :: bmod, phi_elec
      real(8) :: mth, mph  ! Poloidal and toroidal harmonic numbers
      real(8) :: alpha(3)  ! Invariants
      real(8) :: HmReIm(2), taub, delphi, oneout(1)

      integer(4), parameter :: nbox = 101
      real(8) :: sbox(nbox)
      real(8) :: taubox(nbox)
      real(8) :: HmReImbox(2, nbox)

      integer :: npoicut
      real(8) :: rho_pol_max

      ! Inverse relativistic temperature, set >=1e5 to neglect relativistic effects
      rmu = 1.0d5

      ! This will translate to dimensionless units via parmot_mod.ro0
      ! Here we translate an input file in Tesla to computation in Gauss/CGS
      bmod_ref = 1.0d0             ! 1 Gauss in Gauss
      bmod00 = 1.0d0               ! On-axis magnetic field
      tempi1 = 0.17549561306d4     ! ion temperature
      am = 2.0d0                   ! Atomic mass 2 of deuterium ions
      Zb = 1.0d0                   ! Atomic charge 1 of deuterium ions

      ! Initialize common quantities for ions
      mi = am*mu
      qi = Zb*qe

      v0 = sqrt(2.0*tempi1*ev/(am*mu))     ! Reference (thermal) velocity
      ! Reference Larmor radius of thermal particles
      rlarm = v0*am*mu*c/(Zb*qe*bmod_ref)  ! Larmor radius in bmod_ref
      ro0 = rlarm*bmod00                   ! Rescaled Larmor radius

      ! Find Poincare cut:
      npoicut=10000     !number of equidistant points for interpolating the cut
      rho_pol_max=0.8d0 !maximum value of poloidal radius limiting the cut range

      call find_poicut(rho_pol_max,npoicut)

      call init_z(z)

      ! First call for field
      call get_bmod_and_Phi(z(1:3), bmod, phi_elec)
      print *, 'bmod, phielec:'
      print *, bmod, phi_elec

      print *, 'z(t=0)    = ', z
      call timestep(0.0d0, z, zdot)
      print *, 'zdot(t=0) = ', zdot

      ! Testing magnetic field routines
      call test_magfie

      ! Testing bounce average
      call bounce_average(1, z, one, taub, delphi, oneout)
      print *, '<1>_b     = ', oneout

      ! Testing time spent in box
      call init_z(z)
      sbox = linspace(0d0, 1d0, nbox)
      call time_in_box(z, sbox, taub, taubox)
      print *, 'taubox = ', taubox

      ! alpha(z)
      ! perpinv = (1.d0 - z(5)**2)*z(4)**2/bmod(z)
      ! toten = z(4)**2 + phi_elec(z)
      ! p_phi = ro0*z(4)*z(5)*hcovar(2) + psif(z)

      ! call bounce_harmonic(2, z, Hpert, 1, 1, taub, delphi, HmReIm)
      ! print *, 'H_11      = ', HmReIm

    end subroutine main

    subroutine init_z(z)
      ! Initial conditions
      real(8), intent(out) :: z(5)

      z(1) = 170.0d0  ! R
      z(2) = 1.0d0    ! PHI
      z(3) = 20d0     ! Z
      z(4) = 1d0      ! p/p0
      z(5) = 0.9d0    ! p_par/p
    end subroutine init_z


    subroutine one(t, z, res)
    ! Toroidal harmonic of Hamiltonian perturbation: real and imaginary part
      real(8), intent(in)  :: t           ! Orbit time parameter
      real(8), intent(in)  :: z(5)        ! Orbit phase-space variables
      real(8), intent(out) :: res(1)      ! Output of Hamiltonian perturbation

      res(1) = 1d0
    end subroutine one

    subroutine test_magfie()
      use do_magfie_mod, only : inp_swi, s, bfac, do_magfie_init, booz_to_cyl

      real(8) :: x(3), r(3)
      integer :: i, funit

      ! Required by do_magfie_mode to be initialized
      inp_swi = 9
      s = 1.0d0
      bfac = 1.0d0

      call do_magfie_init

      open(newunit=funit, file='test_flux_surface.out')
      write(funit, *) '% s th r z'
      x(1) = s; x(2) = 0.0d0; x(3) = 0.0d0
      do i = 0, 100
        x(3) = i*2d0*pi/100
        call booz_to_cyl(x, r)
        write(funit, *) x(1), x(3), r(1), r(3)
      end do
      close(funit)
    end subroutine

  end program test_orbit
