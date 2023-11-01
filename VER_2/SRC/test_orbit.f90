program test_orbit
  use util
  use parmot_mod, only: rmu, ro0
  use orbit_dim_mod, only : write_orb,iunit1,numbasef
  use get_matrix_mod, only : iclass
  use orbit, only: timestep, neqm, bounce_average, bounce_harmonic, &
    time_in_box, orbit_init
  use transport, only: Hpert
  use poicut_mod, only : npc,rpc_arr,zpc_arr
  use form_classes_doublecount_mod, only : nclasses
  use cc_mod, only : wrbounds,dowrite
  use global_invariants, only : dtau,toten,perpinv,cE_ref,Phi_eff

  implicit none

    call main

    contains

    subroutine main

      real(8) :: z(neqm), zdot(neqm)
      real(8) :: bmod_ref, bmod00, E_alpha, A_alpha, Z_alpha, v0, rlarm
      real(8) :: bmod, phi_elec
      real(8) :: mth, mph  ! Poloidal and toroidal harmonic numbers
      real(8) :: alpha(3)  ! Invariants
      real(8) :: HmReIm(2), taub, delphi, oneout(1)

      real(8) :: Rmax, ntimstep

      integer(4), parameter :: nbox = 101
      real(8) :: sbox(nbox)
      real(8) :: taubox(nbox)
      real(8) :: HmReImbox(2, nbox)

      integer :: npoicut
      real(8) :: rho_pol_max

      integer :: ifdir_type
      double precision :: rho_pol

      integer :: i, ierr
      logical :: plot_poicut = .False.
      logical :: plot_orbits = .True.
      logical :: classes_talk

      integer :: iunit = 71


      double precision, dimension(:,:), allocatable :: resint



      ! Inverse relativistic temperature, set >=1e5 to neglect relativ. effects
      rmu = 1.0d30

      E_alpha=5d3                  ! Normalization energy of alpha-species, eV
      A_alpha=2.d0                 ! Mass number of alpha species: deuterium

      Z_alpha=1.d0                 ! Charge number of alpha species: deuterium
      v0=sqrt(2.d0*E_alpha*eV/(mu*A_alpha)) ! Normalization velocity, cm/s

      ! This will translate to dimensionless units via parmot_mod.ro0
      ! Here we translate an input file in Tesla to computation in Gauss/CGS
      bmod_ref = 1.0d0             ! 1 Gauss in Gauss

      ! Reference Larmor radius of thermal particles
      ro0=v0*mu*A_alpha*c/(qe*Z_alpha*bmod_ref)

      ! Set orbit integration step:
      Rmax=200.d0
      ntimstep=30
      dtau=Rmax/dble(ntimstep)

      call orbit_init(dtau)

      ! Find Poincare cut:
      npoicut=10000     !number of equidistant points for interpolating the cut
      rho_pol_max=0.8d0 !maximum value of poloidal radius limiting the cut range

      call find_poicut(rho_pol_max,npoicut)

      if(plot_poicut) then
        open(iunit,file='poicut.dat')
        do i=0,npc
          write(iunit,*) rpc_arr(i),zpc_arr(i)
        enddo
        close(iunit)
        stop
      endif

      ! Determine mutual direction of poloidal and toroidal fields (not used, for curiosity only):

      call poltor_field_dir(ifdir_type)

      if(ifdir_type.eq.1) then
        print *,'Direction of magnetic field AUG standard: co-passing orbits shifted to the HFS'
      else
        print *,'Direction of magnetic field AUG non-standard: co-passing orbits shifted to the LFS'
      endif

      ! Set outer boundary of the plasma volume where profiles are computed:
        rho_pol=sqrt(0.3d0) !poloidal radius
      !  rho_pol=0.8d0 !sqrt(0.3d0) !poloidal radius

      call rhopol_boundary(rho_pol)

      ! Pre-compute profiles of flux surface labels, safety factor, average nabla psi, flux surface area:

        call load_eqmagprofs

      ! Test the interpolation of flux surface labels, safety factor, average nabla psi, flux surface area:
      !
      !  call test_eqmagprofs  !WARNING: contains "stop" inside, comment this test out to continue
      !

      !.......................................
      !
      ! Plotting the orbits, frequencies, etc
      !
      if(plot_orbits) then
        print *,'Plotting the orbits'

        numbasef=6 !corresponds to unperturbed profile computation (should be > 0 for integrate_class_doublecount)
        allocate(resint(numbasef,2))

        ! example from Fig.1 (electric field ampl=1.12d0):
        toten=1.d0 !1.7521168986242395d0
        perpinv=4.5d-5 !9.9881139234315119d-5

        ! activate writing:
        wrbounds=.true. !write vpar^2 and psi^* curves vs cut parameter R_c, extremum and boundary points
        dowrite=.true. !write canonical frequencies and bounce integrals for adaptive sampling grid and interpolation
        write_orb=.true. !write orbits during adaptive sampling of classes
            call find_bounds_fixpoints(ierr)

        classes_talk=.true.

        call form_classes_doublecount(classes_talk,ierr)

        do iclass=1,nclasses
          ! data is written to fort.* files:
          iunit1=100+iclass

          call sample_class_doublecount(1000+iclass,ierr)

          close(iunit1)      !orbits
          close(1000+iclass) !sampled canonical frequencies and bounce integrals

          call integrate_class_doublecount(2000+iclass,resint)

          close(2000+iclass) !interpolated canonical frequencies and bounce integrals
          close(3000+iclass) !derivatives of interpolated canonical frequencies and bounce integrals
        enddo

        ! deactivate writing:
        wrbounds=.false.
        dowrite=.false.
        write_orb=.false.
      endif

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
      call init_z(z)
      call bounce_average(1, z, one, taub, delphi, oneout)
      print *, '<1>_b     = ', oneout
      call bounce_average(1, z, one, taub, delphi, oneout)
      print *, '<1>_b     = ', oneout

      ! Compute bounce harmonic of a_n = 1
      !call bounce_harmonic(z, cmplxone, 1, 1, taub, delphi, HmReIm)
      !print *, '1_11      = ', HmReIm

      ! Compute bounce harmonic of H - should give some value
      call init_z(z)
      call bounce_harmonic(z, Hpert, 1, 1, taub, delphi, HmReIm)
      print *, 'H_11      = ', HmReIm

      ! Testing time spent in box
      sbox = linspace(0d0, 1d0, nbox)
      call time_in_box(z, sbox, taub, taubox)
      write(99,*) taubox

      ! alpha(z)
      ! perpinv = (1.d0 - z(5)**2)*z(4)**2/bmod(z)
      ! toten = z(4)**2 + phi_elec(z)
      ! p_phi = ro0*z(4)*z(5)*hcovar(2) + psif(z)


    end subroutine main

    subroutine init_z(z)
      ! Initial conditions
      real(8), intent(out) :: z(5)

      z(1) = 130.0d0  ! R
      z(2) = 1.0d0    ! PHI
      z(3) = 20d0     ! Z
      z(4) = 1d0      ! p/p0
      z(5) = 0.8d0    ! p_par/p
    end subroutine init_z


    subroutine one(t, z, res)
    ! Toroidal harmonic of Hamiltonian perturbation: real and imaginary part
      real(8), intent(in)  :: t           ! Orbit time parameter
      real(8), intent(in)  :: z(5)        ! Orbit phase-space variables
      real(8), intent(out) :: res(1)      ! Output of Hamiltonian perturbation

      res(1) = 1d0
    end subroutine one

    subroutine cmplxone(t, z, res)
    ! Toroidal harmonic of Hamiltonian perturbation: real and imaginary part
      real(8), intent(in)  :: t           ! Orbit time parameter
      real(8), intent(in)  :: z(5)        ! Orbit phase-space variables
      real(8), intent(out) :: res(2)      ! Output of Hamiltonian perturbation

      res(1) = 1d0
      res(2) = 0d0
    end subroutine cmplxone

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
