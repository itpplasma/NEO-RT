!
  program classify_orbits
!
  use parmot_mod,                   only : rmu,ro0
  use orbit_dim_mod,                only : write_orb,iunit1,numbasef
  use get_matrix_mod,               only : iclass
  use global_invariants,            only : dtau,toten,perpinv,Phi_eff
  use bounds_fixpoints_mod,         only : nregions
  use form_classes_doublecount_mod, only : nclasses
  use cc_mod,                       only : wrbounds,dowrite
  use resint_mod,                   only : nmodes,marr,narr,delint_mode
!
  implicit none
!
  double precision, parameter :: pi=3.14159265358979d0
  double precision,parameter  :: c=2.9979d10
  double precision,parameter  :: e_charge=4.8032d-10
  double precision,parameter  :: e_mass=9.1094d-28
  double precision,parameter  :: p_mass=1.6726d-24
  double precision,parameter  :: ev=1.6022d-12
!
  logical :: classes_talk,plot_orbits,compute_equibrium_profiles,compute_resonant_flux
!
  integer          :: npoicut,ifdir_type,ierr,ntimstep,itest_type,m
  double precision :: v0,bmod_ref,E_alpha,A_alpha,Z_alpha
  double precision :: rmax,rho_pol
  double precision, dimension(:,:), allocatable :: resint
!
!
  itest_type=3
!
  select case(itest_type)
  case(1)
    plot_orbits=.true.
    compute_equibrium_profiles=.false.
    compute_resonant_flux=.false.
  case(2)
    plot_orbits=.false.
    compute_equibrium_profiles=.true.
    compute_resonant_flux=.false.
  case(3)
    plot_orbits=.false.
    compute_equibrium_profiles=.false.
    compute_resonant_flux=.true.
  case default
    print *,'unknown test'
    stop
  end select
!
!
! inverse relativistic temeprature $m_\alpha c^2 / T_\alpha$:
  rmu=1.d30
! Normalization energy of alpha-speciesi, eV:
  E_alpha=5d3
! Mass number of alpha species:
  A_alpha=2.d0   !deuterium
! Charge number of alpha species:
  Z_alpha=1.d0   !deuterium
!
! Normalization velocity, cm/s:
  v0=sqrt(2.d0*E_alpha*ev/(p_mass*A_alpha))
! Larmor radius:
  bmod_ref=1.d0 !no normalization of magnetic field
  ro0=v0*p_mass*A_alpha*c/(e_charge*Z_alpha*bmod_ref)
! Normalization "potential" Eq.(83):
  Phi_eff=c*E_alpha*ev/(e_charge*Z_alpha*v0)
!
! Set orbit integration step:
  Rmax=200.d0
  ntimstep=30
  dtau=Rmax/dble(ntimstep)
!
! Find Poincare cut:
  npoicut=10000 !number of equidistant points for interpolating the cut
!
  call find_poicut(npoicut)
!
! Determine mutual direction of poloidal and toroidal fields (not used, for curiosity only):
!
  call poltor_field_dir(ifdir_type)
!
  if(ifdir_type.eq.1) then
    print *,'Direction of magnetic field AUG standard: co-passing orbits shifted to the HFS'
  else
    print *,'Direction of magnetic field AUG non-standard: co-passing orbits shifted to the LFS'
  endif
!
! Set outer boundary of the plasma volume where profiles are computed:
  rho_pol=sqrt(0.3d0) !poloidal radius
!
  call rhopol_boundary(rho_pol)
!
! Pre-compute profiles of flux surface labels, safety factor, average nabla psi, flux surface area:
!
  call load_eqmagprofs
!
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
!
    numbasef=6 !corresponds to unperturbed profile computation (should be > 0 for integrate_class_doublecount)
    allocate(resint(numbasef,2))
!
! example from Fig.1 (electric field ampl=1.12d0):
    toten=1.7521168986242395d0
    perpinv=9.9881139234315119d-5
!
! activate writing:
    wrbounds=.true. !write vpar^2 and psi^* curves vs cut parameter R_c, extremum and boundary points
    dowrite=.true. !write canonical frequencies and bounce integrals for adaptive sampling grid and interpolation
    write_orb=.true. !write orbits during adaptive sampling of classes
!
    call find_bounds_fixpoints(ierr)
!
    classes_talk=.true.
!
    call form_classes_doublecount(classes_talk,ierr)
!
    do iclass=1,nclasses
! data is written to fort.* files:
      iunit1=100+iclass
!
      call sample_class_doublecount(1000+iclass,ierr)
!
      close(iunit1)      !orbits
      close(1000+iclass) !sampled canonical frequencies and bounce integrals

      call integrate_class_doublecount(2000+iclass,resint)
!
      close(2000+iclass) !interpolated canonical frequencies and bounce integrals
      close(3000+iclass) !derivatives of interpolated canonical frequencies and bounce integrals
    enddo
!
! deactivate writing:
    wrbounds=.false.
    dowrite=.false.
    write_orb=.false.
  endif
!
! End plotting the orbits, frequencies, etc
!
!.......................................
!
! Compute equilibrium profiles of density and toroidal guiding center flux density
!
  if(compute_equibrium_profiles) then
    print *,'Computing equilibrium profiles'
!
!compute analytical approximation for density and toroidal flow:
!
    call plotrot(Phi_eff,rho_pol)
!
!    wrbounds=.true. !<=uncomment for debugging
!    dowrite=.true.  !<=uncomment for debugging
!    iunit1=101      !<=uncomment for debugging
!
    call fit_profiles
!
  endif
!
! End compute equilibrium profiles of density and toroidal guiding center flux density
!
!.......................................
!
! Compute resonant flux
!
  if(compute_resonant_flux) then
    print *,'Computing the resonant flux'
    nmodes=11
    allocate(marr(nmodes),narr(nmodes),delint_mode(nmodes))
    narr=2           !  toroidal harmonic index  m_3
    nmodes=0
    do m=-5,5
      nmodes=nmodes+1
      marr(nmodes)=m ! "poloidal" harmonic index m_2
    enddo
!
!    dowrite=.true. !write canonical frequencies and bounce integrals for adaptive sampling grid and interpolation
!
    call resonant_flux
!
  endif
!
! End compute resonant flux
!
  stop
!
  end program classify_orbits
!
