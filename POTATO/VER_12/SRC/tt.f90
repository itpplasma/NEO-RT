!
  program classify_orbits
!
  use parmot_mod,                   only : rmu,ro0
  use orbit_dim_mod,                only : write_orb,iunit1,numbasef
  use get_matrix_mod,               only : iclass
  use global_invariants,            only : dtau,toten,perpinv,cE_ref,Phi_eff
  use bounds_fixpoints_mod,         only : nregions
  use form_classes_doublecount_mod, only : nclasses
  use cc_mod,                       only : wrbounds,dowrite
  use resint_mod,                   only : nmodes,marr,narr,delint_mode
  use poicut_mod,                   only : npc,rpc_arr,zpc_arr
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
  logical :: classes_talk,plot_orbits,compute_equibrium_profiles,compute_resonant_torque,plot_poicut
!
  integer          :: npoicut,ifdir_type,ierr,ntimstep,itest_type,m,iunit,i
  double precision :: v0,bmod_ref,E_alpha,A_alpha,Z_alpha
  double precision :: rmax,rho_pol,rho_pol_max
  double precision, dimension(:,:), allocatable :: resint
!
  iunit=71
!
! Uncomment to plot B-mod perturbation Fourier amplitude Bmod_n, stops execution after writing data:
!
!  call test_bmodpert
!
! Type of test: 1 - plot orbits, 2 - equilibrium profiles, 3 - resonant torque
  itest_type=2
!
  select case(itest_type)
  case(1)
    plot_orbits=.true.
    compute_equibrium_profiles=.false.
    compute_resonant_torque=.false.
  case(2)
    plot_orbits=.false.
    compute_equibrium_profiles=.true.
    compute_resonant_torque=.false.
  case(3)
    plot_orbits=.false.
    compute_equibrium_profiles=.false.
    compute_resonant_torque=.true.
  case default
    print *,'unknown test'
    stop
  end select
!
!
! inverse relativistic temeprature $m_\alpha c^2 / T_\alpha$:
  rmu=1.d30
! Normalization energy of alpha-species, eV:
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
! Reference energy $\cE_{ref}$:
  cE_ref=E_alpha*ev
! Normalization "potential" Eq.(83):
  Phi_eff=c*E_alpha*ev/(e_charge*Z_alpha*v0)
!
! Set orbit integration step:
  Rmax=200.d0
  ntimstep=30
  dtau=Rmax/dble(ntimstep)
!
! Find Poincare cut:
  npoicut=10000     !number of equidistant points for interpolating the cut
  rho_pol_max=0.8   !maximum value of poloidal radius limiting the cut range
!  plot_poicut=.true.
  plot_poicut=.false.
!
  call find_poicut(rho_pol_max,npoicut)
!
  if(plot_poicut) then
    open(iunit,file='poicut.dat')
    do i=0,npc
      write(iunit,*) rpc_arr(i),zpc_arr(i)
    enddo
    close(iunit)
    stop
  endif
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
!  rho_pol=0.8d0 !sqrt(0.3d0) !poloidal radius
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
    toten=1.d0 !1.7521168986242395d0
    perpinv=4.5d-5 !9.9881139234315119d-5
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
! Compute resonant torque
!
  if(compute_resonant_torque) then
    print *,'Computing the resonant torque'
    nmodes=11
    allocate(marr(nmodes),narr(nmodes),delint_mode(nmodes))
    narr=1           !  toroidal harmonic index  m_3
    nmodes=0
    do m=-5,5
      nmodes=nmodes+1
      marr(nmodes)=m ! "poloidal" harmonic index m_2
    enddo
!
!    dowrite=.true. !write canonical frequencies and bounce integrals for adaptive sampling grid and interpolation
!
    call resonant_torque
!
  endif
!
! End compute resonant torque
!
  stop
!
  end program classify_orbits
!
