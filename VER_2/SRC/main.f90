program neo_rt

use util
use orbit, only: timestep, bounce_average, bounce_harmonic, neqm
use parmot_mod, only: rmu, ro0

implicit none

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

end program neo_rt
