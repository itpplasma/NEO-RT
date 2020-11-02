program neo_rt
  use common
  use orbit, only: timestep, bounce_harmonic
  use parmot_mod, only: rmu, ro0

  implicit none

  real(8) :: z(neqm), zdot(neqm)
  real(8) :: bmod_ref, bmod00, tempi1, am, Zb, v0, rlarm
  real(8) :: bmod, phi_elec
  real(8) :: mth, mph  ! Poloidal and toroidal harmonic numbers
  real(8) :: alpha(3)  ! Invariants
  real(8) :: HmReIm(2)

  ! Inverse relativistic temperature, set >=1e5 to neglect relativistic effects
  rmu = 1.0d5

  ! This will translate to dimensionless units via parmot_mod.ro0
  ! Here we translate an input file in Tesla to computation in Gauss/CGS
  bmod_ref = 1.0d0             ! 1 Gauss in Gauss
  bmod00 = 1.0d0               ! On-axis magnetic field
  tempi1 = 0.17549561306d4     ! ion temperature
  am = 2.0d0                   ! Atomic mass 2 of deuterium ions
  Zb = 1.0d0                   ! Atomic charge 1 of deuterium ions

  v0 = sqrt(2.0*tempi1*ev/(am*mu))     ! Reference (thermal) velocity
  ! Reference Larmor radius of thermal particles
  rlarm = v0*am*mu*c/(Zb*qe*bmod_ref)  ! Larmor radius in bmod_ref
  ro0 = rlarm*bmod00                   ! Rescaled Larmor radius

  ! Initial conditions
  z(1) = 170.0d0  ! R
  z(2) = 0.0d0    ! PHI
  z(3) = 40d0     ! Z
  z(4) = 1d0      ! v/v0
  z(5) = 0.9d0    ! v_par/v

  ! First call for field
  call get_bmod_and_Phi(z(1:3), bmod, phi_elec)
  print *, 'bmod, phielec:'
  print *, bmod, phi_elec

  print *, z
  call timestep(0.0d0, z, zdot)
  print *, zdot

  ! alpha(z)
  ! perpinv = (1.d0 - z(5)**2)*z(4)**2/bmod(z)
  ! toten = z(4)**2 + phi_elec(z)
  ! p_phi = ro0*z(4)*z(5)*hcovar(2) + psif(z)

  call bounce_harmonic(2, z, 1, 1, HmReIm, Hn)

  print *, 'HmReIm: ', HmReIm


  contains

  subroutine Hn(z2, out)

    real(8), intent(in)  :: z2(5)
    real(8), intent(out) :: out(2)

    out(1) = z2(1)
    out(2) = z2(3)
  end subroutine Hn
end program neo_rt
