!
  program classify_orbits
!
  use parmot_mod,                   only : rmu,ro0
  use orbit_dim_mod,                only : write_orb,iunit1,numbasef,neqm, &
                                          orbit_clip_resonance_classes => clip_resonance_classes
  use get_matrix_mod,               only : iclass
  use global_invariants,            only : dtau,toten,perpinv,cE_ref,Phi_eff
  use bounds_fixpoints_mod,         only : nregions
  use form_classes_doublecount_mod, only : nclasses
  use cc_mod,                       only : wrbounds,dowrite
  use resint_mod,                   only : nmodes,marr,narr,delint_mode
  use poicut_mod,                   only : npc,rpc_arr,zpc_arr
  use phielec_of_psi_mod,           only : npolyphi, polyphi, polydens, polytemp
  use logging_mod,                  only : init_logging, close_logging, get_log_unit, &
                                          is_logging_enabled, tee_message
  use potato_input_mod,             only : read_potato_input, print_potato_input, &
                                          itest_type, E_alpha, A_alpha, Z_alpha, &
                                          rho_pol, rho_pol_max, scalfac_energy, &
                                          scalfac_efield, Rmax_orbit, ntimstep, &
                                          npoicut, m_min, m_max, n_tor, &
                                          toten_plot, perpinv_plot, profile_file, &
                                          edge_extension, &
                                          input_clip_resonance_classes => clip_resonance_classes, &
                                          orbit_Rstart, orbit_Zstart, orbit_lambda, &
                                          freq_Rmin, freq_Rmax, freq_n
  use field_eq_mod,                 only : allow_sol, psi_axis, psi_sep
  use field_sub,                    only : psif
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
  logical :: trace_single_orbit,freq_scan
!
  integer          :: ifdir_type,ierr,m,iunit,i,log_u
  double precision :: v0,bmod_ref
  double precision, dimension(:,:), allocatable :: resint
  double precision, dimension(neqm) :: z_single
  double precision :: taub_single,delphi_single,extraset_single(1)
  external :: find_bounce, velo, magfie
! frequency radial scan (itest_type=5):
  integer          :: ifreq
  double precision :: Rstep,omega_b,omega_phi,s_norm,rhopol
  double precision :: bmod_d,sqrtg_d,bder_d(3),hcovar_d(3),hctrvr_d(3),hcurl_d(3),x_d(3)

!
!
  call init_logging('potato.log')
  log_u = get_log_unit()
!
  call read_potato_input('potato.in')
  call print_potato_input(6)
  if (is_logging_enabled()) call print_potato_input(log_u)
!
! Allow orbits to cross the separatrix into the scrape-off layer when requested:
  allow_sol = edge_extension
  orbit_clip_resonance_classes = input_clip_resonance_classes
!
  iunit=71
!
  trace_single_orbit=.false.
  freq_scan=.false.
  select case(itest_type)
  case(1)
    call tee_message('Plot orbits')
    plot_orbits=.true.
    compute_equibrium_profiles=.false.
    compute_resonant_torque=.false.
  case(2)
    call tee_message('Equilibrium profiles')
    plot_orbits=.false.
    compute_equibrium_profiles=.true.
    compute_resonant_torque=.false.
  case(3)
    call tee_message('Resonant torque')
    plot_orbits=.false.
    compute_equibrium_profiles=.false.
    compute_resonant_torque=.true.
  case(4)
    call tee_message('Trace single orbit')
    plot_orbits=.false.
    compute_equibrium_profiles=.false.
    compute_resonant_torque=.false.
    trace_single_orbit=.true.
  case(5)
    call tee_message('Frequency radial scan')
    plot_orbits=.false.
    compute_equibrium_profiles=.false.
    compute_resonant_torque=.false.
    freq_scan=.true.
  case default
    call tee_message('unknown test')
    call close_logging()
    stop
  end select
!
!
! Apply scaling to input parameters:
  E_alpha=E_alpha/scalfac_energy
!
! inverse relativistic temeprature $m_\alpha c^2 / T_\alpha$:
  rmu=1.d30
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
  dtau=Rmax_orbit/dble(ntimstep)

! Load profiles
  call tee_message('Loading profiles from '//trim(profile_file))
  open(iunit,file=trim(profile_file),status='old',action='read')
  read(iunit,*)  ! dummy
  read(iunit,*)  ! dummy
  read(iunit,*) polydens
  read(iunit,*) polytemp  ! dummy
  read(iunit,*) polytemp
  read(iunit,*) polyphi
  close(iunit)
!
  polytemp = polytemp/scalfac_energy
  polyphi = polyphi/scalfac_efield
!  polyphi = polyphi/scalfac_banana
!
  polytemp = polytemp/E_alpha
  polyphi = polyphi*Z_alpha*e_charge/(E_alpha*ev)
!
! Find Poincare cut:
!  plot_poicut=.true.
  plot_poicut=.false.
!
  call tee_message('Computing Poincare cut')
  call find_poicut(rho_pol_max,npoicut)
  call tee_message('Poincare cut done')
!
  if(plot_poicut) then
    open(iunit,file='poicut.dat',status='replace',action='write')
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
    call tee_message( &
      'Direction of magnetic field AUG standard: ' // &
      'co-passing orbits shifted to the HFS')
  else
    call tee_message( &
      'Direction of magnetic field AUG non-standard: ' // &
      'co-passing orbits shifted to the LFS')
  endif
!
!
  call rhopol_boundary(rho_pol)
!
! Pre-compute profiles of flux surface labels, safety factor, average nabla psi, flux surface area:
!
  call tee_message('Loading equilibrium magnetic profiles')
  call load_eqmagprofs
  call tee_message('Equilibrium magnetic profiles loaded')
!
! Test the interpolation of flux surface labels, safety factor, average nabla psi, flux surface area:
!
!  call test_eqmagprofs  !WARNING: contains "stop" inside, comment this test out to continue
!
!.......................................
!
! Trace a single guiding-center orbit over one bounce period and write it.
! With edge_extension=.true. the orbit may cross the separatrix into the SOL.
!
  if(trace_single_orbit) then
    call tee_message('Tracing single orbit over one bounce')
    z_single(1)=orbit_Rstart   ! R [cm]
    z_single(2)=0.d0           ! phi
    z_single(3)=orbit_Zstart   ! Z [cm]
    z_single(4)=1.d0           ! p, monoenergetic at E_alpha
    z_single(5)=orbit_lambda   ! pitch cosine v_par/v
    write_orb=.true.
    iunit1=100
    extraset_single=0.d0
    call find_bounce(0,velo,dtau,z_single,taub_single,delphi_single, &
                     extraset_single,ierr)
    close(iunit1)
    write_orb=.false.
    write(*,'(A,ES14.6,A,ES14.6,A,I0)') ' single orbit: taub=',taub_single, &
          '  delphi=',delphi_single,'  ierr=',ierr
  endif
!
!.......................................
!
! Frequency radial scan: at each midplane start radius trace one bounce and
! record bounce frequency omega_b = 2 pi v0 / taub and toroidal precession
! frequency omega_phi = delphi v0 / taub against the flux label rho_pol.
!
  if(freq_scan) then
    call tee_message('Frequency radial scan')
    write_orb=.false.
    open(iunit,file='freq_scan.dat',status='replace',action='write')
    write(iunit,'(A)') '# R_start[cm] rho_pol omega_b[1/s] omega_phi[1/s] '// &
                       'taub delphi ierr'
    Rstep=(freq_Rmax-freq_Rmin)/dble(max(freq_n-1,1))
    do ifreq=1,freq_n
      z_single(1)=freq_Rmin+Rstep*dble(ifreq-1)
      z_single(2)=0.d0
      z_single(3)=orbit_Zstart
      z_single(4)=1.d0
      z_single(5)=orbit_lambda
! flux label of the start point: evaluate the field there, then normalize psi.
      x_d=z_single(1:3)
      call magfie(x_d,bmod_d,sqrtg_d,bder_d,hcovar_d,hctrvr_d,hcurl_d)
      s_norm=(psif-psi_axis)/(psi_sep-psi_axis)
      rhopol=sqrt(max(s_norm,0.d0))
      extraset_single=0.d0
      call find_bounce(0,velo,dtau,z_single,taub_single,delphi_single, &
                       extraset_single,ierr)
      if(ierr.eq.0 .and. taub_single.gt.0.d0) then
        omega_b=2.d0*pi*v0/taub_single
        omega_phi=delphi_single*v0/taub_single
        write(iunit,'(6ES16.7,I3)') z_single(1),rhopol,omega_b,omega_phi, &
                                    taub_single,delphi_single,ierr
      else
        write(iunit,'(6ES16.7,I3)') z_single(1),rhopol,0.d0,0.d0, &
                                    taub_single,delphi_single,ierr
      endif
    enddo
    close(iunit)
    call tee_message('Frequency radial scan done -> freq_scan.dat')
  endif
!








!.......................................
!
! Plotting the orbits, frequencies, etc
!
  if(plot_orbits) then
    call tee_message('Plotting the orbits')
!
    numbasef=6 !corresponds to unperturbed profile computation (should be > 0 for integrate_class_doublecount)
    allocate(resint(numbasef,2))
!
    toten=toten_plot
    perpinv=perpinv_plot
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
    call tee_message('Computing equilibrium profiles')
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
call test_prfs
call plot_canfreqs(v0)
!stop
    call tee_message('Computing the resonant torque')
    nmodes=m_max-m_min+1
    allocate(marr(nmodes),narr(nmodes),delint_mode(nmodes))
    narr=n_tor
    do m=m_min,m_max
      marr(m-m_min+1)=m
    enddo
!
!    dowrite=.true. !write canonical frequencies and bounce integrals for adaptive sampling grid and interpolation
!
!    call test_prfs
!    stop
    call resonant_torque
!
  endif
!
! End compute resonant torque
!
  call close_logging()
  stop
!
  end program classify_orbits
!
!
!****************************************************************************************************
!
subroutine test_prfs
 use field_eq_mod,       only : psi_sep,psi_axis
integer :: i,nsurf
double precision :: psi,dens,temp,ddens,dtemp,phi_elec,dPhi_dpsi
nsurf=1000
    open(newunit=iunit,file='profile.out',status='replace',action='write')
        do i=1,nsurf
            psi = psi_axis + i*1d0/nsurf*(psi_sep-psi_axis)
            call phielec_of_psi(psi, phi_elec, dPhi_dpsi)
            call denstemp_of_psi(psi, dens, temp, ddens, dtemp)
            write(iunit, *) psi, phi_elec, dens, temp, dPhi_dpsi, ddens, dtemp
        end do
    close(iunit)
end subroutine test_prfs
!
!****************************************************************************************************
!
  subroutine crude_outer_startpoint(psipol,Rst,bmod,phi_elec,ierr)

! Finds outer point on the Poincare cut with given poloidal flux
! Computes module of B and electrostatic potential there
!
  use field_sub, only : psif
  use poicut_mod, only : npc,rpc_arr,zpc_arr,rmagaxis
  use logging_mod, only : tee_message
!
  implicit none
!
  logical :: initprev
  integer :: ierr,i
  double precision :: psipol,Rst,bmod,phi_elec,psif_prev
  double precision, dimension(3)    :: x
  character(len=256) :: msg
!
  x(2)=0.d0
  ierr=1
  initprev=.false.
!
  do i=1,npc
    if(rpc_arr(i).lt.rmagaxis) cycle
    x(1)=rpc_arr(i)
    x(3)=zpc_arr(i)
!
    call get_bmod_and_Phi(x,bmod,phi_elec)
!
    if(initprev) then
      if((psif-psipol)*(psif_prev-psipol).lt.0.d0) then
        ierr=0
        Rst=x(1)
        write(msg, '(A,3ES14.6)') &
          'R_cut_in < Rst < R_cut_out = ', &
          rpc_arr(1), Rst, rpc_arr(npc)
        call tee_message(trim(msg))
        return
      endif
      psif_prev=psif
    else
      psif_prev=psif
      initprev=.true.
    endif
  enddo
!
  end subroutine crude_outer_startpoint
!
!--------------------------------------------------------------------
!
  subroutine plot_canfreqs(v0)
!
  use orbit_dim_mod,     only : neqm
  use global_invariants, only : dtau
  use field_eq_mod,      only : psi_axis,psi_sep
  use potato_input_mod,  only : enkin_over_temp,rho_pol
!
  implicit none
!
  integer, parameter :: next=0
  double precision, parameter :: pi=3.14159265358979d0
!
  integer :: ierr,i,nperpinv,iunit
  double precision :: psipol,Rst,bmod,phi_elec,sigma
  double precision :: enkin,perpinv_max,toten,perpinv,h_perpinv,v0
  double precision :: dens,temp,ddens,dtemp,psiast,dpsiast_dRst,taub,delphi
  double precision :: eta,omega_b,omega_tor
  double precision, dimension(neqm) :: z
  double precision, dimension(next) :: extraset
!
  external :: velo
!
  psipol=psi_axis+rho_pol**2*(psi_sep-psi_axis)
!
  nperpinv=1000
  iunit=871
!
  call crude_outer_startpoint(psipol,Rst,bmod,phi_elec,ierr)
!
  call denstemp_of_psi(psipol,dens,temp,ddens,dtemp)
!
  enkin=enkin_over_temp*temp           !<=kinetic energy normalized by reference energy
!
  toten=enkin+phi_elec
  perpinv_max=enkin/bmod
  h_perpinv=perpinv_max/dble(nperpinv+1)
!
  open(iunit,file='canonical_freqs_vs_eta_posvpar.dat',status='replace',action='write')
  write(iunit,*) '# eta [1/G], omega_b [rad/s], Omega_tor [rad/s] for starting v_par > 0'
  sigma=1.d0
!
  do i=1,nperpinv
    perpinv=h_perpinv*dble(i)
!
    call starter_doublecount(toten,perpinv,sigma,Rst,   &
                             psiast,dpsiast_dRst,z,ierr)
!
    call find_bounce(next,velo,dtau,z,taub,delphi,extraset,ierr)
    if(ierr.ne.0) cycle
!
    eta=perpinv/enkin
    omega_b=2.d0*pi*v0/taub
    omega_tor=delphi*v0/taub
    write(iunit,*) eta,omega_b,omega_tor
  enddo
!
  close(iunit)
!
  open(iunit,file='canonical_freqs_vs_eta_negvpar.dat',status='replace',action='write')
  write(iunit,*) '# eta [1/G], omega_b [rad/s], Omega_tor [rad/s] for starting v_par < 0'
  sigma=-1.d0
!
  do i=1,nperpinv
    perpinv=h_perpinv*dble(i)
!
    call starter_doublecount(toten,perpinv,sigma,Rst,   &
                             psiast,dpsiast_dRst,z,ierr)
!
    call find_bounce(next,velo,dtau,z,taub,delphi,extraset,ierr)
    if(ierr.ne.0) cycle
!
    eta=perpinv/enkin
    omega_b=2.d0*pi*v0/taub
    omega_tor=delphi*v0/taub
    write(iunit,*) eta,omega_b,omega_tor
  enddo
!
  close(iunit)
!
  end subroutine plot_canfreqs
