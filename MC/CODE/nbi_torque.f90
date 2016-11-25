!
  use parmot_mod, only : rmu,ro0
  use collis_alp, only : swcoll,iswmod,ns,efcolf,velrat,enrat,efcolf_arr,velrat_arr,enrat_arr
!  use odeint_mod, only : adaptive
  USE polylag_3,  ONLY : mp,indef, plag1d
  use neo_input,  only : flux, pertscale, bscale
  USE probstart_mod, ONLY : calc_probstart
  use elefie_mod, only: Mtprofile, rbig, plasma, amb,am1,am2,Zb,Z1,Z2,densi1,densi2,tempi1,&
       tempi2,tempe,v0,escale
  use constants, only: pi,c,e_charge,e_mass,p_mass,ev
  use spline_settings, only: flux_surf_dist
!
  implicit none
!
  integer          :: ierr,npoiper,i,ntestpart
  integer          :: ipart,icpu,iskip,istep
  real             :: zzg,gauss
  double precision :: dphi,bmod00,rlarm
  double precision :: dtau,xi,bmod_ref,E_beam
  double precision, dimension(5) :: z
  integer,          dimension(:),   allocatable :: ibinsrc_x
  double precision, dimension(:),   allocatable :: s_ionized
  double precision, dimension(:,:), allocatable :: zstart     
  double precision, dimension(:), allocatable :: savg, ds2avg
  integer :: iout
!
  double precision :: bmod,sqrtg
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
  integer :: nplasma
  integer,          dimension(mp) :: indu
  double precision, dimension(mp) :: xp,fp
  double precision :: s0, s, der, dxm1
  integer :: ks0 ! index of starting flux surface
  
!  double precision :: x_nrl
!
!
!  safety:
  integer, parameter :: nstepmax=100000000 ! TODO: set to original: 100000000
!
!  thermalization:
  double precision, parameter :: vmin_therm=sqrt(1.5d0)
!
  double precision :: dummy, taumax, runtime

  character(len=20) :: buffer
  
  ! inverse relativistic temperature
  rmu=1d5
  
! process command line arguments

  call get_command_argument(1, buffer)
  read (buffer, '(I3)') icpu
  print *, "CPU: ", icpu
  
!
  open(1,file='mc_torque.inp')
  read (1,*) npoiper           !number of points per period 
  read (1,*) ntestpart         !number of test particles
  read (1,*,end=1) swcoll      !collision switch: .true. - collisions on
  read (1,*) iswmod        !switch of the collision mode:
!                              1 - full operator (pitch-angle and energy scattering and drag)
!                              2 - energy scattering and drag only
!                              3 - drag only
!                              4 - pitch-angle scattering only
 read(1,*) pertscale
 read(1,*) bscale
 read(1,*) escale
 read(1,*) s0
 read(1,*) runtime
 read(1,*) flux_surf_dist
1 continue
  close(1)
!
!

  
  open(1,file='plasma.dat')
  read (1,*)
  read (1,*) nplasma,am1,am2,Z1,Z2
  read (1,*)
  allocate(plasma(nplasma,6))
  do i=1,nplasma
    read (1,*) plasma(i,:)
  enddo
  dxm1=1.d0/(plasma(2,1)-plasma(1,1))
  close(1)

  amb=am1
  Zb=Z1
  
  allocate(Mtprofile(nplasma,3))
  open(1,file='Mtprofile.dat')
  do i=1,nplasma
     read (1,*) Mtprofile(i,:)
  enddo
  close(1)
!
  call indef(s0,plasma(1,1),dxm1,nplasma,indu)

  xp=plasma(indu,1)
  fp=plasma(indu,4)

  call plag1d(s0,fp,dxm1,xp,tempi1,der)

  v0 = sqrt(2.d0*tempi1*ev/(am1*p_mass))
  E_beam = amb*p_mass*v0**2/(2d0*ev)
print *, v0, tempi1, ev, amb, p_mass
  
  if(swcoll) then
!
    do i=1,ns
      s=dfloat(i-1)/dfloat(ns-1)
!
      call indef(s,plasma(1,1),dxm1,nplasma,indu)
!
      xp=plasma(indu,1)
      fp=plasma(indu,2)
!
      call plag1d(s,fp,dxm1,xp,densi1,der)
!
      fp=plasma(indu,3)
!
      call plag1d(s,fp,dxm1,xp,densi2,der)
!
      fp=plasma(indu,4)
!
      call plag1d(s,fp,dxm1,xp,tempi1,der)
!
      fp=plasma(indu,5)
!
      call plag1d(s,fp,dxm1,xp,tempi2,der)
!
      fp=plasma(indu,6)
!
      call plag1d(s,fp,dxm1,xp,tempe,der)
      if (abs(s0-s)<0.5d0/dfloat(ns-1)) then
         ks0 = i
         print *, v0, tempi1, ev, amb, p_mass
      end if
!if (i==ns/2) print *, amb,am1,am2,Zb,Z1,Z2
!if (i==ns/2) print *, densi1,densi2,tempi1,tempi2,tempe,E_beam
      call loacol_nbi(amb,am1,am2,Zb,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,E_beam,v0)
!if (i==ns/2) print *, v0
!if (i==ns/2) x_nrl=densi1/tempi1**1.5d0
!   
      efcolf_arr(:,i)=efcolf
      velrat_arr(:,i)=velrat
      enrat_arr(:,i)=enrat
   enddo
   !
endif
!
s = s0

!do i=1,300
!p = 1.d-2*i
!call coleff(s,p,dpp,dhh,fpeff)
!print *, v0*dpp, v0*dhh
!if(p.lt.1.d0) then
!write (1,*) p,dhh*v0,1.4d-7*x_nrl/sqrt(2.d0)*18.d0/p**2 ! 18.d0 - Coulomb logarithm
!else
!write (1,*) p,dhh*v0,1.8d-7*x_nrl/sqrt(2.d0)*18.d0/p**3 ! 18.d0 - Coulomb logarithm
!endif
!enddo
!
  x=0.d0
  x(1)=1.d-8
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
!  major radius:
  rbig=hcovar(2)
!  reference field - one Tesla:
  bmod_ref=1.d4
  bmod00=1.d0
print *, v0, tempi1, ev, amb, p_mass
  v0=sqrt(2.d0*E_beam*ev/(amb*p_mass))
print *, v0, tempi1, ev, amb, p_mass
!
!  Larmor radius:
  print *, v0,amb,Zb,bmod_ref,bmod00
  rlarm=v0*amb*p_mass*c/(Zb*e_charge*bmod_ref)
  ro0=rlarm*bmod00
  print *, rlarm
!

!  orbit integration time step:
  dphi=2.d0*pi/npoiper
  dtau=dphi*rbig
!
!
  
  if (icpu==0) then
     open(1,file='nbi_torque.log', recl=1024)
     write (1,*) 'npoiper = ',npoiper
     write (1,*) 'ntestpart = ',ntestpart
     write (1,*) 'dphi = ',dphi
     write (1,*) 'v0 = ',v0
     write (1,*) 'rlarm = ',rlarm
     write (1,*) 'dtau = ',dtau
     write (1,*) 'E_beam = ',E_beam
     write (1,*) 'swcoll = ',swcoll
     write (1,*) 'rbig = ',rbig
     write (1,*) 'flux = ',flux
     write (1,*) 'psi_pr = ',1.0d8*flux/(2*pi)
     if(swcoll) then
        write (1,*) 'iswmod = ',iswmod
        write (1,*) 'am1 = ',am1
        write (1,*) 'am2 = ',am2
        write (1,*) 'Z1 = ',Z1
        write (1,*) 'Z2 = ',Z2
        write (1,*) 'tau_c = ', 1.d0/efcolf_arr(1,ns/2)
        write (1,*) 'tau_max = ', 5.d0/efcolf_arr(1,ns/2)
        write (1,*) 'factor = ', v0*efcolf_arr(1,ns/2)*0.1
     endif
     close(1)
  endif
!
  do iskip=1,icpu+1
    do ipart=1,ntestpart
      xi=zzg()
      xi=zzg()
      xi=gauss()
      xi=gauss()
      xi=gauss()
      xi=zzg()
    enddo
  enddo
!
  allocate(zstart(5,ntestpart),ibinsrc_x(ntestpart),s_ionized(ntestpart))

!
  do ipart=1,ntestpart
     ! determine the starting point:
     !
     !call binsrc(probstart,1,nbeam,xi,i)
     !
     !ibinsrc_x(ipart)=i
     ! coordinates: z(1) = s, z(2) = phi, z(3) = theta
     zstart(1,ipart)=s0
     xi=zzg()
     zstart(2,ipart)=xi*2.d0*pi
     xi=zzg()
     zstart(3,ipart)=xi*2.d0*pi
     ! normalized velocity module z(4) = v / v_0:
!     zstart(4,ipart)=1.d0
     zstart(4,ipart)=sqrt((gauss()**2+gauss()**2+gauss()**2)/3.d0)
     !
     x=zstart(1:3,ipart)
     !
     call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
     !
     !  pitch:
     ! starting pitch z(5)=v_\parallel / v:
     xi=zzg()
     zstart(5,ipart) = 2.0d0*xi-1.0d0
     print *, ipart, zstart(:,ipart)
  enddo
!
!

!  open(2000+icpu,recl=1024)

  ! TODO: con't hardcode 200 steps
  allocate(savg(200),ds2avg(200))

  savg = 0.0
  ds2avg = 0.0
  
  do ipart=1,ntestpart
     print *,ipart,' / ',ntestpart
!
    z=zstart(:,ipart)
! TODO: disable next lines
!    z(1) = s0
!    z(2) = 0.0
!    z(3) = 0.0
!    z(4) = 1.0
!    z(5) = 0.1
!
    taumax = 1.0d0
    if(swcoll) taumax = runtime/efcolf_arr(1,ks0)
    
    dummy = 0.d0
    iout = 1
    do istep=1,nstepmax
       if(swcoll) then
          if (istep*dtau > taumax) exit

          if (istep*dtau > dummy) then
             !write(2000+icpu,*) ipart, istep, istep*dtau/efcolf_arr(1,ks0), z(1),&
             !     (z(1)-s0)**2, z(2), z(3)

             savg(iout) = savg(iout) + z(1)/ntestpart
             ds2avg(iout) = ds2avg(iout) + (z(1)-s0)**2/ntestpart
             iout = iout + 1             
             dummy = dummy + .005d0*taumax
          endif
       endif     

       if (z(1)<0.d0 .OR. z(1)>1.d0) then
          print *, "particle lost:"
          print *, ipart, z
          exit
       endif

       call regst(z,dtau,ierr)
       !
       if(swcoll) then
          call stost(z,dtau,iswmod,ierr)
       endif
       !if(mod(istep,nstepmax/40000)==0) then
       !   write(3000,*) istep*dtau, z
       !endif
       
   enddo
   print *,istep,'  steps'
!
  
enddo
!close(2000+icpu)

open(4000+icpu,recl=1024)
do iout=1,200
   write(4000+icpu,*) savg(iout), ds2avg(iout)
enddo
close(4000+icpu)

! 
!
!
end program
