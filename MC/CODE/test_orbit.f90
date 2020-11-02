!
  use parmot_mod, only : rmu,ro0
  use collis_alp, only : swcoll,iswmod,ns,efcolf,velrat,enrat,efcolf_arr,velrat_arr,enrat_arr
!  use odeint_mod, only : adaptive
  USE polylag_3,  ONLY : mp,indef, plag1d
  use neo_input,  only : flux
  USE probstart_mod, ONLY : calc_probstart
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
  integer          :: ierr,npoiper,i,ntestpart
  integer          :: ipart,icpu,iskip,istep
  real             :: zzg,gauss
  double precision :: dphi,bmod00,rlarm
  double precision :: dtau,xi,v0,bmod_ref,E_beam
  double precision, dimension(5) :: z
  integer,          dimension(:),   allocatable :: ibinsrc_x
  double precision, dimension(:),   allocatable :: s_ionized
  double precision, dimension(:,:), allocatable :: zstart
  double precision :: amb,am1,am2,Zb,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe             
!
  double precision :: bmod,sqrtg,rbig
  double precision, dimension(3) :: x,bder,hcovar,hctrvr,hcurl
  integer :: nplasma
  double precision, dimension(:,:), allocatable :: plasma
  integer,          dimension(mp) :: indu
  double precision, dimension(mp) :: xp,fp
  double precision :: s, der, dxm1
!
  double precision :: x_nrl
!
!
!  safety:
  integer, parameter :: nstepmax=50000 ! original: 1e8
!
!  thermalization:
  double precision, parameter :: vmin_therm=sqrt(1.5d0)
!
  double precision :: dummy, debugstep
  double precision :: p,dpp,dhh,fpeff

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
  if(swcoll) then
    read (1,*) iswmod        !switch of the collision mode:
!                              1 - full operator (pitch-angle and energy scattering and drag)
!                              2 - energy scattering and drag only
!                              3 - drag only
!                              4 - pitch-angle scattering only
  endif
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
!
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

      if (i==ns/2) then
         v0 = sqrt(2.d0*tempi1*ev/(am1*p_mass))
         amb   = am1
         Zb    = Z1
         E_beam = amb*p_mass*v0**2/(2d0*ev)
      endif

   enddo
   
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
      if (i==ns/2) print *, v0, tempi1, ev, amb, p_mass
      if (i==ns/2) print *, amb,am1,am2,Zb,Z1,Z2
      if (i==ns/2) print *, densi1,densi2,tempi1,tempi2,tempe,E_beam
      call loacol_nbi(amb,am1,am2,Zb,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,E_beam,v0)
      if (i==ns/2) print *, v0
      if (i==ns/2) x_nrl=densi1/tempi1**1.5d0

      !   
      efcolf_arr(:,i)=efcolf
      velrat_arr(:,i)=velrat
      enrat_arr(:,i)=enrat
   enddo
   !
endif
!
s = .5
x = 0.d0
x(1) = 1.d-8
!
  call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
!  major radius:
  rbig=hcovar(2)
!  reference field - one Tesla:
  bmod_ref=1.d4
  bmod00=1.d0
!
!
  v0=sqrt(2.d0*E_beam*ev/(amb*p_mass))
!
!  Larmor radius:
  print *, v0,amb,Zb,bmod_ref,bmod00
  rlarm=v0*amb*p_mass*c/(Zb*e_charge*bmod_ref)
  ro0=rlarm*bmod00
  print *, rlarm

!  orbit integration time step:
  dphi=2.d0*pi/npoiper
  dtau=dphi*rbig
!
!
  if (icpu==0) then
     open(1,file='nbi_torque.log')
     write (1,*) 'npoiper = ',npoiper
     write (1,*) 'ntestpart = ',ntestpart
     write (1,*) 'dphi = ',dphi
     write (1,*) 'v0 = ',v0
     write (1,*) 'rlarm = ',rlarm
     write (1,*) 'dtau = ',dtau
     write (1,*) 'E_beam = ',E_beam
     write (1,*) 'swcoll = ',swcoll
     if(swcoll) then
        write (1,*) 'iswmod = ',iswmod
        write (1,*) 'am1 = ',am1
        write (1,*) 'am2 = ',am2
        write (1,*) 'Z1 = ',Z1
        write (1,*) 'Z2 = ',Z2
        write (1,*) 'tau_max = ',5.d0/efcolf_arr(1,ns/2)
        write (1,*) 'factor = ',v0*efcolf_arr(1,ns/2)*0.1
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
     zstart(1,ipart)=0.5
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
     zstart(5,ipart)=xi
     print *, ipart, zstart(:,ipart)
  enddo
!
!
  do ipart=1,ntestpart
     print *,ipart,' / ',ntestpart
!
    z=zstart(:,ipart)
!
    debugstep = .01d0*5.d0/efcolf_arr(1,ns/2)
    dummy = 0.d0
    do istep=1,nstepmax
       !if (istep*dtau > 5.d0/efcolf_arr(1,ns/2)) exit

       ! TODO: TESTING at half radius
       if (istep*dtau > dummy) then
          write(2000+icpu,*) ipart, istep, istep*dtau/(5.d0/efcolf_arr(1,ns/2)), z(1),&
               (z(1)-0.5)**2, z(2), z(3)
          dummy = dummy + debugstep
       endif

       if (z(1)<0.d0 .OR. z(1)>1.d0) then
          print *, "particle lost:"
          print *, ipart, z
          exit
       endif

       call regst(z,dtau,ierr)
!
       call stost(z,dtau,iswmod,ierr)
       write(3000,*) istep*dtau, z
   enddo
   print *,istep,'  steps'
!
  
enddo
! 
!
!
  end
