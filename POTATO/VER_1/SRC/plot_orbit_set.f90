!
  subroutine plot_orbit_set(wrorb,iun1,iun2,npoi,npoi_pot,dtau,eps,toten,perpinv, &
                            trapped,copass,ctrpass,triplet,rst,zst,sigma_x,Z_tp,  &
                            rxpoi,zxpoi,raxis_ctr,zaxis_ctr,raxis_co,zaxis_co)
!
  use orbit_dim_mod,     only : neqm,next,write_orb,iunit1
  use vparzero_line_mod, only : Z_up
!
  implicit none
!
  logical :: wrorb,trapped,copass,ctrpass,triplet
  integer :: i,npoi,npoi_pot,iun1,iun2
  double precision :: eps,toten,perpinv,sigma_x,Z_tp
  double precision :: rxpoi,zxpoi,raxis_ctr,zaxis_ctr,raxis_co,zaxis_co
  double precision :: bmod,phi_elec,p2,rst,zst,alam2,hr,hz,dR_dZ,p_phi
  double precision :: dtau,taub,delphi,rbeg,zbeg,rend,zend,sigma,potato
  double precision, dimension(neqm) :: z,z_start,vz
  double precision, dimension(next) :: extraset
!
  write_orb=wrorb
  iunit1=iun1
!
! Plot orbits:
!
! Co-passing particles:
!
  if(copass) then
    open(1,file='co-passing.dat')
    open(2,file='cut_taub_dephi_co_passing.dat')
!
    if(triplet) then
      if(sigma_x.gt.0.d0) then
! between separatrix and axis only
        rbeg=rxpoi
        zbeg=zxpoi
      else
! connects to potatoes
        rbeg=rxpoi
        zbeg=zxpoi
        rend=rst
        zend=zst
        potato=1.d0
        sigma=-1.d0
!
        call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato,eps)
!
        rbeg=rst
        zbeg=zst
      endif
    else
      rbeg=rst
      zbeg=zst
    endif
!
    rend=raxis_co
    zend=zaxis_co
    potato=0.d0
    sigma=1.d0
!
    call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato,eps)
!
    close(1)
    close(2)
  endif
!
! End co-passing particles
!
! Counter-passing particles:
!
  if(ctrpass) then
    open(1,file='ctr-passing.dat')
    open(2,file='cut_taub_dephi_ctr_passing.dat')
!
    if(triplet) then
      if(sigma_x.lt.0.d0) then
! between separatrix and axis only
        rbeg=rxpoi
        zbeg=zxpoi
      else
! connects to potatoes
        rbeg=rxpoi
        zbeg=zxpoi
        rend=rst
        zend=zst
        potato=1.d0
        sigma=1.d0
!
        call passing_orbits(npoi_pot,rbeg,zbeg,rend,zend,sigma,potato,eps)
!
        rbeg=rst
        zbeg=zst
      endif
    else
      rbeg=rst
      zbeg=zst
    endif
!
    rend=raxis_ctr
    zend=zaxis_ctr
    potato=0.d0
    sigma=-1.d0
!
    call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato,eps)
!
    close(1)
    close(2)
  endif
!
! End counter-passing particles
!
! Trapped particles:
!
  if(trapped) then
    open(1,file='trapped.dat')
    open(2,file='cut_taub_dephi_trapped.dat')
!
    hZ=(Z_up-z_tp)*(1.d0-eps)/(dble(npoi)+eps)
!
    do i=0,npoi
      z(2)=0.d0
      z(3)=z_tp+hZ*(dble(i)+eps)
!
      call vparzero_line(z(3),z(1),dR_dZ)
!
      call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
      p2=toten-phi_elec
      if(p2.le.0.d0) then
        print *,'tt: negative kinetic energy '
        exit
       endif
      z(4)=sqrt(p2)
!
      z(5)=0.d0
!
      extraset=0.d0
!
      call find_bounce(dtau,z,taub,delphi,extraset)
!
      call get_tormom(z,p_phi)
!
      write(iun2,*) z(1),z(3),taub,delphi,p_phi
      write(iun1,*) 'NaN NaN NaN NaN NaN NaN'
    enddo
!
    close(1)
    close(2)
  endif
!
! End trapped particles
!
!---------------------
!
  contains
!
!---------------------
!
  subroutine passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato,eps)
!
  implicit none
!
  integer :: npoi
  double precision :: rbeg,zbeg,rend,zend,sigma,potato,eps,hR,hZ
!
!
  hR=(rend-rbeg)/(dble(npoi)+2.d0*eps)
  hZ=(zend-zbeg)/(dble(npoi)+2.d0*eps)
!
  do i=0,npoi
    z(1)=rbeg+hR*(dble(i)+eps)
    z(2)=0.d0
    z(3)=zbeg+hZ*(dble(i)+eps)
!
    call get_bmod_and_Phi(z(1:3),bmod,phi_elec)
!
    p2=toten-phi_elec
    if(p2.le.0.d0) then
      print *,'tt passing_orbits: negative kinetic energy '
      exit
    endif
    z(4)=sqrt(p2)
!
    alam2=1.d0-perpinv*bmod/p2
    if(alam2.lt.0.d0) then
      print *,'tt passing_orbits: negative lambda^2 '
      exit
    endif
    z(5)=sign(sqrt(alam2),sigma)
!
    extraset=0.d0
!
    call find_bounce(dtau,z,taub,delphi,extraset)
!
    call get_tormom(z,p_phi)
!
    write(iun2,*) z(1),z(3),taub,delphi,p_phi,potato
    write(iun1,*) 'NaN NaN NaN NaN NaN NaN'
  enddo
!
  end subroutine passing_orbits
!
!---------------------
!
  end subroutine plot_orbit_set
