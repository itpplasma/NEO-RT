!
  subroutine plot_orbit_set(wrorb,iun1,iun2,npoi,npoi_pot,dtau,icase,toten,perpinv,        &
                            trapped,copass,ctrpass,triplet,sst,rst,zst,sst2,rst2,zst2,     &
                            sigma_x,s_tp,rxpoi,zxpoi,raxis_ctr,zaxis_ctr,raxis_co,zaxis_co)
!
  use orbit_dim_mod,     only : neqm,next,write_orb,iunit1
  use vparzero_line_mod, only : nvpline,s_line,psi_line,rinb,zinb
!
  implicit none
!
  double precision, parameter :: eps_fixp=1.d-5
!
  logical :: wrorb,trapped,copass,ctrpass,triplet
  integer :: icase,i,npoi,npoi_pot,iun1,iun2
  double precision :: eps,toten,perpinv,sigma_x,s_tp,delL_fixp
  double precision :: rxpoi,zxpoi,raxis_ctr,zaxis_ctr,raxis_co,zaxis_co
  double precision :: bmod,phi_elec,p2,sst,rst,zst,sst2,rst2,zst2
  double precision :: alam2,hs,hr,hz,s,dR_ds,dZ_ds,p_phi
  double precision :: dtau,taub,delphi,rbeg,zbeg,rend,zend,sigma,potato
  double precision :: b_marg,e_marg,s_beg,s_end
  double precision, dimension(neqm) :: z,z_start,vz
  double precision, dimension(next) :: extraset
!
  eps=0.d0
  delL_fixp=psi_line(1,0)*eps_fixp
!
  write_orb=wrorb
  iunit1=iun1
!
! Plot orbits:
!
! Co-passing particles:
!
  if(icase.eq.4) then
    open(1,file='co-passing.dat')
    open(2,file='cut_taub_dephi_co_passing.dat')
    rbeg=rinb
    zbeg=zinb
    b_marg=0.d0
    rend=raxis_co
    zend=zaxis_co
    e_marg=1.d0
!
    call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
    close(1)
    close(2)
!
    open(1,file='ctr-passing.dat')
    open(2,file='cut_taub_dephi_ctr_passing.dat')
    rbeg=rinb
    zbeg=zinb
    b_marg=0.d0
    rend=raxis_ctr
    zend=zaxis_ctr
    e_marg=1.d0
!
    call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
    close(1)
    close(2)
!
    return
  endif
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
        b_marg=1.d0
      else
! connects to potatoes
        rbeg=rxpoi
        zbeg=zxpoi
        b_marg=1.d0
        rend=rst
        zend=zst
        e_marg=1.d0
        potato=1.d0
        sigma=-1.d0
!
        call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
        rbeg=rst
        zbeg=zst
      endif
    else
      rbeg=rst
      zbeg=zst
      b_marg=1.d0
    endif
!
    rend=raxis_co
    zend=zaxis_co
    e_marg=1.d0
    potato=0.d0
    sigma=1.d0
!
    call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
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
        b_marg=1.d0
      else
! connects to potatoes
        rbeg=rxpoi
        zbeg=zxpoi
        b_marg=1.d0
        rend=rst
        zend=zst
        e_marg=1.d0
        potato=1.d0
        sigma=1.d0
!
        call passing_orbits(npoi_pot,rbeg,zbeg,rend,zend,sigma,potato)
!
        rbeg=rst
        zbeg=zst
      endif
    else
      rbeg=rst
      zbeg=zst
      b_marg=1.d0
    endif
!
    rend=raxis_ctr
    zend=zaxis_ctr
    e_marg=1.d0
    potato=0.d0
    sigma=-1.d0
!
    call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
    close(1)
    close(2)
  endif
!
! End counter-passing particles
!
! Trapped particles:
!
  if(icase.eq.2) then
    trapped=.true.
    s_beg=s_tp
    s_end=s_line(nvpline)
  elseif(icase.eq.3) then
    trapped=.true.
    s_beg=sst
    s_end=s_tp
  else
    trapped=.false.
  endif
!
  if(trapped) then
    open(1,file='trapped.dat')
    open(2,file='cut_taub_dephi_trapped.dat')
!
    hs=(s_end-s_beg)/dble(npoi)
!
    do i=0,npoi
      z(2)=0.d0
      s=s_beg+hs*(dble(i)+eps)
!
      call vparzero_line(s,z(1),z(3),dR_ds,dZ_ds)
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
      call find_bounce(next,velo_ext,dtau,z,taub,delphi,extraset)
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
  if(icase.eq.3) then
    if(sigma_x.lt.0) then
      sigma=1.d0
      open(1,file='co-passing2.dat')
      open(2,file='cut_taub_dephi_co_passing2.dat')
      rbeg=rinb
      zbeg=zinb
      b_marg=0.d0
      rend=rst2
      zend=zst2
      e_marg=1.d0
!
      call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
      close(1)
      close(2)
      sigma=-1.d0
      open(1,file='ctr-passing2.dat')
      open(2,file='cut_taub_dephi_ctr_passing2.dat')
      rend=rxpoi
      zend=zxpoi
      e_marg=1.d0
!
      call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
      write(2,*) ' '
      rbeg=rst2
      zbeg=zst2
      b_marg=1.d0
!
      call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
      close(1)
      close(2)
    else
      sigma=-1.d0
      open(1,file='ctr-passing2.dat')
      open(2,file='cut_taub_dephi_ctr_passing2.dat')
      rbeg=rinb
      zbeg=zinb
      b_marg=0.d0
      rend=rst2
      zend=zst2
      e_marg=1.d0
!
      call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
      close(1)
      close(2)
      sigma=1.d0
      open(1,file='co-passing2.dat')
      open(2,file='cut_taub_dephi_co_passing2.dat')
      rend=rxpoi
      zend=zxpoi
      e_marg=1.d0
!
      call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
      write(2,*) ' '
      rbeg=rst2
      zbeg=zst2
      b_marg=1.d0
!
      call passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
      close(1)
      close(2)
    endif
  endif
!
!---------------------
!
  contains
!
!---------------------
!
  subroutine passing_orbits(npoi,rbeg,zbeg,rend,zend,sigma,potato)
!
  implicit none
!
  integer :: npoi
  double precision :: rbeg,zbeg,rend,zend,sigma,potato,hR,hZ,del_R,del_Z,fac_abs
!
  del_R=rend-rbeg
  del_Z=zend-zbeg
  fac_abs=delL_fixp/sqrt(del_R**2+del_Z**2)
  if(fac_abs*(b_marg+e_marg).gt.0.5d0) return
  del_R=del_R*fac_abs
  del_Z=del_Z*fac_abs
  hR=(rend-rbeg-del_R*(b_marg+e_marg))/dble(npoi)
  hZ=(zend-zbeg-del_Z*(b_marg+e_marg))/dble(npoi)
!
  do i=0,npoi
    z(1)=rbeg+hR*dble(i)+del_R*b_marg
    z(2)=0.d0
    z(3)=zbeg+hZ*dble(i)+del_Z*e_marg
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
    call find_bounce(next,velo_ext,dtau,z,taub,delphi,extraset)
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
  subroutine velo_ext(dtau2,z2,vz2)
!
  use orbit_dim_mod, only : neqm,next,ndim
!
  double precision :: dtau2
  double precision, dimension(ndim) :: z2,vz2
!
  call velo(dtau2,z2(1:neqm),vz2(1:neqm))
!
  vz2(neqm+1)=z2(1)
  vz2(neqm+2)=z2(3)
  vz2(neqm+3)=z2(4)*z2(5)
!
  end subroutine velo_ext
!
!---------------------
!
  end subroutine plot_orbit_set
