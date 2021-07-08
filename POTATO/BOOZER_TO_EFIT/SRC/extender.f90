!
  module extender_mod
    logical :: prop=.true.
    integer, parameter :: npoi_ext=10000
    double precision :: hdec=3.d-2
    double precision :: s_bou,R_0,Z_0,twopi
    double precision, dimension(0:npoi_ext) :: rho_w,tht_w,ds_drho,dt_drho,thtb
  end module extender_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine extender(R_in,Z_in,outside,s_out,theta_out)
!
  use extender_mod, only : prop,npoi_ext,R_0,Z_0,rho_w,tht_w,ds_drho,dt_drho, &
                           twopi,s_bou,hdec,thtb
!
  implicit none
!
  logical :: outside
  integer :: i
  double precision :: R_in,Z_in,s_out,theta_out,R,Z,rho,tht,w,wm,rho_c
  double precision :: htheta,theta,phi,dR_ds,dR_dt,dZ_ds,dZ_dt,det_J
  double precision :: dsdrho,dtdrho,drho,theta_b
  integer, dimension(1) :: ibt
!
!----------- 1st call --------------------------------------------------------
  if(prop) then
    prop=.false.
    twopi = atan(1.d0)*8.d0
    htheta=twopi/dble(npoi_ext)
    phi=0.d0
    theta=0.d0
!
    call boozer_data_in_symfluxcoord(s_bou,theta,phi,R,Z,dR_ds,dR_dt,dZ_ds,dZ_dt)
!
    if(dR_ds*dZ_dt - dR_dt*dZ_ds.lt.0.d0) htheta=-htheta
    R_0=0.d0
!
    do i=1,npoi_ext
      theta=dble(i)*htheta
!
      call boozer_data_in_symfluxcoord(s_bou,theta,phi,R,Z,dR_ds,dR_dt,dZ_ds,dZ_dt)
!
      R_0=R_0+R
      Z_0=Z_0+Z
    enddo
!
    R_0=R_0/dble(npoi_ext)
    Z_0=Z_0/dble(npoi_ext)
!
    do i=0,npoi_ext
      theta=dble(i)*htheta
!
      call boozer_data_in_symfluxcoord(s_bou,theta,phi,R,Z,dR_ds,dR_dt,dZ_ds,dZ_dt)
!
      rho_w(i)=sqrt((R-R_0)**2+(Z-Z_0)**2)
      tht_w(i)=atan2(Z-Z_0,R-R_0)
      if(tht_w(i).lt.tht_w(0)) tht_w(i)=tht_w(i)+twopi
      det_J = dR_ds*dZ_dt - dR_dt*dZ_ds
      ds_drho(i) = (dZ_dt*cos(tht_w(i))-dR_dt*sin(tht_w(i)))/det_J
      dt_drho(i) = (dR_ds*sin(tht_w(i))-dZ_ds*cos(tht_w(i)))/det_J
      thtb(i) = theta
    enddo
!
  endif
!----------- end of the 1st call --------------------------------------------
!
  if(abs(R_in-R_0).lt.epsilon(1.d0)) then
    outside=.false.
    return
  endif
!
  rho = sqrt((R_in-R_0)**2 + (Z_in-Z_0)**2)
  tht = atan2(Z_in-Z_0,R_in-R_0)
  if(tht.lt.tht_w(0)) tht=tht+twopi
!
  call binsrc(tht_w,0,npoi_ext,tht,i)
!
  wm=(tht_w(i)-tht)/(tht_w(i)-tht_w(i-1))
  w=1.d0-wm
  rho_c = rho_w(i-1)*wm + rho_w(i)*w
!
  if(rho .le. rho_c) then
    outside=.false.
    return
  endif
!
  outside=.true.
  dsdrho = ds_drho(i-1)*wm + ds_drho(i)*w
  dtdrho = dt_drho(i-1)*wm + dt_drho(i)*w
  theta_b = thtb(i-1)*wm + thtb(i)*w
  drho = hdec*atan2((rho-rho_c), hdec)
  s_out = s_bou + dsdrho*drho
  theta_out = theta_b + dtdrho*drho
!
  return
  end subroutine extender
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
