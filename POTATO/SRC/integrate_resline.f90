!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  module integrate_resline_mod
! Contains output of the routine integrate_resline
! np          - number of integration points:
    integer :: np
! w_orbs(np)   - weights of integration points:
    double precision, dimension(:),   allocatable :: w_orbs
! z_orbs(2,np) - coordinates of integration points:
    double precision, dimension(:,:), allocatable :: z_orbs
  end module integrate_resline_mod
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine integrate_resline(fun2d)
!
! Finds all 2D resonant lines crossing the boundary of the rectangular 
! region [0:1][0:1]. Resonant line is zero line of the 2D function computed
! by the external subroutine fun2d.
! Prepares the integration along those lines specified by np points with
! coordinates z_orb(2,np) and computes weights w_orb(np). Integral of any 
! function f(z) along the line will be sum over index i of f(z_orbs(:,i))*w_orbs(i)
!
  use integrate_resline_mod, only : np,w_orbs,z_orbs
  use find_all_roots_mod,    only : customgrid,nsearch_min,nroots,roots
!
  implicit none
!
  double precision, parameter :: f_lev=0.d0, err_dist=1.d-10 !parameters for Newton adjustment
!
  integer :: ifsw,i,ierr,nstp,n,j,n_min,n_max,i_exit,ibou
  double precision :: x1in,x2in,f,df,h,delth_rec,h_in,h_in_loc,bou,w,h_prev
  double precision, dimension(2) :: z,df_dz
  double precision, dimension(:),   allocatable :: t_orb,w_orb,gradm1
  double precision, dimension(:,:), allocatable :: zbeg,dummy2d,z_orb
!
  external :: fun2d
!
! settings for root finder:
  customgrid=.false.
  nsearch_min=100
  x1in=0.d0
  x2in=1.d0
!
  h_in=3.d-2 !0.01d0         !maximum resline integration step
  delth_rec=3.d0*h_in !maximum angle between steps
  n_min=10            !minimum number of points per resline
!
! Search for roots - resline starting points at the domain boundary:
!
! First run - determine the number of roots:
  nstp=0
!
  do ifsw=1,4
!
    call find_all_roots(fun1d,x1in,x2in,ierr)
!
    do i=1,nroots
!
      call fun1d(roots(i),f,df)
!
! keep only those with positive derivative over boundary parameter:
      select case(ifsw)
      case(1,4)
        if(df.gt.0.d0) then
          nstp=nstp+1
        endif
      case(2,3)
        if(df.lt.0.d0) then
          nstp=nstp+1
        endif
      end select
    enddo
  enddo
!
  allocate(zbeg(2,nstp))
  nstp=0
!
! Second run - store the roots:
  do ifsw=1,4
!
    call find_all_roots(fun1d,x1in,x2in,ierr)
!
    do i=1,nroots
!
      call fun1d(roots(i),f,df)
!
      select case(ifsw)
      case(1)
        if(df.gt.0.d0) then
          nstp=nstp+1
          zbeg(1,nstp)=roots(i)
          zbeg(2,nstp)=0.d0
        endif
      case(2)
        if(df.lt.0.d0) then
          nstp=nstp+1
          zbeg(1,nstp)=0.d0
          zbeg(2,nstp)=roots(i)
        endif
      case(3)
        if(df.lt.0.d0) then
          nstp=nstp+1
          zbeg(1,nstp)=roots(i)
          zbeg(2,nstp)=1.d0
        endif
      case(4)
        if(df.gt.0.d0) then
          nstp=nstp+1
          zbeg(1,nstp)=1.d0
          zbeg(2,nstp)=roots(i)
        endif
      end select
    enddo
  enddo
!
! End search for roots - resline starting points at the domain boundary
!
!........
!
! Integration of reslines
!
  np=0 !initialize the number of output points
!
  do i=1,nstp         !loop over reslines
!
! first integration of reslines - determine integration step and array sizes:
!
    h_in_loc=h_in
    do 
      n=1
      z=zbeg(:,i)
      h_prev=h_in_loc
!
      do
!
        call choose_step(fun2d,z,delth_rec,h_in_loc,h)
!
        h=min(h,h_prev*2.d0)
        h_prev=h
!
        call level_set_step_2D(fun2d, h, z)
!
! determine exit boundary:
        if(z(1).lt.0.d0) then
          i_exit=1
          exit
        elseif(z(1).gt.1.d0) then
          i_exit=2
          exit
        elseif(z(2).lt.0.d0) then
          i_exit=3
          exit
        elseif(z(2).gt.1.d0) then
          i_exit=4
          exit
        endif
        n=n+1
!
! use Newton adjustment only within the domain:
!
        call newt_adjust_root_2d(fun2d,f_lev,err_dist,z)
!
      enddo
!
      if(n.ge.n_min) exit
      h_in_loc=h*dble(n)/dble(n_min)
    enddo
!
    if(allocated(z_orb)) deallocate(z_orb,t_orb,w_orb,gradm1)
    allocate(z_orb(2,0:n),t_orb(0:n),w_orb(0:n),gradm1(0:n))
!
! second integration of reslines - use Newton adjustment and store the points:
!
    z=zbeg(:,i)
    n=0
    z_orb(:,n)=z
    t_orb(n)=0.d0
!
    call fun2d(z,f,df_dz)
!
    gradm1(n)=1.d0/sqrt(df_dz(1)**2+df_dz(2)**2)
    h_prev=h_in_loc
!
    do
!
      call choose_step(fun2d,z,delth_rec,h_in_loc,h)
!
      h=min(h,h_prev*2.d0)
      h_prev=h
!
      call level_set_step_2D(fun2d, h, z)
!
      n=n+1
      z_orb(:,n)=z
      t_orb(n)=t_orb(n-1)+h
!
! determine exit boundary:
      if(z(1).lt.0.d0) then
        i_exit=1
        exit
      elseif(z(1).gt.1.d0) then
        i_exit=2
        exit
      elseif(z(2).lt.0.d0) then
        i_exit=3
        exit
      elseif(z(2).gt.1.d0) then
        i_exit=4
        exit
      endif
!
! use Newton adjustment only within the domain:
!
      call newt_adjust_root_2d(fun2d,f_lev,err_dist,z)
!
      z_orb(:,n)=z
!
      call fun2d(z,f,df_dz)
!
      gradm1(n)=1.d0/sqrt(df_dz(1)**2+df_dz(2)**2)
    enddo
!
! cut last integration step to stop at the boundary (linear order):
!
    select case(i_exit)
    case(1)
      ibou=1
      bou=0.d0
    case(2)
      ibou=1
      bou=1.d0
    case(3)
      ibou=2
      bou=0.d0
    case(4)
      ibou=2
      bou=1.d0
    end select
!
    t_orb(n)=t_orb(n-1)+(t_orb(n)-t_orb(n-1)) &
            *(bou-z_orb(ibou,n-1))/(z_orb(ibou,n)-z_orb(ibou,n-1))
!
! determine integration weights of inner orbit points (trapez fomula):
!
    w_orb=0.d0
! intergation over inner intervals:
    do j=2,n-1
      w=0.5d0*(t_orb(j)-t_orb(j-1))
      w_orb(j-1)=w_orb(j-1)+w
      w_orb(j)=w_orb(j)+w
    enddo
! intergation over first interval:
    w=0.5d0*(t_orb(1)-t_orb(0))/(t_orb(2)-t_orb(1))
    w_orb(1)=w_orb(1)+(1.d0+w)*(t_orb(1)-t_orb(0))
    w_orb(2)=w_orb(2)-w*(t_orb(1)-t_orb(0))
! intergation over last  interval:
    w=0.5d0*(t_orb(n)-t_orb(n-1))/(t_orb(n-1)-t_orb(n-2))
    w_orb(n-1)=w_orb(n-1)+(1.d0+w)*(t_orb(n)-t_orb(n-1))
    w_orb(n-2)=w_orb(n-2)-w*(t_orb(n)-t_orb(n-1))
!
! store weights and points in output arrays
!
    if(np.eq.0) then
      np=n-1
      if(allocated(w_orbs)) deallocate(w_orbs,z_orbs)
      allocate(w_orbs(np),z_orbs(2,np))
      z_orbs(:,:)=z_orb(:,1:n-1)
      w_orbs=w_orb(1:n-1)*gradm1(1:n-1)
    else
      allocate(dummy2d(3,np))
      dummy2d(1:2,:)=z_orbs
      dummy2d(3,:)=w_orbs
      deallocate(z_orbs,w_orbs)
      allocate(w_orbs(np+n-1),z_orbs(2,np+n-1))
      z_orbs(:,1:np)=dummy2d(1:2,:)
      w_orbs(1:np)=dummy2d(3,:)
      z_orbs(:,np+1:np+n-1)=z_orb(:,1:n-1)
      w_orbs(np+1:np+n-1)=w_orb(1:n-1)*gradm1(1:n-1)
      np=np+n-1
      deallocate(dummy2d)
    endif
!
  enddo
!
! End integration of reslines
!
!--------------------------------------------------
  contains
!--------------------------------------------------
!
  subroutine fun1d(x,f,df)
!
  implicit none
!
  double precision :: x,f,df
  double precision, dimension(2) :: z,df_dz
!
  select case(ifsw)
  case(1)
    z(1)=x
    z(2)=0.d0
!
    call fun2d(z,f,df_dz)
!
    df=df_dz(1)
  case(2)
    z(1)=0.d0
    z(2)=x
!
    call fun2d(z,f,df_dz)
!
    df=df_dz(2)
  case(3)
    z(1)=x
    z(2)=1.d0
!
    call fun2d(z,f,df_dz)
!
    df=df_dz(1)
  case(4)
    z(1)=1.d0
    z(2)=x
!
    call fun2d(z,f,df_dz)
!
    df=df_dz(2)
  end select
!
  end subroutine fun1d
!--------------------------------------------------
!
  end subroutine integrate_resline
!
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
