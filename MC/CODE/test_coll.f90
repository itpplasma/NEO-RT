!
  implicit none
!
  integer, parameter :: nl=50,nv=1000
  integer :: i,j,k,n,iswmode,ierr
  double precision :: am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ealpha, &
                      v0,dchichi,slowrate,dchichi_norm,slowrate_norm,dtauc,hv,hl
!
  double precision, dimension(5) :: z
  double precision, dimension(nv,-nl:nl) :: fun
!
  am1=2.d0
  am2=3.d0
  Z1=1.d0
  Z2=1.d0
  densi1=1.d14
  densi2=0.d13
  tempi1=1.d4
  tempi2=1.d4
  tempe=1.d4
  ealpha=3.5d6
!
  call loacol_alpha(am1,am2,Z1,Z2,densi1,densi2,tempi1,tempi2,tempe,ealpha, &
                    v0,dchichi,slowrate,dchichi_norm,slowrate_norm)
!
  print *,'v0 = ',v0
  print *,'dchichi = ',dchichi, &
  (23.d0-log(max(epsilon(1.d0), &
   sqrt(densi1*Z1**2/tempi1)*2.d0*Z1*(4.d0+am1)/(4.d0*tempi1+am1*ealpha)))) &
  *1.8d-7*densi1*4.d0/2.d0/3.5d6**1.5d0  /4.d0
  print *,'slowrate = ',slowrate, &
  (24.d0-log(sqrt(densi1)/tempe))*1.6d-9*densi1*4.d0/4.d0/tempe**1.5d0
  print *,'dchichi_norm = ',dchichi_norm
  print *,'slowrate_norm = ',slowrate_norm
!
  fun=0.d0
  dtauc=0.0001d0/slowrate_norm
  iswmode=1
  hl=1.d0/nl
  hv=1.1d0/nv
  n=10000000
!
  z=0.d0
  z(4)=1.d0
!
  do i=1,n
    call stost(z,dtauc,iswmode,ierr)
    k=int(z(4)/hv)
    k=max(1,min(nv,k))
    j=int(z(5)/hl)
    j=max(-nl,min(nl,j))
    fun(k,j)=fun(k,j)+1.d0
  enddo
!
  fun=fun/(n*hl*hv)
!
  do i=1,nv
    write (1,*) i*hv,sum(fun(i,:))*hl
  enddo
!
  stop
  end
