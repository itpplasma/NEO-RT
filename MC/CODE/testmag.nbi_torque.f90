!
!
      use magfield_mod, only : ierrfield
!
      implicit none
!
      integer :: i,n
      double precision :: phi
      double precision, dimension(5) :: y,dery
      double precision x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl
      dimension x(3),bder(3),hcovar(3),hctrvr(3),hcurl(3)
!
      n=1000
!
      do i=1,n
        x(1)=dfloat(i)/dfloat(n)
        x(2:3)=0.d0
!
        call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
        write (1,*) x(1),bmod
!
     enddo
!
      do i=1,n
        x(1)=0.5d0
        x(2)=0.d0
        x(3)=2.d0*3.1415927d0*dfloat(i)/dfloat(n)
!
        call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
        write (3,*) x(3),bmod
!
     enddo
!
      do i=1,n
        x(1)=0.5d0
        x(3)=0.d0
        x(2)=2.d0*3.1415927d0*dfloat(i)/dfloat(n)
!
        call magfie(x,bmod,sqrtg,bder,hcovar,hctrvr,hcurl)
!
        write (2,*) x(2),bmod
!
     enddo
!
     end


