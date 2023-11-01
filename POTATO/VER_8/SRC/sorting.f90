  subroutine sortin(a,ipoi,n)
!
! Permutates elements of array a to the increasing sequence
! a(i),       i=1,..,n - original sequence
! a(ipoi(i)), i=1,..,n - increasing sequence
!
  implicit none
!
  integer :: n,i,j,isave
  integer,          dimension(n) :: ipoi
  double precision, dimension(n) :: a
!
  do i=1,n
    ipoi(i)=i
  enddo
!
  do i=1,n-1
    do j=i+1,n
      if(a(ipoi(j)).lt.a(ipoi(i))) then
        isave=ipoi(i)
        ipoi(i)=ipoi(j)
        ipoi(j)=isave
      endif
    enddo
  enddo
!   
  return
  end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
