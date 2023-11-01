!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine window_filter(n,nw,arr_in,arr_out)
!
  implicit none
!
  integer, intent(in) :: n,nw
  double precision, dimension(n), intent(in) :: arr_in
  double precision, dimension(n), intent(out) :: arr_out
  integer :: nwa,i
!
  do i=1,n
    nwa=min(nw,i-1,n-i)
    arr_out(i)=sum(arr_in(i-nwa:i+nwa))/(2*nwa+1)
  enddo
!
  return
  end subroutine window_filter
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
!
  subroutine gauss_filter(n,sigma,arr_in,arr_out)
!
  implicit none
!
  double precision, parameter :: sigma_min=1.d-10
!
  integer, intent(in) :: n
  double precision, intent(in) :: sigma
  double precision, dimension(n), intent(in)  :: arr_in
  double precision, dimension(n), intent(out) :: arr_out
  double precision, dimension(:), allocatable :: gwin
  double precision :: dummy,sigma_loc
  integer :: nwa,i,nwa1,nwa2
!
  nwa=int(sigma*7.d0)
  allocate(gwin(-nwa:nwa))
  sigma_loc=max(sigma,sigma_min)
!
  do i=-nwa,nwa
    gwin(i)=exp(-0.5d0*(dble(i)/sigma_loc)**2)
  enddo
!
  dummy=sum(gwin)
  gwin=gwin/dummy
!
  do i=nwa+1,n-nwa
    arr_out(i)=sum(arr_in(i-nwa:i+nwa)*gwin)
  enddo
!
  nwa1=nwa+1
  nwa2=nwa1+1
  do i=1,nwa
    arr_out(i)=arr_out(nwa2)*dble(i-nwa1)+arr_out(nwa1)*dble(nwa2-i)
  enddo
!
  nwa2=n-nwa
  nwa1=nwa2-1
  do i=nwa2+1,n
    arr_out(i)=arr_out(nwa2)*dble(i-nwa1)+arr_out(nwa1)*dble(nwa2-i)
  enddo
!
  end subroutine gauss_filter
!
!cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
