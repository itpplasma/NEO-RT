!
  subroutine intsrc(N,npoi,xx,fg_point,step,ibeg,iend)
!
! input : xx,fg_point,step
! output : n_index,i_part 
! fg_point - first grid poing
! n_index - index of the node
! xx - point 
!
  implicit none
! 
  integer :: i_part,N,npoi,ibeg,iend
  double precision :: n_index,step,xx,fg_point
!
  n_index=(xx-fg_point)/step 
  i_part=nint(n_index)
!
  ibeg=max(1,i_part-npoi/2)
  iend=ibeg+npoi-1
  if(iend.gt.N) then
    iend=N
    ibeg=iend+1-npoi
  endif
!  
  return
  end
!ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
