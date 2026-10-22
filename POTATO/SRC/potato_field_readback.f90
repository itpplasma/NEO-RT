program potato_field_readback
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    integer :: n,i,u
    real(dp), allocatable :: r(:),z(:),br(:),bi(:)
    complex(dp) :: a
    open(newunit=u,file='query.bin',access='stream',form='unformatted',status='old')
    read(u) n
    allocate(r(n),z(n),br(n),bi(n))
    read(u) r,z
    close(u)
    !$omp parallel do private(a)
    do i=1,n
        call bmod_pert(r(i),z(i),a)
        br(i)=real(a,dp)*1.e-4_dp
        bi(i)=aimag(a)*1.e-4_dp
    end do
    !$omp end parallel do
    open(newunit=u,file='readback.bin',access='stream',form='unformatted',status='replace')
    write(u) br,bi
    close(u)
end program potato_field_readback
