program test_class_return_segments
    use class_return_segments_mod, only : partition_valid_samples
    implicit none

    double precision :: x(6),begins(6),ends(6)
    logical :: valid(6)
    integer :: n

    x=(/0.d0,1.d0,2.d0,3.d0,4.d0,5.d0/)
    valid=(/.true.,.true.,.false.,.false.,.true.,.true./)
    call partition_valid_samples(x,valid,begins,ends,n)
    if(n.ne.2) error stop 'valid samples were not split into two segments'
    if(abs(begins(1)-0.d0).gt.1.d-14 .or. abs(ends(1)-1.d0).gt.1.d-14) then
        error stop 'left valid segment is wrong'
    endif
    if(abs(begins(2)-4.d0).gt.1.d-14 .or. abs(ends(2)-5.d0).gt.1.d-14) then
        error stop 'right valid segment is wrong'
    endif

    valid=.true.
    call partition_valid_samples(x,valid,begins,ends,n)
    if(n.ne.1 .or. abs(begins(1)-x(1)).gt.1.d-14 .or. &
       abs(ends(1)-x(6)).gt.1.d-14) then
        error stop 'connected valid samples were changed'
    endif

    valid=(/.false.,.true.,.true.,.true.,.true.,.false./)
    call partition_valid_samples(x,valid,begins,ends,n)
    if(n.ne.1 .or. abs(begins(1)-x(2)).gt.1.d-14 .or. &
       abs(ends(1)-x(5)).gt.1.d-14) then
        error stop 'endpoint-invalid segment was not retained'
    endif
end program test_class_return_segments
