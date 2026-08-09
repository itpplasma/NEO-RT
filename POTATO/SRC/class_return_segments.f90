module class_return_segments_mod
    implicit none

    integer :: nsegments=0
    integer :: current_segment=0
    double precision, allocatable :: segment_beg(:),segment_end(:)
    !$omp threadprivate(nsegments,current_segment,segment_beg,segment_end)

contains

    subroutine clear_class_segments()
        nsegments=0
        current_segment=0
        if(allocated(segment_beg)) deallocate(segment_beg)
        if(allocated(segment_end)) deallocate(segment_end)
    end subroutine clear_class_segments

    subroutine set_class_segments(n,begins,ends)
        integer, intent(in) :: n
        double precision, intent(in) :: begins(:),ends(:)

        call clear_class_segments()
        if(n.le.0) return
        allocate(segment_beg(n),segment_end(n))
        segment_beg=begins(1:n)
        segment_end=ends(1:n)
        nsegments=n
        current_segment=0
    end subroutine set_class_segments

    subroutine partition_valid_samples(x,valid,begins,ends,n)
        double precision, intent(in) :: x(:)
        logical, intent(in) :: valid(:)
        double precision, intent(out) :: begins(:),ends(:)
        integer, intent(out) :: n
        integer :: i,j,npoints

        n=0
        npoints=size(x)
        if(size(valid).ne.npoints) error stop 'sample topology arrays differ in size'
        if(size(begins).lt.npoints) error stop 'sample topology begin array is too small'
        if(size(ends).lt.npoints) error stop 'sample topology end array is too small'

        i=1
        do while(i.le.npoints)
            if(.not.valid(i)) then
                i=i+1
                cycle
            endif
            j=i
            do while(j.lt.npoints)
                if(.not.valid(j+1)) exit
                j=j+1
            enddo
            n=n+1
            begins(n)=x(i)
            ends(n)=x(j)
            i=j+1
        enddo
    end subroutine partition_valid_samples

end module class_return_segments_mod
