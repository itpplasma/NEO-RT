module topology_gap_quadrature_mod
    !! Evaluate a topology bracket with the generated square-root endpoint map.
    !!
    !! The physical matrix callback has already been sampled on each open
    !! branch.  The returned value is an average contribution compatible with
    !! the existing sampler ABI; the sampler multiplies it by the geometric
    !! gap width.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use potato_input_mod, only : topology_gap_fit_points
    use potato_symbolic_kernel_mod, only : potato_gap_square_map, &
        potato_gap_sqrt_coefficient,potato_gap_sqrt_integral
    use sample_matrix_out_mod, only : sample_matrix_success, &
        sample_matrix_contribution_unresolved,xbeg,xend, &
        npoi,xarr,amat_arr,topology_arr
    implicit none
    private
    public :: estimate_topology_gap_from_samples

    integer, parameter :: minimum_gap_fit_points=2

contains

    subroutine estimate_topology_gap_from_samples(gap_lo,gap_hi,average,ierr)
        !! Estimate the omitted contribution from already sampled branch data.
        !!
        !! On a reflecting boundary the matrix contribution has the leading
        !! form C/sqrt(|x-x_b|).  The generated coefficient and integral
        !! kernels integrate that limiting expression exactly.  The nearest
        !! branch samples provide C, with the largest of the configured local
        !! fit points used as the one-sided envelope.  No new orbit or
        !! resonance solve is performed here.
        real(real64), intent(in) :: gap_lo,gap_hi
        real(real64), intent(out) :: average
        integer, intent(out) :: ierr

        integer :: first_side,last_side,side
        real(real64) :: width,total,side_contribution

        average=0.0_real64
        ierr=sample_matrix_success
        width=gap_hi-gap_lo
        if(width.le.0.0_real64 .or. .not.allocated(xarr) .or. &
           .not.allocated(amat_arr) .or. npoi.le.0) then
            ierr=sample_matrix_contribution_unresolved
            return
        endif

        if(gap_lo.eq.xbeg) then
            first_side=1
            last_side=1
        elseif(gap_hi.eq.xend) then
            first_side=2
            last_side=2
        else
            first_side=1
            last_side=2
        endif

        total=0.0_real64
        do side=first_side,last_side
            call estimate_side_from_samples(gap_lo,gap_hi,side, &
                side_contribution,ierr)
            if(ierr.ne.sample_matrix_success) then
                print *,'topology gap side failed: gap,side,npoi,xbeg,xend = ', &
                    gap_lo,gap_hi,side,npoi,xbeg,xend
                return
            endif
            total=total+side_contribution
        enddo

        average=total/width
        if(.not.ieee_is_finite(average) .or. average.lt.0.0_real64) then
            ierr=sample_matrix_contribution_unresolved
        endif
    end subroutine estimate_topology_gap_from_samples

    subroutine estimate_side_from_samples(gap_lo,gap_hi,side,contribution,ierr)
        real(real64), intent(in) :: gap_lo,gap_hi
        integer, intent(in) :: side
        real(real64), intent(out) :: contribution
        integer, intent(out) :: ierr

        integer :: j,k,count,fit_count,side_signature
        real(real64) :: boundary,direction,width,coordinate,jacobian
        real(real64) :: distance,integrand,coefficient,limit_coefficient
        real(real64), allocatable :: distances(:),coefficients(:)
        integer, allocatable :: signatures(:)
        logical :: belongs

        contribution=0.0_real64
        ierr=sample_matrix_success
        call side_geometry(gap_lo,gap_hi,side,boundary,direction,width)
        call potato_gap_square_map(boundary,direction,width,1.0_real64, &
            coordinate,jacobian)
        if(.not.ieee_is_finite(coordinate) .or. &
           .not.ieee_is_finite(jacobian) .or. jacobian.le.0.0_real64 .or. &
           width.le.0.0_real64) then
            ierr=sample_matrix_contribution_unresolved
            return
        endif

        allocate(distances(npoi),coefficients(npoi),signatures(npoi))
        count=0
        do j=1,npoi
            if(side.eq.1) then
                if(gap_lo.eq.xbeg) then
                    belongs=xarr(j).ge.gap_hi
                else
                    belongs=xarr(j).le.gap_lo
                endif
            else
                if(gap_hi.eq.xend) then
                    belongs=xarr(j).le.gap_lo
                else
                    belongs=xarr(j).ge.gap_hi
                endif
            endif
            if(.not.belongs) cycle
            distance=abs(xarr(j)-boundary)
            if(distance.le.0.0_real64) cycle
            integrand=sum(abs(amat_arr(:,:,j)))
            if(.not.ieee_is_finite(integrand) .or. integrand.lt.0.0_real64) then
                ierr=sample_matrix_contribution_unresolved
                deallocate(distances,coefficients,signatures)
                return
            endif
            call potato_gap_sqrt_coefficient(integrand,distance,coefficient)
            if(.not.ieee_is_finite(coefficient) .or. coefficient.lt.0.0_real64) then
                ierr=sample_matrix_contribution_unresolved
                deallocate(distances,coefficients,signatures)
                return
            endif
            count=count+1
            distances(count)=distance
            coefficients(count)=coefficient
            signatures(count)=topology_arr(j)
        enddo
        if(count.lt.minimum_gap_fit_points) then
            ierr=sample_matrix_contribution_unresolved
            deallocate(distances,coefficients,signatures)
            return
        endif

        do j=2,count
            distance=distances(j)
            coefficient=coefficients(j)
            side_signature=signatures(j)
            k=j
            do while(k.gt.1)
                if(distances(k-1).le.distance) exit
                distances(k)=distances(k-1)
                coefficients(k)=coefficients(k-1)
                signatures(k)=signatures(k-1)
                k=k-1
            enddo
            distances(k)=distance
            coefficients(k)=coefficient
            signatures(k)=side_signature
        enddo

        side_signature=signatures(1)
        fit_count=0
        do j=1,count
            if(signatures(j).ne.side_signature) exit
            fit_count=fit_count+1
            if(fit_count.ge.topology_gap_fit_points) exit
        enddo
        if(fit_count.lt.minimum_gap_fit_points) then
            ierr=sample_matrix_contribution_unresolved
            deallocate(distances,coefficients,signatures)
            return
        endif
        limit_coefficient=maxval(coefficients(1:fit_count))
        call potato_gap_sqrt_integral(limit_coefficient,width,contribution)
        if(.not.ieee_is_finite(contribution) .or. contribution.lt.0.0_real64) then
            ierr=sample_matrix_contribution_unresolved
        endif
        deallocate(distances,coefficients,signatures)
    end subroutine estimate_side_from_samples

    subroutine side_geometry(gap_lo,gap_hi,side,boundary,direction,gap_width)
        real(real64), intent(in) :: gap_lo,gap_hi
        integer, intent(in) :: side
        real(real64), intent(out) :: boundary,direction,gap_width

        if(side.eq.1) then
            if(gap_lo.eq.xbeg) then
                boundary=gap_lo
                direction=1.0_real64
                gap_width=gap_hi-gap_lo
            else
                boundary=0.5_real64*(gap_lo+gap_hi)
                direction=-1.0_real64
                gap_width=boundary-gap_lo
            endif
        else
            if(gap_hi.eq.xend) then
                boundary=gap_hi
                direction=-1.0_real64
                gap_width=gap_hi-gap_lo
            else
                boundary=0.5_real64*(gap_lo+gap_hi)
                direction=1.0_real64
                gap_width=gap_hi-boundary
            endif
        endif
    end subroutine side_geometry

end module topology_gap_quadrature_mod
