module topology_gap_quadrature_mod
    !! Evaluate a topology bracket with the generated square-root endpoint map.
    !!
    !! The physical matrix callback is sampled only at open-side points.  A
    !! coarse/fine pair is used because the transformed integrand is regular
    !! at a reflecting U-turn, while the original x-coordinate integrand can
    !! have an integrable inverse-square-root singularity.  The returned value
    !! is an average contribution compatible with the existing sampler ABI;
    !! the sampler multiplies it by the geometric gap width.
    use, intrinsic :: iso_fortran_env, only : real64
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use potato_input_mod, only : topology_gap_quadrature_coarse, &
        topology_gap_quadrature_fine
    use potato_symbolic_kernel_mod, only : potato_gap_square_map
    use sample_matrix_out_mod, only : n1,n2,x,amat,topology_error, &
        topology_signature,topology_probe_only,topology_contribution_probe, &
        sample_matrix_success,sample_matrix_contribution_unresolved,xbeg,xend
    implicit none
    private
    public :: integrate_topology_gap

    abstract interface
        subroutine matrix_callback()
        end subroutine matrix_callback
    end interface

contains

    subroutine integrate_topology_gap(get_matrix,gap_lo,gap_hi,average,ierr)
        procedure(matrix_callback) :: get_matrix
        real(real64), intent(in) :: gap_lo,gap_hi
        real(real64), intent(out) :: average
        integer, intent(out) :: ierr

        real(real64) :: coarse,fine,contribution,width
        logical :: old_probe_only,old_contribution_probe

        average=0.0_real64
        ierr=sample_matrix_success
        width=gap_hi-gap_lo
        if(.not.ieee_is_finite(gap_lo) .or. .not.ieee_is_finite(gap_hi) .or. &
           width.le.0.0_real64) then
            ierr=sample_matrix_contribution_unresolved
            return
        endif

        if(.not.allocated(amat)) allocate(amat(n1,n2))
        if(size(amat,1).ne.n1 .or. size(amat,2).ne.n2) then
            deallocate(amat)
            allocate(amat(n1,n2))
        endif

        old_probe_only=topology_probe_only
        old_contribution_probe=topology_contribution_probe
        topology_probe_only=.false.
        topology_contribution_probe=.true.

        call integrate_with_nodes(get_matrix,gap_lo,gap_hi, &
            topology_gap_quadrature_coarse,coarse,ierr)
        if(ierr.eq.sample_matrix_success) then
            call integrate_with_nodes(get_matrix,gap_lo,gap_hi, &
                topology_gap_quadrature_fine,fine,ierr)
        endif

        topology_probe_only=old_probe_only
        topology_contribution_probe=old_contribution_probe
        if(ierr.ne.sample_matrix_success) return

        ! The refinement difference is retained as a conservative numerical
        ! enclosure of the unresolved bracket contribution.  It is deliberately
        ! added before division so the callback remains compatible with the
        ! sampler's existing average-times-width accounting.
        contribution=max(coarse,fine)+abs(fine-coarse)
        if(.not.ieee_is_finite(contribution) .or. contribution.lt.0.0_real64) then
            ierr=sample_matrix_contribution_unresolved
            return
        endif
        average=contribution/width
        if(.not.ieee_is_finite(average)) then
            ierr=sample_matrix_contribution_unresolved
        endif
    end subroutine integrate_topology_gap

    subroutine integrate_with_nodes(get_matrix,gap_lo,gap_hi,node_count,integral,ierr)
        procedure(matrix_callback) :: get_matrix
        real(real64), intent(in) :: gap_lo,gap_hi
        integer, intent(in) :: node_count
        real(real64), intent(out) :: integral
        integer, intent(out) :: ierr

        integer :: node,side,first_side,last_side,side_signature
        real(real64) :: boundary,direction,gap_width,parameter
        real(real64) :: coordinate,jacobian,integrand

        integral=0.0_real64
        ierr=sample_matrix_success
        if(node_count.le.0) then
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

        do side=first_side,last_side
            call side_geometry(gap_lo,gap_hi,side,boundary,direction,gap_width)
            side_signature=0
            do node=1,node_count
                parameter=(real(node,real64)-0.5_real64)/real(node_count,real64)
                call potato_gap_square_map(boundary,direction,gap_width,parameter, &
                    coordinate,jacobian)
                x=coordinate
                call get_matrix()
                if(topology_error.ne.0 .or. topology_signature.eq.0) then
                    ierr=sample_matrix_contribution_unresolved
                    return
                endif
                if(side_signature.eq.0) then
                    side_signature=topology_signature
                elseif(topology_signature.ne.side_signature) then
                    ierr=sample_matrix_contribution_unresolved
                    return
                endif
                integrand=sum(abs(amat))
                if(.not.ieee_is_finite(integrand) .or. integrand.lt.0.0_real64 .or. &
                   .not.ieee_is_finite(jacobian) .or. jacobian.lt.0.0_real64) then
                    ierr=sample_matrix_contribution_unresolved
                    return
                endif
                integral=integral+integrand*jacobian/real(node_count,real64)
            enddo
        enddo
    end subroutine integrate_with_nodes

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
