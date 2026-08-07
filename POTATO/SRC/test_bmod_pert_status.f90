program test_bmod_pert_status
    use, intrinsic :: ieee_arithmetic, only : ieee_is_finite
    use bmod_pert_mod, only : nrad,nzet,rad,zet,bmod_pert_success, &
        bmod_pert_outside_grid
    implicit none

    integer :: ierr
    double precision :: rinside,zinside,rlo,rhi,zlo,zhi
    complex(8) :: strict_value,legacy_value
    external :: bmod_pert_status,bmod_pert

    ! The test working directory contains the public circular bmod_n.dat
    ! fixture.  The first strict call also exercises one-time spline loading.
    call bmod_pert_status(0.d0,0.d0,strict_value,ierr)
    call require(ierr.eq.bmod_pert_outside_grid, &
        'strict perturbation evaluator accepted an outside point')
    call require(nrad.ge.2 .and. nzet.ge.2,'perturbation grid was not loaded')

    rlo=min(rad(1),rad(nrad))
    rhi=max(rad(1),rad(nrad))
    zlo=min(zet(1),zet(nzet))
    zhi=max(zet(1),zet(nzet))
    rinside=0.5d0*(rlo+rhi)
    zinside=0.5d0*(zlo+zhi)
    call bmod_pert_status(rinside,zinside,strict_value,ierr)
    call require(ierr.eq.bmod_pert_success, &
        'strict perturbation evaluator rejected an interior grid point')
    call require(ieee_is_finite(real(strict_value)) .and. &
        ieee_is_finite(aimag(strict_value)), &
        'strict perturbation evaluator returned a nonfinite interior value')

    call bmod_pert_status(rlo,zlo,strict_value,ierr)
    call require(ierr.eq.bmod_pert_success, &
        'strict evaluator rejected the lower grid edge')
    call bmod_pert_status(rhi,zhi,strict_value,ierr)
    call require(ierr.eq.bmod_pert_success, &
        'strict evaluator rejected the upper grid edge')

    call bmod_pert_status(rlo-1.d-8,zinside,strict_value,ierr)
    call require(ierr.eq.bmod_pert_outside_grid, &
        'strict evaluator accepted a point below the R grid')
    call bmod_pert_status(rhi+1.d-8,zinside,strict_value,ierr)
    call require(ierr.eq.bmod_pert_outside_grid, &
        'strict evaluator accepted a point above the R grid')
    call bmod_pert_status(rinside,zlo-1.d-8,strict_value,ierr)
    call require(ierr.eq.bmod_pert_outside_grid, &
        'strict evaluator accepted a point below the Z grid')
    call bmod_pert_status(rinside,zhi+1.d-8,strict_value,ierr)
    call require(ierr.eq.bmod_pert_outside_grid, &
        'strict evaluator accepted a point above the Z grid')

    ! The compatibility wrapper still clamps and therefore returns a finite
    ! legacy edge value for the same outside point.
    call bmod_pert(rhi+1.d0,zinside,legacy_value)
    call require(ieee_is_finite(real(legacy_value)) .and. &
        ieee_is_finite(aimag(legacy_value)), &
        'legacy perturbation wrapper no longer preserves clamped semantics')

    print *, 'bmod_pert strict grid-domain oracle: passed'

contains

    subroutine require(ok,message)
        logical, intent(in) :: ok
        character(*), intent(in) :: message
        if(.not.ok) error stop trim(message)
    end subroutine require

end program test_bmod_pert_status
