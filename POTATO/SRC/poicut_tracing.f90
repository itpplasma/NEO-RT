module poicut_tracing_mod
    ! Tracing of the Poincare cut, the line grad(psi) x grad(B) = 0. The field
    ! enters as a procedure argument so the tracing can be checked against an
    ! analytic field without the equilibrium provider.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none

    abstract interface
        ! Field state at one point: poloidal flux psi and its first derivatives,
        ! and the phi component of grad(psi) x grad(B) with its derivatives.
        subroutine poicut_field(r, z, psi, dpsi_dr, dpsi_dz, g, gr, gz)
            import :: dp
            implicit none
            real(dp), intent(in) :: r, z
            real(dp), intent(out) :: psi, dpsi_dr, dpsi_dz, g, gr, gz
        end subroutine poicut_field
    end interface

contains

    subroutine poicut_z_of_R(r, z, err_dist, field, psi)
        ! Newton iterations in Z onto grad(psi) x grad(B)=0. The starting guess
        ! must remain on the magnetic-axis branch.
        implicit none
        real(dp), intent(in) :: r, err_dist
        real(dp), intent(inout) :: z
        procedure(poicut_field) :: field
        real(dp), intent(out), optional :: psi
        real(dp) :: psi_here
        integer :: ierr

        call try_poicut_z_of_R(r, z, err_dist, field, psi_here, ierr)
        if (ierr /= 0) error stop 'poicut_z_of_R: Newton iterations did not converge'
        if (present(psi)) psi = psi_here
    end subroutine poicut_z_of_R

    subroutine try_poicut_z_of_R(r, z, err_dist, field, psi, ierr)
        ! A continuation step can leave the field provider's valid rectangle
        ! before crossing a near-edge target. Report the failed trial for retry.
        implicit none
        real(dp), intent(in) :: r, err_dist
        real(dp), intent(inout) :: z
        procedure(poicut_field) :: field
        real(dp), intent(out) :: psi
        integer, intent(out) :: ierr
        integer, parameter :: niter = 100
        integer :: i
        real(dp) :: dz, dpsi_dr, dpsi_dz, g, gr, gz

        ierr = 1
        do i = 1, niter
            call field(r, z, psi, dpsi_dr, dpsi_dz, g, gr, gz)
            if (.not. all(ieee_is_finite([psi, g, gr, gz]))) return
            if (gz == 0.0_dp) return
            dz = g/gz
            if (.not. ieee_is_finite(dz)) return
            z = z - dz
            if (abs(dz) < err_dist) then
                call field(r, z, psi, dpsi_dr, dpsi_dz, g, gr, gz)
                if (.not. all(ieee_is_finite([psi, g, gr, gz]))) return
                if (gz == 0.0_dp) return
                ierr = 0
                return
            end if
        end do
    end subroutine try_poicut_z_of_R

    subroutine trace_poicut_end(r0, z0, h_r, psils, sigpsi, err_dist, field, rend, zend)
        ! Trace grad(psi) x grad(B)=0 from (r0,z0) until crossing psils, then
        ! refine both conditions. Halve a failed continuation step so a valid
        ! target near the wall is not lost by overshooting into zero field.
        implicit none
        real(dp), intent(in) :: r0, z0, h_r, psils, sigpsi, err_dist
        procedure(poicut_field) :: field
        real(dp), intent(out) :: rend, zend
        integer, parameter :: nsplit = 100, niter = 100
        integer :: i, ierr
        real(dp) :: r, z, psi, dpsi_dr, dpsi_dz, det, dr, dz, stepfac, stepmod
        real(dp) :: step_r, trial_r, trial_z, g, gr, gz

        r = r0
        z = z0
        step_r = h_r
        ierr = 0
        do i = 1, nsplit
            trial_r = r + step_r
            trial_z = z
            call try_poicut_z_of_R(trial_r, trial_z, err_dist, field, psi, ierr)
            if (ierr /= 0) then
                step_r = 0.5_dp*step_r
                if (abs(step_r) < err_dist) exit
                cycle
            end if
            r = trial_r
            z = trial_z
            if ((psi - psils)*sigpsi > 0.0_dp) exit
        end do
        if (i > nsplit .or. ierr /= 0) &
            error stop 'trace_poicut_end: flux surface rho_pol_max not '// &
            'reached - check convex wall coverage and equilibrium box size'

        do i = 1, niter
            call field(r, z, psi, dpsi_dr, dpsi_dz, g, gr, gz)
            det = dpsi_dr*gz - dpsi_dz*gr
            dr = ((psi - psils)*gz - dpsi_dz*g)/det
            dz = (dpsi_dr*g - (psi - psils)*gr)/det
            stepmod = sqrt(dr**2 + dz**2)
            stepfac = abs(h_r)/max(abs(h_r), stepmod)
            r = r - dr*stepfac
            z = z - dz*stepfac
            if (stepmod < err_dist) exit
        end do
        if (i > niter) &
            error stop 'trace_poicut_end: Newton iterations did not converge'
        call field(r, z, psi, dpsi_dr, dpsi_dz, g, gr, gz)
        rend = r
        zend = z
    end subroutine trace_poicut_end

end module poicut_tracing_mod
