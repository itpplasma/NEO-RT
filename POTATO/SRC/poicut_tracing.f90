subroutine poicut_z_of_R(r, z, err_dist)
    ! Newton iterations in Z onto grad(psi) x grad(B)=0. The starting guess
    ! must remain on the magnetic-axis branch. Return the converged field state.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    implicit none
    real(dp), intent(in) :: r, err_dist
    real(dp), intent(inout) :: z
    integer :: ierr

    call try_poicut_z_of_R(r, z, err_dist, ierr)
    if (ierr /= 0) error stop 'poicut_z_of_R: Newton iterations did not converge'
end subroutine poicut_z_of_R

subroutine try_poicut_z_of_R(r, z, err_dist, ierr)
    ! A continuation step can leave the field provider's valid rectangle
    ! before crossing a near-edge target. Report the failed trial for retry.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    real(dp), intent(in) :: r, err_dist
    real(dp), intent(inout) :: z
    integer, intent(out) :: ierr
    integer, parameter :: niter = 100
    integer :: i
    real(dp) :: dz, g, gr, gz

    ierr = 1
    do i = 1, niter
        call gpsigb_and_ders(r, z, g, gr, gz)
        if (.not. all(ieee_is_finite([g, gz]))) return
        if (gz == 0.0_dp) return
        dz = g/gz
        if (.not. ieee_is_finite(dz)) return
        z = z - dz
        if (abs(dz) < err_dist) then
            call gpsigb_and_ders(r, z, g, gr, gz)
            if (.not. all(ieee_is_finite([g, gz]))) return
            if (gz == 0.0_dp) return
            ierr = 0
            return
        end if
    end do
end subroutine try_poicut_z_of_R

subroutine trace_poicut_end(h_r, psils, sigpsi, err_dist, rend, zend)
    ! Trace grad(psi) x grad(B)=0 from the axis until crossing psils, then
    ! refine both conditions. Halve a failed continuation step so a valid
    ! target near the wall is not lost by overshooting into zero field.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_sub, only: psif, dpsidr, dpsidz
    use poicut_mod, only: rmagaxis, zmagaxis
    implicit none
    real(dp), intent(in) :: h_r, psils, sigpsi, err_dist
    real(dp), intent(out) :: rend, zend
    integer, parameter :: nsplit = 100, niter = 100
    integer :: i, ierr
    real(dp) :: r, z, det, dr, dz, stepfac, stepmod
    real(dp) :: step_r, trial_r, trial_z, g, gr, gz

    r = rmagaxis
    z = zmagaxis
    step_r = h_r
    do i = 1, nsplit
        trial_r = r + step_r
        trial_z = z
        call try_poicut_z_of_R(trial_r, trial_z, err_dist, ierr)
        if (ierr /= 0) then
            step_r = 0.5_dp*step_r
            if (abs(step_r) < err_dist) exit
            cycle
        end if
        r = trial_r
        z = trial_z
        if ((psif - psils)*sigpsi > 0.0_dp) exit
    end do
    if (i > nsplit .or. ierr /= 0) &
        error stop 'trace_poicut_end: flux surface rho_pol_max not '// &
        'reached - check convex wall coverage and equilibrium box size'

    do i = 1, niter
        call gpsigb_and_ders(r, z, g, gr, gz)
        det = dpsidr*gz - dpsidz*gr
        dr = ((psif - psils)*gz - dpsidz*g)/det
        dz = (dpsidr*g - (psif - psils)*gr)/det
        stepmod = sqrt(dr**2 + dz**2)
        stepfac = abs(h_r)/max(abs(h_r), stepmod)
        r = r - dr*stepfac
        z = z - dz*stepfac
        if (stepmod < err_dist) exit
    end do
    if (i > niter) &
        error stop 'trace_poicut_end: Newton iterations did not converge'
    call gpsigb_and_ders(r, z, g, gr, gz)
    rend = r
    zend = z
end subroutine trace_poicut_end
