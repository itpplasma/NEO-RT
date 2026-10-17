module neort_lorentz_controls
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_species, only: mixture_enabled
    use neort_mixture_collisions, only: selected_deflection_rate
    use neort_krook, only: krook_collisions
    implicit none
    private
    public :: lorentz_pitch_collisions, lorentz_nu, lorentz_rate_source
    public :: lorentz_rate_scale, lorentz_pitch_steps, reset_lorentz
    public :: validate_lorentz_controls, lorentz_deflection_rate
    public :: lorentz_susceptibility_valid

    logical :: lorentz_pitch_collisions = .false.
    real(dp) :: lorentz_nu = 0.0_dp
    integer :: lorentz_rate_source = 0
    real(dp) :: lorentz_rate_scale = 1.0_dp
    integer :: lorentz_pitch_steps = 32
    logical :: lorentz_susceptibility_valid = .true.
    !$omp threadprivate(lorentz_susceptibility_valid)

contains

    subroutine reset_lorentz()
        lorentz_pitch_collisions = .false.
        lorentz_nu = 0.0_dp
        lorentz_rate_source = 0
        lorentz_rate_scale = 1.0_dp
        lorentz_pitch_steps = 32
    end subroutine reset_lorentz

    subroutine validate_lorentz_controls(require_background)
        use driftorbit, only: nonlin, drift_phase, fow_radial, supban, nopassing, &
            magdrift, magdrift_passing
        logical, intent(in), optional :: require_background
        logical :: require_mixture

        require_mixture = .true.
        if (present(require_background)) require_mixture = require_background
        if (.not. ieee_is_finite(lorentz_nu) .or. lorentz_nu < 0.0_dp) &
            error stop 'lorentz_nu must be finite and nonnegative [s^-1]'
        if (.not. ieee_is_finite(lorentz_rate_scale) .or. lorentz_rate_scale < 0.0_dp) &
            error stop 'lorentz_rate_scale must be finite and nonnegative'
        if (lorentz_rate_source < 0 .or. lorentz_rate_source > 1) &
            error stop 'lorentz_rate_source must be 0 (supplied) or 1 (all-background)'
        if (lorentz_pitch_steps < 2) error stop 'lorentz_pitch_steps must be at least 2'
        if (.not. lorentz_pitch_collisions) return
        if (krook_collisions) &
            error stop 'select either Lorentz diffusion or Krook damping'
        if (require_mixture .and. lorentz_rate_source == 1 .and. &
                .not. mixture_enabled) &
            error stop 'all-background Lorentz rates require '// &
                'an explicit species mixture'
        if (nonlin .or. drift_phase .or. fow_radial .or. supban .or. nopassing) &
            error stop 'Lorentz common-phase response requires '// &
                'linear thin orbits and all branches'
        if (magdrift .and. magdrift_passing /= 1) &
            error stop 'Lorentz common-phase response needs '// &
                'the same drift model on all branches'
    end subroutine validate_lorentz_controls

    real(dp) function lorentz_deflection_rate(speed_cm_s) result(rate)
        real(dp), intent(in) :: speed_cm_s
        if (lorentz_rate_source == 0) then
            rate = lorentz_rate_scale*lorentz_nu
        else
            rate = lorentz_rate_scale*selected_deflection_rate(speed_cm_s/100.0_dp)
        end if
    end function lorentz_deflection_rate

end module neort_lorentz_controls
