module neort_fp_controls
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_fp_mode, only: tp_fp_collisions
    implicit none
    private
    public :: tp_fp_collisions, fp_rate_scale, fp_speed_steps, fp_pitch_steps
    public :: fp_umax, reset_fp, validate_fp_controls, validate_fp_background

    real(dp) :: fp_rate_scale = 1.0_dp
    integer :: fp_speed_steps = 16, fp_pitch_steps = 16
    real(dp) :: fp_umax = 4.0_dp

contains

    subroutine reset_fp()
        tp_fp_collisions = .false.
        fp_rate_scale = 1.0_dp
        fp_speed_steps = 16
        fp_pitch_steps = 16
        fp_umax = 4.0_dp
    end subroutine reset_fp

    subroutine validate_fp_controls(require_background)
        use neort_species, only: mixture_enabled
        use neort_krook, only: krook_collisions
        use neort_lorentz_controls, only: lorentz_pitch_collisions
        use driftorbit, only: nonlin, drift_phase, fow_radial, supban, nopassing, &
            magdrift, magdrift_passing
        logical, intent(in), optional :: require_background
        logical :: require_mixture

        if (.not. ieee_is_finite(fp_rate_scale) .or. fp_rate_scale <= 0.0_dp) &
            error stop 'fp_rate_scale must be finite and strictly positive'
        if (.not. ieee_is_finite(fp_umax) .or. fp_umax <= 0.0_dp) &
            error stop 'fp_umax must be finite and positive'
        if (fp_speed_steps < 2 .or. fp_pitch_steps < 2) &
            error stop 'FP speed and pitch steps must be at least 2'
        if (.not. tp_fp_collisions) return
        if (krook_collisions .or. lorentz_pitch_collisions) &
            error stop 'select only one of TP FP, Lorentz or Krook collisions'
        if (nonlin .or. drift_phase .or. fow_radial .or. supban .or. nopassing) &
            error stop 'TP FP requires linear common-phase thin orbits and all branches'
        if (magdrift .and. magdrift_passing /= 1) &
            error stop 'TP FP requires the same drift model on all branches'
        require_mixture = .true.
        if (present(require_background)) require_mixture = require_background
        if (require_mixture .and. .not. mixture_enabled) &
            error stop 'test-particle FP requires an explicit species mixture'
    end subroutine validate_fp_controls

    subroutine validate_fp_background()
        use neort_species, only: mixture_at_s, selected_species
        if (.not. tp_fp_collisions) return
        if (.not. allocated(mixture_at_s)) &
            error stop 'evaluate the species mixture before TP FP'
        if (any(mixture_at_s%temperature_ev /= &
                mixture_at_s(selected_species)%temperature_ev)) &
            error stop 'test-particle FP currently requires exactly equal temperatures'
    end subroutine validate_fp_background
end module neort_fp_controls
