module neort_fp_controls
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_fp_mode, only: tp_fp_collisions, conserving_fp_collisions, &
        fp_physical_mode
    implicit none
    private
    public :: tp_fp_collisions, fp_rate_scale, fp_speed_steps, fp_pitch_steps
    public :: fp_umax, reset_fp, validate_fp_controls, validate_fp_background
    public :: conserving_fp_collisions, validate_conserving_controls
    public :: conserving_speed_steps, conserving_pitch_steps, conserving_umax
    public :: conserving_rate_scale
    public :: conserving_theta_order, conserving_velocity_order, conserving_cm_order
    public :: conserving_rel_order, conserving_polar_order, conserving_azimuth_order
    public :: conserving_rank_rtol, conserving_compatibility_rtol

    real(dp) :: fp_rate_scale = 1.0_dp
    integer :: fp_speed_steps = 16, fp_pitch_steps = 16
    real(dp) :: fp_umax = 4.0_dp
    integer :: conserving_speed_steps = 2, conserving_pitch_steps = 2
    real(dp) :: conserving_umax = 4.0_dp
    real(dp) :: conserving_rate_scale = 1.0_dp
    integer :: conserving_theta_order = 8, conserving_velocity_order = 4
    integer :: conserving_cm_order = 4, conserving_rel_order = 4
    integer :: conserving_polar_order = 4, conserving_azimuth_order = 8
    real(dp) :: conserving_rank_rtol = 1.0e-11_dp
    real(dp) :: conserving_compatibility_rtol = 1.0e-9_dp

contains

    subroutine reset_fp()
        tp_fp_collisions = .false.
        conserving_fp_collisions = .false.
        conserving_speed_steps = 2
        conserving_pitch_steps = 2
        conserving_umax = 4.0_dp
        conserving_rate_scale = 1.0_dp
        conserving_theta_order = 8
        conserving_velocity_order = 4
        conserving_cm_order = 4
        conserving_rel_order = 4
        conserving_polar_order = 4
        conserving_azimuth_order = 8
        conserving_rank_rtol = 1.0e-11_dp
        conserving_compatibility_rtol = 1.0e-9_dp
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
        call validate_conserving_controls(require_background)
        if (.not. tp_fp_collisions) return
        if (krook_collisions .or. lorentz_pitch_collisions .or. &
                conserving_fp_collisions) &
            error stop 'select only one of TP FP, conserving FP, '// &
                'Lorentz or Krook collisions'
        if (nonlin .or. drift_phase .or. fow_radial .or. supban .or. nopassing) &
            error stop 'TP FP requires linear common-phase thin orbits and all branches'
        if (magdrift .and. magdrift_passing /= 1) &
            error stop 'TP FP requires the same drift model on all branches'
        require_mixture = .true.
        if (present(require_background)) require_mixture = require_background
        if (require_mixture .and. .not. mixture_enabled) &
            error stop 'test-particle FP requires an explicit species mixture'
    end subroutine validate_fp_controls

    subroutine validate_conserving_controls(require_background)
        use neort_species, only: mixture_enabled
        use neort_krook, only: krook_collisions
        use neort_lorentz_controls, only: lorentz_pitch_collisions
        use driftorbit, only: nonlin, drift_phase, fow_radial, supban, nopassing, &
            magdrift, magdrift_passing
        logical, intent(in), optional :: require_background
        logical :: require_mixture

        if (.not. ieee_is_finite(conserving_rate_scale) .or. &
                conserving_rate_scale <= 0.0_dp) &
            error stop 'conserving_rate_scale must be finite and strictly positive'
        if (conserving_speed_steps < 2 .or. conserving_pitch_steps < 2) &
            error stop 'conserving FP speed and pitch steps must be at least 2'
        if (.not. ieee_is_finite(conserving_umax) .or. conserving_umax <= 0.0_dp) &
            error stop 'conserving_umax must be finite and positive'
        if (min(conserving_theta_order, conserving_velocity_order, &
                conserving_cm_order, conserving_rel_order, conserving_polar_order, &
                conserving_azimuth_order) < 2) &
            error stop 'conserving FP quadrature orders must be at least 2'
        if (conserving_azimuth_order < 4) &
            error stop 'conserving FP pair azimuth order must be at least 4'
        if (.not. ieee_is_finite(conserving_rank_rtol)) &
            error stop 'conserving_rank_rtol must be finite and in (0,1)'
        if (conserving_rank_rtol <= 0.0_dp .or. conserving_rank_rtol >= 1.0_dp) &
            error stop 'conserving_rank_rtol must be finite and in (0,1)'
        if (.not. ieee_is_finite(conserving_compatibility_rtol)) &
            error stop 'conserving_compatibility_rtol must be finite and in (0,1)'
        if (conserving_compatibility_rtol <= 0.0_dp .or. &
                conserving_compatibility_rtol >= 1.0_dp) &
            error stop 'conserving_compatibility_rtol must be finite and in (0,1)'
        if (.not. conserving_fp_collisions) return
        if (tp_fp_collisions .or. krook_collisions .or. lorentz_pitch_collisions) &
            error stop 'select only one of TP FP, conserving FP, '// &
                'Lorentz or Krook collisions'
        if (nonlin .or. drift_phase .or. fow_radial .or. supban .or. nopassing) &
            error stop 'conserving FP requires linear common-phase '// &
                'thin orbits and all branches'
        if (magdrift .and. magdrift_passing /= 1) &
            error stop 'conserving FP requires the same drift model on all branches'
        require_mixture = .true.
        if (present(require_background)) require_mixture = require_background
        if (require_mixture .and. .not. mixture_enabled) &
            error stop 'conserving FP requires an explicit species mixture'
    end subroutine validate_conserving_controls

    subroutine validate_fp_background()
        use neort_species, only: mixture_at_s, selected_species
        if (.not. fp_physical_mode()) return
        if (.not. allocated(mixture_at_s)) &
            error stop 'evaluate the species mixture before FP'
        if (any(mixture_at_s%temperature_ev /= &
                mixture_at_s(selected_species)%temperature_ev)) &
            error stop 'FP currently requires exactly equal temperatures'
    end subroutine validate_fp_background
end module neort_fp_controls
