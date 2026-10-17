module neort_mixture_collisions
    use iso_fortran_env, only: dp => real64
    use ieee_arithmetic, only: ieee_is_finite
    use neort_species, only: mixture_enabled, species_mixture, mixture_at_s, selected_species
    use neort_collision_rates, only: collision_rates_t, initialize_collision_rates, &
        deflection_frequency, reference_frequency
    use neort_krook, only: krook_collisions, krook_nu, krook_model, krook_rate
    implicit none
    private
    public :: krook_rate_source, krook_rate_scale, validate_collision_background
    public :: initialize_collision_background, finite_krook_active, effective_krook_rate
    public :: selected_thermal_deflection_rate, selected_reference_frequency
    public :: selected_deflection_rate

    ! 0 supplied legacy width; 1 all-background thermal deflection amplitude;
    ! 2 all-background deflection rate at the actual test-particle speed.
    ! These remain diagonal Krook models, not Lorentz or coupled FP operators.
    integer :: krook_rate_source = 0
    real(dp) :: krook_rate_scale = 1.0_dp
    type(collision_rates_t) :: background
    logical :: background_ready = .false.
    !$omp threadprivate(background, background_ready)

contains

    subroutine validate_collision_background(require_mixture)
        logical, intent(in), optional :: require_mixture
        if (krook_rate_source < 0 .or. krook_rate_source > 2) &
            error stop 'krook_rate_source must be 0, 1 or 2'
        if (.not. ieee_is_finite(krook_rate_scale) .or. krook_rate_scale < 0.0_dp) &
            error stop 'krook_rate_scale must be finite and nonnegative'
        if (present(require_mixture)) then
            if (require_mixture .and. krook_rate_source /= 0 .and. .not. mixture_enabled) &
                error stop 'all-background rates require an explicit species mixture'
        end if
        if (krook_rate_source == 2 .and. krook_model /= 0) &
            error stop 'velocity-dependent deflection rates require krook_model=0'
    end subroutine validate_collision_background

    subroutine initialize_collision_background()
        background_ready = .false.
        call validate_collision_background(.true.)
        if (.not. mixture_enabled) return
        if (.not. allocated(mixture_at_s)) error stop 'evaluate mixture before collision setup'
        call initialize_collision_rates(mixture_at_s%mass_kg, mixture_at_s%charge_e, &
            mixture_at_s%density_m3, mixture_at_s%temperature_ev, &
            species_mixture%coulomb_log, background)
        background_ready = .true.
    end subroutine initialize_collision_background

    logical function finite_krook_active()
        finite_krook_active = .false.
        if (.not. krook_collisions) return
        if (krook_rate_source == 0) then
            finite_krook_active = krook_nu > 0.0_dp
        else
            finite_krook_active = krook_rate_scale > 0.0_dp
        end if
    end function finite_krook_active

    real(dp) function selected_thermal_deflection_rate() result(rate)
        if (.not. background_ready) error stop 'collision background is not initialized'
        rate = deflection_frequency(background, selected_species, &
            background%vth_m_s(selected_species))
    end function selected_thermal_deflection_rate

    real(dp) function selected_reference_frequency() result(rate)
        if (.not. background_ready) error stop 'collision background is not initialized'
        rate = reference_frequency(background, selected_species)
    end function selected_reference_frequency

    real(dp) function selected_deflection_rate(speed_m_s) result(rate)
        real(dp), intent(in) :: speed_m_s
        if (.not. background_ready) error stop 'collision background is not initialized'
        rate = deflection_frequency(background, selected_species, speed_m_s)
    end function selected_deflection_rate

    real(dp) function effective_krook_rate(speed_over_vth, harmonic, geometry, speed_m_s) result(rate)
        real(dp), intent(in) :: speed_over_vth, harmonic, geometry
        real(dp), intent(in), optional :: speed_m_s
        real(dp) :: base, physical_speed
        select case (krook_rate_source)
        case (0)
            rate = krook_rate(krook_nu, speed_over_vth, harmonic, krook_model, geometry)
        case (1)
            base = krook_rate_scale*selected_thermal_deflection_rate()
            rate = krook_rate(base, speed_over_vth, harmonic, krook_model, geometry)
        case (2)
            if (.not. background_ready) error stop 'collision background is not initialized'
            physical_speed = speed_over_vth*background%vth_m_s(selected_species)
            if (present(speed_m_s)) physical_speed = speed_m_s
            rate = krook_rate_scale*geometry*deflection_frequency(background, &
                selected_species, physical_speed)
        case default
            error stop 'invalid collision rate source'
        end select
    end function effective_krook_rate
end module neort_mixture_collisions
