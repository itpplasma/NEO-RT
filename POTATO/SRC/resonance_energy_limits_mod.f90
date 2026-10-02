module resonance_energy_limits_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    implicit none
    private
    public :: resonance_energy_limits
contains
    subroutine resonance_energy_limits(phi_min, phi_max, temperature, &
            kinetic_max, cutoff_over_temp, total_min, total_max, first_bin, valid)
        real(dp), intent(in) :: phi_min, phi_max, temperature, kinetic_max
        real(dp), intent(in) :: cutoff_over_temp
        real(dp), intent(out) :: total_min, total_max
        integer, intent(out) :: first_bin
        logical, intent(out) :: valid

        valid = .false.
        total_min = 0.0_dp
        total_max = 0.0_dp
        first_bin = 0
        if (.not. all(ieee_is_finite([phi_min, phi_max, temperature, &
            kinetic_max, cutoff_over_temp]))) return
        if (phi_max < phi_min) return
        if (temperature <= 0.0_dp) return
        if (kinetic_max <= 0.0_dp) return
        if (cutoff_over_temp < 0.0_dp) return
        if (phi_max > huge(1.0_dp) - kinetic_max) return
        total_max = phi_max + kinetic_max
        total_min = phi_min
        first_bin = 2
        if (cutoff_over_temp > 0.0_dp) then
            if (cutoff_over_temp > 1.0_dp) then
                if (temperature > huge(1.0_dp)/cutoff_over_temp) return
            end if
            total_min = cutoff_over_temp*temperature
            if (phi_max > huge(1.0_dp) - total_min) return
            total_min = phi_max + total_min
            first_bin = 1
        end if
        if (.not. all(ieee_is_finite([total_min, total_max]))) return
        if (total_min < 0.0_dp) then
            if (total_max > huge(1.0_dp) + total_min) return
        end if
        valid = total_max > total_min
    end subroutine resonance_energy_limits
end module resonance_energy_limits_mod
