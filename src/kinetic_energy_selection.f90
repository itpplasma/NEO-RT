module kinetic_energy_selection
    use iso_fortran_env, only: dp => real64

    implicit none

    private
    public :: maxwellian_speed_density, speed_from_energy_ratio, validate_energy_ratio
    public :: energy_sample_count, energy_sample_speed, energy_sample_weight
    public :: total_energy_sample, total_energy_measure

    real(dp), parameter :: pi = acos(-1.0_dp)

contains

    pure real(dp) function maxwellian_speed_density(speed)
        real(dp), intent(in) :: speed

        maxwellian_speed_density = 4.0_dp/sqrt(pi)*speed**2*exp(-speed**2)
    end function maxwellian_speed_density

    pure real(dp) function speed_from_energy_ratio(energy_ratio)
        real(dp), intent(in) :: energy_ratio

        speed_from_energy_ratio = sqrt(energy_ratio)
    end function speed_from_energy_ratio

    subroutine validate_energy_ratio(energy_ratio)
        real(dp), intent(in) :: energy_ratio

        if (energy_ratio /= -1.0_dp .and. energy_ratio <= 0.0_dp) &
            error stop "monoenergetic_x must be -1 (thermal) or positive"
    end subroutine validate_energy_ratio

    pure integer function energy_sample_count(energy_ratio, thermal_steps)
        real(dp), intent(in) :: energy_ratio
        integer, intent(in) :: thermal_steps

        if (energy_ratio > 0.0_dp) then
            energy_sample_count = 1
        else
            energy_sample_count = thermal_steps
        end if
    end function energy_sample_count

    pure real(dp) function energy_sample_speed(energy_ratio, sample_index, &
            minimum_speed, maximum_speed, thermal_steps)
        real(dp), intent(in) :: energy_ratio, minimum_speed, maximum_speed
        integer, intent(in) :: sample_index, thermal_steps
        real(dp) :: step

        if (energy_ratio > 0.0_dp) then
            energy_sample_speed = speed_from_energy_ratio(energy_ratio)
        else
            step = (maximum_speed - minimum_speed)/real(thermal_steps, dp)
            energy_sample_speed = minimum_speed + (real(sample_index, dp) - 0.5_dp)*step
        end if
    end function energy_sample_speed

    pure real(dp) function energy_sample_weight(energy_ratio, speed, thermal_step)
        real(dp), intent(in) :: energy_ratio, speed, thermal_step

        if (energy_ratio > 0.0_dp) then
            energy_sample_weight = 1.0_dp
        else
            energy_sample_weight = thermal_step*maxwellian_speed_density(speed)
        end if
    end function energy_sample_weight

    pure real(dp) function total_energy_sample(energy_ratio, sample_index, thermal_steps, &
            minimum_energy, energy_range, reference_temperature, reference_potential)
        real(dp), intent(in) :: energy_ratio, minimum_energy, energy_range
        real(dp), intent(in) :: reference_temperature, reference_potential
        integer, intent(in) :: sample_index, thermal_steps

        if (energy_ratio > 0.0_dp) then
            total_energy_sample = reference_potential &
                + energy_ratio*reference_temperature
        else
            total_energy_sample = minimum_energy &
                + (real(sample_index, dp) - 0.5_dp)*energy_range/thermal_steps
        end if
    end function total_energy_sample

    pure real(dp) function total_energy_measure(energy_ratio, thermal_steps, &
            energy_range, reference_temperature)
        real(dp), intent(in) :: energy_ratio, energy_range, reference_temperature
        integer, intent(in) :: thermal_steps

        if (energy_ratio > 0.0_dp) then
            total_energy_measure = reference_temperature
        else
            total_energy_measure = energy_range/thermal_steps
        end if
    end function total_energy_measure

end module kinetic_energy_selection
