program test_energy_selection
    use iso_fortran_env, only: dp => real64
    use kinetic_energy_selection, only: total_energy_sample, total_energy_measure
    use potato_input_mod, only: monoenergetic_x

    implicit none

    real(dp), parameter :: tolerance = 20.0_dp*epsilon(1.0_dp)
    real(dp) :: actual

    if (monoenergetic_x /= -1.0_dp) error stop "thermal mode must remain default"

    actual = total_energy_sample(-1.0_dp, 3, 4, 2.0_dp, 10.0_dp, 0.0_dp, 0.0_dp)
    if (abs(actual - 8.25_dp) > tolerance) &
        error stop "thermal total-energy midpoint changed"

    actual = total_energy_measure(-1.0_dp, 4, 10.0_dp, 2.0_dp)
    if (abs(actual - 2.5_dp) > tolerance) &
        error stop "thermal total-energy measure changed"

    actual = total_energy_sample(2.25_dp, 1, 60, 0.0_dp, 0.0_dp, 2.0_dp, 0.5_dp)
    if (abs(actual - 5.0_dp) > tolerance) &
        error stop "fixed x must select Phi_ref + x*T_ref"

    actual = total_energy_measure(2.25_dp, 60, 0.0_dp, 2.0_dp)
    if (abs(actual - 2.0_dp) > tolerance) &
        error stop "delta(x-x0) measure must transform as dH=T_ref*dx"
end program test_energy_selection
