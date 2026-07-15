program test_nonlinear_contract
    use iso_fortran_env, only: dp => real64
    use driftorbit, only: nonlin
    use neort_nonlin, only: nonlinear_attenuation
    use neort_orbit, only: nvar
    implicit none

    real(dp) :: bounceavg(nvar), theta, dnorm

    bounceavg = 0.0_dp
    nonlin = .false.
    theta = nonlinear_attenuation(1.0_dp, 0.5_dp, bounceavg, 1.0_dp, &
        0.0_dp, 0.0_dp, 1.0_dp, dnorm)

    if (theta /= 1.0_dp) error stop "nonlin=false must be unattenuated"
    if (dnorm /= -1.0_dp) error stop "D_nl must be marked unevaluated"

end program test_nonlinear_contract
