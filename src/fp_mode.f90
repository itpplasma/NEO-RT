module neort_fp_mode
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use neort_collision_rates, only: collision_elementary_charge
    implicit none
    private
    public :: tp_fp_collisions, fp_ev_cgs, fp_charge_over_c

    ! Dependency-light switch shared with profile setup; ordinary paths keep
    ! their historical constants and outputs when this mode is disabled.
    logical :: tp_fp_collisions = .false.
    real(dp), parameter :: fp_ev_cgs = collision_elementary_charge*1.0e7_dp
    real(dp), parameter :: fp_charge_over_c = collision_elementary_charge/10.0_dp
end module neort_fp_mode
