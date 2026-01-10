! Resonant transport regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, since 2015

module driftorbit
    use iso_fortran_env, only: dp => real64
    use util
    use do_magfie_mod
    use do_magfie_pert_mod, only: do_magfie_pert_amp, mph
    use spline
    use neort_profiles, only: vth, M_t, ni1, Om_tE, A1, A2, dM_tds, dni1ds, dOm_tEds
    use collis_alp, only: coleff

    implicit none

    ! Normalization factor for radial electric field
    real(dp) :: efac = 1.0_dp

    ! Harmonics TODO: make dynamic, multiple harmonics
    ! Default values are overridden by config file in driftorbit_test:read_control
    real(dp) :: epsmn = 1.0_dp            ! perturbation amplitude B1/B0
    integer :: m0 = 1                 ! Boozer poloidal perturbation mode
    integer :: mth = 1                ! canonical poloidal mode
    logical :: magdrift = .true.      ! consider magnetic drift
    logical :: nopassing = .false.    ! neglect passing particles
    logical :: pertfile = .false.     ! read perturbation from file with neo_magfie_pert
    logical :: comptorque = .true.    ! compute torque

    ! Flux surface TODO: make a dynamic, multiple flux surfaces support
    real(dp) :: dVds = 0.0_dp, etadt = 0.0_dp, etatp = 0.0_dp
    real(dp) :: etamin = 0.0_dp, etamax = 0.0_dp

    ! TODO: better B0 calculation (average magnetic field on flux surface)
    real(dp) :: B0 = 0.0_dp
    real(dp) :: Bmin = 0.0_dp, Bmax = 0.0_dp

    real(dp), parameter :: epst_spl = 1.0e-6_dp, epsp_spl = 1.0e-6_dp   ! dist to tpb for spline
    real(dp), parameter :: epsst_spl = 1.0e-3_dp, epssp_spl = 1.0e-3_dp ! dist to deep for spline
    real(dp), parameter :: epst = 1.0e-8_dp, epsp = 1.0e-8_dp ! smallest eta distance to tp bound



    ! Number of levels for coarse root finding
    integer, parameter :: nlev = 100
    real(dp) :: sign_vpar = 1.0_dp, sign_vpar_htheta = 1.0_dp

    ! Nonlinear calculation switch
    logical :: nonlin = .false.

    ! Flux-surface dependent and working variables (computed/modified per-thread)
    !$omp threadprivate (mth, dVds, etadt, etatp, etamin, etamax)
    !$omp threadprivate (B0, Bmin, Bmax, sign_vpar, sign_vpar_htheta)

    ! Shared read-only configuration (NOT threadprivate): efac, epsmn, m0,
    ! magdrift, nopassing, pertfile, comptorque, nonlin

end module driftorbit
