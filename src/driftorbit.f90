! Resonant transport regimes in tokamaks
! in the action-angle formalism
! Christopher Albert, since 2015

module driftorbit
    use util
    use do_magfie_mod
    use do_magfie_pert_mod, only: do_magfie_pert_amp, mph
    use spline
    use neort_profiles, only: vth, M_t, ni1, Om_tE, A1, A2, dM_tds, dni1ds, dOm_tEds
    use collis_alp, only: coleff

    implicit none

    ! Normalization factor for radial electric field
    real(8) :: efac = 1d0

    ! Harmonics TODO: make dynamic, multiple harmonics
    ! Default values are overridden by config file in driftorbit_test:read_control
    real(8) :: epsmn = 1d0            ! perturbation amplitude B1/B0
    integer :: m0 = 1                 ! Boozer poloidal perturbation mode
    integer :: mth = 1                ! canonical poloidal mode
    logical :: magdrift = .true.      ! consider magnetic drift
    logical :: nopassing = .false.    ! neglect passing particles
    logical :: pertfile = .false.     ! read perturbation from file with neo_magfie_pert
    logical :: comptorque = .true.    ! compute torque

    ! Flux surface TODO: make a dynamic, multiple flux surfaces support
    real(8) :: dVds, etadt, etatp
    real(8) :: etamin, etamax

    ! TODO: better B0 calculation (average magnetic field on flux surface)
    real(8) :: B0
    real(8) :: Bmin, Bmax

    real(8), parameter :: epst_spl = 1d-6, epsp_spl = 1d-6   ! dist to tpb for spline
    real(8), parameter :: epsst_spl = 1d-3, epssp_spl = 1d-3 ! dist to deep for spline
    real(8), parameter :: epst = 1d-8, epsp = 1d-8 ! smallest eta distance to tp bound



    ! Number of levels for coarse root finding
    integer, parameter :: nlev = 100
    real(8) :: sign_vpar = 1d0, sign_vpar_htheta = 1d0

    ! Nonlinear calculation switch
    logical :: nonlin = .false.

    !$omp threadprivate (efac, epsmn, m0, mth, magdrift, nopassing, pertfile, comptorque)
    !$omp threadprivate (dVds, etadt, etatp, etamin, etamax)
    !$omp threadprivate (B0, Bmin, Bmax, sign_vpar, sign_vpar_htheta, nonlin)

end module driftorbit
