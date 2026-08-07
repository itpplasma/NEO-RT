module neort_gc_nonlocal_resonance_types
    !! Data contracts for the nonlocal finite-orbit-width resonance integral.
    !!
    !! The class coordinate x is a disconnected Poincare-orbit coordinate
    !! (normally R_c).  H0 and J_perp remain fixed arguments of the local
    !! class integral; their outer quadrature is deliberately owned by the
    !! transport layer.  No Boozer eta coordinate occurs in this contract.
    use, intrinsic :: iso_fortran_env, only: dp => real64

    implicit none
    private

    integer, parameter, public :: GC_NONLOCAL_MAX_FORCE_VALUES = 8
    integer, parameter, public :: GC_NONLOCAL_MAX_COMPONENTS = 256
    integer, parameter, public :: GC_NONLOCAL_MAX_ROOTS = 4096
    integer, parameter, public :: GC_NONLOCAL_MAX_SCAN_INTERVALS = 32768
    integer, parameter, public :: GC_NONLOCAL_MAX_ROOT_ITERATIONS = 1024

    integer, parameter, public :: GC_NONLOCAL_SUCCESS = 0
    integer, parameter, public :: GC_NONLOCAL_INVALID_INPUT = 1
    integer, parameter, public :: GC_NONLOCAL_PARTIAL = 2
    integer, parameter, public :: GC_NONLOCAL_CALLBACK_FAILURE = 3
    integer, parameter, public :: GC_NONLOCAL_ROOT_NOT_CONVERGED = 4
    integer, parameter, public :: GC_NONLOCAL_SINGULAR_RESONANCE = 5
    integer, parameter, public :: GC_NONLOCAL_COMPONENT_IDENTITY = 6
    integer, parameter, public :: GC_NONLOCAL_FORCE_CONTRACT = 7
    integer, parameter, public :: GC_NONLOCAL_NONFINITE = 8
    integer, parameter, public :: GC_NONLOCAL_CAPACITY = 9
    integer, parameter, public :: GC_NONLOCAL_DERIVATIVE_MISSING = 10

    integer, parameter, public :: GC_NONLOCAL_SAMPLE_VALID = 0
    integer, parameter, public :: GC_NONLOCAL_SAMPLE_WALL = 1
    integer, parameter, public :: GC_NONLOCAL_SAMPLE_UNRESOLVED = 2
    integer, parameter, public :: GC_NONLOCAL_SAMPLE_INVALID = 3

    type, public :: gc_nonlocal_resonance_options_t
        integer :: scan_intervals = 128
        integer :: max_root_iterations = 100
        integer :: max_roots = 256
        integer :: force_count = 1
        real(dp) :: residual_tolerance = 1.0e-11_dp
        real(dp) :: x_tolerance = 1.0e-11_dp
        real(dp) :: derivative_tolerance = 1.0e-13_dp
    end type gc_nonlocal_resonance_options_t

    type, public :: gc_nonlocal_component_t
        !! Identity is the tuple (sigma, component_id).  Components with
        !! equal sigma represent disconnected x intervals and must not
        !! overlap; opposite-sigma branches may share the same x interval.
        integer :: component_id = 0
        integer :: sigma = 0
        integer :: class_kind = 0
        real(dp) :: x_min = 0.0_dp
        real(dp) :: x_max = 0.0_dp
    end type gc_nonlocal_component_t

    type, public :: gc_nonlocal_orbit_sample_t
        integer :: status = GC_NONLOCAL_SAMPLE_INVALID
        integer :: component_id = 0
        integer :: sigma = 0
        integer :: class_kind = 0
        real(dp) :: psi_star = 0.0_dp
        real(dp) :: dpsi_star_dx = 0.0_dp
        real(dp) :: tau_b = 0.0_dp
        real(dp) :: omega_b = 0.0_dp
        real(dp) :: omega_phi = 0.0_dp
        real(dp) :: domega_b_dx = 0.0_dp
        real(dp) :: domega_phi_dx = 0.0_dp
        complex(dp) :: h_m = cmplx(0.0_dp, 0.0_dp, kind=dp)
        real(dp) :: thermodynamic_force(GC_NONLOCAL_MAX_FORCE_VALUES) = 0.0_dp
        integer :: nforce = 0
        logical :: derivatives_available = .false.
        logical :: class_behavior_certified = .false.
    end type gc_nonlocal_orbit_sample_t

    type, public :: gc_nonlocal_resonance_result_t
        !! When certified is false, contribution and component_contribution
        !! are diagnostic partial sums and are not accepted torque results.
        integer :: status = GC_NONLOCAL_INVALID_INPUT
        integer :: harmonic_m = 0
        integer :: harmonic_n = 0
        integer :: nforce = 0
        integer :: ncomponents = 0
        integer :: nroots = 0
        integer :: valid_samples = 0
        integer :: unresolved_samples = 0
        integer :: wall_samples = 0
        real(dp) :: h0 = 0.0_dp
        real(dp) :: jperp = 0.0_dp
        real(dp) :: contribution(GC_NONLOCAL_MAX_FORCE_VALUES) = 0.0_dp
        logical :: certified = .false.
        integer, allocatable :: component_id(:)
        integer, allocatable :: sigma(:)
        real(dp), allocatable :: component_contribution(:, :)
        real(dp), allocatable :: root_x(:)
        real(dp), allocatable :: root_residual_derivative(:)
        integer, allocatable :: root_component_id(:)
        integer, allocatable :: root_sigma(:)
        integer, allocatable :: root_class_kind(:)
        real(dp), allocatable :: root_force_contribution(:, :)
    end type gc_nonlocal_resonance_result_t

end module neort_gc_nonlocal_resonance_types
