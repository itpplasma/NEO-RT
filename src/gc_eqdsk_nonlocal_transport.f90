module neort_gc_eqdsk_nonlocal_transport
    !! Standalone direct-EQDSK full-FOW transport factory.
    !!
    !! This is the production seam for a real-space nonlocal calculation.  It
    !! owns a concrete direct-EQDSK field, the native NEO-RT profiles, the
    !! direct R-Z perturbation, and the validated wall.  It instantiates the
    !! generic fixed-(H0,Jperp) cylindrical transport provider; it does not
    !! manufacture a field, a potential, a perturbation, or a frequency.
    !!
    !! The class coordinate is physical cylindrical R on the complete
    !! certified Eq. 13 cut from its HFS endpoint to its LFS endpoint.  No
    !! outboard-only branch or flux-coordinate hole is a valid production
    !! chart.
    !!
    !!     (grad B cross grad psi) dot grad phi = 0.
    !!
    !! The cut is located in direct R-Z space for every surface.  H0 and Jperp
    !! are held fixed.  A complete same-oriented return is integrated in
    !! (R,Z,phi), with Phi already in the Hamiltonian.  No Boozer angle or
    !! local flux-surface reduction is used in this module.
    !!
    !! The provider's single canonical conversion is
    !! psi_star=(c/q)*p_phi.  W_outer contains the explicit Eq. 17 profile
    !! factor and the shell residence average, but no tau_b, n^2, or canonical
    !! Jacobian.  The latter three factors are owned by the generic transport
    !! and resonance layers.
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use fortnum_ode_vode, only: vode_init, vode_integrate_to, vode_state_t
    use fortnum_status, only: FORTNUM_OK, fortnum_status_t
    use do_magfie_mod, only: R0, a, bfac, inp_swi, psi_pr, sample_eqdsk_field, s
    use field_eq_mod, only: btf, hfpol, nrad, nzet, psi_sep, rad, rtf, &
        splfpol, use_fpol, zet
    use neort_gc_callback_context, only: gc_callback_context_t
    use do_magfie_pert_mod, only: inp_swi_pert, MAGFIE_PERT_OK, &
        do_magfie_pert_amp_cylindrical, read_boozer_pert_file, rz_nrad, &
        rz_nzet, set_mph
    use geoflux_coordinates, only: geoflux_get_axis, &
        geoflux_get_flux_profiles
    use neort_gc_cylindrical_class_adapter, only: &
        GC_CYL_CLASS_CUT_ERROR, GC_CYL_CLASS_INTERIOR_INVALID, &
        GC_CYL_CLASS_SPLITTER_FAILURE, &
        GC_CYL_CLASS_SUCCESS, &
        gc_cylindrical_class_adapter_t, &
        gc_cylindrical_class_interval_t, gc_cylindrical_class_options_t, &
        gc_cylindrical_class_launch_t, gc_cylindrical_class_point_t, &
        gc_cylindrical_class_result_t, &
        enumerate_gc_cylindrical_classes, &
        evaluate_gc_cylindrical_class_point, &
        initialize_gc_cylindrical_class_adapter, &
        launch_gc_cylindrical_class
    use neort_gc_cylindrical_topology, only: &
        gc_cylindrical_allowed_region_set_t
    use neort_gc_cylindrical_dynamics, only: gc_cylindrical_rhs
    use neort_gc_cylindrical_model, only: &
        GC_CYL_EQUILIBRIUM_DOMAIN, GC_CYL_FIELD_ERROR, GC_CYL_INTEGRATOR_ERROR, &
        GC_CYL_INVARIANT_ERROR, GC_CYL_NO_RETURN, GC_CYL_PERTURBATION_ERROR, &
        GC_CYL_POTENTIAL_ERROR, GC_CYL_START_ERROR, GC_CYL_STATE_ERROR, &
        GC_CYL_SUCCESS, GC_CYL_WALL_ERROR, GC_CYL_WALL_LOSS, &
        gc_cylindrical_allowed_component_t, &
        gc_cylindrical_field_sample_t, gc_cylindrical_field_t, &
        gc_cylindrical_invariants_t, gc_cylindrical_potential_t, &
        gc_cylindrical_polygon_wall_t, gc_cylindrical_state_t, &
        gc_cylindrical_invariant_residuals, invariants_from_cylindrical_state
    use neort_gc_cylindrical_nonlocal_provider, only: &
        GC_CYL_NONLOCAL_CALLBACK_FAILURE, GC_CYL_NONLOCAL_ORBIT_ERROR_STATUS, &
        GC_CYL_NONLOCAL_ORBIT_UNRESOLVED, GC_CYL_NONLOCAL_ORBIT_VALID, &
        GC_CYL_NONLOCAL_ORBIT_WALL, GC_CYL_NONLOCAL_SUCCESS, &
        GC_CYL_NONLOCAL_WALL_CLEAR, GC_CYL_NONLOCAL_WALL_HIT, &
        gc_cylindrical_nonlocal_context_t, gc_cylindrical_nonlocal_orbit_t, &
        initialize_gc_cylindrical_nonlocal_provider
    use neort_gc_cylindrical_orbit, only: gc_cylindrical_orbit_options_t
    use neort_full_fow_cycle_frequency_symbolic, only: &
        evaluate_neort_full_fow_cycle_frequency
    use neort_full_fow_cycle_average_symbolic, only: &
        evaluate_neort_full_fow_cycle_average
    use neort_gc_cylindrical_physical_return, only: &
        GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR, &
        GC_CYL_PHYSICAL_EVENT_RADIAL_DOMAIN, GC_CYL_PHYSICAL_EVENT_RETURN, &
        GC_CYL_PHYSICAL_EVENT_WALL, gc_cylindrical_physical_return_options_t, &
        attach_gc_cylindrical_physical_return_certificate, &
        gc_cylindrical_physical_return_certificate_t, &
        gc_cylindrical_physical_return_t, &
        compute_gc_cylindrical_physical_return
    use neort_gc_cylindrical_transport_provider, only: &
        GC_CYL_TRANSPORT_SUCCESS, gc_cylindrical_transport_node_factory_i, &
        gc_cylindrical_transport_class_kind_i, &
        gc_cylindrical_transport_evidence_i, gc_cylindrical_transport_reset_i, &
        gc_cylindrical_transport_outer_factor_i, &
        gc_cylindrical_transport_quadrature_builder_i, &
        gc_cylindrical_transport_provider_t, &
        clear_gc_cylindrical_transport_provider, &
        initialize_gc_cylindrical_transport_provider
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        eqdsk_cylindrical_field_t, initialize_eqdsk_cylindrical_field, &
        map_eqdsk_flux_position
    use neort_gc_eqdsk_composite_cut_atlas, only: &
        EQDSK_COMPOSITE_ATLAS_SUCCESS, EQDSK_COMPOSITE_CUT_INBOARD, &
        EQDSK_COMPOSITE_CUT_OUTBOARD, build_eqdsk_composite_cut_atlas, &
        eqdsk_composite_cut_atlas_options_t, &
        eqdsk_composite_cut_atlas_t, get_eqdsk_composite_cut_radius_bounds, &
        map_eqdsk_composite_cut_atlas_radius, &
        map_eqdsk_composite_cut_atlas_rho
    use neort_gc_eqdsk_composite_r_ownership, only: &
        EQDSK_R_OWNERSHIP_SUCCESS, build_eqdsk_composite_r_partition, &
        eqdsk_composite_r_partition_t
    use neort_gc_eqdsk_allowed_region_cut_box, only: &
        EQDSK_CUT_BOX_SUCCESS, eqdsk_potential_profile_nodes_t, &
        validate_eqdsk_potential_profile_nodes
    use neort_gc_eqdsk_certified_allowed_provider, only: &
        GC_EQDSK_ALLOWED_PROVIDER_SUCCESS, &
        GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_ID, &
        build_gc_eqdsk_certified_allowed_regions, &
        gc_eqdsk_certified_allowed_provider_context_t, &
        initialize_gc_eqdsk_certified_allowed_provider, &
        verify_gc_eqdsk_certified_allowed_regions
    use neort_gc_eqdsk_cut_endpoint_certificate, only: &
        EQDSK_ENDPOINT_CERT_SUCCESS, build_eqdsk_cut_endpoint_certificate, &
        eqdsk_cut_endpoint_certificate_t
    use neort_gc_eqdsk_cut_jet, only: &
        EQDSK_CUT_JET_NONFINITE, EQDSK_CUT_JET_OUT_OF_DOMAIN, &
        EQDSK_CUT_JET_SUCCESS, eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
    use neort_gc_eqdsk_flux_profile_map, only: &
        EQDSK_FLUX_MAP_NONFINITE, EQDSK_FLUX_MAP_OUT_OF_RANGE, &
        EQDSK_FLUX_MAP_SUCCESS, clear_eqdsk_flux_profile_map, &
        eqdsk_flux_profile_map_t, initialize_eqdsk_flux_profile_map, &
        map_eqdsk_rho_tor_to_psihat, map_eqdsk_scaled_psi_to_s_tor
    use neort_gc_full_fow_normalization_runtime, only: &
        GC_FULL_FOW_NORMALIZATION_SUCCESS, &
        evaluate_gc_full_fow_canonical_flux, &
        evaluate_gc_full_fow_degree5_cell_enclosure, &
        evaluate_gc_full_fow_eq17, evaluate_gc_full_fow_jk_envelope, &
        evaluate_gc_full_fow_scaled_magnitude, &
        gc_full_fow_degree5_enclosure_t, gc_full_fow_energy_quadrature_t, &
        gc_full_fow_eq17_result_t, gc_full_fow_jk_envelope_t, &
        gc_full_fow_paired_quadrature_t, &
        gc_full_fow_phase_space_bound_certificate_t, &
        gc_full_fow_reference_scales_t, initialize_gc_full_fow_reference, &
        map_gc_full_fow_energy_quadrature, &
        map_gc_full_fow_paired_quadrature
    use neort_gc_unit_quadrature, only: GC_UNIT_QUADRATURE_SUCCESS, &
        build_gc_unit_gauss_legendre
    use neort_gc_nonlocal_resonance_types, only: &
        GC_NONLOCAL_MAX_FORCE_VALUES, GC_NONLOCAL_SAMPLE_VALID, &
        gc_nonlocal_component_t, gc_nonlocal_orbit_sample_t
    use neort_gc_nonlocal_transport_types, only: &
        GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES, &
        GC_NONLOCAL_CLASS_COUNTERPASSING, GC_NONLOCAL_CLASS_COPASSING, &
        GC_NONLOCAL_CLASS_TRAPPED, gc_nonlocal_transport_observed_evidence_t, &
        gc_nonlocal_transport_options_t, gc_nonlocal_transport_quadrature_t, &
        gc_nonlocal_transport_reference_t
    use neort_gc_perpendicular_invariant, only: &
        gc_buchholz_jk_from_mu_phys, gc_mu_phys_from_buchholz_jk, &
        gc_mu_phys_from_vperp
    use neort_profile_endpoint_symbolic, only: &
        evaluate_neort_profile_endpoints
    use neort_profile_potential_segment_symbolic, only: &
        evaluate_neort_profile_potential_segment
    use neort_profile_potential_map_symbolic, only: &
        evaluate_neort_profile_potential_map
    use neort_full_fow_harmonic_symbolic, only: &
        evaluate_neort_full_fow_harmonic_integrand
    use neort_profiles, only: A1, A2, M_t, Om_tE, Ti1, ni1, vth, &
        am1_global, Z1_global, init_plasma_at_s, init_profile_at_s, &
        init_thermodynamic_forces, prepare_plasma_splines, &
        prepare_profile_splines, read_plasma_input
    use neort_wall_io, only: WALL_IO_OK, load_wall_polygon, wall_polygon_t
    use util, only: c, eV, mu, pi, qe, readdata

    implicit none
    private

    integer, parameter, public :: GC_EQDSK_NONLOCAL_SUCCESS = 0
    integer, parameter, public :: GC_EQDSK_NONLOCAL_INVALID_INPUT = 1
    integer, parameter, public :: GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE = 1001
    integer, parameter, public :: GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE = 1002
    integer, parameter, public :: &
        GC_EQDSK_NONLOCAL_PERTURBATION_UNAVAILABLE = 1003
    integer, parameter, public :: GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE = 1004
    integer, parameter, public :: GC_EQDSK_NONLOCAL_WALL_UNAVAILABLE = 1005
    integer, parameter, public :: GC_EQDSK_NONLOCAL_DERIVATIVE_UNAVAILABLE = 1006
    integer, parameter, public :: GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED = 1007
    integer, parameter, public :: GC_EQDSK_NONLOCAL_SPECIES_MISMATCH = 1008
    integer, parameter, public :: GC_EQDSK_NONLOCAL_DOMAIN_ERROR = 1009
    integer, parameter, public :: GC_EQDSK_NONLOCAL_ORBIT_ERROR = 1010
    integer, parameter, public :: GC_EQDSK_NONLOCAL_NONFINITE = 1011
    integer, parameter, public :: GC_EQDSK_TWO_CUT_MULTIPLICITY_CERTIFICATE_ID = &
        130041

    integer, parameter :: FACTORY_FORCE_COUNT = 3
    integer, parameter, public :: GC_EQDSK_PROFILE_DENSITY = 1
    integer, parameter, public :: GC_EQDSK_PROFILE_TEMPERATURE = 2
    integer, parameter, public :: GC_EQDSK_PROFILE_PHI = 3
    integer, parameter, public :: GC_EQDSK_PROFILE_A1STAR = 4
    integer, parameter, public :: GC_EQDSK_PROFILE_A2STAR = 5
    integer, parameter, public :: GC_EQDSK_PROFILE_OMEGA_E = 6
    integer, parameter :: PROFILE_DENSITY = GC_EQDSK_PROFILE_DENSITY
    integer, parameter :: PROFILE_TEMPERATURE = GC_EQDSK_PROFILE_TEMPERATURE
    integer, parameter :: PROFILE_PHI = GC_EQDSK_PROFILE_PHI
    integer, parameter :: PROFILE_A1STAR = GC_EQDSK_PROFILE_A1STAR
    integer, parameter :: PROFILE_A2STAR = GC_EQDSK_PROFILE_A2STAR
    integer, parameter :: PROFILE_OMEGA_E = GC_EQDSK_PROFILE_OMEGA_E
    integer, parameter :: AUGMENTED_STATE_SIZE = 8
    integer, parameter :: MAX_CUT_ROOTS = 64

    type, public :: gc_eqdsk_nonlocal_species_t
        !! All quantities are in the direct NEO-RT CGS convention.
        character(len=32) :: name = ''
        real(dp) :: mass_g = 0.0_dp
        real(dp) :: charge_esu = 0.0_dp
        real(dp) :: reference_energy_erg = 0.0_dp
        real(dp) :: reference_velocity_cm_s = 0.0_dp
        character(len=16) :: mass_units = 'g'
        character(len=16) :: charge_units = 'statC'
        character(len=16) :: energy_units = 'erg'
        character(len=16) :: velocity_units = 'cm/s'
    end type gc_eqdsk_nonlocal_species_t

    type, public :: gc_eqdsk_nonlocal_options_t
        real(dp) :: field_scale = 1.0_dp
        real(dp) :: profile_electric_factor = 1.0_dp
        real(dp) :: profile_bfactor = 1.0_dp
        !! These are rho_tor bounds of the physical volume, not a profile
        !! interpolation cutoff.  The default includes the axis; profile
        !! inputs must cover the complete requested domain or initialization
        !! fails closed.  Every conversion to s_tor or psi_hat is explicit.
        real(dp) :: surface_min = 0.0_dp
        real(dp) :: surface_max = 1.0_dp
        real(dp) :: reference_surface = 0.50_dp
        integer :: cut_theta_points = 257
        real(dp) :: cut_root_tolerance = 1.0e-10_dp
        real(dp) :: cut_derivative_step = 1.0e-5_dp
        real(dp) :: orbit_derivative_fraction = 2.0e-4_dp
        real(dp) :: orbit_derivative_tolerance = 2.0e-3_dp
        real(dp) :: invariant_relative_tolerance = 3.0e-8_dp
        real(dp) :: orbit_maximum_step = 0.0_dp
        logical :: require_step_refinement = .true.
        real(dp) :: step_refinement_relative_tolerance = 1.0e-7_dp
        real(dp) :: step_refinement_absolute_tolerance = 1.0e-12_dp
        real(dp) :: topology_probe_fraction = 0.125_dp
        integer :: topology_probe_count = 5
        real(dp) :: endpoint_initial_grid_fraction = 0.20_dp
        integer :: endpoint_max_box_expansions = 7
        ! Negative selects launch rho_tor, converted to s_tor by the generated
        ! flux map, as the Eq. 17 target.  Otherwise this is an explicit s_tor
        ! threshold, never a geometric minor-radius approximation.
        real(dp) :: residence_target_s_tor = -1.0_dp
        !! Certified global lower bounds used by the paired (H0,J_K)
        !! domain.  A finite scan is retained only as a lower-bound witness;
        !! it is never promoted to a certificate.
        logical :: phase_space_bound_certified = .false.
        real(dp) :: qphi_min_certificate = 0.0_dp
        real(dp) :: bmod_min_certificate = 0.0_dp
        real(dp) :: phase_space_bound_absolute_tolerance = 1.0e-12_dp
        real(dp) :: phase_space_bound_relative_tolerance = 1.0e-10_dp
        character(len=32) :: phase_space_bound_method = ''
        type(gc_cylindrical_class_options_t) :: class_options = &
            gc_cylindrical_class_options_t()
        type(gc_cylindrical_orbit_options_t) :: orbit_options = &
            gc_cylindrical_orbit_options_t()
        type(eqdsk_composite_cut_atlas_options_t) :: cut_atlas_options
    end type gc_eqdsk_nonlocal_options_t

    type, extends(gc_cylindrical_potential_t), public :: &
            gc_eqdsk_profile_potential_t
        !! Phi is integrated from native Omega_E(s_tor) once.
        real(dp), allocatable :: s_tor(:)
        real(dp), allocatable :: psi_pol(:)
        real(dp), allocatable :: phi(:)
        real(dp), allocatable :: omega_e(:)
        real(dp) :: c_light = 0.0_dp
        logical :: initialized = .false.
        character(len=64) :: derivative_convention = &
            'dPhi/dpsi=+Omega_E/c'
        character(len=64) :: gauge_convention = 'Phi(psi_axis)=0'
        character(len=64) :: segment_integral_method = ''
        character(len=64) :: endpoint_reconstruction_method = ''
    contains
        procedure :: evaluate => evaluate_eqdsk_profile_potential
    end type gc_eqdsk_profile_potential_t

    type, public :: gc_eqdsk_orbit_result_t
        integer :: status = GC_EQDSK_NONLOCAL_ORBIT_ERROR
        real(dp) :: period = 0.0_dp
        real(dp) :: delta_phi = 0.0_dp
        real(dp) :: shell_average = 0.0_dp
        real(dp) :: omega_b = 0.0_dp
        real(dp) :: omega_phi = 0.0_dp
        complex(dp) :: h_m = cmplx(0.0_dp, 0.0_dp, kind=dp)
        logical :: wall_hit = .false.
        real(dp) :: energy_error = 0.0_dp
        real(dp) :: magnetic_moment_error = 0.0_dp
        real(dp) :: canonical_momentum_error = 0.0_dp
        integer :: return_orientation = 0
        integer :: event_kind = 0
        integer :: intersection_count = 0
        integer :: intersection_orientations(2) = 0
        real(dp) :: intersection_times(2) = 0.0_dp
        real(dp) :: intersection_rates(2) = 0.0_dp
        logical :: intersection_multiplicity_certified = .false.
        logical :: physical_return_found = .false.
        logical :: radial_domain_hit = .false.
        integer :: parallel_sign_changes = 0
        integer :: winding_number = 0
        logical :: complete_cycle_certified = .false.
        logical :: step_refinement_certified = .false.
        real(dp) :: base_step = 0.0_dp
        real(dp) :: refined_step = 0.0_dp
        real(dp) :: period_refinement_error = 0.0_dp
        real(dp) :: delta_phi_refinement_error = 0.0_dp
        real(dp) :: omega_b_refinement_error = 0.0_dp
        real(dp) :: omega_phi_refinement_error = 0.0_dp
        real(dp) :: h_m_refinement_error = 0.0_dp
        real(dp) :: shell_refinement_error = 0.0_dp
        type(gc_cylindrical_state_t) :: state_at_return
    end type gc_eqdsk_orbit_result_t

    type, public :: gc_eqdsk_cut_branch_t
        integer :: branch_id = 0
        real(dp) :: rc_min = 0.0_dp
        real(dp) :: rc_max = 0.0_dp
        real(dp) :: arclength_min = 0.0_dp
        real(dp) :: arclength_max = 0.0_dp
        integer :: orientation = 0
        logical :: continuation_certified = .false.
        logical :: root_isolation_certified = .false.
        character(len=24) :: lower_boundary_kind = 'unresolved'
        character(len=24) :: upper_boundary_kind = 'unresolved'
        character(len=32) :: limiting_chart = 'unresolved'
    end type gc_eqdsk_cut_branch_t

    type, public :: gc_eqdsk_cut_atlas_t
        integer :: nbranches = 0
        logical :: certified = .false.
        logical :: two_cut_multiplicity_certified = .false.
        character(len=64) :: method = ''
        type(gc_eqdsk_cut_branch_t), allocatable :: branches(:)
    end type gc_eqdsk_cut_atlas_t

    type, public, extends(gc_callback_context_t) :: gc_eqdsk_nonlocal_factory_t
        !! The concrete field is owned here.  It is never a pointer to a
        !! temporary local object, so cached class adapters retain a valid
        !! target for their whole factory lifetime.
        type(eqdsk_cylindrical_field_t) :: field
        type(gc_cylindrical_transport_provider_t) :: provider
        type(gc_nonlocal_transport_quadrature_t) :: quadrature
        type(gc_nonlocal_transport_options_t) :: transport_options
        type(gc_eqdsk_nonlocal_species_t) :: species
        type(gc_full_fow_reference_scales_t) :: normalization
        type(gc_eqdsk_nonlocal_options_t) :: options
        type(gc_eqdsk_profile_potential_t) :: potential
        type(gc_cylindrical_polygon_wall_t) :: wall
        type(gc_eqdsk_cut_atlas_t) :: cut_atlas
        type(eqdsk_composite_cut_atlas_t) :: certified_cut_atlas
        type(eqdsk_composite_r_partition_t) :: r_partition
        type(eqdsk_potential_profile_nodes_t) :: certified_profile
        type(gc_eqdsk_certified_allowed_provider_context_t) :: &
            allowed_region_context
        type(eqdsk_flux_profile_map_t) :: flux_map
        character(len=1024) :: eqdsk_path = ''
        character(len=1024) :: wall_path = ''
        character(len=1024) :: perturbation_path = ''
        character(len=1024) :: plasma_path = ''
        character(len=1024) :: profile_path = ''
        character(len=16) :: wall_units = ''
        character(len=16) :: wall_backend_units = 'cm'
        character(len=64) :: wall_hash = ''
        character(len=64) :: section_reference_id = ''
        character(len=32) :: section_coordinate = &
            'rho_tor:certified_Eq13_outboard'
        character(len=32) :: section_units = '1'
        real(dp) :: section_position(3) = 0.0_dp
        real(dp) :: section_flux = 0.0_dp
        integer :: section_orientation = 0
        real(dp) :: minor_radius_cm = 0.0_dp
        real(dp) :: psi_pr_effective = 0.0_dp
        !! Buchholz Eq. 17 Phi_eff=c*Eref/(|q|*v0), not a transport fit.
        real(dp) :: phi_eff_normalization = 0.0_dp
        real(dp) :: refinement_length_scale = 0.0_dp
        real(dp) :: refinement_time_scale = 0.0_dp
        real(dp) :: refinement_frequency_scale = 0.0_dp
        real(dp) :: refinement_momentum_scale = 0.0_dp
        real(dp) :: refinement_mu_scale = 0.0_dp
        character(len=64) :: perturbation_amplitude_convention = ''
        character(len=64) :: perturbation_coordinate_dependency = ''
        logical :: perturbation_interpolation_certified = .false.
        logical :: perturbation_provenance_certified = .false.
        real(dp), allocatable :: profile_s(:)
        real(dp), allocatable :: profile_psi(:)
        real(dp), allocatable :: profile_dpsi_ds(:)
        real(dp), allocatable :: profile_values(:, :)
        logical :: field_ready = .false.
        logical :: profile_ready = .false.
        logical :: perturbation_ready = .false.
        logical :: wall_ready = .false.
        logical :: cut_ready = .false.
        logical :: allowed_region_ready = .false.
        logical :: topology_ready = .false.
        logical :: cut_atlas_certified = .false.
        character(len=64) :: cut_atlas_method = ''
        character(len=64) :: profile_endpoint_reconstruction_method = &
            'fortsym-affine-endpoint-required'
        character(len=64) :: profile_potential_derivative_convention = &
            'dPhi/dpsi=+Omega_E/c'
        character(len=64) :: profile_potential_gauge = 'Phi(psi_axis)=0'
        integer :: cut_atlas_branch_count = 0
        integer :: active_cut_branch_id = 0
        character(len=64) :: flux_orientation_convention = &
            'B_R=-psi_Z/R; B_Z=psi_R/R'
        character(len=64) :: cut_parameter_orientation = ''
        character(len=64) :: canonical_one_form_convention = ''
        character(len=64) :: fourier_phase_convention = ''
        logical :: native_orientation_certified = .false.
        logical :: initialized = .false.
        integer :: physical_return_attempts = 0
        integer :: physical_return_successes = 0
        integer :: wall_return_count = 0
        integer :: radial_return_count = 0
        integer :: no_return_count = 0
        integer :: invariant_rejection_count = 0
        integer :: field_error_count = 0
        integer :: potential_error_count = 0
        integer :: state_error_count = 0
        integer :: start_error_count = 0
        integer :: integrator_error_count = 0
        integer :: other_return_error_count = 0
        integer :: invariant_success_count = 0
        integer :: invariant_failure_count = 0
        integer :: return_error_count = 0
        integer :: wall_not_hit_count = 0
        integer :: wall_error_count = 0
        integer :: last_return_status = GC_EQDSK_NONLOCAL_ORBIT_ERROR
        integer :: last_return_event_kind = 0
        real(dp) :: last_launch_event_value = 0.0_dp
        real(dp) :: last_return_event_value = 0.0_dp
        real(dp) :: last_invariant_scaled_drift = 0.0_dp
        real(dp) :: last_energy_error = 0.0_dp
        real(dp) :: last_magnetic_moment_error = 0.0_dp
        real(dp) :: last_canonical_momentum_error = 0.0_dp
        integer :: topology_certification_attempts = 0
        integer :: topology_certification_successes = 0
        integer :: two_cut_multiplicity_certificate_id = 0
        integer :: harmonic_average_successes = 0
        integer :: residence_average_successes = 0
        logical :: step_refinement_certified = .false.
        real(dp) :: base_step = 0.0_dp
        real(dp) :: refined_step = 0.0_dp
        real(dp) :: period_refinement_error = 0.0_dp
        real(dp) :: delta_phi_refinement_error = 0.0_dp
        real(dp) :: omega_b_refinement_error = 0.0_dp
        real(dp) :: omega_phi_refinement_error = 0.0_dp
        real(dp) :: h_m_refinement_error = 0.0_dp
        real(dp) :: shell_refinement_error = 0.0_dp
    end type gc_eqdsk_nonlocal_factory_t

    type, public :: gc_eqdsk_nonlocal_diagnostics_t
        integer :: physical_return_attempts = 0
        integer :: physical_return_successes = 0
        integer :: wall_returns = 0
        integer :: radial_domain_returns = 0
        integer :: no_returns = 0
        integer :: invariant_rejections = 0
        integer :: field_errors = 0
        integer :: potential_errors = 0
        integer :: state_errors = 0
        integer :: start_errors = 0
        integer :: integrator_errors = 0
        integer :: other_return_errors = 0
        integer :: invariant_successes = 0
        integer :: invariant_failures = 0
        integer :: return_errors = 0
        integer :: wall_not_hit = 0
        integer :: wall_loss = 0
        integer :: wall_errors = 0
        integer :: categorized_returns = 0
        logical :: return_accounting_complete = .false.
        integer :: last_return_status = GC_EQDSK_NONLOCAL_ORBIT_ERROR
        integer :: last_return_event_kind = 0
        real(dp) :: last_launch_event_value = 0.0_dp
        real(dp) :: last_return_event_value = 0.0_dp
        real(dp) :: last_invariant_scaled_drift = 0.0_dp
        real(dp) :: last_energy_error = 0.0_dp
        real(dp) :: last_magnetic_moment_error = 0.0_dp
        real(dp) :: last_canonical_momentum_error = 0.0_dp
        integer :: topology_certification_attempts = 0
        integer :: topology_certification_successes = 0
        integer :: harmonic_average_successes = 0
        integer :: residence_average_successes = 0
        logical :: step_refinement_certified = .false.
        real(dp) :: base_step = 0.0_dp
        real(dp) :: refined_step = 0.0_dp
        real(dp) :: period_refinement_error = 0.0_dp
        real(dp) :: delta_phi_refinement_error = 0.0_dp
        real(dp) :: omega_b_refinement_error = 0.0_dp
        real(dp) :: omega_phi_refinement_error = 0.0_dp
        real(dp) :: h_m_refinement_error = 0.0_dp
        real(dp) :: shell_refinement_error = 0.0_dp
        logical :: field_ready = .false.
        logical :: profile_ready = .false.
        logical :: perturbation_ready = .false.
        logical :: wall_ready = .false.
        logical :: cut_ready = .false.
        logical :: topology_ready = .false.
        logical :: initialized = .false.
    end type gc_eqdsk_nonlocal_diagnostics_t

    public :: initialize_gc_eqdsk_nonlocal_transport
    public :: clear_gc_eqdsk_nonlocal_transport
    public :: evaluate_gc_eqdsk_native_profile
    public :: gc_eqdsk_nonlocal_status_message
    public :: get_gc_eqdsk_nonlocal_diagnostics
    public :: poincare_cut_value_position

contains

    pure subroutine get_gc_eqdsk_nonlocal_diagnostics(factory, diagnostics)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        type(gc_eqdsk_nonlocal_diagnostics_t), intent(out) :: diagnostics

        diagnostics = gc_eqdsk_nonlocal_diagnostics_t()
        diagnostics%physical_return_attempts = factory%physical_return_attempts
        diagnostics%physical_return_successes = &
            factory%physical_return_successes
        diagnostics%wall_returns = factory%wall_return_count
        diagnostics%radial_domain_returns = factory%radial_return_count
        diagnostics%no_returns = factory%no_return_count
        diagnostics%invariant_rejections = factory%invariant_rejection_count
        diagnostics%field_errors = factory%field_error_count
        diagnostics%potential_errors = factory%potential_error_count
        diagnostics%state_errors = factory%state_error_count
        diagnostics%start_errors = factory%start_error_count
        diagnostics%integrator_errors = factory%integrator_error_count
        diagnostics%other_return_errors = factory%other_return_error_count
        diagnostics%invariant_successes = factory%invariant_success_count
        diagnostics%invariant_failures = factory%invariant_failure_count
        diagnostics%return_errors = factory%return_error_count
        diagnostics%wall_not_hit = factory%wall_not_hit_count
        diagnostics%wall_loss = factory%wall_return_count
        diagnostics%wall_errors = factory%wall_error_count
        diagnostics%categorized_returns = diagnostics%physical_return_successes &
            +diagnostics%wall_returns+diagnostics%radial_domain_returns &
            +diagnostics%no_returns+diagnostics%invariant_rejections &
            +diagnostics%field_errors+diagnostics%potential_errors &
            +diagnostics%state_errors+diagnostics%start_errors &
            +diagnostics%integrator_errors+diagnostics%other_return_errors
        diagnostics%return_accounting_complete = &
            diagnostics%invariant_successes+diagnostics%invariant_failures == &
            diagnostics%physical_return_attempts .and. &
            diagnostics%physical_return_successes+diagnostics%no_returns+ &
            diagnostics%radial_domain_returns+diagnostics%wall_returns+ &
            diagnostics%return_errors == diagnostics%physical_return_attempts .and. &
            diagnostics%wall_not_hit+diagnostics%wall_loss+diagnostics%wall_errors == &
            diagnostics%physical_return_attempts
        diagnostics%last_return_status = factory%last_return_status
        diagnostics%last_return_event_kind = factory%last_return_event_kind
        diagnostics%last_launch_event_value = factory%last_launch_event_value
        diagnostics%last_return_event_value = factory%last_return_event_value
        diagnostics%last_invariant_scaled_drift = &
            factory%last_invariant_scaled_drift
        diagnostics%last_energy_error = factory%last_energy_error
        diagnostics%last_magnetic_moment_error = &
            factory%last_magnetic_moment_error
        diagnostics%last_canonical_momentum_error = &
            factory%last_canonical_momentum_error
        diagnostics%topology_certification_attempts = &
            factory%topology_certification_attempts
        diagnostics%topology_certification_successes = &
            factory%topology_certification_successes
        diagnostics%harmonic_average_successes = &
            factory%harmonic_average_successes
        diagnostics%residence_average_successes = &
            factory%residence_average_successes
        diagnostics%step_refinement_certified = factory%step_refinement_certified
        diagnostics%base_step = factory%base_step
        diagnostics%refined_step = factory%refined_step
        diagnostics%period_refinement_error = factory%period_refinement_error
        diagnostics%delta_phi_refinement_error = factory%delta_phi_refinement_error
        diagnostics%omega_b_refinement_error = factory%omega_b_refinement_error
        diagnostics%omega_phi_refinement_error = factory%omega_phi_refinement_error
        diagnostics%h_m_refinement_error = factory%h_m_refinement_error
        diagnostics%shell_refinement_error = factory%shell_refinement_error
        diagnostics%field_ready = factory%field_ready
        diagnostics%profile_ready = factory%profile_ready
        diagnostics%perturbation_ready = factory%perturbation_ready
        diagnostics%wall_ready = factory%wall_ready
        diagnostics%cut_ready = factory%cut_ready
        diagnostics%topology_ready = factory%topology_ready
        diagnostics%initialized = factory%initialized
    end subroutine get_gc_eqdsk_nonlocal_diagnostics

    subroutine initialize_gc_eqdsk_nonlocal_transport(eqdsk_path, wall_path, &
            wall_units, perturbation_path, plasma_path, profile_path, species, &
            harmonic_n, factory, status, options, transport_options)
        character(len=*), intent(in) :: eqdsk_path, wall_path, wall_units
        character(len=*), intent(in) :: perturbation_path, plasma_path
        character(len=*), intent(in) :: profile_path
        type(gc_eqdsk_nonlocal_species_t), intent(in) :: species
        integer, intent(in) :: harmonic_n
        type(gc_eqdsk_nonlocal_factory_t), target, intent(out) :: factory
        integer, intent(out) :: status
        type(gc_eqdsk_nonlocal_options_t), intent(in), optional :: options
        type(gc_nonlocal_transport_options_t), intent(in), optional :: &
            transport_options

        type(wall_polygon_t) :: input_wall
        type(gc_nonlocal_transport_reference_t) :: reference
        character(len=128) :: wall_message
        integer :: local_status, wall_status

        factory = gc_eqdsk_nonlocal_factory_t()
        status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        if (present(options)) factory%options = options
        if (present(transport_options)) then
            factory%transport_options = transport_options
        end if
        factory%species = species
        factory%eqdsk_path = trim(eqdsk_path)
        factory%wall_path = trim(wall_path)
        factory%wall_units = trim(wall_units)
        factory%perturbation_path = trim(perturbation_path)
        factory%plasma_path = trim(plasma_path)
        factory%profile_path = trim(profile_path)

        local_status = validate_factory_inputs(factory, harmonic_n)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        if (.not. readable_file(factory%eqdsk_path)) then
            status = GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE
            return
        end if
        if (.not. readable_file(factory%wall_path)) then
            status = GC_EQDSK_NONLOCAL_WALL_UNAVAILABLE
            return
        end if
        if (.not. readable_file(factory%perturbation_path)) then
            status = GC_EQDSK_NONLOCAL_PERTURBATION_UNAVAILABLE
            return
        end if
        if (.not. readable_file(factory%plasma_path)) then
            status = GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE
            return
        end if
        if (.not. readable_file(factory%profile_path)) then
            status = GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE
            return
        end if

        call initialize_direct_field(factory, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        call load_direct_wall(factory, input_wall, wall_message, wall_status)
        if (wall_status /= WALL_IO_OK) then
            status = GC_EQDSK_NONLOCAL_WALL_UNAVAILABLE
            return
        end if
        call factory%wall%set_vertices(input_wall%vertices, local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            status = GC_EQDSK_NONLOCAL_WALL_UNAVAILABLE
            return
        end if
        factory%wall_units = input_wall%input_units
        factory%wall_hash = input_wall%hash
        factory%wall_ready = .true.

        call load_native_profiles(factory, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        call load_direct_perturbation(factory, harmonic_n, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        call build_poincare_reference(factory, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if

        call build_paired_phase_space_quadrature(factory, &
            factory%transport_options%h0_order, factory%transport_options%jk_order, &
            factory%quadrature, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if

        reference = make_transport_reference(factory)
        call initialize_gc_cylindrical_transport_provider(factory%provider, &
            factory%quadrature, reference, 0, harmonic_n, FACTORY_FORCE_COUNT, &
            factory_node_factory, factory_outer_factor, factory_class_kind, &
            factory_reset_evidence, factory_get_evidence, local_status, &
            user_data=factory, &
            section_reference_id=factory%section_reference_id, &
            quadrature_builder=factory_build_quadrature)
        if (local_status /= GC_CYL_TRANSPORT_SUCCESS) then
            status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
            return
        end if
        factory%initialized = .true.
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine initialize_gc_eqdsk_nonlocal_transport

    subroutine clear_gc_eqdsk_nonlocal_transport(factory)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory

        call clear_gc_cylindrical_transport_provider(factory%provider)
        factory%normalization = gc_full_fow_reference_scales_t()
        factory%quadrature = gc_nonlocal_transport_quadrature_t()
        factory%transport_options = gc_nonlocal_transport_options_t()
        if (allocated(factory%profile_s)) deallocate(factory%profile_s)
        if (allocated(factory%profile_psi)) deallocate(factory%profile_psi)
        if (allocated(factory%profile_dpsi_ds)) deallocate(factory%profile_dpsi_ds)
        if (allocated(factory%profile_values)) deallocate(factory%profile_values)
        call clear_eqdsk_flux_profile_map(factory%flux_map)
        call clear_profile_potential(factory%potential)
        if (allocated(factory%wall%vertices)) deallocate(factory%wall%vertices)
        if (allocated(factory%cut_atlas%branches)) deallocate(factory%cut_atlas%branches)
        factory%cut_atlas = gc_eqdsk_cut_atlas_t()
        factory%certified_cut_atlas = eqdsk_composite_cut_atlas_t()
        factory%r_partition = eqdsk_composite_r_partition_t()
        factory%certified_profile = eqdsk_potential_profile_nodes_t()
        factory%allowed_region_context = &
            gc_eqdsk_certified_allowed_provider_context_t()
        factory%wall%initialized = .false.
        factory%field = eqdsk_cylindrical_field_t()
        factory%field_ready = .false.
        factory%profile_ready = .false.
        factory%perturbation_ready = .false.
        factory%wall_ready = .false.
        factory%cut_ready = .false.
        factory%allowed_region_ready = .false.
        factory%topology_ready = .false.
        factory%cut_atlas_certified = .false.
        factory%cut_atlas_method = ''
        factory%cut_atlas_branch_count = 0
        factory%active_cut_branch_id = 0
        factory%cut_parameter_orientation = ''
        factory%canonical_one_form_convention = ''
        factory%fourier_phase_convention = ''
        factory%native_orientation_certified = .false.
        factory%initialized = .false.
        factory%physical_return_attempts = 0
        factory%physical_return_successes = 0
        factory%wall_return_count = 0
        factory%radial_return_count = 0
        factory%no_return_count = 0
        factory%invariant_rejection_count = 0
        factory%field_error_count = 0
        factory%potential_error_count = 0
        factory%state_error_count = 0
        factory%start_error_count = 0
        factory%integrator_error_count = 0
        factory%other_return_error_count = 0
        factory%invariant_success_count = 0
        factory%invariant_failure_count = 0
        factory%return_error_count = 0
        factory%wall_not_hit_count = 0
        factory%wall_error_count = 0
        factory%last_return_status = GC_EQDSK_NONLOCAL_ORBIT_ERROR
        factory%last_return_event_kind = 0
        factory%last_launch_event_value = 0.0_dp
        factory%last_return_event_value = 0.0_dp
        factory%last_invariant_scaled_drift = 0.0_dp
        factory%last_energy_error = 0.0_dp
        factory%last_magnetic_moment_error = 0.0_dp
        factory%last_canonical_momentum_error = 0.0_dp
        factory%topology_certification_attempts = 0
        factory%topology_certification_successes = 0
        factory%harmonic_average_successes = 0
        factory%residence_average_successes = 0
    end subroutine clear_gc_eqdsk_nonlocal_transport

    subroutine evaluate_gc_eqdsk_native_profile(factory, psi_star, density, &
            temperature_erg, potential, a1_star, a2_star, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        real(dp), intent(in) :: psi_star
        real(dp), intent(out) :: density, temperature_erg, potential
        real(dp), intent(out) :: a1_star, a2_star
        integer, intent(out) :: status

        real(dp) :: values(FACTORY_FORCE_COUNT + 3)

        density = 0.0_dp
        temperature_erg = 0.0_dp
        potential = 0.0_dp
        a1_star = 0.0_dp
        a2_star = 0.0_dp
        status = GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE
        if (.not. factory%profile_ready) return
        call interpolate_profile(factory, psi_star, values, status)
        if (status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        density = values(PROFILE_DENSITY)
        temperature_erg = values(PROFILE_TEMPERATURE)
        potential = values(PROFILE_PHI)
        a1_star = values(PROFILE_A1STAR)
        a2_star = values(PROFILE_A2STAR)
        if (.not. all(ieee_is_finite([density, temperature_erg, potential, &
            a1_star, a2_star]))) then
            density = 0.0_dp
            temperature_erg = 0.0_dp
            potential = 0.0_dp
            a1_star = 0.0_dp
            a2_star = 0.0_dp
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        if (density <= 0.0_dp .or. temperature_erg <= 0.0_dp) then
            status = GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE
            return
        end if
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine evaluate_gc_eqdsk_native_profile

    function gc_eqdsk_nonlocal_status_message(status) result(message)
        integer, intent(in) :: status
        character(len=64) :: message

        select case (status)
        case (GC_EQDSK_NONLOCAL_SUCCESS)
            message = 'success'
        case (GC_EQDSK_NONLOCAL_INVALID_INPUT)
            message = 'invalid input'
        case (GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE)
            message = 'direct EQDSK field unavailable'
        case (GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE)
            message = 'native profile unavailable'
        case (GC_EQDSK_NONLOCAL_PERTURBATION_UNAVAILABLE)
            message = 'direct R-Z perturbation unavailable'
        case (GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE)
            message = 'physical topology unavailable'
        case (GC_EQDSK_NONLOCAL_WALL_UNAVAILABLE)
            message = 'wall certification unavailable'
        case (GC_EQDSK_NONLOCAL_DERIVATIVE_UNAVAILABLE)
            message = 'orbit derivative unavailable'
        case (GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED)
            message = 'production certification failed'
        case (GC_EQDSK_NONLOCAL_SPECIES_MISMATCH)
            message = 'requested species does not match native profile'
        case (GC_EQDSK_NONLOCAL_DOMAIN_ERROR)
            message = 'direct EQDSK domain error'
        case (GC_EQDSK_NONLOCAL_ORBIT_ERROR)
            message = 'full-FOW orbit return failed'
        case (GC_EQDSK_NONLOCAL_NONFINITE)
            message = 'nonfinite physical quantity'
        case default
            message = 'unknown direct EQDSK nonlocal status'
        end select
    end function gc_eqdsk_nonlocal_status_message

    subroutine evaluate_eqdsk_profile_potential(self, position, field, potential, &
            gradient, status)
        class(gc_eqdsk_profile_potential_t), intent(in) :: self
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        real(dp), intent(out) :: potential, gradient(3)
        integer, intent(out) :: status

        real(dp) :: values(2), psi

        associate (unused_position => position)
        end associate
        potential = 0.0_dp
        gradient = 0.0_dp
        status = GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE
        if (.not. self%initialized) return
        if (self%c_light <= 0.0_dp) return
        psi = field%psi
        call interpolate_potential(self, psi, values, status)
        if (status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            write (*, '(a,1x,i0,3(1x,es24.16))') &
                'potential interpolation diagnostic=', status, psi, &
                self%psi_pol(1), self%psi_pol(size(self%psi_pol))
            return
        end if
        potential = values(1)
        gradient = values(2)*field%grad_psi
        if (.not. all(ieee_is_finite([potential, gradient]))) then
            potential = 0.0_dp
            gradient = 0.0_dp
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        status = GC_CYL_SUCCESS
    end subroutine evaluate_eqdsk_profile_potential

    subroutine initialize_direct_field(factory, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        integer, intent(out) :: status
        integer :: local_status

        factory%field_ready = .false.
        call initialize_eqdsk_cylindrical_field(trim(factory%eqdsk_path), &
            factory%options%field_scale, factory%field, local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            status = GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE
            return
        end if
        if (.not. factory%field%domain_initialized) then
            status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
            return
        end if
        if (.not. all(ieee_is_finite([R0, a, psi_pr]))) then
            status = GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE
            return
        end if
        if (a <= 0.0_dp) then
            status = GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE
            return
        end if
        factory%minor_radius_cm = a
        factory%psi_pr_effective = factory%options%field_scale*psi_pr
        factory%field_ready = .true.
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine initialize_direct_field

    subroutine load_direct_wall(factory, input_wall, message, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        type(wall_polygon_t), intent(out) :: input_wall
        character(len=*), intent(out) :: message
        integer, intent(out) :: status

        call load_wall_polygon(trim(factory%wall_path), input_wall, status, &
            message, trim(factory%wall_units))
    end subroutine load_direct_wall

    subroutine load_native_profiles(factory, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        integer, intent(out) :: status

        integer :: nplasma, nrotation, i, local_status
        real(dp) :: am1, am2, z1, z2, q_value, dqds_value
        real(dp) :: psi_value, dpsi_ds, psi_edge
        real(dp), allocatable :: plasma(:, :), rotation(:, :)
        real(dp) :: s_value
        character(len=256) :: normalization_message

        status = GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE
        factory%profile_ready = .false.
        call read_plasma_input(trim(factory%plasma_path), nplasma, am1, am2, &
            z1, z2, plasma)
        if (.not. allocated(plasma)) return
        if (nplasma < 3 .or. size(plasma, 2) /= 6) then
            deallocate(plasma)
            return
        end if
        if (.not. all(ieee_is_finite(plasma))) then
            deallocate(plasma)
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        if (abs(am1*mu - factory%species%mass_g) > &
            3.0e-10_dp*max(abs(factory%species%mass_g), tiny(1.0_dp))) then
            deallocate(plasma)
            status = GC_EQDSK_NONLOCAL_SPECIES_MISMATCH
            return
        end if
        if (abs(z1*qe - factory%species%charge_esu) > &
            3.0e-10_dp*max(abs(factory%species%charge_esu), tiny(1.0_dp))) then
            deallocate(plasma)
            status = GC_EQDSK_NONLOCAL_SPECIES_MISMATCH
            return
        end if
        do i = 2, nplasma
            if (plasma(i, 1) <= plasma(i - 1, 1)) then
                deallocate(plasma)
                return
            end if
        end do
        call extend_profile_to_closed_flux_domain(plasma, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            deallocate(plasma)
            status = local_status
            return
        end if
        nplasma = size(plasma, 1)
        if (plasma(1, 1) > factory%options%surface_min .or. &
                plasma(nplasma, 1) < factory%options%surface_max) then
            deallocate(plasma)
            status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
            return
        end if

        call prepare_plasma_splines(nplasma, am1, am2, z1, z2, plasma)
        call readdata(trim(factory%profile_path), 2, rotation)
        if (.not. allocated(rotation)) then
            deallocate(plasma)
            return
        end if
        nrotation = size(rotation, 1)
        if (nrotation < 2 .or. size(rotation, 2) /= 2) then
            deallocate(plasma, rotation)
            return
        end if
        if (.not. all(ieee_is_finite(rotation))) then
            deallocate(plasma, rotation)
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        do i = 2, nrotation
            if (rotation(i, 1) <= rotation(i - 1, 1)) then
                deallocate(plasma, rotation)
                return
            end if
        end do
        call extend_profile_to_closed_flux_domain(rotation, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            deallocate(plasma, rotation)
            status = local_status
            return
        end if
        nrotation = size(rotation, 1)
        if (rotation(1, 1) > factory%options%surface_min .or. &
                rotation(nrotation, 1) < factory%options%surface_max) then
            deallocate(plasma, rotation)
            status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
            return
        end if
        call prepare_profile_splines(rotation)

        if (allocated(factory%profile_s)) deallocate(factory%profile_s)
        if (allocated(factory%profile_psi)) deallocate(factory%profile_psi)
        if (allocated(factory%profile_dpsi_ds)) then
            deallocate(factory%profile_dpsi_ds)
        end if
        if (allocated(factory%profile_values)) deallocate(factory%profile_values)
        allocate(factory%profile_s(nplasma), factory%profile_psi(nplasma), &
            factory%profile_dpsi_ds(nplasma), &
            factory%profile_values(nplasma, 6))
        factory%profile_s = plasma(:, 1)
        factory%profile_values = 0.0_dp

        do i = 1, nplasma
            s_value = factory%profile_s(i)
            call geoflux_get_flux_profiles(s_value, q_value, dqds_value, &
                psi_value, dpsi_ds, psi_edge)
            call sample_one_profile_surface(factory, s_value, q_value, &
                psi_value, dpsi_ds, local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                deallocate(plasma, rotation)
                status = local_status
                return
            end if
            if (Ti1 <= 0.0_dp .or. ni1 <= 0.0_dp .or. &
                abs(dpsi_ds) <= tiny(dpsi_ds)) then
                deallocate(plasma, rotation)
                return
            end if
            factory%profile_psi(i) = factory%options%field_scale*psi_value
            factory%profile_dpsi_ds(i) = factory%options%field_scale*dpsi_ds
            factory%profile_values(i, PROFILE_DENSITY) = ni1
            factory%profile_values(i, PROFILE_TEMPERATURE) = Ti1*eV
            factory%profile_values(i, PROFILE_OMEGA_E) = Om_tE
            factory%profile_values(i, PROFILE_A1STAR) = &
                A1/factory%profile_dpsi_ds(i)
            factory%profile_values(i, PROFILE_A2STAR) = &
                A2/factory%profile_dpsi_ds(i)
        end do
        deallocate(plasma, rotation)

        call initialize_eqdsk_flux_profile_map(factory%profile_s, &
            factory%profile_psi, factory%options%field_scale, psi_sep, &
            factory%flux_map, local_status)
        if (local_status /= EQDSK_FLUX_MAP_SUCCESS) then
            status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
            return
        end if
        call initialize_profile_potential(factory, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        call initialize_gc_full_fow_reference(factory%species%mass_g, &
            factory%species%charge_esu, c, &
            factory%species%reference_energy_erg, &
            factory%species%reference_velocity_cm_s, factory%normalization, &
            local_status, normalization_message)
        if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) then
            status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
            return
        end if
        factory%phi_eff_normalization = factory%normalization%phi_eff
        factory%refinement_length_scale = factory%minor_radius_cm
        factory%refinement_time_scale = factory%minor_radius_cm/ &
            factory%normalization%reference_velocity
        factory%refinement_frequency_scale = &
            factory%normalization%reference_velocity/factory%minor_radius_cm
        factory%refinement_momentum_scale = factory%normalization%mass* &
            factory%normalization%reference_velocity
        if (.not. all(ieee_is_finite([factory%refinement_length_scale, &
                factory%refinement_time_scale, &
                factory%refinement_frequency_scale, &
                factory%refinement_momentum_scale]))) then
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        if (min(factory%refinement_length_scale, factory%refinement_time_scale, &
                factory%refinement_frequency_scale, &
                factory%refinement_momentum_scale) <= 0.0_dp) then
            status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
            return
        end if
        factory%profile_ready = .true.
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine load_native_profiles

    subroutine extend_profile_to_closed_flux_domain(values, status)
        !! Extend cell-centred data to s=0 and s=1 using only the generated
        !! affine endpoint kernel.  The input nodes remain untouched.
        real(dp), allocatable, intent(inout) :: values(:, :)
        integer, intent(out) :: status

        real(dp), allocatable :: extended(:, :)
        real(dp) :: left_value, right_value, unused(6)
        integer :: old_count, new_count, first_row, column
        logical :: add_axis, add_edge

        status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
        if (.not. allocated(values)) return
        old_count = size(values, 1)
        if (old_count < 2 .or. size(values, 2) < 2) return
        if (.not. all(ieee_is_finite(values))) return
        if (values(1, 1) < 0.0_dp .or. values(old_count, 1) > 1.0_dp) return
        add_axis = values(1, 1) > 0.0_dp
        add_edge = values(old_count, 1) < 1.0_dp
        if (.not. add_axis .and. .not. add_edge) then
            status = GC_EQDSK_NONLOCAL_SUCCESS
            return
        end if

        new_count = old_count+merge(1, 0, add_axis)+merge(1, 0, add_edge)
        allocate(extended(new_count, size(values, 2)))
        first_row = merge(2, 1, add_axis)
        extended(first_row:first_row+old_count-1, :) = values
        if (add_axis) extended(1, 1) = 0.0_dp
        if (add_edge) extended(new_count, 1) = 1.0_dp
        do column = 2, size(values, 2)
            call evaluate_neort_profile_endpoints(values(1, 1), values(2, 1), &
                values(1, column), values(2, column), 0.0_dp, left_value, &
                unused(1), unused(2), unused(3), unused(4), unused(5), &
                unused(6), right_value)
            if (add_axis) extended(1, column) = left_value
            call evaluate_neort_profile_endpoints(values(old_count-1, 1), &
                values(old_count, 1), values(old_count-1, column), &
                values(old_count, column), 1.0_dp, unused(1), right_value, &
                unused(2), unused(3), unused(4), unused(5), unused(6), &
                left_value)
            if (add_edge) extended(new_count, column) = right_value
        end do
        if (.not. all(ieee_is_finite(extended))) return
        call move_alloc(extended, values)
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine extend_profile_to_closed_flux_domain

    subroutine sample_one_profile_surface(factory, s_value, q_value, psi_value, &
            dpsi_ds, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: s_value
        real(dp), intent(inout) :: q_value, psi_value, dpsi_ds
        integer, intent(out) :: status

        real(dp) :: dqds_unused, psi_edge_unused

        call geoflux_get_flux_profiles(s_value, q_value, dqds_unused, &
            psi_value, dpsi_ds, psi_edge_unused)
        if (.not. all(ieee_is_finite([q_value, psi_value, dpsi_ds]))) then
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        if (abs(q_value) <= tiny(q_value)) then
            status = GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE
            return
        end if
        s = s_value
        call init_plasma_at_s()
        call init_profile_at_s(R0, factory%options%profile_electric_factor, &
            factory%options%profile_bfactor)
        call init_thermodynamic_forces(factory%psi_pr_effective, q_value)
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine sample_one_profile_surface

    subroutine initialize_profile_potential(factory, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        integer, intent(out) :: status

        real(dp) :: delta_phi, reversed_delta, constant_limit, scale
        integer :: i, node_count

        status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
        call clear_profile_potential(factory%potential)
        factory%potential%derivative_convention = &
            factory%profile_potential_derivative_convention
        factory%potential%gauge_convention = factory%profile_potential_gauge
        factory%potential%segment_integral_method = &
            'fortsym-affine-segment-integral-required'
        factory%potential%endpoint_reconstruction_method = &
            factory%profile_endpoint_reconstruction_method
        if (.not. allocated(factory%profile_s)) return
        if (.not. allocated(factory%profile_psi)) return
        if (.not. allocated(factory%profile_values)) return
        node_count = size(factory%profile_s)
        if (node_count < 2) return
        if (size(factory%profile_psi) /= node_count) return
        if (size(factory%profile_values, 1) /= node_count) return
        if (.not. all(ieee_is_finite(factory%profile_s))) return
        if (.not. all(ieee_is_finite(factory%profile_psi))) return
        if (.not. all(ieee_is_finite( &
                factory%profile_values(:, PROFILE_OMEGA_E)))) return
        if (factory%profile_s(1) /= 0.0_dp) return
        if (factory%profile_s(node_count) /= 1.0_dp) return
        allocate(factory%potential%s_tor(node_count), &
            factory%potential%psi_pol(node_count), &
            factory%potential%phi(node_count), &
            factory%potential%omega_e(node_count))
        factory%potential%s_tor = factory%profile_s
        factory%potential%psi_pol = factory%profile_psi
        factory%potential%omega_e = &
            factory%profile_values(:, PROFILE_OMEGA_E)
        factory%potential%phi = 0.0_dp
        do i = 1, node_count-1
            call evaluate_neort_profile_potential_segment( &
                factory%potential%psi_pol(i), factory%potential%psi_pol(i+1), &
                factory%potential%omega_e(i), factory%potential%omega_e(i+1), &
                factory%potential%omega_e(i), c, delta_phi, reversed_delta, &
                constant_limit)
            if (.not. all(ieee_is_finite([delta_phi, reversed_delta, &
                    constant_limit]))) return
            scale = max(1.0_dp, max(abs(delta_phi), abs(reversed_delta)))
            if (abs(delta_phi+reversed_delta) > 2.0e-12_dp*scale) return
            factory%potential%phi(i+1) = &
                factory%potential%phi(i)+delta_phi
        end do
        if (.not. all(ieee_is_finite(factory%potential%phi))) return
        factory%potential%c_light = c
        factory%potential%segment_integral_method = &
            'fortsym-affine-segment-integral'
        factory%potential%endpoint_reconstruction_method = &
            'fortsym-affine-endpoints-s=0,1'
        factory%potential%initialized = .true.
        call validate_eqdsk_potential_profile_nodes( &
            factory%potential%psi_pol, factory%potential%phi, &
            factory%potential%omega_e, factory%certified_profile, status)
        if (status /= EQDSK_CUT_BOX_SUCCESS) then
            status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
            return
        end if
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine initialize_profile_potential

    subroutine clear_profile_potential(potential)
        type(gc_eqdsk_profile_potential_t), intent(inout) :: potential

        if (allocated(potential%s_tor)) deallocate(potential%s_tor)
        if (allocated(potential%psi_pol)) deallocate(potential%psi_pol)
        if (allocated(potential%phi)) deallocate(potential%phi)
        if (allocated(potential%omega_e)) deallocate(potential%omega_e)
        potential%c_light = 0.0_dp
        potential%initialized = .false.
        potential%derivative_convention = 'dPhi/dpsi=+Omega_E/c'
        potential%gauge_convention = 'Phi(psi_axis)=0'
        potential%segment_integral_method = ''
        potential%endpoint_reconstruction_method = ''
    end subroutine clear_profile_potential

    subroutine load_direct_perturbation(factory, harmonic_n, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        integer, intent(in) :: harmonic_n
        integer, intent(out) :: status

        factory%perturbation_ready = .false.
        factory%perturbation_provenance_certified = .false.
        factory%perturbation_amplitude_convention = ''
        factory%perturbation_coordinate_dependency = ''
        factory%perturbation_interpolation_certified = .false.
        if (harmonic_n == 0) then
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
            return
        end if
        if (factory%options%profile_bfactor <= 0.0_dp) then
            status = GC_EQDSK_NONLOCAL_PERTURBATION_UNAVAILABLE
            return
        end if
        if (inp_swi /= 11) then
            status = GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE
            return
        end if
        ! The direct R-Z reader uses the shared NEO-RT perturbation b-factor.
        ! This is a serialized standalone factory boundary; the legacy global
        ! state has no instance API, so concurrent direct factories are not
        ! certified by this seam.
        bfac = factory%options%profile_bfactor
        inp_swi_pert = 11
        call read_boozer_pert_file(trim(factory%perturbation_path))
        call set_mph(harmonic_n)
        if (inp_swi_pert /= 11 .or. rz_nrad < 2 .or. rz_nzet < 2) then
            status = GC_EQDSK_NONLOCAL_PERTURBATION_UNAVAILABLE
            return
        end if
        ! do_magfie_pert_amp returns the full real-field single-n amplitude
        ! A in delta B=Re[A exp(i*n*phi)].  It is not a one-sided Fourier
        ! coefficient A/2; the implicit conjugate is handled by Eq.17 /4.
        factory%perturbation_amplitude_convention = &
            'real_field_amplitude_one_signed_n'
        factory%perturbation_coordinate_dependency = 'direct_cylindrical_RZ'
        factory%perturbation_interpolation_certified = .true.
        factory%perturbation_provenance_certified = &
            trim(factory%perturbation_path) /= ''
        factory%perturbation_ready = &
            factory%perturbation_interpolation_certified .and. &
            factory%perturbation_provenance_certified
        if (.not. factory%perturbation_ready) then
            status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
            return
        end if
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine load_direct_perturbation

    subroutine build_poincare_reference(factory, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        integer, intent(out) :: status

        type(gc_cylindrical_field_sample_t) :: section_field
        real(dp) :: dposition(3), distance
        integer :: wall_status, local_status

        factory%cut_ready = .false.
        status = GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE
        call build_cut_atlas(factory, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        if (.not. factory%cut_atlas%certified) return
        if (.not. factory%certified_cut_atlas% &
                surface_intersection_pair_certified) return
        call physical_cut_map(factory, factory%options%reference_surface, &
            factory%section_position, dposition, status)
        if (status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        call factory%wall%evaluate(factory%section_position, distance, wall_status)
        if (wall_status /= GC_CYL_SUCCESS .or. distance <= 0.0_dp) then
            status = GC_EQDSK_NONLOCAL_WALL_UNAVAILABLE
            return
        end if
        call factory%field%evaluate(factory%section_position, section_field, &
            local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            status = GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE
            return
        end if
        factory%section_flux = section_field%psi
        if (.not. ieee_is_finite(factory%section_flux)) then
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        if (abs(factory%section_orientation) /= 1 .or. &
                .not. factory%native_orientation_certified) then
            status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
            return
        end if
        ! This reference certifies only the geometric section.  Finite-width
        ! orbit multiplicity remains an orbit-evidence gate and is not needed
        ! to construct or lock the section itself.
        factory%section_reference_id = 'direct-eqdsk-C0-cut-atlas-v3'
        factory%cut_ready = .true.
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine build_poincare_reference

    subroutine build_cut_atlas(factory, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        integer, intent(out) :: status

        real(dp) :: axis_R, axis_Z, target_s_tor, target_psihat
        real(dp) :: unused_derivative
        real(dp) :: inboard_seed(2), outboard_seed(2)
        real(dp) :: axis_box(4), inboard_box(4), outboard_box(4)
        real(dp) :: inboard_theta, outboard_theta, position(3)
        integer :: axis_R_index, axis_Z_index, local_status

        if (allocated(factory%cut_atlas%branches)) then
            deallocate(factory%cut_atlas%branches)
        end if
        factory%cut_atlas = gc_eqdsk_cut_atlas_t()
        factory%certified_cut_atlas = eqdsk_composite_cut_atlas_t()
        factory%r_partition = eqdsk_composite_r_partition_t()
        factory%allowed_region_context = &
            gc_eqdsk_certified_allowed_provider_context_t()
        factory%allowed_region_ready = .false.
        factory%cut_atlas%method = 'fortsym-interval-cut-atlas-required'
        factory%cut_atlas_method = trim(factory%cut_atlas%method)
        factory%cut_atlas_certified = .false.
        factory%cut_atlas_branch_count = 0
        factory%active_cut_branch_id = 0
        factory%section_orientation = 0
        factory%cut_parameter_orientation = ''
        factory%canonical_one_form_convention = ''
        factory%fourier_phase_convention = ''
        factory%native_orientation_certified = .false.
        status = GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE
        if (.not. factory%field_ready .or. .not. factory%flux_map%initialized) &
            return
        if (nrad < 3 .or. nzet < 3) return
        call map_eqdsk_rho_tor_to_psihat(factory%flux_map, &
            factory%options%surface_max, target_s_tor, target_psihat, &
            unused_derivative, local_status)
        if (local_status /= EQDSK_FLUX_MAP_SUCCESS) return

        call geoflux_get_axis(axis_R, axis_Z)
        if (.not. all(ieee_is_finite([axis_R, axis_Z]))) return
        axis_R_index = minloc(abs(rad-axis_R), dim=1)
        axis_Z_index = minloc(abs(zet-axis_Z), dim=1)
        if (axis_R_index <= 1 .or. axis_R_index >= nrad) return
        if (axis_Z_index <= 1 .or. axis_Z_index >= nzet) return
        axis_box = [rad(axis_R_index-1), rad(axis_R_index+1), &
            zet(axis_Z_index-1), zet(axis_Z_index+1)]

        ! The geoflux scan is only a nonlinear-solver preconditioner.  It
        ! cannot certify a root or multiplicity; both seeds are subsequently
        ! enclosed by generated Krawczyk maps and the complete interval atlas.
        call find_cut_endpoint_thetas(factory, target_s_tor, inboard_theta, &
            outboard_theta, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        call map_eqdsk_flux_position(target_s_tor, inboard_theta, 0.0_dp, &
            position, local_status)
        if (local_status /= GC_CYL_SUCCESS) return
        inboard_seed = position(1:2)
        call map_eqdsk_flux_position(target_s_tor, outboard_theta, 0.0_dp, &
            position, local_status)
        if (local_status /= GC_CYL_SUCCESS) return
        outboard_seed = position(1:2)
        call certify_endpoint_preconditioner(factory, target_psihat, &
            axis_box, EQDSK_COMPOSITE_CUT_INBOARD, inboard_seed, &
            inboard_box, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        call certify_endpoint_preconditioner(factory, target_psihat, &
            axis_box, EQDSK_COMPOSITE_CUT_OUTBOARD, outboard_seed, &
            outboard_box, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        call build_eqdsk_composite_cut_atlas(axis_box, inboard_box, &
            outboard_box, inboard_seed, outboard_seed, target_psihat, &
            factory%options%field_scale, zet(1), zet(nzet), &
            factory%options%cut_atlas_options, &
            factory%certified_cut_atlas, local_status)
        if (local_status /= EQDSK_COMPOSITE_ATLAS_SUCCESS) then
            write (*, '(a,1x,i0,8(1x,i0),4(1x,es24.16))') 'cut atlas diagnostic=', local_status, &
                factory%certified_cut_atlas%inboard_graph%failure_stage, &
                factory%certified_cut_atlas%inboard_graph%failure_r_depth, &
                factory%certified_cut_atlas%inboard_graph%failure_z_depth, &
                factory%certified_cut_atlas%axis_graph%failure_stage, &
                factory%certified_cut_atlas%axis_graph%failure_r_depth, &
                factory%certified_cut_atlas%axis_graph%failure_z_depth, &
                factory%certified_cut_atlas%outboard_graph%failure_stage, &
                factory%certified_cut_atlas%outboard_graph%failure_r_depth, &
                factory%certified_cut_atlas%inboard_graph%failure_r_lo, &
                factory%certified_cut_atlas%inboard_graph%failure_r_hi, &
                factory%certified_cut_atlas%inboard_graph%failure_z_lo, &
                factory%certified_cut_atlas%inboard_graph%failure_z_hi
            return
        end if
        call build_eqdsk_composite_r_partition(factory%certified_cut_atlas, &
            factory%r_partition, local_status)
        if (local_status /= EQDSK_R_OWNERSHIP_SUCCESS) return
        call initialize_gc_eqdsk_certified_allowed_provider( &
            factory%allowed_region_context, factory%certified_cut_atlas, &
            factory%r_partition, factory%certified_profile, &
            factory%options%field_scale, &
            factory%certified_cut_atlas%inboard_graph%raw_psi_sep, &
            factory%species%mass_g, factory%species%charge_esu, c, local_status)
        if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) return
        factory%allowed_region_ready = .true.

        allocate(factory%cut_atlas%branches(2))
        factory%cut_atlas%nbranches = 2
        factory%cut_atlas%branches(1)%branch_id = 1
        factory%cut_atlas%branches(1)%rc_min = factory%options%surface_min
        factory%cut_atlas%branches(1)%rc_max = factory%options%surface_max
        factory%cut_atlas%branches(1)%orientation = &
            EQDSK_COMPOSITE_CUT_INBOARD
        factory%cut_atlas%branches(1)%continuation_certified = .true.
        factory%cut_atlas%branches(1)%root_isolation_certified = .true.
        factory%cut_atlas%branches(1)%lower_boundary_kind = 'axis-limit'
        factory%cut_atlas%branches(1)%upper_boundary_kind = 'regular'
        factory%cut_atlas%branches(1)%limiting_chart = &
            'fortsym-rho-tor-axis'
        factory%cut_atlas%branches(2) = factory%cut_atlas%branches(1)
        factory%cut_atlas%branches(2)%branch_id = 2
        factory%cut_atlas%branches(2)%orientation = &
            EQDSK_COMPOSITE_CUT_OUTBOARD
        factory%cut_atlas%certified = &
            factory%certified_cut_atlas%geometric_completeness_certified &
            .and. factory%certified_cut_atlas%branch_connectivity_certified &
            .and. factory%certified_cut_atlas% &
                surface_intersection_pair_certified
        ! Physical-orbit multiplicity remains a separate theorem gate.  The
        ! geometric certificate must never promote this flag by itself.
        factory%cut_atlas%two_cut_multiplicity_certified = .false.
        factory%cut_atlas%method = 'fortsym-interval-composite-Eq13-atlas'
        factory%cut_atlas_method = trim(factory%cut_atlas%method)
        factory%cut_atlas_certified = factory%cut_atlas%certified
        factory%cut_atlas_branch_count = factory%cut_atlas%nbranches
        factory%active_cut_branch_id = 2
        factory%section_orientation = 1
        factory%cut_parameter_orientation = &
            'rho_tor:axis-to-outboard;C_native=+1'
        factory%canonical_one_form_convention = &
            'P_phi=(q/c)psi+m*v_parallel*R*b_phi'
        factory%fourier_phase_convention = 'Re[A_n*exp(+i*n*phi)]'
        factory%native_orientation_certified = factory%cut_atlas%certified
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine build_cut_atlas

    subroutine certify_endpoint_preconditioner(factory, target_psihat, &
            axis_box, branch_sign, seed, box, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        real(dp), intent(in) :: target_psihat, axis_box(4)
        integer, intent(in) :: branch_sign
        real(dp), intent(inout) :: seed(2)
        real(dp), intent(out) :: box(4)
        integer, intent(out) :: status

        type(eqdsk_cut_endpoint_certificate_t) :: endpoint
        real(dp) :: base_R, base_Z, half_R, half_Z, scale, guard
        integer :: expansion, local_status

        box = 0.0_dp
        status = GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE
        if (abs(branch_sign) /= 1) return
        if (.not. all(ieee_is_finite([target_psihat, axis_box, seed]))) return
        base_R = minval(rad(2:nrad)-rad(1:nrad-1))
        base_Z = minval(zet(2:nzet)-zet(1:nzet-1))
        if (min(base_R, base_Z) <= 0.0_dp) return
        guard = 1024.0_dp*epsilon(max(1.0_dp, maxval(abs(axis_box))))
        do expansion = 0, factory%options%endpoint_max_box_expansions-1
            scale = factory%options%endpoint_initial_grid_fraction &
                *2.0_dp**expansion
            half_R = scale*base_R
            half_Z = scale*base_Z
            box = [max(rad(1), seed(1)-half_R), &
                min(rad(nrad), seed(1)+half_R), &
                max(zet(1), seed(2)-half_Z), &
                min(zet(nzet), seed(2)+half_Z)]
            if (branch_sign == EQDSK_COMPOSITE_CUT_INBOARD) then
                box(2) = min(box(2), axis_box(1)-guard)
            else
                box(1) = max(box(1), axis_box(2)+guard)
            end if
            if (box(1) >= seed(1) .or. box(2) <= seed(1) .or. &
                    box(3) >= seed(2) .or. box(4) <= seed(2)) cycle
            call build_eqdsk_cut_endpoint_certificate(box(1), box(2), &
                box(3), box(4), target_psihat, factory%options%field_scale, &
                seed(1), seed(2), factory%options%cut_atlas_options%endpoint, &
                endpoint, local_status)
            if (local_status == EQDSK_ENDPOINT_CERT_SUCCESS) then
                seed = [endpoint%newton_point_R, endpoint%newton_point_Z]
                status = GC_EQDSK_NONLOCAL_SUCCESS
                return
            end if
        end do
    end subroutine certify_endpoint_preconditioner

    subroutine interpolate_profile(factory, psi, values, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        real(dp), intent(in) :: psi
        real(dp), intent(out) :: values(:)
        integer, intent(out) :: status

        integer :: left, right, n, i
        real(dp) :: weight

        values = 0.0_dp
        status = GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE
        if (.not. factory%profile_ready) return
        n = size(factory%profile_psi)
        if (size(values) < size(factory%profile_values, 2)) return
        call locate_monotonic(factory%profile_psi, psi, left, right, weight, &
            status)
        if (status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        do i = 1, size(factory%profile_values, 2)
            values(i) = (1.0_dp - weight)*factory%profile_values(left, i) &
                +weight*factory%profile_values(right, i)
        end do
        if (n < 2 .or. .not. all(ieee_is_finite(values))) then
            values = 0.0_dp
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine interpolate_profile

    subroutine interpolate_potential(potential, psi, values, status, &
            second_derivative)
        class(gc_eqdsk_profile_potential_t), intent(in) :: potential
        real(dp), intent(in) :: psi
        real(dp), intent(out) :: values(2)
        integer, intent(out) :: status
        real(dp), intent(out), optional :: second_derivative

        integer :: left, right
        real(dp) :: weight, dpsi, local_second_derivative

        values = 0.0_dp
        if (present(second_derivative)) second_derivative = 0.0_dp
        status = GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE
        if (.not. potential%initialized) return
        call locate_monotonic(potential%psi_pol, psi, left, right, weight, &
            status)
        if (status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        dpsi = potential%psi_pol(right) - potential%psi_pol(left)
        if (abs(dpsi) <= tiny(dpsi)) then
            status = GC_EQDSK_NONLOCAL_DERIVATIVE_UNAVAILABLE
            return
        end if
        call evaluate_neort_profile_potential_map(psi, &
            potential%psi_pol(left), potential%psi_pol(right), &
            potential%phi(left), potential%omega_e(left), &
            potential%omega_e(right), potential%c_light, values(1), &
            values(2), local_second_derivative)
        if (.not. all(ieee_is_finite([values, local_second_derivative]))) then
            values = 0.0_dp
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        if (present(second_derivative)) &
            second_derivative = local_second_derivative
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine interpolate_potential

    subroutine locate_monotonic(axis, value, left, right, weight, status)
        real(dp), intent(in) :: axis(:), value
        integer, intent(out) :: left, right
        real(dp), intent(out) :: weight
        integer, intent(out) :: status

        integer :: low, high, mid, n
        logical :: ascending
        real(dp) :: tolerance, denominator, interpolation_value, axis_scale

        left = 0
        right = 0
        weight = 0.0_dp
        status = GC_EQDSK_NONLOCAL_PROFILE_UNAVAILABLE
        n = size(axis)
        if (n < 2) return
        if (.not. all(ieee_is_finite([axis, value]))) then
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        ascending = axis(n) > axis(1)
        axis_scale = max(1.0_dp, maxval(abs([axis(1), axis(n), value])))
        tolerance = 100.0_dp*epsilon(axis_scale)*axis_scale
        if (ascending) then
            if (value < axis(1)-tolerance .or. value > axis(n)+tolerance) return
            interpolation_value = min(axis(n), max(axis(1), value))
        else
            if (value > axis(1)+tolerance .or. value < axis(n)-tolerance) return
            interpolation_value = max(axis(n), min(axis(1), value))
        end if
        low = 1
        high = n
        do while (high - low > 1)
            mid = (low + high)/2
            if (ascending) then
                if (axis(mid) <= interpolation_value) then
                    low = mid
                else
                    high = mid
                end if
            else
                if (axis(mid) >= interpolation_value) then
                    low = mid
                else
                    high = mid
                end if
            end if
        end do
        if (interpolation_value == axis(n)) then
            low = n - 1
            high = n
        end if
        denominator = axis(high) - axis(low)
        tolerance = 100.0_dp*epsilon(max(1.0_dp, abs(value)))
        if (abs(denominator) <= tolerance) then
            status = GC_EQDSK_NONLOCAL_DERIVATIVE_UNAVAILABLE
            return
        end if
        left = low
        right = high
        weight = (interpolation_value - axis(left))/denominator
        weight = max(0.0_dp, min(1.0_dp, weight))
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine locate_monotonic

    subroutine physical_cut_map(factory, surface, position, dposition_ds, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: surface
        real(dp), intent(out) :: position(3), dposition_ds(3)
        integer, intent(out) :: status

        real(dp) :: theta, theta_minus, theta_plus, step, surface_step
        real(dp) :: position_minus(3), position_plus(3)
        logical :: has_minus, has_plus
        integer :: local_status

        position = 0.0_dp
        dposition_ds = 0.0_dp
        status = GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE
        ! Production mapping is owned by the certified cut atlas.  The old
        ! largest-R finite-theta root is intentionally not a fallback: it
        ! omits the inboard intersection and has no axis-safe branch
        ! continuation or multiplicity proof.
        if (.not. factory%cut_atlas_certified) return
        call map_certified_cut_atlas(factory, surface, position, dposition_ds, &
            status)
        return

        if (.not. factory%field_ready) return
        if (surface < factory%options%surface_min .or. &
            surface > factory%options%surface_max) return
        call find_outboard_cut_theta(factory, surface, theta, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        call map_eqdsk_flux_position(surface, theta, 0.0_dp, position, &
            local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
            return
        end if
        surface_step = factory%options%cut_derivative_step &
            *max(1.0_dp, abs(surface))
        step = min(surface_step, 0.25_dp &
            *(factory%options%surface_max-factory%options%surface_min))
        if (step <= 0.0_dp) return
        has_minus = surface - step >= factory%options%surface_min
        has_plus = surface + step <= factory%options%surface_max
        if (has_minus .and. has_plus) then
            call map_cut_without_derivative(factory, surface-step, position_minus, &
                local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                status = local_status
                return
            end if
            call map_cut_without_derivative(factory, surface+step, position_plus, &
                local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                status = local_status
                return
            end if
            dposition_ds = (position_plus-position_minus)/(2.0_dp*step)
        else if (has_plus) then
            call map_cut_without_derivative(factory, surface+step, position_plus, &
                local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                status = local_status
                return
            end if
            dposition_ds = (position_plus-position)/step
        else if (has_minus) then
            call map_cut_without_derivative(factory, surface-step, position_minus, &
                local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                status = local_status
                return
            end if
            dposition_ds = (position-position_minus)/step
        else
            return
        end if
        if (.not. all(ieee_is_finite([position, dposition_ds]))) then
            position = 0.0_dp
            dposition_ds = 0.0_dp
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        if (position(1) <= 0.0_dp) then
            status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
            return
        end if
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine physical_cut_map

    subroutine map_certified_cut_atlas(factory, rc, position, dposition_drc, &
            status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: rc
        real(dp), intent(out) :: position(3), dposition_drc(3)
        integer, intent(out) :: status

        integer :: branch_sign, local_status

        position = 0.0_dp
        dposition_drc = 0.0_dp
        status = GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE
        if (.not. factory%cut_atlas_certified) return
        select case (factory%active_cut_branch_id)
        case (1)
            branch_sign = EQDSK_COMPOSITE_CUT_INBOARD
        case (2)
            branch_sign = EQDSK_COMPOSITE_CUT_OUTBOARD
        case default
            return
        end select
        call map_eqdsk_composite_cut_atlas_rho( &
            factory%certified_cut_atlas, factory%flux_map, rc, branch_sign, &
            position, dposition_drc, local_status)
        if (local_status /= EQDSK_COMPOSITE_ATLAS_SUCCESS) return
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine map_certified_cut_atlas

    subroutine map_cut_without_derivative(factory, surface, position, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: surface
        real(dp), intent(out) :: position(3)
        integer, intent(out) :: status
        real(dp) :: theta

        position = 0.0_dp
        call find_outboard_cut_theta(factory, surface, theta, status)
        if (status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        call map_eqdsk_flux_position(surface, theta, 0.0_dp, position, status)
        if (status /= GC_CYL_SUCCESS) status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
    end subroutine map_cut_without_derivative

    subroutine find_outboard_cut_theta(factory, surface, theta, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: surface
        real(dp), intent(out) :: theta
        integer, intent(out) :: status

        integer :: n, i, nroots, local_status
        real(dp) :: theta_left, theta_right, f_left, f_right, root
        real(dp) :: roots(MAX_CUT_ROOTS), positions(3, MAX_CUT_ROOTS)
        real(dp) :: theta_scale, best_radius

        theta = 0.0_dp
        status = GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE
        n = max(32, factory%options%cut_theta_points)
        nroots = 0
        theta_scale = 2.0_dp*pi/real(n, dp)
        theta_left = -pi
        call poincare_cut_value(factory, surface, theta_left, f_left, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        do i = 1, n
            theta_right = -pi + real(i, dp)*theta_scale
            call poincare_cut_value(factory, surface, theta_right, f_right, &
                local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                status = local_status
                return
            end if
            if (f_left == 0.0_dp) then
                root = theta_left
                call append_cut_root(root, roots, nroots)
            end if
            if (f_left*f_right < 0.0_dp) then
                call bisect_cut_root(factory, surface, theta_left, theta_right, &
                    f_left, f_right, root, local_status)
                if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                    status = local_status
                    return
                end if
                call append_cut_root(root, roots, nroots)
            end if
            theta_left = theta_right
            f_left = f_right
        end do
        if (abs(theta_left-pi) <= 10.0_dp*epsilon(pi) .and. &
            f_left == 0.0_dp) then
            call append_cut_root(-pi, roots, nroots)
        end if
        if (nroots < 1) return
        do i = 1, nroots
            call map_eqdsk_flux_position(surface, roots(i), 0.0_dp, &
                positions(:, i), local_status)
            if (local_status /= GC_CYL_SUCCESS) then
                status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
                return
            end if
        end do
        best_radius = -huge(1.0_dp)
        do i = 1, nroots
            if (positions(1, i) > best_radius) then
                best_radius = positions(1, i)
                theta = roots(i)
            end if
        end do
        if (.not. ieee_is_finite(theta)) then
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine find_outboard_cut_theta

    subroutine find_cut_endpoint_thetas(factory, surface, inboard_theta, &
            outboard_theta, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: surface
        real(dp), intent(out) :: inboard_theta, outboard_theta
        integer, intent(out) :: status

        real(dp) :: theta_left, theta_right, f_left, f_right, root
        real(dp) :: roots(MAX_CUT_ROOTS), positions(3, MAX_CUT_ROOTS)
        real(dp) :: theta_scale
        integer :: n, i, nroots, local_status, inboard_index, outboard_index

        inboard_theta = 0.0_dp
        outboard_theta = 0.0_dp
        status = GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE
        n = max(32, factory%options%cut_theta_points)
        nroots = 0
        theta_scale = 2.0_dp*pi/real(n, dp)
        theta_left = -pi
        call poincare_cut_value(factory, surface, theta_left, f_left, &
            local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        do i = 1, n
            theta_right = -pi+real(i, dp)*theta_scale
            call poincare_cut_value(factory, surface, theta_right, f_right, &
                local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                status = local_status
                return
            end if
            if (f_left == 0.0_dp) then
                call append_cut_root(theta_left, roots, nroots)
            end if
            if (f_left*f_right < 0.0_dp) then
                call bisect_cut_root(factory, surface, theta_left, &
                    theta_right, f_left, f_right, root, local_status)
                if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                    status = local_status
                    return
                end if
                call append_cut_root(root, roots, nroots)
            end if
            theta_left = theta_right
            f_left = f_right
        end do
        if (abs(theta_left-pi) <= 10.0_dp*epsilon(pi) .and. &
                f_left == 0.0_dp) then
            call append_cut_root(-pi, roots, nroots)
        end if
        ! This is only a seed gate, but accepting a third resolved extremum
        ! here would contradict the requested regular two-endpoint surface.
        if (nroots /= 2) return
        do i = 1, nroots
            call map_eqdsk_flux_position(surface, roots(i), 0.0_dp, &
                positions(:, i), local_status)
            if (local_status /= GC_CYL_SUCCESS) then
                status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
                return
            end if
        end do
        inboard_index = minloc(positions(1, 1:nroots), dim=1)
        outboard_index = maxloc(positions(1, 1:nroots), dim=1)
        if (inboard_index == outboard_index) return
        inboard_theta = roots(inboard_index)
        outboard_theta = roots(outboard_index)
        if (.not. all(ieee_is_finite([inboard_theta, outboard_theta]))) then
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine find_cut_endpoint_thetas

    subroutine append_cut_root(root, roots, nroots)
        real(dp), intent(in) :: root
        real(dp), intent(inout) :: roots(:)
        integer, intent(inout) :: nroots
        integer :: i

        do i = 1, nroots
            if (abs(root-roots(i)) <= 1.0e-8_dp) return
        end do
        if (nroots < size(roots)) then
            nroots = nroots + 1
            roots(nroots) = root
        end if
    end subroutine append_cut_root

    subroutine bisect_cut_root(factory, surface, left, right, fleft, fright, &
            root, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: surface, left, right, fleft, fright
        real(dp), intent(out) :: root
        integer, intent(out) :: status

        integer :: iteration
        real(dp) :: a_left, a_right, fa, fb, midpoint, fm

        a_left = left
        a_right = right
        fa = fleft
        fb = fright
        root = 0.5_dp*(left+right)
        status = GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE
        do iteration = 1, 100
            midpoint = 0.5_dp*(a_left+a_right)
            call poincare_cut_value(factory, surface, midpoint, fm, status)
            if (status /= GC_EQDSK_NONLOCAL_SUCCESS) return
            if (abs(fm) <= factory%options%cut_root_tolerance .or. &
                abs(a_right-a_left) <= 1.0e-12_dp) then
                root = midpoint
                status = GC_EQDSK_NONLOCAL_SUCCESS
                return
            end if
            if (fa*fm <= 0.0_dp) then
                a_right = midpoint
                fb = fm
            else
                a_left = midpoint
                fa = fm
            end if
        end do
        root = 0.5_dp*(a_left+a_right)
        call poincare_cut_value(factory, surface, root, fm, status)
        if (status == GC_EQDSK_NONLOCAL_SUCCESS) then
            if (abs(fm) > factory%options%cut_root_tolerance) then
                status = GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE
            end if
        end if
        associate (unused_fb => fb)
        end associate
    end subroutine bisect_cut_root

    subroutine poincare_cut_value(factory, surface, theta, value, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: surface, theta
        real(dp), intent(out) :: value
        integer, intent(out) :: status
        real(dp) :: position(3)

        value = 0.0_dp
        call map_eqdsk_flux_position(surface, theta, 0.0_dp, position, status)
        if (status /= GC_CYL_SUCCESS) then
            status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
            return
        end if
        call poincare_cut_value_position(factory, position, value, status)
    end subroutine poincare_cut_value

    subroutine poincare_cut_value_position(factory, position, value, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: value
        integer, intent(out) :: status
        type(eqdsk_cut_jet_t) :: cut
        integer :: cut_status

        value = 0.0_dp
        call evaluate_eqdsk_cut_jet(position, factory%field%field_scale, 1, &
            [0.0_dp, 0.0_dp, 0.0_dp], cut, cut_status)
        if (cut_status == EQDSK_CUT_JET_OUT_OF_DOMAIN) then
            status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
            return
        end if
        if (cut_status == EQDSK_CUT_JET_NONFINITE) then
            status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        if (cut_status /= EQDSK_CUT_JET_SUCCESS) then
            status = GC_EQDSK_NONLOCAL_FIELD_UNAVAILABLE
            return
        end if
        value = cut%cut_value
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine poincare_cut_value_position

    subroutine initialize_factory_class_adapter(factory, h0, jperp, adapter, &
            status)
        !! Own the complete adapter contract for every factory callback seam.
        !!
        !! Keep future allowed-region provider, verifier, and certificate
        !! associations here so the factory call sites cannot diverge.  The
        !! splitter is borrowed by the adapter just like the factory context;
        !! associating it does not participate in point evaluation.
        type(gc_eqdsk_nonlocal_factory_t), target, intent(inout) :: factory
        real(dp), intent(in) :: h0, jperp
        type(gc_cylindrical_class_adapter_t), intent(out) :: adapter
        integer, intent(out) :: status

        real(dp) :: radius_lo, radius_hi
        integer :: local_status

        adapter = gc_cylindrical_class_adapter_t()
        call get_eqdsk_composite_cut_radius_bounds( &
            factory%certified_cut_atlas, radius_lo, radius_hi, local_status)
        if (local_status /= EQDSK_COMPOSITE_ATLAS_SUCCESS) then
            status = GC_CYL_CLASS_CUT_ERROR
            return
        end if
        call initialize_gc_cylindrical_class_adapter(factory%field, &
            factory%potential, h0, jperp, factory%species%mass_g, &
            factory%species%charge_esu, c, radius_lo, radius_hi, &
            physical_class_cut_map_callback, adapter, &
            status, options=factory%options%class_options, &
            splitter=certified_splitter_callback, user_data=factory, &
            allowed_region_provider=factory_allowed_region_provider, &
            allowed_region_verifier=factory_allowed_region_verifier, &
            allowed_region_certificate_id= &
                GC_EQDSK_ALLOWED_PROVIDER_CERTIFICATE_ID)
    end subroutine initialize_factory_class_adapter

    subroutine factory_allowed_region_provider(h0, jperp, sigma, user_data, &
            regions, status)
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: sigma
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        type(gc_cylindrical_allowed_region_set_t), intent(out) :: regions
        integer, intent(out) :: status

        integer :: local_status

        regions = gc_cylindrical_allowed_region_set_t()
        status = GC_CYL_CLASS_INTERIOR_INVALID
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
                factory%last_return_status = 2500 + sigma
                if (.not. factory%allowed_region_ready) return
                call build_gc_eqdsk_certified_allowed_regions( &
                    factory%allowed_region_context, h0, jperp, sigma, regions, &
                    local_status)
                if (local_status == GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) then
                    status = GC_CYL_CLASS_SUCCESS
                else
                    ! Preserve the typed provider failure for the factory
                    ! diagnostic; the adapter still exposes only its public
                    ! class-contract status to the caller.
                    factory%last_return_status = 2000 + local_status
                end if
            class default
                return
        end select
    end subroutine factory_allowed_region_provider

    subroutine factory_allowed_region_verifier(h0, jperp, sigma, rc_min, &
            rc_max, regions, user_data, certificate_id, status)
        real(dp), intent(in) :: h0, jperp, rc_min, rc_max
        integer, intent(in) :: sigma
        type(gc_cylindrical_allowed_region_set_t), intent(in) :: regions
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        integer, intent(out) :: certificate_id, status

        integer :: local_status

        certificate_id = 0
        status = GC_CYL_CLASS_SPLITTER_FAILURE
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
                factory%last_return_status = 3000 + sigma
                if (.not. factory%allowed_region_ready) return
                call verify_gc_eqdsk_certified_allowed_regions( &
                    factory%allowed_region_context, h0, jperp, sigma, rc_min, &
                    rc_max, regions, certificate_id, local_status)
                if (local_status == GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) then
                    status = GC_CYL_CLASS_SUCCESS
                else
                    certificate_id = 0
                    factory%last_return_status = 3000 + local_status
                end if
            class default
                return
        end select
    end subroutine factory_allowed_region_verifier

    subroutine factory_node_factory(h0, jperp, user_data, adapter, context, &
            status)
        real(dp), intent(in) :: h0, jperp
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        type(gc_cylindrical_class_adapter_t), intent(out) :: adapter
        type(gc_cylindrical_nonlocal_context_t), intent(out) :: context
        integer, intent(out) :: status

        integer :: local_status

        adapter = gc_cylindrical_class_adapter_t()
        context = gc_cylindrical_nonlocal_context_t()
        status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            if (.not. factory%field_ready .or. .not. factory%profile_ready) return
            if (.not. factory%wall_ready .or. .not. factory%cut_ready) return
            call initialize_factory_class_adapter(factory, h0, jperp, adapter, &
                local_status)
            if (local_status /= GC_CYL_CLASS_SUCCESS) then
                factory%last_return_status = 2100 + local_status
                status = GC_EQDSK_NONLOCAL_TOPOLOGY_UNAVAILABLE
                return
            end if
            call initialize_gc_cylindrical_nonlocal_provider(h0, jperp, context, &
                local_status, particle_charge=factory%species%charge_esu, &
                c_light=c, component_provider=factory_component_provider, &
                orbit_provider=factory_orbit_provider, &
                harmonic_provider=factory_harmonic_provider, &
                force_provider=factory_force_provider, &
                section_coordinate=factory%section_coordinate, &
                section_reference=factory%section_position, &
                section_reference_id=factory%section_reference_id, &
                required_return_crossings=2, &
                user_data=factory)
            if (local_status /= GC_CYL_NONLOCAL_SUCCESS) then
                factory%last_return_status = 2200 + local_status
                status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
                return
            end if
            factory%last_return_status = 2300
            status = GC_EQDSK_NONLOCAL_SUCCESS
        class default
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_node_factory

    subroutine factory_outer_factor(reference, h0, jperp, x, sigma, component_id, &
            launch, sample, user_data, outer_factor, status)
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_cylindrical_class_launch_t), intent(in) :: launch
        type(gc_nonlocal_orbit_sample_t), intent(in) :: sample
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: outer_factor
        integer, intent(out) :: status

        real(dp) :: density, temperature, potential, a1_star, a2_star
        type(gc_eqdsk_orbit_result_t) :: shell_result
        type(gc_full_fow_eq17_result_t) :: eq17
        integer :: local_status
        character(len=256) :: normalization_message

        outer_factor = 0.0_dp
        status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
        associate (unused_reference => reference, unused_jperp => jperp, &
                unused_x => x, unused_sigma => sigma, &
                unused_component_id => component_id)
        end associate
        if (.not. associated(user_data)) return
        if (sample%status /= GC_NONLOCAL_SAMPLE_VALID) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            if (.not. factory%initialized) return
            if (.not. factory%normalization%initialized) return
            call evaluate_gc_eqdsk_native_profile(factory, launch%psi_star, &
                density, temperature, potential, a1_star, a2_star, local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) return
            call integrate_physical_cycle(factory, launch, 0, 0, .false., .true., &
                0.0_dp, 0.0_dp, shell_result)
            if (shell_result%status /= GC_EQDSK_NONLOCAL_SUCCESS) return
            call evaluate_gc_full_fow_eq17(factory%normalization, h0, potential, &
                temperature, a1_star, a2_star, density, &
                shell_result%shell_average, eq17, local_status, &
                normalization_message)
            if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
            outer_factor = eq17%outer_factor
            status = GC_CYL_TRANSPORT_SUCCESS
        class default
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_outer_factor

    subroutine factory_class_kind(reference, h0, jperp, x, sigma, component_id, &
            user_data, class_kind, status)
        type(gc_nonlocal_transport_reference_t), intent(in) :: reference
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        integer, intent(out) :: class_kind, status

        type(gc_cylindrical_class_adapter_t) :: adapter
        type(gc_cylindrical_class_result_t) :: classes
        type(gc_cylindrical_class_launch_t) :: launch
        type(gc_eqdsk_orbit_result_t) :: orbit_result
        integer :: local_status, i
        logical :: found

        class_kind = 0
        status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            call initialize_factory_class_adapter(factory, h0, jperp, adapter, &
                local_status)
            if (local_status /= GC_CYL_CLASS_SUCCESS) return
            call enumerate_gc_cylindrical_classes(adapter, classes, local_status)
            if (local_status /= GC_CYL_CLASS_SUCCESS) return
            found = .false.
            do i = 1, classes%nclasses
                if (classes%classes(i)%component_id /= component_id) cycle
                if (classes%classes(i)%sigma /= sigma) cycle
                if (x < classes%classes(i)%rc_min .or. &
                        x > classes%classes(i)%rc_max) cycle
                found = .true.
                call launch_gc_cylindrical_class(adapter, x, sigma, &
                    component_id, launch, local_status)
                if (local_status /= GC_CYL_CLASS_SUCCESS) return
                call integrate_physical_cycle(factory, launch, 0, 0, .false., &
                    .false., 0.0_dp, 0.0_dp, orbit_result, &
                    track_behavior=.true.)
                if (orbit_result%status /= GC_EQDSK_NONLOCAL_SUCCESS .or. &
                        .not. orbit_result%complete_cycle_certified) return
                if (orbit_result%parallel_sign_changes > 0) then
                    class_kind = GC_NONLOCAL_CLASS_TRAPPED
                else if (launch%vparallel > 0.0_dp) then
                    class_kind = GC_NONLOCAL_CLASS_COPASSING
                else if (launch%vparallel < 0.0_dp) then
                    class_kind = GC_NONLOCAL_CLASS_COUNTERPASSING
                else
                    return
                end if
                exit
            end do
            if (.not. found) return
            status = GC_EQDSK_NONLOCAL_SUCCESS
        class default
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_class_kind

    subroutine factory_reset_evidence(user_data, status)
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        integer, intent(out) :: status

        status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            call reset_execution_evidence(factory)
            status = GC_EQDSK_NONLOCAL_SUCCESS
        class default
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_reset_evidence

    subroutine factory_get_evidence(user_data, evidence, status)
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        type(gc_nonlocal_transport_observed_evidence_t), intent(out) :: evidence
        integer, intent(out) :: status

        evidence = gc_nonlocal_transport_observed_evidence_t()
        status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            evidence%physical_return_attempts = factory%physical_return_attempts
            evidence%invariant_successes = factory%invariant_success_count
            evidence%invariant_failures = factory%invariant_failure_count
            evidence%return_successes = factory%physical_return_successes
            evidence%return_no_return = factory%no_return_count
            evidence%return_radial_domain = factory%radial_return_count
            evidence%return_wall_loss = factory%wall_return_count
            evidence%return_errors = factory%return_error_count
            evidence%wall_not_hit = factory%wall_not_hit_count
            evidence%wall_loss = factory%wall_return_count
            evidence%wall_errors = factory%wall_error_count
            evidence%topology_certification_attempts = &
                factory%topology_certification_attempts
            evidence%topology_certification_successes = &
                factory%topology_certification_successes
            evidence%harmonic_average_successes = factory%harmonic_average_successes
            evidence%residence_average_successes = factory%residence_average_successes
            evidence%return_accounting_complete = &
                evidence%physical_return_attempts == &
                evidence%invariant_successes+evidence%invariant_failures .and. &
                evidence%physical_return_attempts == &
                evidence%return_successes+evidence%return_no_return+ &
                evidence%return_radial_domain+evidence%return_wall_loss+ &
                evidence%return_errors .and. evidence%physical_return_attempts == &
                evidence%wall_not_hit+evidence%wall_loss+evidence%wall_errors
            status = GC_EQDSK_NONLOCAL_SUCCESS
        class default
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_get_evidence

    subroutine reset_execution_evidence(factory)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory

        factory%physical_return_attempts = 0
        factory%physical_return_successes = 0
        factory%wall_return_count = 0
        factory%radial_return_count = 0
        factory%no_return_count = 0
        factory%invariant_rejection_count = 0
        factory%field_error_count = 0
        factory%potential_error_count = 0
        factory%state_error_count = 0
        factory%start_error_count = 0
        factory%integrator_error_count = 0
        factory%other_return_error_count = 0
        factory%invariant_success_count = 0
        factory%invariant_failure_count = 0
        factory%return_error_count = 0
        factory%wall_not_hit_count = 0
        factory%wall_error_count = 0
        factory%topology_certification_attempts = 0
        factory%topology_certification_successes = 0
        factory%harmonic_average_successes = 0
        factory%residence_average_successes = 0
    end subroutine reset_execution_evidence

    subroutine factory_component_provider(h0, jperp, user_data, components, &
            status)
        real(dp), intent(in) :: h0, jperp
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        type(gc_nonlocal_component_t), allocatable, intent(out) :: components(:)
        integer, intent(out) :: status

        type(gc_cylindrical_class_adapter_t) :: adapter
        type(gc_cylindrical_class_result_t) :: classes
        integer :: local_status, i

        if (allocated(components)) deallocate(components)
        status = GC_CYL_NONLOCAL_CALLBACK_FAILURE
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            call initialize_factory_class_adapter(factory, h0, jperp, adapter, &
                local_status)
            if (local_status /= GC_CYL_CLASS_SUCCESS) then
                factory%last_return_status = 2400 + local_status
                return
            end if
            call enumerate_gc_cylindrical_classes(adapter, classes, local_status)
            if (local_status /= GC_CYL_CLASS_SUCCESS) then
                factory%last_return_status = 2500 + local_status
                return
            end if
            if (.not. classes%class_complete) then
                factory%last_return_status = 2600
                return
            end if
            if (classes%nclasses < 0) then
                factory%last_return_status = 2601
                return
            end if
            if (classes%nclasses == 0) then
                allocate(components(0))
                status = GC_CYL_NONLOCAL_SUCCESS
                return
            end if
            if (.not. allocated(classes%classes)) return
            allocate(components(classes%nclasses))
            do i = 1, classes%nclasses
                components(i)%component_id = classes%classes(i)%component_id
                components(i)%sigma = classes%classes(i)%sigma
                components(i)%x_min = classes%classes(i)%rc_min
                components(i)%x_max = classes%classes(i)%rc_max
            end do
            status = GC_CYL_NONLOCAL_SUCCESS
        class default
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_component_provider

    subroutine factory_orbit_provider(h0, jperp, x, sigma, component_id, &
            user_data, orbit, status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        type(gc_cylindrical_nonlocal_orbit_t), intent(out) :: orbit
        integer, intent(out) :: status

        type(gc_cylindrical_class_launch_t) :: launch
        type(gc_eqdsk_orbit_result_t) :: result
        integer :: local_status

        orbit = gc_cylindrical_nonlocal_orbit_t()
        status = GC_CYL_NONLOCAL_CALLBACK_FAILURE
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            call make_raw_launch(factory, h0, jperp, x, sigma, component_id, &
                launch, local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                status = GC_CYL_NONLOCAL_CALLBACK_FAILURE
                return
            end if
            call integrate_physical_cycle(factory, launch, 0, 0, .false., .false., &
                0.0_dp, 0.0_dp, result, track_behavior=.true.)
            call fill_orbit_metadata(factory, h0, jperp, x, sigma, component_id, &
                launch, result, orbit, local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                status = GC_CYL_NONLOCAL_SUCCESS
                return
            end if
            if (orbit%status == GC_CYL_NONLOCAL_ORBIT_WALL .or. &
                orbit%status == GC_CYL_NONLOCAL_ORBIT_UNRESOLVED) then
                status = GC_CYL_NONLOCAL_SUCCESS
                return
            end if
            call compute_orbit_derivatives(factory, h0, jperp, x, sigma, &
                component_id, result, orbit, local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                orbit%status = GC_CYL_NONLOCAL_ORBIT_UNRESOLVED
                status = GC_CYL_NONLOCAL_SUCCESS
                return
            end if
            status = GC_CYL_NONLOCAL_SUCCESS
        class default
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_orbit_provider

    subroutine factory_harmonic_provider(h0, jperp, x, sigma, component_id, &
            harmonic_m, harmonic_n, orbit, user_data, h_m, status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id, harmonic_m, harmonic_n
        type(gc_cylindrical_nonlocal_orbit_t), intent(in) :: orbit
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        complex(dp), intent(out) :: h_m
        integer, intent(out) :: status

        type(gc_cylindrical_class_launch_t) :: launch
        type(gc_eqdsk_orbit_result_t) :: result
        integer :: local_status

        h_m = cmplx(0.0_dp, 0.0_dp, kind=dp)
        status = GC_CYL_NONLOCAL_CALLBACK_FAILURE
        if (.not. associated(user_data)) return
        if (orbit%status /= GC_CYL_NONLOCAL_ORBIT_VALID) return
        if (orbit%tau_b <= 0.0_dp) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            call make_raw_launch(factory, h0, jperp, x, sigma, component_id, &
                launch, local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) return
            call integrate_physical_cycle(factory, launch, harmonic_m, harmonic_n, &
                .true., .false., orbit%omega_b, &
                orbit%omega_phi, result)
            if (result%status /= GC_EQDSK_NONLOCAL_SUCCESS) return
            if (abs(result%period-orbit%tau_b) > &
                5.0e-8_dp*max(1.0_dp, orbit%tau_b)) return
            h_m = result%h_m
            if (.not. all(ieee_is_finite([real(h_m, dp), aimag(h_m)]))) then
                h_m = cmplx(0.0_dp, 0.0_dp, kind=dp)
                status = GC_CYL_NONLOCAL_CALLBACK_FAILURE
                return
            end if
            status = GC_CYL_NONLOCAL_SUCCESS
        class default
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_harmonic_provider

    subroutine factory_force_provider(h0, jperp, x, sigma, component_id, orbit, &
            user_data, force, status)
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_cylindrical_nonlocal_orbit_t), intent(in) :: orbit
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: force(:)
        integer, intent(out) :: status

        real(dp) :: density, temperature, potential, a1_star, a2_star
        real(dp) :: psi_star, dpsi_star_dp_phi
        type(gc_full_fow_eq17_result_t) :: eq17
        integer :: local_status
        character(len=256) :: normalization_message

        force = 0.0_dp
        status = GC_CYL_NONLOCAL_CALLBACK_FAILURE
        associate (unused_jperp => jperp, unused_x => x, &
                unused_sigma => sigma, unused_component_id => component_id)
        end associate
        if (size(force) < FACTORY_FORCE_COUNT) return
        if (.not. associated(user_data)) return
        if (orbit%status /= GC_CYL_NONLOCAL_ORBIT_VALID) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            if (.not. factory%profile_ready) return
            if (.not. factory%normalization%initialized) return
            call evaluate_gc_full_fow_canonical_flux(factory%normalization, &
                orbit%p_phi, psi_star, dpsi_star_dp_phi, local_status, &
                normalization_message)
            if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
            call evaluate_gc_eqdsk_native_profile(factory, psi_star, density, &
                temperature, potential, a1_star, a2_star, local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) return
            call evaluate_gc_full_fow_eq17(factory%normalization, h0, potential, &
                temperature, a1_star, a2_star, density, 0.0_dp, eq17, &
                local_status, normalization_message)
            if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
            force(1) = a1_star
            force(2) = a2_star
            force(3) = eq17%force_bracket
            status = GC_CYL_NONLOCAL_SUCCESS
        class default
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_force_provider

    subroutine physical_class_cut_map_callback(radius, user_data, position, &
            dposition_dradius, status)
        real(dp), intent(in) :: radius
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: position(3), dposition_dradius(3)
        integer, intent(out) :: status

        integer :: local_status

        position = 0.0_dp
        dposition_dradius = 0.0_dp
        status = GC_CYL_EQUILIBRIUM_DOMAIN
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            call map_eqdsk_composite_cut_atlas_radius( &
                factory%certified_cut_atlas, radius, position, &
                dposition_dradius, local_status)
            if (local_status /= EQDSK_COMPOSITE_ATLAS_SUCCESS) then
                write (*, '(a,1x,i0,1x,es24.16)') &
                    'composite cut diagnostic=', local_status, radius
                return
            end if
            status = GC_CYL_SUCCESS
        class default
            status = GC_CYL_EQUILIBRIUM_DOMAIN
        end select
    end subroutine physical_class_cut_map_callback

    subroutine certified_splitter_callback(h0, jperp, sigma, candidate, user_data, &
            split_classes, certified, status)
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: sigma
        type(gc_cylindrical_class_interval_t), intent(in) :: candidate
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        type(gc_cylindrical_class_interval_t), allocatable, intent(out) :: &
            split_classes(:)
        logical, intent(out) :: certified
        integer, intent(out) :: status

        if (allocated(split_classes)) deallocate(split_classes)
        certified = .false.
        status = GC_CYL_CLASS_SUCCESS
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            call certify_homoclinic_component(factory, h0, jperp, sigma, candidate, &
                certified, status)
            if (status /= GC_CYL_CLASS_SUCCESS) return
            if (certified) then
                allocate(split_classes(1))
                split_classes(1) = candidate
                split_classes(1)%topology_certified = .true.
                split_classes(1)%orbit_return_certified = .true.
            end if
        class default
            status = GC_CYL_CLASS_SUCCESS
        end select
    end subroutine certified_splitter_callback

    subroutine certify_homoclinic_component(factory, h0, jperp, sigma, candidate, &
            certified, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: h0, jperp
        integer, intent(in) :: sigma
        type(gc_cylindrical_class_interval_t), intent(in) :: candidate
        logical, intent(out) :: certified
        integer, intent(out) :: status

        type(gc_cylindrical_allowed_region_set_t) :: same_sign_regions
        type(gc_cylindrical_allowed_region_set_t) :: opposite_sign_regions
        integer :: local_status

        certified = .false.
        status = GC_CYL_CLASS_SPLITTER_FAILURE
        factory%last_return_status = 4000 + sigma
        factory%topology_certification_attempts = &
            factory%topology_certification_attempts + 1
        factory%topology_ready = .false.

        ! The regular chart has a different theorem from a turning chart.  A
        ! component ending at v_parallel=0 needs the one-sided Puiseux
        ! certificate; it must not enter this path through a finite sample.
        if (candidate%lower_root .or. candidate%upper_root) return
        if (abs(sigma) /= 1) return

        ! For fixed (H0,J_K), P_phi=(q/c) psi_star is constant on an orbit.
        ! The generated interval evaluator has already excluded zero from
        ! dpsi_star/dR on every regular component.  The provider call below
        ! supplies the complete root-isolated decomposition for both signs;
        ! requiring pairwise disjoint canonical ranges then proves that a
        ! fixed P_phi has at most one cut point for each sign of v_parallel.
        ! The physical integrator still has to observe the two transverse
        ! crossings and the return-state checks before a sample is usable.
        call build_gc_eqdsk_certified_allowed_regions( &
            factory%allowed_region_context, h0, jperp, sigma, &
            same_sign_regions, local_status)
        if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) return
        call build_gc_eqdsk_certified_allowed_regions( &
            factory%allowed_region_context, h0, jperp, -sigma, &
            opposite_sign_regions, local_status)
        if (local_status /= GC_EQDSK_ALLOWED_PROVIDER_SUCCESS) return
        if (.not. regular_canonical_ranges_certified(same_sign_regions)) return
        if (.not. regular_canonical_ranges_certified(opposite_sign_regions)) return
        if (.not. candidate_is_provider_component(candidate, same_sign_regions)) return

        ! This is an independent topology certificate, not a claim made by
        ! the ODE event sequence.  Its ID is carried into the physical-return
        ! options so that the integrator cannot manufacture multiplicity.
        factory%cut_atlas%two_cut_multiplicity_certified = .true.
        factory%two_cut_multiplicity_certificate_id = &
            GC_EQDSK_TWO_CUT_MULTIPLICITY_CERTIFICATE_ID
        factory%topology_certification_successes = &
            factory%topology_certification_successes + 1
        factory%topology_ready = .true.
        certified = .true.
        status = GC_CYL_CLASS_SUCCESS
    end subroutine certify_homoclinic_component

    logical function regular_canonical_ranges_certified(regions)
        type(gc_cylindrical_allowed_region_set_t), intent(in) :: regions

        integer :: i, j
        real(dp) :: lower_i, upper_i, lower_j, upper_j
        type(gc_cylindrical_allowed_component_t) :: component_i
        type(gc_cylindrical_allowed_component_t) :: component_j

        regular_canonical_ranges_certified = .false.
        if (.not. regions%topology_certified) return
        if (regions%nroots /= 0) return
        if (regions%ncomponents < 1) return
        if (.not. allocated(regions%components)) return
        if (size(regions%components) /= regions%ncomponents) return
        do i = 1, regions%ncomponents
            component_i = regions%components(i)
            if (component_i%lower_root .or. component_i%upper_root) return
            if (.not. component_i%canonical_measure_certified) return
            if (component_i%canonical_measure <= 0.0_dp) return
            if (.not. all(ieee_is_finite([component_i%canonical_begin, &
                    component_i%canonical_end, &
                    component_i%canonical_measure]))) return
            lower_i = min(component_i%canonical_begin, component_i%canonical_end)
            upper_i = max(component_i%canonical_begin, component_i%canonical_end)
            if (upper_i <= lower_i) return
            do j = 1, i - 1
                component_j = regions%components(j)
                lower_j = min(component_j%canonical_begin, &
                    component_j%canonical_end)
                upper_j = max(component_j%canonical_begin, &
                    component_j%canonical_end)
                if (max(lower_i, lower_j) < min(upper_i, upper_j)) return
            end do
        end do
        regular_canonical_ranges_certified = .true.
    end function regular_canonical_ranges_certified

    logical function candidate_is_provider_component(candidate, regions)
        type(gc_cylindrical_class_interval_t), intent(in) :: candidate
        type(gc_cylindrical_allowed_region_set_t), intent(in) :: regions

        integer :: i

        candidate_is_provider_component = .false.
        if (.not. allocated(regions%components)) return
        do i = 1, size(regions%components)
            if (regions%components(i)%component_id /= candidate%component_id) cycle
            if (regions%components(i)%sigma /= candidate%sigma) cycle
            if (regions%components(i)%x_begin /= candidate%rc_min) cycle
            if (regions%components(i)%x_end /= candidate%rc_max) cycle
            candidate_is_provider_component = .true.
            return
        end do
    end function candidate_is_provider_component

    subroutine make_raw_launch(factory, h0, jperp, x, sigma, component_id, &
            launch, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_cylindrical_class_launch_t), intent(out) :: launch
        integer, intent(out) :: status

        type(gc_cylindrical_class_adapter_t) :: adapter
        type(gc_cylindrical_class_point_t) :: point
        integer :: local_status
        real(dp) :: magnetic_moment, energy, p_phi, canonical

        launch = gc_cylindrical_class_launch_t()
        launch%h0 = h0
        launch%jperp = jperp
        launch%rc = x
        launch%sigma = sigma
        launch%component_id = component_id
        status = GC_EQDSK_NONLOCAL_ORBIT_ERROR
        call initialize_factory_class_adapter(factory, h0, jperp, adapter, &
            local_status)
        if (local_status /= GC_CYL_CLASS_SUCCESS) return
        call evaluate_gc_cylindrical_class_point(adapter, x, sigma, point, &
            local_status)
        if (local_status /= GC_CYL_CLASS_SUCCESS .or. .not. point%allowed) return
        magnetic_moment = gc_mu_phys_from_buchholz_jk(jperp, &
            factory%species%mass_g, factory%species%charge_esu, c)
        launch%position = point%position
        launch%vparallel_squared = max(0.0_dp, point%vparallel_squared)
        launch%vparallel = point%vparallel
        launch%omega_c = point%omega_c
        launch%psi_star = point%psi_star
        launch%dpsi_star_drc = point%dpsi_star_drc
        launch%derivative_available = point%derivative_available
        launch%endpoint_tangent = point%at_turning_point
        launch%state%R = point%position(1)
        launch%state%Z = point%position(2)
        launch%state%phi = point%position(3)
        launch%state%p_parallel = factory%species%mass_g*point%vparallel
        launch%state%mu = magnetic_moment
        p_phi = factory%species%charge_esu/c*point%psi_star
        canonical = factory%species%charge_esu/c*point%field%psi &
            +launch%state%p_parallel*point%position(1)*point%field%bhat(2)
        energy = launch%state%p_parallel**2/(2.0_dp*factory%species%mass_g) &
            +magnetic_moment*point%field%bmod &
            +factory%species%charge_esu*point%potential
        launch%p_phi = p_phi
        launch%energy_residual = energy-h0
        launch%jperp_residual = gc_buchholz_jk_from_mu_phys(magnetic_moment, &
            factory%species%mass_g, factory%species%charge_esu, c)-jperp
        launch%canonical_residual = canonical-p_phi
        if (.not. all(ieee_is_finite([launch%position, launch%state%p_parallel, &
            launch%state%mu, launch%p_phi, launch%energy_residual, &
            launch%jperp_residual, launch%canonical_residual]))) return
        launch%class_certified = .false.
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine make_raw_launch

    subroutine fill_orbit_metadata(factory, h0, jperp, x, sigma, component_id, &
            launch, result, orbit, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_cylindrical_class_launch_t), intent(in) :: launch
        type(gc_eqdsk_orbit_result_t), intent(in) :: result
        type(gc_cylindrical_nonlocal_orbit_t), intent(out) :: orbit
        integer, intent(out) :: status
        real(dp) :: unused_omega_prec, unused_domega_b, unused_domega_phi
        real(dp) :: unused_domega_prec

        orbit = gc_cylindrical_nonlocal_orbit_t()
        orbit%component_id = component_id
        orbit%sigma = sigma
        orbit%section%coordinate = factory%section_coordinate
        orbit%section%reference = factory%section_position
        orbit%section%reference_id = factory%section_reference_id
        orbit%section%locked = factory%cut_ready
        orbit%section%required_return_crossings = 2
        orbit%section%return_crossings = 2
        orbit%p_phi = launch%p_phi
        orbit%dp_phi_dx = factory%species%charge_esu/c*launch%dpsi_star_drc
        orbit%wall_checked = .true.
        orbit%wall_status = GC_CYL_NONLOCAL_WALL_CLEAR
        orbit%winding = 1
        orbit%section_return_crossings = result%intersection_count
        orbit%intersection_orientations = result%intersection_orientations
        orbit%intersection_times = result%intersection_times
        orbit%intersection_rates = result%intersection_rates
        orbit%intersection_multiplicity_certified = &
            result%intersection_multiplicity_certified
        orbit%winding_available = .true.
        orbit%section_return_available = .true.
        orbit%p_phi_mapping_certified = .false.
        orbit%mapping_orientation = sign_of_real(orbit%dp_phi_dx)
        orbit%electric_potential_included = factory%potential%initialized
        orbit%delta_phi_b = result%delta_phi
        orbit%parallel_sign_changes = result%parallel_sign_changes
        orbit%class_behavior_certified = result%complete_cycle_certified
        if (result%parallel_sign_changes > 0) then
            orbit%class_kind = GC_NONLOCAL_CLASS_TRAPPED
        else if (launch%vparallel > 0.0_dp) then
            orbit%class_kind = GC_NONLOCAL_CLASS_COPASSING
        else if (launch%vparallel < 0.0_dp) then
            orbit%class_kind = GC_NONLOCAL_CLASS_COUNTERPASSING
        else
            orbit%class_behavior_certified = .false.
        end if
        orbit%phase_advance_available = .false.
        orbit%period_derivative_available = .false.
        orbit%derivatives_available = .false.
        status = GC_EQDSK_NONLOCAL_ORBIT_ERROR
        associate (unused_h0 => h0, unused_jperp => jperp, unused_x => x)
        end associate
        if (result%wall_hit) then
            orbit%status = GC_CYL_NONLOCAL_ORBIT_WALL
            orbit%wall_status = GC_CYL_NONLOCAL_WALL_HIT
            status = GC_EQDSK_NONLOCAL_SUCCESS
            return
        end if
        if (result%status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            orbit%status = GC_CYL_NONLOCAL_ORBIT_UNRESOLVED
            status = GC_EQDSK_NONLOCAL_SUCCESS
            return
        end if
        if (.not. physical_return_has_identity(factory, launch, result)) then
            orbit%status = GC_CYL_NONLOCAL_ORBIT_UNRESOLVED
            status = GC_EQDSK_NONLOCAL_SUCCESS
            return
        end if
        if (.not. ieee_is_finite(result%period) .or. result%period <= 0.0_dp) then
            orbit%status = GC_CYL_NONLOCAL_ORBIT_UNRESOLVED
            status = GC_EQDSK_NONLOCAL_SUCCESS
            return
        end if
        orbit%tau_b = result%period
        ! Direct real-space FOW has no Boozer q reduction.  The explicit +1
        ! first-return winding gives the positive bounce frequency and keeps
        ! delta_phi/period signed; precession is intentionally not stored in
        ! this nonlocal sample.
        call evaluate_neort_full_fow_cycle_frequency(result%period, 1.0_dp, &
            result%delta_phi, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            orbit%omega_b, orbit%omega_phi, unused_omega_prec, &
            unused_domega_b, unused_domega_phi, unused_domega_prec)
        orbit%delta_phi_b = result%delta_phi
        orbit%complete_cycle_return = .true.
        orbit%p_phi_mapping_certified = launch%derivative_available
        orbit%phase_advance_available = .true.
        orbit%status = GC_CYL_NONLOCAL_ORBIT_VALID
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine fill_orbit_metadata

    subroutine compute_orbit_derivatives(factory, h0, jperp, x, sigma, &
            component_id, base, orbit, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: h0, jperp, x
        integer, intent(in) :: sigma, component_id
        type(gc_eqdsk_orbit_result_t), intent(in) :: base
        type(gc_cylindrical_nonlocal_orbit_t), intent(inout) :: orbit
        integer, intent(out) :: status

        type(gc_cylindrical_class_launch_t) :: minus_launch, plus_launch
        type(gc_eqdsk_orbit_result_t) :: minus_result, plus_result
        real(dp) :: lower, upper, step, x_minus, x_plus, denominator
        integer :: local_status

        status = GC_EQDSK_NONLOCAL_DERIVATIVE_UNAVAILABLE
        call find_cached_component_bounds(factory, sigma, component_id, lower, &
            upper, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        step = factory%options%orbit_derivative_fraction*max(upper-lower, &
            max(1.0_dp, abs(x)))
        step = min(step, 0.25_dp*(upper-lower))
        if (step <= 0.0_dp) return
        x_minus = max(lower, x-step)
        x_plus = min(upper, x+step)
        denominator = x_plus-x_minus
        if (denominator <= 0.0_dp) return
        call make_raw_launch(factory, h0, jperp, x_minus, sigma, component_id, &
            minus_launch, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        call integrate_physical_cycle(factory, minus_launch, 0, 0, .false., &
            .false., 0.0_dp, 0.0_dp, minus_result)
        if (.not. physical_return_has_identity(factory, minus_launch, &
            minus_result)) return
        call make_raw_launch(factory, h0, jperp, x_plus, sigma, component_id, &
            plus_launch, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        call integrate_physical_cycle(factory, plus_launch, 0, 0, .false., &
            .false., 0.0_dp, 0.0_dp, plus_result)
        if (.not. physical_return_has_identity(factory, plus_launch, &
            plus_result)) return
        orbit%dtau_b_dx = (plus_result%period-minus_result%period)/denominator
        orbit%ddelta_phi_b_dx = &
            (plus_result%delta_phi-minus_result%delta_phi)/denominator
        orbit%domega_b_dx = -2.0_dp*pi*orbit%dtau_b_dx/base%period**2
        orbit%domega_phi_dx = (orbit%ddelta_phi_b_dx*base%period &
            -base%delta_phi*orbit%dtau_b_dx)/base%period**2
        if (.not. all(ieee_is_finite([orbit%dtau_b_dx, &
            orbit%ddelta_phi_b_dx, orbit%domega_b_dx, &
            orbit%domega_phi_dx]))) return
        orbit%period_derivative_available = .true.
        orbit%derivatives_available = .true.
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine compute_orbit_derivatives

    subroutine find_cached_component_bounds(factory, sigma, component_id, lower, &
            upper, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        integer, intent(in) :: sigma, component_id
        real(dp), intent(out) :: lower, upper
        integer, intent(out) :: status

        integer :: i

        lower = 0.0_dp
        upper = 0.0_dp
        status = GC_EQDSK_NONLOCAL_DERIVATIVE_UNAVAILABLE
        if (.not. allocated(factory%provider%node_class_result%classes)) return
        do i = 1, size(factory%provider%node_class_result%classes)
            if (factory%provider%node_class_result%classes(i)%sigma /= sigma) cycle
            if (factory%provider%node_class_result%classes(i)%component_id /= &
                component_id) cycle
            lower = factory%provider%node_class_result%classes(i)%rc_min
            upper = factory%provider%node_class_result%classes(i)%rc_max
            if (upper <= lower) return
            status = GC_EQDSK_NONLOCAL_SUCCESS
            return
        end do
    end subroutine find_cached_component_bounds

    subroutine integrate_physical_cycle(factory, launch, harmonic_m, harmonic_n, &
            include_harmonic, include_shell, omega_b, omega_phi, result, &
            track_behavior, refinement_pass, maximum_step_override)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        type(gc_cylindrical_class_launch_t), intent(in) :: launch
        integer, intent(in) :: harmonic_m, harmonic_n
        logical, intent(in) :: include_harmonic, include_shell
        real(dp), intent(in) :: omega_b, omega_phi
        type(gc_eqdsk_orbit_result_t), intent(out) :: result
        logical, intent(in), optional :: track_behavior
        logical, intent(in), optional :: refinement_pass
        real(dp), intent(in), optional :: maximum_step_override

        type(gc_cylindrical_invariants_t) :: invariants
        type(gc_cylindrical_physical_return_options_t) :: return_options
        type(gc_cylindrical_physical_return_t) :: physical_return
        type(gc_cylindrical_physical_return_certificate_t) :: &
            multiplicity_certificate
        integer :: local_status
        character(len=256) :: certificate_message
        logical :: behavior_requested, is_refinement_pass
        real(dp) :: base_step, refined_step, tolerance_factor
        real(dp) :: pass_omega_b, pass_omega_phi
        type(gc_eqdsk_orbit_result_t) :: refined_result
        real(dp) :: saved_base_step, saved_refined_step
        real(dp) :: saved_period_error, saved_delta_phi_error
        real(dp) :: saved_omega_b_error, saved_omega_phi_error
        real(dp) :: saved_h_m_error, saved_shell_error
        real(dp) :: unused_omega_prec, unused_domega_b, unused_domega_phi
        real(dp) :: unused_domega_prec

        result = gc_eqdsk_orbit_result_t()
        behavior_requested = present(track_behavior) .and. track_behavior
        is_refinement_pass = present(refinement_pass) .and. refinement_pass
        if (present(maximum_step_override)) then
            base_step = maximum_step_override
        else if (factory%options%orbit_maximum_step > 0.0_dp) then
            base_step = factory%options%orbit_maximum_step
        else
            base_step = 2.0_dp*pi*factory%minor_radius_cm/(256.0_dp* &
                factory%species%reference_velocity_cm_s)
        end if
        if (.not. ieee_is_finite(base_step) .or. base_step <= 0.0_dp) then
            result%status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
            return
        end if
        tolerance_factor = 1.0_dp
        if (is_refinement_pass) tolerance_factor = 0.25_dp
        if (.not. factory%field_ready) return
        if (.not. factory%profile_ready) return
        if (.not. factory%wall_ready) return
        if (.not. factory%cut_ready) return
        factory%physical_return_attempts = factory%physical_return_attempts + 1
        if (.not. launch_state_is_consistent(factory, launch)) then
            factory%invariant_failure_count = factory%invariant_failure_count + 1
            factory%return_error_count = factory%return_error_count + 1
            factory%wall_error_count = factory%wall_error_count + 1
            return
        end if
        factory%invariant_success_count = factory%invariant_success_count + 1
        invariants%energy = launch%h0
        invariants%magnetic_moment = gc_mu_phys_from_buchholz_jk(launch%jperp, &
            factory%species%mass_g, factory%species%charge_esu, c)
        invariants%canonical_toroidal_momentum = launch%p_phi
        call make_physical_return_options(factory, invariants, return_options, &
            base_step, tolerance_factor)
        call compute_gc_cylindrical_physical_return(factory%field, &
            factory%potential, launch%state, invariants, &
            factory%species%mass_g, factory%species%charge_esu, c, &
            factory_physical_return_event, return_options, physical_return, &
            wall_model=factory%wall, radial_domain=factory_radial_domain_event, &
            user_data=factory, return_event_rate=factory_physical_return_rate)
        if (factory%cut_atlas%two_cut_multiplicity_certified .and. &
                factory%two_cut_multiplicity_certificate_id == &
                    GC_EQDSK_TWO_CUT_MULTIPLICITY_CERTIFICATE_ID) then
            multiplicity_certificate = &
                gc_cylindrical_physical_return_certificate_t()
            multiplicity_certificate%certificate_id = &
                factory%two_cut_multiplicity_certificate_id
            multiplicity_certificate%crossing_count = 2
            multiplicity_certificate%exactly_two_proved = .true.
            call attach_gc_cylindrical_physical_return_certificate( &
                physical_return, multiplicity_certificate, local_status, &
                certificate_message)
            if (local_status /= GC_CYL_SUCCESS) then
                physical_return%status = GC_CYL_PHYSICAL_EVENT_CALLBACK_ERROR
            end if
        end if
        call record_physical_return(factory, physical_return)
        result%event_kind = physical_return%event_kind
        result%intersection_count = physical_return%intersection_count
        result%intersection_orientations = &
            physical_return%intersection_orientations
        result%intersection_times = physical_return%intersection_times
        result%intersection_rates = physical_return%intersection_rates
        result%intersection_multiplicity_certified = &
            physical_return%intersection_multiplicity_certified
        result%return_orientation = physical_return%return_orientation
        result%wall_hit = physical_return%wall_hit
        result%radial_domain_hit = physical_return%radial_domain_hit
        result%physical_return_found = physical_return%physical_return_found
        result%state_at_return = physical_return%state_at_event
        result%energy_error = physical_return%energy_error
        result%magnetic_moment_error = physical_return%magnetic_moment_error
        result%canonical_momentum_error = physical_return%canonical_momentum_error
        if (physical_return%event_kind == GC_CYL_PHYSICAL_EVENT_WALL) then
            result%status = GC_CYL_WALL_LOSS
            return
        end if
        if (physical_return%event_kind == GC_CYL_PHYSICAL_EVENT_RADIAL_DOMAIN) then
            return
        end if
        if (physical_return%status == GC_CYL_NO_RETURN) then
            return
        end if
        if (physical_return%status /= GC_CYL_SUCCESS) return
        if (physical_return%event_kind /= GC_CYL_PHYSICAL_EVENT_RETURN) return
        if (.not. physical_return%physical_return_found) return
        if (return_options%require_opposite_intersection) then
            if (physical_return%intersection_count /= 2 .or. &
                    physical_return%intersection_orientations(1) == &
                    physical_return%intersection_orientations(2)) then
                result%status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
                return
            end if
        end if
        if (.not. ieee_is_finite(physical_return%period) .or. &
                physical_return%period <= 0.0_dp) return
        result%period = physical_return%period
        result%delta_phi = physical_return%delta_phi
        call evaluate_neort_full_fow_cycle_frequency(result%period, 1.0_dp, &
            result%delta_phi, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            result%omega_b, result%omega_phi, unused_omega_prec, &
            unused_domega_b, unused_domega_phi, unused_domega_prec)
        pass_omega_b = result%omega_b
        pass_omega_phi = result%omega_phi
        if (.not. all(ieee_is_finite([pass_omega_b, pass_omega_phi]))) then
            result%status = GC_EQDSK_NONLOCAL_NONFINITE
            return
        end if
        if (abs(omega_b) > tiny(omega_b)) then
            if (.not. refinement_observable_close(omega_b, pass_omega_b, &
                    1.0e-10_dp, 1.0e-10_dp, &
                    factory%refinement_frequency_scale)) then
                result%status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
                return
            end if
        end if
        if (abs(omega_phi) > tiny(omega_phi)) then
            if (.not. refinement_observable_close(omega_phi, pass_omega_phi, &
                    1.0e-10_dp, 1.0e-10_dp, &
                    factory%refinement_frequency_scale)) then
                result%status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
                return
            end if
        end if
        if (include_harmonic .or. include_shell .or. behavior_requested) then
            call integrate_cycle_averages(factory, launch, result%period, &
                harmonic_m, harmonic_n, include_harmonic, include_shell, &
                pass_omega_b, pass_omega_phi, result, local_status, &
                behavior_requested, base_step, tolerance_factor)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                result%status = local_status
                return
            end if
        end if
        if (behavior_requested) result%complete_cycle_certified = .true.
        if (factory%options%require_step_refinement .and. .not. is_refinement_pass) then
            refined_step = 0.5_dp*base_step
            if (.not. ieee_is_finite(base_step) .or. base_step <= 0.0_dp .or. &
                    .not. ieee_is_finite(refined_step) .or. refined_step <= 0.0_dp) then
                result%status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
                return
            end if
            call integrate_physical_cycle(factory, launch, harmonic_m, harmonic_n, &
                include_harmonic, include_shell, omega_b, omega_phi, refined_result, &
                track_behavior, refinement_pass=.true., &
                maximum_step_override=refined_step)
            result%base_step = base_step
            result%refined_step = refined_step
            call certify_step_refinement(result, refined_result, factory, &
                local_status)
            if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                result%status = local_status
                result%step_refinement_certified = .false.
                return
            end if
            saved_base_step = result%base_step
            saved_refined_step = result%refined_step
            saved_period_error = result%period_refinement_error
            saved_delta_phi_error = result%delta_phi_refinement_error
            saved_omega_b_error = result%omega_b_refinement_error
            saved_omega_phi_error = result%omega_phi_refinement_error
            saved_h_m_error = result%h_m_refinement_error
            saved_shell_error = result%shell_refinement_error
            result = refined_result
            result%step_refinement_certified = .true.
            result%base_step = saved_base_step
            result%refined_step = saved_refined_step
            result%period_refinement_error = saved_period_error
            result%delta_phi_refinement_error = saved_delta_phi_error
            result%omega_b_refinement_error = saved_omega_b_error
            result%omega_phi_refinement_error = saved_omega_phi_error
            result%h_m_refinement_error = saved_h_m_error
            result%shell_refinement_error = saved_shell_error
            factory%step_refinement_certified = .true.
            factory%base_step = base_step
            factory%refined_step = refined_step
            factory%period_refinement_error = result%period_refinement_error
            factory%delta_phi_refinement_error = result%delta_phi_refinement_error
            factory%omega_b_refinement_error = result%omega_b_refinement_error
            factory%omega_phi_refinement_error = result%omega_phi_refinement_error
            factory%h_m_refinement_error = result%h_m_refinement_error
            factory%shell_refinement_error = result%shell_refinement_error
        end if
        result%status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine integrate_physical_cycle

    subroutine certify_step_refinement(coarse, refined, factory, status)
        type(gc_eqdsk_orbit_result_t), intent(inout) :: coarse
        type(gc_eqdsk_orbit_result_t), intent(in) :: refined
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        integer, intent(out) :: status

        status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
        coarse%step_refinement_certified = .false.
        if (refined%status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        if (.not. refined%physical_return_found) return
        if (refined%event_kind /= coarse%event_kind) return
        if (refined%return_orientation /= coarse%return_orientation) return
        if (refined%intersection_count /= coarse%intersection_count) return
        if (refined%intersection_multiplicity_certified .neqv. &
                coarse%intersection_multiplicity_certified) return
        if (any(refined%intersection_orientations /= &
                coarse%intersection_orientations)) return
        if (refined%wall_hit .neqv. coarse%wall_hit) return
        if (refined%radial_domain_hit .neqv. coarse%radial_domain_hit) return
        if (refined%complete_cycle_certified .neqv. &
                coarse%complete_cycle_certified) return
        if (.not. all(ieee_is_finite([factory%refinement_length_scale, &
                factory%refinement_time_scale, &
                factory%refinement_frequency_scale, &
                factory%refinement_momentum_scale, factory%refinement_mu_scale]))) &
            return
        if (min(factory%refinement_length_scale, factory%refinement_time_scale, &
                factory%refinement_frequency_scale, &
                factory%refinement_momentum_scale, factory%refinement_mu_scale) &
                <= 0.0_dp) return

        coarse%period_refinement_error = abs(refined%period-coarse%period)
        coarse%delta_phi_refinement_error = &
            abs(refined%delta_phi-coarse%delta_phi)
        coarse%omega_b_refinement_error = abs(refined%omega_b-coarse%omega_b)
        coarse%omega_phi_refinement_error = &
            abs(refined%omega_phi-coarse%omega_phi)
        coarse%h_m_refinement_error = max( &
            abs(real(refined%h_m-coarse%h_m, dp)), &
            abs(aimag(refined%h_m-coarse%h_m)))
        coarse%shell_refinement_error = &
            abs(refined%shell_average-coarse%shell_average)

        if (.not. refinement_observable_close(coarse%period, refined%period, &
                factory%options%step_refinement_absolute_tolerance, &
                factory%options%step_refinement_relative_tolerance, &
                factory%refinement_time_scale)) return
        if (.not. refinement_observable_close(coarse%delta_phi, &
                refined%delta_phi, &
                factory%options%step_refinement_absolute_tolerance, &
                factory%options%step_refinement_relative_tolerance, &
                2.0_dp*pi)) return
        if (.not. refinement_observable_close(coarse%omega_b, refined%omega_b, &
                factory%options%step_refinement_absolute_tolerance, &
                factory%options%step_refinement_relative_tolerance, &
                factory%refinement_frequency_scale)) return
        if (.not. refinement_observable_close(coarse%omega_phi, &
                refined%omega_phi, &
                factory%options%step_refinement_absolute_tolerance, &
                factory%options%step_refinement_relative_tolerance, &
                factory%refinement_frequency_scale)) return
        if (.not. refinement_observable_close(real(coarse%h_m, dp), &
                real(refined%h_m, dp), &
                factory%options%step_refinement_absolute_tolerance, &
                factory%options%step_refinement_relative_tolerance, &
                factory%normalization%e_ref)) return
        if (.not. refinement_observable_close(aimag(coarse%h_m), &
                aimag(refined%h_m), &
                factory%options%step_refinement_absolute_tolerance, &
                factory%options%step_refinement_relative_tolerance, &
                factory%normalization%e_ref)) return
        if (.not. refinement_observable_close(coarse%shell_average, &
                refined%shell_average, &
                factory%options%step_refinement_absolute_tolerance, &
                factory%options%step_refinement_relative_tolerance, 1.0_dp)) return
        if (.not. refinement_state_close(coarse%state_at_return, &
                refined%state_at_return, factory)) return
        coarse%step_refinement_certified = .true.
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine certify_step_refinement

    logical function refinement_observable_close(first, second, atol, rtol, &
            reference_scale)
        real(dp), intent(in) :: first, second, atol, rtol, reference_scale

        real(dp) :: first_scaled, second_scaled, error_scaled, tolerance
        integer :: local_status
        character(len=256) :: normalization_message

        refinement_observable_close = .false.
        if (.not. all(ieee_is_finite([first, second, atol, rtol, &
                reference_scale]))) return
        if (atol < 0.0_dp .or. rtol < 0.0_dp .or. reference_scale <= 0.0_dp) return
        call evaluate_gc_full_fow_scaled_magnitude(first, reference_scale, &
            first_scaled, local_status, normalization_message)
        if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
        call evaluate_gc_full_fow_scaled_magnitude(second, reference_scale, &
            second_scaled, local_status, normalization_message)
        if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
        call evaluate_gc_full_fow_scaled_magnitude(first-second, &
            reference_scale, error_scaled, local_status, normalization_message)
        if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
        tolerance = atol+rtol*max(first_scaled, second_scaled)
        refinement_observable_close = error_scaled <= tolerance
    end function refinement_observable_close

    logical function refinement_state_close(first, second, factory)
        type(gc_cylindrical_state_t), intent(in) :: first, second
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory

        real(dp) :: atol, rtol

        refinement_state_close = .false.
        atol = factory%options%step_refinement_absolute_tolerance
        rtol = factory%options%step_refinement_relative_tolerance
        if (.not. refinement_observable_close(first%R, second%R, atol, rtol, &
                factory%refinement_length_scale)) return
        if (.not. refinement_observable_close(first%Z, second%Z, atol, rtol, &
                factory%refinement_length_scale)) return
        if (.not. refinement_observable_close(first%phi, second%phi, atol, rtol, &
                2.0_dp*pi)) return
        if (.not. refinement_observable_close(first%p_parallel, &
                second%p_parallel, atol, rtol, &
                factory%refinement_momentum_scale)) return
        if (.not. refinement_observable_close(first%mu, second%mu, atol, rtol, &
                factory%refinement_mu_scale)) return
        refinement_state_close = .true.
    end function refinement_state_close

    subroutine record_physical_return(factory, result)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        type(gc_cylindrical_physical_return_t), intent(in) :: result

        factory%last_return_status = result%status
        factory%last_return_event_kind = result%event_kind
        factory%last_launch_event_value = result%launch_event_value
        factory%last_return_event_value = result%return_event_value
        factory%last_invariant_scaled_drift = &
            result%maximum_invariant_scaled_drift
        factory%last_energy_error = result%energy_error
        factory%last_magnetic_moment_error = result%magnetic_moment_error
        factory%last_canonical_momentum_error = &
            result%canonical_momentum_error
        if (result%status == GC_CYL_SUCCESS .and. &
            result%event_kind == GC_CYL_PHYSICAL_EVENT_RETURN .and. &
            result%physical_return_found) then
            factory%physical_return_successes = &
                factory%physical_return_successes + 1
            factory%wall_not_hit_count = factory%wall_not_hit_count + 1
            return
        end if
        if (result%event_kind == GC_CYL_PHYSICAL_EVENT_WALL .or. &
            result%status == GC_CYL_WALL_LOSS) then
            factory%wall_return_count = factory%wall_return_count + 1
            return
        end if
        if (result%event_kind == GC_CYL_PHYSICAL_EVENT_RADIAL_DOMAIN .or. &
            result%status == GC_CYL_EQUILIBRIUM_DOMAIN) then
            factory%radial_return_count = factory%radial_return_count + 1
            factory%wall_not_hit_count = factory%wall_not_hit_count + 1
            return
        end if
        select case (result%status)
        case (GC_CYL_NO_RETURN)
            factory%no_return_count = factory%no_return_count + 1
            factory%wall_not_hit_count = factory%wall_not_hit_count + 1
        case (GC_CYL_INVARIANT_ERROR)
            factory%invariant_rejection_count = &
                factory%invariant_rejection_count + 1
        case (GC_CYL_FIELD_ERROR)
            factory%field_error_count = factory%field_error_count + 1
        case (GC_CYL_POTENTIAL_ERROR)
            factory%potential_error_count = factory%potential_error_count + 1
        case (GC_CYL_STATE_ERROR)
            factory%state_error_count = factory%state_error_count + 1
        case (GC_CYL_START_ERROR)
            factory%start_error_count = factory%start_error_count + 1
        case (GC_CYL_INTEGRATOR_ERROR)
            factory%integrator_error_count = factory%integrator_error_count + 1
        case default
            factory%other_return_error_count = &
                factory%other_return_error_count + 1
        end select
        if (result%event_kind /= GC_CYL_PHYSICAL_EVENT_WALL .and. &
                result%event_kind /= GC_CYL_PHYSICAL_EVENT_RADIAL_DOMAIN .and. &
                result%status /= GC_CYL_NO_RETURN) then
            factory%return_error_count = factory%return_error_count + 1
            factory%wall_error_count = factory%wall_error_count + 1
        else if (result%event_kind == GC_CYL_PHYSICAL_EVENT_RADIAL_DOMAIN) then
            factory%wall_not_hit_count = factory%wall_not_hit_count + 1
        end if
    end subroutine record_physical_return

    subroutine make_physical_return_options(factory, invariants, options, &
            maximum_step_override, tolerance_factor)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        type(gc_cylindrical_invariants_t), intent(in) :: invariants
        type(gc_cylindrical_physical_return_options_t), intent(out) :: options
        real(dp), intent(in) :: maximum_step_override
        real(dp), intent(in), optional :: tolerance_factor

        real(dp) :: tolerance, local_tolerance_factor

        local_tolerance_factor = 1.0_dp
        if (present(tolerance_factor)) local_tolerance_factor = tolerance_factor
        if (.not. ieee_is_finite(local_tolerance_factor) .or. &
                local_tolerance_factor <= 0.0_dp) then
            options = gc_cylindrical_physical_return_options_t()
            return
        end if

        options = gc_cylindrical_physical_return_options_t()
        options%relative_tolerance = local_tolerance_factor* &
            factory%options%orbit_options%relative_tolerance
        options%absolute_tolerance = &
            local_tolerance_factor*factory%options%orbit_options%absolute_tolerance
        options%invariant_relative_tolerance = &
            factory%options%invariant_relative_tolerance
        tolerance = local_tolerance_factor*factory%options%invariant_relative_tolerance
        options%invariant_absolute_tolerance(1) = &
            max(tiny(1.0_dp), tolerance*abs(invariants%energy))
        options%invariant_absolute_tolerance(2) = &
            max(tiny(1.0_dp), tolerance*abs(invariants%magnetic_moment))
        options%invariant_absolute_tolerance(3) = &
            max(tiny(1.0_dp), &
            tolerance*abs(invariants%canonical_toroidal_momentum))
        options%event_time_tolerance = local_tolerance_factor* &
            factory%options%orbit_options%event_time_tolerance
        options%event_value_tolerance = local_tolerance_factor* &
            factory%options%cut_root_tolerance
        options%minimum_return_time = &
            factory%options%orbit_options%minimum_event_time
        options%maximum_time = factory%options%orbit_options%maximum_time
        if (.not. ieee_is_finite(maximum_step_override) .or. &
                maximum_step_override <= 0.0_dp) then
            options%maximum_step = 0.0_dp
        else
            options%maximum_step = maximum_step_override
        end if
        options%radial_distance_scale = 1.0_dp
        options%wall_distance_scale = max(1.0_dp, factory%minor_radius_cm)
        options%return_orientation = 0
        ! Buchholz Eq. 13--17 requires the opposite-oriented and the
        ! subsequent same-oriented intersections of the complete cut atlas.
        ! A same-oriented return alone is not a certified full-bounce period.
        options%require_opposite_intersection = .true.
        options%require_transverse_intersection = .true.
        options%complete_atlas_multiplicity_certified = &
            factory%cut_atlas%two_cut_multiplicity_certified .and. &
            factory%two_cut_multiplicity_certificate_id == &
                GC_EQDSK_TWO_CUT_MULTIPLICITY_CERTIFICATE_ID
        options%cut_rate_tolerance = factory%options%cut_root_tolerance
    end subroutine make_physical_return_options

    subroutine factory_physical_return_rate(position, state, field, user_data, &
            rate, status)
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: rate
        integer, intent(out) :: status

        type(eqdsk_cut_jet_t) :: cut
        real(dp) :: derivative(5), potential, potential_gradient(3)
        integer :: cut_status, dynamics_status, potential_status

        rate = 0.0_dp
        status = GC_CYL_FIELD_ERROR
        if (.not. associated(user_data)) return
        select type (factory => user_data)
        type is (gc_eqdsk_nonlocal_factory_t)
            if (abs(factory%section_orientation) /= 1) return
            call factory%potential%evaluate(position, field, potential, &
                potential_gradient, potential_status)
            if (potential_status /= GC_CYL_SUCCESS) then
                status = GC_CYL_POTENTIAL_ERROR
                return
            end if
            call gc_cylindrical_rhs(field, potential, potential_gradient, &
                factory%species%mass_g, factory%species%charge_esu, c, state, &
                derivative, dynamics_status)
            if (dynamics_status /= GC_CYL_SUCCESS) then
                status = dynamics_status
                return
            end if
            ! coordinate_velocity follows the state order (R,Z,phi).  The
            ! generated adapter owns the conversion from dot(phi) to physical
            ! arc velocity used by grad(C).Xdot.
            call evaluate_eqdsk_cut_jet(position, factory%field%field_scale, &
                factory%section_orientation, &
                [derivative(1), derivative(2), derivative(3)], cut, cut_status)
            status = map_cut_jet_gc_status(cut_status)
            if (status /= GC_CYL_SUCCESS) return
            rate = cut%cut_rate
        class default
            status = GC_CYL_FIELD_ERROR
        end select
    end subroutine factory_physical_return_rate

    subroutine factory_physical_return_event(position, state, field, user_data, &
            value, status)
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: value
        integer, intent(out) :: status

        type(eqdsk_cut_jet_t) :: cut
        integer :: cut_status

        value = 0.0_dp
        status = GC_CYL_FIELD_ERROR
        associate (unused_state => state, unused_field => field)
        end associate
        if (.not. associated(user_data)) return
        select type (factory => user_data)
        type is (gc_eqdsk_nonlocal_factory_t)
            if (abs(factory%section_orientation) /= 1) return
            call evaluate_eqdsk_cut_jet(position, factory%field%field_scale, &
                factory%section_orientation, [0.0_dp, 0.0_dp, 0.0_dp], cut, &
                cut_status)
            status = map_cut_jet_gc_status(cut_status)
            if (status /= GC_CYL_SUCCESS) return
            value = cut%cut_value
        class default
            status = GC_CYL_FIELD_ERROR
        end select
    end subroutine factory_physical_return_event

    integer function map_cut_jet_gc_status(cut_status) result(status)
        integer, intent(in) :: cut_status

        select case (cut_status)
        case (EQDSK_CUT_JET_SUCCESS)
            status = GC_CYL_SUCCESS
        case (EQDSK_CUT_JET_OUT_OF_DOMAIN)
            status = GC_CYL_EQUILIBRIUM_DOMAIN
        case default
            status = GC_CYL_FIELD_ERROR
        end select
    end function map_cut_jet_gc_status

    subroutine factory_radial_domain_event(position, state, field, user_data, &
            margin, status)
        real(dp), intent(in) :: position(3)
        type(gc_cylindrical_state_t), intent(in) :: state
        type(gc_cylindrical_field_sample_t), intent(in) :: field
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        real(dp), intent(out) :: margin
        integer, intent(out) :: status

        real(dp) :: lower, upper, scale

        associate (unused_position => position, unused_state => state)
        end associate
        margin = -1.0_dp
        status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
            ! Profile nodes are cell-centred data support, not a physical
            ! loss surface.  The radial gate is the accepted equilibrium
            ! domain; wall_model is the actual material-loss boundary.
            if (.not. factory%field%domain_initialized) then
                status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
                return
            end if
            lower = factory%field%domain_R_min
            upper = factory%field%domain_R_max
            scale = upper-lower
            if (scale <= 0.0_dp) return
            margin = min(position(1)-lower, upper-position(1))/scale
            if (.not. ieee_is_finite(margin)) return
            status = GC_CYL_SUCCESS
        class default
            status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_radial_domain_event

    subroutine integrate_cycle_averages(factory, launch, period, harmonic_m, &
            harmonic_n, include_harmonic, include_shell, omega_b, omega_phi, &
            result, status, track_behavior, maximum_step_override, &
            tolerance_factor)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        type(gc_cylindrical_class_launch_t), intent(in) :: launch
        real(dp), intent(in) :: period, omega_b, omega_phi
        integer, intent(in) :: harmonic_m, harmonic_n
        logical, intent(in) :: include_harmonic, include_shell
        type(gc_eqdsk_orbit_result_t), intent(inout) :: result
        integer, intent(out) :: status
        logical, intent(in), optional :: track_behavior
        real(dp), intent(in) :: maximum_step_override
        real(dp), intent(in), optional :: tolerance_factor

        type(vode_state_t) :: integrator
        type(fortnum_status_t) :: integration_status
        real(dp) :: initial(AUGMENTED_STATE_SIZE), atol(AUGMENTED_STATE_SIZE)
        real(dp), allocatable :: final(:)
        integer :: callback_status
        logical :: behavior_requested
        integer :: last_parallel_sign, parallel_sign_changes
        real(dp) :: parallel_sign_tolerance, local_tolerance_factor
        real(dp) :: residence_target_s_tor, unused_psihat, unused_derivative
        integer :: flux_map_status
        real(dp) :: average_real, average_imag, residence_average, field_average

        status = GC_EQDSK_NONLOCAL_ORBIT_ERROR
        callback_status = GC_CYL_SUCCESS
        behavior_requested = present(track_behavior) .and. track_behavior
        if (.not. ieee_is_finite(period) .or. period <= 0.0_dp) return
        local_tolerance_factor = 1.0_dp
        if (present(tolerance_factor)) local_tolerance_factor = tolerance_factor
        if (.not. ieee_is_finite(local_tolerance_factor) .or. &
                local_tolerance_factor <= 0.0_dp) return
        last_parallel_sign = 0
        parallel_sign_changes = 0
        parallel_sign_tolerance = max(tiny(1.0_dp), 1.0e-12_dp* &
            factory%species%mass_g*factory%species%reference_velocity_cm_s)
        residence_target_s_tor = factory%options%residence_target_s_tor
        if (include_shell .and. residence_target_s_tor < 0.0_dp) then
            call map_eqdsk_rho_tor_to_psihat(factory%flux_map, launch%rc, &
                residence_target_s_tor, unused_psihat, unused_derivative, &
                flux_map_status)
            if (flux_map_status /= EQDSK_FLUX_MAP_SUCCESS) return
        end if
        initial = 0.0_dp
        initial(1) = launch%state%R
        initial(2) = launch%state%Z
        initial(3) = launch%state%phi
        initial(4) = launch%state%p_parallel
        initial(5) = launch%state%mu
        atol(1:5) = local_tolerance_factor* &
            factory%options%orbit_options%absolute_tolerance
        atol(6:7) = max(tiny(1.0_dp), &
            local_tolerance_factor*factory%options%orbit_options%relative_tolerance &
            *factory%species%reference_energy_erg*period)
        atol(8) = max(tiny(1.0_dp), &
            local_tolerance_factor*factory%options%orbit_options%relative_tolerance*period)
        call vode_init(integrator, AUGMENTED_STATE_SIZE, 0.0_dp, initial)
        if (.not. ieee_is_finite(maximum_step_override) .or. &
                maximum_step_override <= 0.0_dp) return
        integrator%hmax = maximum_step_override
        integrator%max_steps = 500000
        call vode_integrate_to(average_rhs, integrator, period, &
            local_tolerance_factor*factory%options%orbit_options%relative_tolerance, &
            atol, final, &
            integration_status)
        if (callback_status /= GC_CYL_SUCCESS) return
        if (integration_status%code /= FORTNUM_OK) return
        if (size(final) /= AUGMENTED_STATE_SIZE) return
        if (.not. all(ieee_is_finite(final))) return
        if (include_harmonic) then
            call evaluate_neort_full_fow_cycle_average(final(6), final(7), final(8), &
                final(9), period, average_real, average_imag, residence_average, &
                field_average)
            result%h_m = cmplx(average_real, average_imag, kind=dp)
            if (.not. all(ieee_is_finite([real(result%h_m, dp), &
                aimag(result%h_m)]))) return
            factory%harmonic_average_successes = &
                factory%harmonic_average_successes + 1
        end if
        if (include_shell) then
            if (.not. include_harmonic) then
                call evaluate_neort_full_fow_cycle_average(final(6), final(7), &
                    final(8), final(9), period, average_real, average_imag, &
                    residence_average, field_average)
            end if
            result%shell_average = residence_average
            if (.not. ieee_is_finite(result%shell_average)) return
            if (result%shell_average < -1.0e-8_dp) return
            if (result%shell_average > 1.0_dp+1.0e-8_dp) return
            result%shell_average = max(0.0_dp, min(1.0_dp, &
                result%shell_average))
            factory%residence_average_successes = &
                factory%residence_average_successes + 1
        end if
        if (behavior_requested) then
            result%parallel_sign_changes = parallel_sign_changes
            result%winding_number = nint(result%delta_phi/(2.0_dp*pi))
            if (.not. ieee_is_finite(result%delta_phi)) return
        end if
        status = GC_EQDSK_NONLOCAL_SUCCESS

    contains

        subroutine average_rhs(time, state_array, derivative, context)
            real(dp), intent(in) :: time
            real(dp), intent(in) :: state_array(:)
            real(dp), intent(out) :: derivative(:)
            class(*), intent(in), optional :: context

            type(gc_cylindrical_state_t) :: state
            type(gc_cylindrical_field_sample_t) :: field
            real(dp) :: potential, gradient(3)
            real(dp) :: orbit_surface
            real(dp) :: target_surface
            real(dp) :: harmonic_values(8)
            complex(dp) :: perturbation_amplitude
            integer :: field_status, potential_status, dynamics_status
            integer :: surface_status, perturbation_status

            associate (unused_context => context)
            end associate
            derivative = 0.0_dp
            if (size(state_array) /= AUGMENTED_STATE_SIZE) then
                callback_status = GC_CYL_INTEGRATOR_ERROR
                return
            end if
            state%R = state_array(1)
            state%Z = state_array(2)
            state%phi = state_array(3)
            state%p_parallel = state_array(4)
            state%mu = state_array(5)
            if (behavior_requested .and. &
                    abs(state%p_parallel) > parallel_sign_tolerance) then
                if (state%p_parallel > 0.0_dp) then
                    if (last_parallel_sign < 0) parallel_sign_changes = &
                        parallel_sign_changes + 1
                    last_parallel_sign = 1
                else
                    if (last_parallel_sign > 0) parallel_sign_changes = &
                        parallel_sign_changes + 1
                    last_parallel_sign = -1
                end if
            end if
            call factory%field%evaluate(state_array(1:3), field, field_status)
            if (field_status /= GC_CYL_SUCCESS) then
                callback_status = field_status
                return
            end if
            call factory%potential%evaluate(state_array(1:3), field, potential, &
                gradient, potential_status)
            if (potential_status /= GC_CYL_SUCCESS) then
                callback_status = potential_status
                return
            end if
            call gc_cylindrical_rhs(field, potential, gradient, &
                factory%species%mass_g, factory%species%charge_esu, c, state, &
                derivative(1:5), dynamics_status)
            if (dynamics_status /= GC_CYL_SUCCESS) then
                callback_status = dynamics_status
                return
            end if
            if (include_harmonic) then
                if (.not. factory%perturbation_ready) then
                    callback_status = GC_CYL_PERTURBATION_ERROR
                    return
                end if
                call do_magfie_pert_amp_cylindrical(state_array(1:3), &
                    perturbation_amplitude, perturbation_status)
                if (perturbation_status /= MAGFIE_PERT_OK) then
                    callback_status = GC_CYL_PERTURBATION_ERROR
                    return
                end if
                call evaluate_neort_full_fow_harmonic_integrand( &
                    factory%species%mass_g, factory%species%charge_esu, &
                    state%mu, field%bmod, launch%h0, potential, &
                    state%p_parallel, real(perturbation_amplitude, dp), &
                    aimag(perturbation_amplitude), real(harmonic_m, dp), &
                    real(harmonic_n, dp), state%phi, launch%state%phi, &
                    time, omega_b, omega_phi, harmonic_values(1), &
                    harmonic_values(2), harmonic_values(3), &
                    harmonic_values(4), harmonic_values(5), &
                    harmonic_values(6), harmonic_values(7), &
                    harmonic_values(8))
                if (.not. all(ieee_is_finite(harmonic_values))) then
                    callback_status = GC_CYL_PERTURBATION_ERROR
                    return
                end if
                derivative(6) = harmonic_values(7)
                derivative(7) = harmonic_values(8)
            end if
            if (include_shell) then
                call surface_from_profile_psi(factory, field%psi, orbit_surface, &
                    surface_status)
                if (surface_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
                    callback_status = GC_CYL_EQUILIBRIUM_DOMAIN
                    return
                end if
                target_surface = residence_target_s_tor
                if (orbit_surface <= target_surface) then
                    derivative(8) = 1.0_dp
                end if
            end if
        end subroutine average_rhs

    end subroutine integrate_cycle_averages

    subroutine surface_from_profile_psi(factory, psi, surface, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        real(dp), intent(in) :: psi
        real(dp), intent(out) :: surface
        integer, intent(out) :: status

        integer :: flux_map_status

        surface = 0.0_dp
        status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
        if (.not. factory%profile_ready) return
        call map_eqdsk_scaled_psi_to_s_tor(factory%flux_map, psi, surface, &
            flux_map_status)
        select case (flux_map_status)
        case (EQDSK_FLUX_MAP_SUCCESS)
            status = GC_EQDSK_NONLOCAL_SUCCESS
        case (EQDSK_FLUX_MAP_NONFINITE)
            status = GC_EQDSK_NONLOCAL_NONFINITE
        case (EQDSK_FLUX_MAP_OUT_OF_RANGE)
            status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
        case default
            status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
        end select
    end subroutine surface_from_profile_psi

    logical function launch_state_is_consistent(factory, launch) result(valid)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        type(gc_cylindrical_class_launch_t), intent(in) :: launch

        real(dp) :: scale

        valid = .false.
        if (.not. all(ieee_is_finite([launch%state%R, launch%state%Z, &
            launch%state%phi, launch%state%p_parallel, launch%state%mu, &
            launch%h0, launch%jperp, launch%p_phi]))) return
        if (launch%state%R <= 0.0_dp) return
        if (launch%state%mu < 0.0_dp) return
        scale = max(tiny(1.0_dp), abs(launch%h0))
        if (abs(launch%energy_residual) > &
            factory%options%invariant_relative_tolerance*scale) return
        scale = max(tiny(1.0_dp), abs(launch%jperp))
        if (abs(launch%jperp_residual) > &
            factory%options%invariant_relative_tolerance*scale) return
        scale = max(tiny(1.0_dp), abs(launch%p_phi))
        if (abs(launch%canonical_residual) > &
            factory%options%invariant_relative_tolerance*scale) return
        valid = .true.
    end function launch_state_is_consistent

    logical function physical_return_has_identity(factory, launch, result) &
            result(certified)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        type(gc_cylindrical_class_launch_t), intent(in) :: launch
        type(gc_eqdsk_orbit_result_t), intent(in) :: result

        real(dp) :: position_tolerance

        certified = .false.
        if (result%status /= GC_EQDSK_NONLOCAL_SUCCESS) return
        if (.not. result%physical_return_found) return
        if (result%event_kind /= GC_CYL_PHYSICAL_EVENT_RETURN) return
        if (result%intersection_count /= 2 .or. &
                .not. result%intersection_multiplicity_certified) return
        if (abs(result%intersection_orientations(1)) /= 1 .or. &
                abs(result%intersection_orientations(2)) /= 1 .or. &
                result%intersection_orientations(1) == &
                result%intersection_orientations(2) .or. &
                result%intersection_times(2) <= result%intersection_times(1)) return
        if (.not. all(ieee_is_finite(result%intersection_rates))) return
        if (any(abs(result%intersection_rates) <= tiny(1.0_dp))) return
        if (result%return_orientation == 0) return
        position_tolerance = factory%options%orbit_derivative_tolerance &
            *max(1.0_dp, factory%minor_radius_cm)
        if (abs(result%period) <= tiny(result%period)) return
        if (.not. ieee_is_finite(result%delta_phi)) return
        if (position_tolerance <= 0.0_dp) return
        if (sqrt((result%state_at_return%R-launch%state%R)**2 &
            +(result%state_at_return%Z-launch%state%Z)**2) > &
            position_tolerance) return
        if (result%state_at_return%p_parallel*launch%state%p_parallel < 0.0_dp) &
            return
        certified = .true.
    end function physical_return_has_identity

    function make_transport_reference(factory) result(reference)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        type(gc_nonlocal_transport_reference_t) :: reference

        reference = gc_nonlocal_transport_reference_t()
        reference%section_id = 1
        reference%section_coordinate = factory%section_coordinate
        reference%section_units = factory%section_units
        reference%section_position = factory%section_position
        reference%section_flux = factory%section_flux
        reference%p_zeta_orientation = factory%section_orientation
        reference%frequency_semantics = GC_NONLOCAL_COMPLETE_CYCLE_FREQUENCIES
        reference%energy_scale = factory%species%reference_energy_erg
        reference%velocity_scale = factory%species%reference_velocity_cm_s
        reference%psi_star_scale = factory%phi_eff_normalization
        reference%hamiltonian_includes_phi = factory%potential%initialized
        reference%fixed = factory%cut_ready
    end function make_transport_reference

    subroutine factory_build_quadrature(h0_order, jk_order, user_data, quadrature, &
            status)
        integer, intent(in) :: h0_order, jk_order
        class(gc_callback_context_t), pointer, intent(inout) :: user_data
        type(gc_nonlocal_transport_quadrature_t), intent(out) :: quadrature
        integer, intent(out) :: status

        quadrature = gc_nonlocal_transport_quadrature_t()
        status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
        if (.not. associated(user_data)) return
        select type (factory => user_data)
            type is (gc_eqdsk_nonlocal_factory_t)
                call build_paired_phase_space_quadrature(factory, h0_order, &
                    jk_order, quadrature, status)
            class default
                status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        end select
    end subroutine factory_build_quadrature

    subroutine build_paired_phase_space_quadrature(factory, h0_order, jk_order, &
            quadrature, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        integer, intent(in) :: h0_order, jk_order
        type(gc_nonlocal_transport_quadrature_t), intent(out) :: quadrature
        integer, intent(out) :: status

        real(dp), allocatable :: h_nodes(:), h_weights(:), j_nodes(:), j_weights(:)
        type(gc_full_fow_phase_space_bound_certificate_t) :: certificate
        type(gc_full_fow_energy_quadrature_t) :: energy_mapped
        type(gc_full_fow_jk_envelope_t) :: envelope
        type(gc_full_fow_paired_quadrature_t) :: mapped
        real(dp) :: phi_min
        integer :: i, j, nh, nj, node, local_status
        character(len=256) :: normalization_message

        quadrature = gc_nonlocal_transport_quadrature_t()
        status = GC_EQDSK_NONLOCAL_DOMAIN_ERROR
        if (.not. factory%field_ready .or. .not. factory%profile_ready) return
        call certify_phase_space_domain(factory, local_status)
        if (local_status /= GC_EQDSK_NONLOCAL_SUCCESS) then
            status = local_status
            return
        end if
        nh = h0_order
        nj = jk_order
        if (nh < 2 .or. nj < 2) return
        if (nh*nj > factory%transport_options%max_total_nodes) return
        call build_gc_unit_gauss_legendre(nh, h_nodes, h_weights, local_status)
        if (local_status /= GC_UNIT_QUADRATURE_SUCCESS) return
        call build_gc_unit_gauss_legendre(nj, j_nodes, j_weights, local_status)
        if (local_status /= GC_UNIT_QUADRATURE_SUCCESS) return
        if (.not. factory%options%phase_space_bound_certified) return
        phi_min = factory%options%qphi_min_certificate
        if (.not. ieee_is_finite(phi_min)) return
        if (.not. factory%normalization%initialized) return
        certificate%qphi_min_certified = .true.
        certificate%bmod_min_certified = .true.
        certificate%qphi_min = phi_min
        certificate%bmod_min = factory%options%bmod_min_certificate
        certificate%method = factory%options%phase_space_bound_method
        certificate%certificate_id = 'direct-eqdsk-global-phase-space-bound'
        allocate(quadrature%h0(nh*nj), quadrature%j_k(nh*nj), &
            quadrature%weight(nh*nj), quadrature%j_k_upper_bound(nh*nj))
        node = 0
        do i = 1, nh
            call map_gc_full_fow_energy_quadrature(factory%normalization, &
                phi_min, h_nodes(i), h_weights(i), energy_mapped, local_status, &
                normalization_message)
            if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
            call evaluate_gc_full_fow_jk_envelope(factory%normalization, &
                energy_mapped%h_physical, certificate, envelope, local_status, &
                normalization_message)
            if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
            do j = 1, nj
                node = node + 1
                call map_gc_full_fow_paired_quadrature(factory%normalization, &
                    envelope, h_nodes(i), h_weights(i), j_nodes(j), &
                    j_weights(j), mapped, local_status, normalization_message)
                if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
                quadrature%h0(node) = mapped%h_physical
                quadrature%j_k(node) = mapped%jk_physical
                quadrature%weight(node) = mapped%paired_weight
                quadrature%j_k_upper_bound(node) = envelope%jk_max_physical
            end do
        end do
        quadrature%h0_order = nh
        quadrature%jk_order = nj
        quadrature%n_nodes = node
        quadrature%paired_domain = .true.
        quadrature%domain_certified = .true.
        quadrature%upper_bound_is_envelope = .true.
        quadrature%domain_certificate = &
            'axisymmetric_geqdsk_fpol_fluxfunction_zero_nodes'
        !! A rule certifies only the paired physical domain.  Integrand
        !! convergence is certified by the outer N-vs-2N comparison.
        quadrature%converged = .false.
        quadrature%h0_min = phi_min
        quadrature%h0_scale = factory%species%reference_energy_erg
        if (.not. all(ieee_is_finite(quadrature%h0))) return
        if (.not. all(ieee_is_finite(quadrature%j_k))) return
        if (.not. all(ieee_is_finite(quadrature%weight))) return
        if (.not. all(ieee_is_finite(quadrature%j_k_upper_bound))) return
        if (any(quadrature%weight < 0.0_dp)) return
        if (any(quadrature%j_k_upper_bound < 0.0_dp)) return
        if (.not. any(quadrature%weight > 0.0_dp)) return
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine build_paired_phase_space_quadrature

    subroutine certify_phase_space_domain(factory, status)
        !! Certify the only envelope used by the paired physical domain.
        !!
        !! Phi is a flux function and interpolate_profile is piecewise
        !! linear in the monotone physical psi_star axis, hence the exact
        !! minimum of q*Phi on the covered interval is attained at an input
        !! node.  For the axisymmetric GEQDSK field, B >= |F(psi)|/R.  The
        !! quintic libneo fpol coefficients are enclosed on every accepted
        !! psi cell by the coefficient absolute-value bound.  A cell whose
        !! F enclosure crosses zero is rejected; no sampled lower estimate is
        !! accepted.  This certificate is intentionally not a 3-D Phi claim.
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        integer, intent(out) :: status

        real(dp) :: qphi_min, x_lo, x_hi
        real(dp) :: cell_lo, cell_hi, f_lo, f_hi
        real(dp) :: f_lower, r_upper, b_lower, field_scale_local
        type(gc_full_fow_degree5_enclosure_t) :: enclosure
        integer :: i, cell, first_cell, last_cell, local_status
        character(len=256) :: normalization_message

        status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
        factory%options%phase_space_bound_certified = .false.
        factory%refinement_mu_scale = 0.0_dp
        if (.not. factory%field_ready .or. .not. factory%profile_ready) return
        if (.not. allocated(factory%profile_psi) .or. &
                .not. allocated(factory%profile_values)) return
        if (size(factory%profile_psi) < 2) return
        if (.not. all(ieee_is_finite(factory%profile_psi)) .or. &
                .not. all(ieee_is_finite(factory%profile_values))) return
        do i = 2, size(factory%profile_psi)
            if (factory%profile_psi(i) <= factory%profile_psi(i-1)) return
        end do
        qphi_min = huge(1.0_dp)
        do i = 1, size(factory%profile_psi)
            qphi_min = min(qphi_min, factory%species%charge_esu * &
                factory%profile_values(i, PROFILE_PHI))
        end do
        if (.not. ieee_is_finite(qphi_min)) return
        if (.not. factory%flux_map%initialized .or. &
                .not. factory%flux_map%endpoints_certified) return
        r_upper = factory%field%domain_R_max
        field_scale_local = factory%field%field_scale
        if (.not. ieee_is_finite(r_upper) .or. r_upper <= 0.0_dp .or. &
                .not. ieee_is_finite(field_scale_local) .or. &
                field_scale_local <= 0.0_dp) return

        if (use_fpol) then
            if (.not. allocated(splfpol)) return
            if (size(splfpol, 1) < 6 .or. size(splfpol, 2) < 2) return
            if (.not. ieee_is_finite(hfpol) .or. hfpol <= 0.0_dp) return
            if (.not. ieee_is_finite(psi_sep) .or. psi_sep <= 0.0_dp) return
            ! The profile map has already certified its labelled s_tor=0 and
            ! s_tor=1 endpoints against the Fortsym normalization.  Use those
            ! exact domain labels for the global F enclosure; recomputing the
            ! endpoint ratios would turn harmless axis subtraction roundoff
            ! into a false out-of-domain result.
            x_lo = 0.0_dp
            x_hi = 1.0_dp
            first_cell = max(0, min(size(splfpol, 2)-2, int(x_lo/hfpol)))
            last_cell = max(0, min(size(splfpol, 2)-2, int(x_hi/hfpol)))
            f_lower = huge(1.0_dp)
            do cell = first_cell, last_cell
                cell_lo = real(cell, dp)*hfpol
                cell_hi = min(1.0_dp, cell_lo+hfpol)
                ! A full-cell enclosure is conservative for a partial
                ! accepted interval and therefore cannot under-bound |F|.
                call evaluate_gc_full_fow_degree5_cell_enclosure( &
                    splfpol(0:5, cell+1), hfpol, enclosure, local_status, &
                    normalization_message)
                if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
                f_lo = enclosure%lower_bound
                f_hi = enclosure%upper_bound
                if (.not. all(ieee_is_finite([cell_lo, cell_hi, f_lo, f_hi]))) return
                if (f_lo <= 0.0_dp .and. f_hi >= 0.0_dp) return
                f_lower = min(f_lower, min(abs(f_lo), abs(f_hi)))
            end do
        else
            if (.not. all(ieee_is_finite([btf, rtf]))) return
            f_lower = abs(btf*rtf)
            if (f_lower <= 0.0_dp) return
        end if
        b_lower = field_scale_local*f_lower/r_upper
        if (.not. ieee_is_finite(b_lower) .or. b_lower <= 0.0_dp) return
        factory%options%qphi_min_certificate = qphi_min
        factory%options%bmod_min_certificate = b_lower
        factory%refinement_mu_scale = gc_mu_phys_from_vperp( &
            factory%normalization%reference_velocity, &
            factory%normalization%mass, b_lower)
        if (.not. ieee_is_finite(factory%refinement_mu_scale) .or. &
                factory%refinement_mu_scale <= 0.0_dp) return
        factory%options%phase_space_bound_method = &
            'axisymmetric_geqdsk_fpol_fluxfunction'
        factory%options%phase_space_bound_certified = .true.
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine certify_phase_space_domain

    subroutine compute_jk_upper_bound(factory, h0, upper_bound, status)
        type(gc_eqdsk_nonlocal_factory_t), intent(inout) :: factory
        real(dp), intent(in) :: h0
        real(dp), intent(out) :: upper_bound
        integer, intent(out) :: status

        type(gc_full_fow_phase_space_bound_certificate_t) :: certificate
        type(gc_full_fow_jk_envelope_t) :: envelope
        integer :: local_status
        character(len=256) :: normalization_message

        upper_bound = 0.0_dp
        status = GC_EQDSK_NONLOCAL_CERTIFICATION_FAILED
        if (.not. factory%normalization%initialized) return
        if (.not. factory%options%phase_space_bound_certified) return
        certificate%qphi_min_certified = .true.
        certificate%bmod_min_certified = .true.
        certificate%qphi_min = factory%options%qphi_min_certificate
        certificate%bmod_min = factory%options%bmod_min_certificate
        certificate%method = factory%options%phase_space_bound_method
        certificate%certificate_id = 'direct-eqdsk-global-phase-space-bound'
        call evaluate_gc_full_fow_jk_envelope(factory%normalization, h0, &
            certificate, envelope, local_status, normalization_message)
        if (local_status /= GC_FULL_FOW_NORMALIZATION_SUCCESS) return
        upper_bound = envelope%jk_max_physical
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end subroutine compute_jk_upper_bound

    integer function validate_factory_inputs(factory, harmonic_n) result(status)
        type(gc_eqdsk_nonlocal_factory_t), intent(in) :: factory
        integer, intent(in) :: harmonic_n

        status = GC_EQDSK_NONLOCAL_INVALID_INPUT
        if (harmonic_n == 0) return
        if (.not. all(ieee_is_finite([factory%species%mass_g, &
            factory%species%charge_esu, &
            factory%species%reference_energy_erg, &
            factory%species%reference_velocity_cm_s]))) return
        if (factory%species%mass_g <= 0.0_dp) return
        if (abs(factory%species%charge_esu) <= tiny(1.0_dp)) return
        if (factory%species%reference_energy_erg <= 0.0_dp) return
        if (factory%species%reference_velocity_cm_s <= 0.0_dp) return
        if (factory%transport_options%h0_order < 2) return
        if (factory%transport_options%jk_order < 2) return
        if (factory%transport_options%max_total_nodes < &
            factory%transport_options%h0_order*factory%transport_options%jk_order) return
        if (factory%options%field_scale <= 0.0_dp) return
        ! The physical equilibrium domain is the complete normalized flux
        ! interval.  Profile support is checked separately while loading the
        ! data; a cell-centred profile edge must never become a loss surface.
        if (abs(factory%options%surface_min) > 100.0_dp*epsilon(1.0_dp)) return
        if (abs(factory%options%surface_max-1.0_dp) > &
                100.0_dp*epsilon(1.0_dp)) return
        if (factory%options%reference_surface < factory%options%surface_min .or. &
            factory%options%reference_surface > factory%options%surface_max) return
        if (factory%options%cut_theta_points < 32) return
        if (factory%options%cut_root_tolerance <= 0.0_dp) return
        if (factory%options%cut_derivative_step <= 0.0_dp) return
        if (factory%options%orbit_derivative_fraction <= 0.0_dp) return
        if (factory%options%orbit_derivative_tolerance <= 0.0_dp) return
        if (factory%options%invariant_relative_tolerance <= 0.0_dp) return
        if (factory%options%orbit_maximum_step < 0.0_dp) return
        if (factory%options%topology_probe_fraction <= 0.0_dp) return
        if (factory%options%topology_probe_fraction >= 0.5_dp) return
        if (factory%options%topology_probe_count < 3) return
        if (factory%options%endpoint_initial_grid_fraction <= 0.0_dp .or. &
                factory%options%endpoint_initial_grid_fraction >= 0.5_dp) &
            return
        if (factory%options%endpoint_max_box_expansions < 1) return
        if (factory%options%residence_target_s_tor >= 0.0_dp) then
            if (factory%options%residence_target_s_tor < &
                factory%options%surface_min) return
            if (factory%options%residence_target_s_tor > &
                factory%options%surface_max) return
        end if
        if (len_trim(factory%wall_units) == 0) return
        status = GC_EQDSK_NONLOCAL_SUCCESS
    end function validate_factory_inputs

    logical function readable_file(path)
        character(len=*), intent(in) :: path
        logical :: exists
        integer :: ios

        readable_file = .false.
        if (len_trim(path) == 0) return
        inquire(file=trim(path), exist=exists, iostat=ios)
        if (ios == 0) readable_file = exists
    end function readable_file

    integer function sign_of_real(value) result(sign_value)
        real(dp), intent(in) :: value

        sign_value = 0
        if (value > 0.0_dp) sign_value = 1
        if (value < 0.0_dp) sign_value = -1
    end function sign_of_real

end module neort_gc_eqdsk_nonlocal_transport
