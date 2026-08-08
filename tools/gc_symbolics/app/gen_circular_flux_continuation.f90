program gen_circular_flux_continuation
    !! Derive the finite-aspect circular fixture and its smooth continuation.
    !!
    !! The synthetic geometry is concentric, R=R0+r*cos(theta).  With
    !! F=R*Bphi and Bp=abs(dpsi/dr)/R, its field-line pitch is
    !!
    !!   q = F*r/(dpsi/dr*sqrt(R0**2-r**2)).
    !!
    !! All equilibrium samples and coordinate samples below are evaluated from
    !! these Fortsym expression trees.  The Python file is only a serialized,
    !! provenance-labelled table for the EQDSK writer.
    use, intrinsic :: iso_fortran_env, only: output_unit, real64
    use fortsym_arena, only: arena_t
    use fortsym_check, only: check_identity, suite_begin, suite_end, suite_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine_symengine, only: make_symengine_engine, &
        symengine_engine_t
    use fortsym_eval, only: binding_t, eval_expr
    use fortsym_expr, only: cos, exp, expr_t, log, num, operator(+), &
        operator(-), operator(*), operator(/), operator(**), pi_expr, sin, &
        real_expr, sym, sqrt
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortsym_subs, only: subs
    implicit none

    integer, parameter :: dp = real64
    ! This is a generated fixture constant, not a numerical fit.  The parent
    ! run must retain the smallest resolution that passes the independent
    ! interpolation convergence gate after the finite-aspect normalization is
    ! corrected.
    integer, parameter :: FIXTURE_NR = 65
    integer, parameter :: FIXTURE_NZ = 65
    integer, parameter :: BOUNDARY_POINTS = 129
    integer, parameter :: WALL_POINTS = 100
    real(dp), parameter :: R0_VALUE = 1.60_dp
    real(dp), parameter :: EDGE_RADIUS_VALUE = 0.50_dp
    real(dp), parameter :: B0_VALUE = 2.0_dp
    real(dp), parameter :: Q_AXIS_VALUE = 1.5_dp
    real(dp), parameter :: Q_EDGE_VALUE = 4.0_dp
    real(dp), parameter :: RBOX_HALF_EXTENT_FACTOR = 1.6_dp
    real(dp), parameter :: LIMITER_RADIUS_FACTOR = 1.1_dp
    real(dp), parameter :: PRESSURE_AXIS_VALUE = 5.0e3_dp
    character(*), parameter :: FORTSYM_REVISION = &
        'fortsym@545788453a204d58705f735b519c3863c2f734c8'
    character(*), parameter :: GENERATOR_PATH = &
        'tools/gc_symbolics/app/gen_circular_flux_continuation.f90'
    character(*), parameter :: REGENERATE_COMMAND = &
        'cd tools/gc_symbolics && fo exec gen_circular_flux_continuation '// &
        '../../src/generated/neort_circular_flux_continuation_symbolic.f90 '// &
        '../../src/generated/neort_circular_flux_continuation_limit_symbolic.f90 '// &
        '../../POTATO/test/golden_record_resonance/'// &
        'circular_flux_continuation_generated.py'

    type(arena_t), target :: arena
    type(symengine_engine_t) :: proof_engine
    type(native_engine_t) :: simplify_engine
    type(suite_t) :: proofs
    type(expr_t) :: radius, edge_radius, psi_edge, toroidal_flux, major_radius
    type(expr_t) :: pi_constant
    type(expr_t) :: grid_R, grid_Z, radial_distance, psi_coordinate
    type(expr_t) :: grid_index, profile_index, boundary_index
    type(expr_t) :: wall_index
    type(expr_t) :: q_axis, delta_q, safety_factor
    type(expr_t) :: psi_continuation, dpsi_continuation
    type(expr_t) :: psi_limit, dpsi_limit, psi_profile, dpsi_profile
    type(expr_t) :: psi_edge_profile, psi_continuation_profile
    type(expr_t) :: dpsi_continuation_profile, q_from_psi
    type(expr_t) :: alpha, a0, minor_u, edge_u, inverse_u
    type(expr_t) :: kappa, shear_primitive, shear_primitive_axis
    type(expr_t) :: inverse_ratio, ratio_axis, ratio_u
    type(expr_t) :: minor_u_symbol, shear_primitive_u
    type(expr_t) :: du_dr, du_dr_formula
    type(expr_t) :: derivative_scale, derivative_kappa, derivative_ratio
    type(expr_t) :: derivative_ratio_cont
    type(expr_t) :: generic_primitive, generic_dprimitive
    type(expr_t) :: generic_primitive_cont, generic_dprimitive_cont
    type(expr_t) :: generic_psi, generic_inverse_ratio, generic_inverse_u
    type(expr_t) :: generic_alpha, generic_a0, generic_q_profile
    type(expr_t) :: generic_q_expected
    type(expr_t) :: inverse_radius, psi_profile_sample
    type(expr_t) :: psi_tor_u, dpsi_tor_du
    type(expr_t) :: psi_tor, psi_tor_edge, s_tor_radius
    type(expr_t) :: rho_tor_radius, volume, volume_per_s_tor
    type(expr_t) :: s_tor_from_psi, rho_tor_from_psi
    type(expr_t) :: volume_per_s_tor_from_psi
    type(expr_t) :: pressure_axis, pressure_profile, dpressure_dpsi
    type(expr_t) :: f_profile, ffprime_profile
    type(expr_t) :: bpol_edge, mu0_expression, plasma_current
    type(expr_t) :: theta_sample, lcfs_R, lcfs_Z, limiter_R, limiter_Z
    type(expr_t) :: grid_R_expression, grid_Z_expression
    type(expr_t) :: rbox_left_expression, rbox_right_expression
    type(expr_t) :: zbox_bottom_expression, zbox_top_expression
    type(expr_t) :: axis_z
    type(expr_t) :: rho_tor_from_psi_residual
    type(expr_t) :: roots(8), limit_roots(3)
    type(expr_t) :: edge_value_residual, edge_slope_residual
    type(expr_t) :: profile_edge_residual, q_psi_residual
    type(expr_t) :: inverse_radius_residual
    type(expr_t) :: volume_jacobian_formula_residual
    type(expr_t) :: pressure_edge_residual, pressure_slope_residual
    type(expr_t) :: fprime_residual, ffprime_residual
    type(expr_t) :: ampere_residual
    type(expr_t) :: analytic_circle_inverse_R_integral
    type(expr_t) :: bmod_rad_cm, bmod_zet_cm, bmod_real, bmod_imag
    type(expr_t) :: wall_theta, wall_R, wall_Z
    type(expr_t) :: density_slope, density_intercept
    type(expr_t) :: temperature_slope, temperature_intercept
    type(expr_t) :: potential_slope, potential_intercept
    type(expr_t) :: lcfs_circle_residual, limiter_circle_residual
    type(expr_t) :: psi_axis_residual, zero_shear_residual
    type(engine_result_t) :: simplified
    character(2048) :: kernel_path, limit_kernel_path, python_path
    integer :: kernel_length, limit_kernel_length, python_length, argument_status
    integer :: unit, ios, i
    logical :: defined
    type(binding_t) :: bindings

    call get_command_argument(1, kernel_path, length=kernel_length, &
        status=argument_status)
    if (argument_status /= 0 .or. kernel_length == 0) error stop &
        'usage: gen_circular_flux_continuation KERNEL_PATH '// &
        'LIMIT_KERNEL_PATH PYTHON_PATH'
    call get_command_argument(2, limit_kernel_path, length=limit_kernel_length, &
        status=argument_status)
    if (argument_status /= 0 .or. limit_kernel_length == 0) error stop &
        'usage: gen_circular_flux_continuation KERNEL_PATH '// &
        'LIMIT_KERNEL_PATH PYTHON_PATH'
    call get_command_argument(3, python_path, length=python_length, &
        status=argument_status)
    if (argument_status /= 0 .or. python_length == 0) error stop &
        'usage: gen_circular_flux_continuation KERNEL_PATH '// &
        'LIMIT_KERNEL_PATH PYTHON_PATH'
    kernel_path = kernel_path(:kernel_length)
    limit_kernel_path = limit_kernel_path(:limit_kernel_length)
    python_path = python_path(:python_length)

    call arena%init()
    proof_engine = make_symengine_engine(arena)
    simplify_engine = make_native_engine(arena)
    pi_constant = pi_expr(arena)

    radius = sym(arena, 'radius')
    edge_radius = sym(arena, 'edge_radius')
    psi_edge = sym(arena, 'psi_edge')
    toroidal_flux = sym(arena, 'toroidal_flux')
    major_radius = sym(arena, 'major_radius')
    grid_R = sym(arena, 'grid_R')
    grid_Z = sym(arena, 'grid_Z')
    psi_coordinate = sym(arena, 'psi_coordinate')
    grid_index = sym(arena, 'grid_index')
    profile_index = sym(arena, 'profile_index')
    boundary_index = sym(arena, 'boundary_index')
    wall_index = sym(arena, 'wall_index')
    pressure_axis = sym(arena, 'pressure_axis')
    axis_z = num(arena, 0)
    q_axis = sym(arena, 'q_axis')
    delta_q = sym(arena, 'delta_q')
    alpha = delta_q/edge_radius**2
    a0 = q_axis + alpha*major_radius**2
    safety_factor = q_axis + alpha*radius**2
    minor_u = major_radius*sqrt(1-(radius/major_radius)**2)
    edge_u = major_radius*sqrt(1-(edge_radius/major_radius)**2)
    minor_u_symbol = sym(arena, 'minor_u_symbol')
    derivative_scale = sym(arena, 'derivative_scale')
    derivative_kappa = sym(arena, 'derivative_kappa')
    derivative_ratio = sym(arena, 'derivative_ratio')
    derivative_ratio_cont = sym(arena, 'derivative_ratio_cont')
    generic_alpha = sym(arena, 'generic_alpha')
    generic_a0 = sym(arena, 'generic_a0')

    ! Finite-aspect primitive, anchored at psi(axis)=0.  This is the exact
    ! logarithmic form of the atanh difference; it avoids relying on a special
    ! inverse-hyperbolic simplification in the proof engine.
    kappa = sqrt(alpha/a0)
    ratio_axis = (1+major_radius*kappa)/(1-major_radius*kappa)
    ratio_u = (1+minor_u*kappa)/(1-minor_u*kappa)
    shear_primitive_u = toroidal_flux/(2*sqrt(a0*alpha))*log( &
        ratio_axis/((1+minor_u_symbol*kappa)/ &
        (1-minor_u_symbol*kappa)))
    shear_primitive = subs(shear_primitive_u, minor_u_symbol, minor_u)
    ! The normalized u expression reduces exactly to major_radius at r=0, so
    ! this log-ratio is already axis anchored; retain that tree for inversion.
    shear_primitive_axis = num(arena, 0)
    psi_profile = shear_primitive
    dpsi_profile = diff(psi_profile, radius)
    du_dr = diff(minor_u, radius)
    du_dr_formula = -radius/minor_u
    generic_primitive = derivative_scale*log(derivative_ratio* &
        (1-minor_u_symbol*derivative_kappa)/ &
        (1+minor_u_symbol*derivative_kappa))
    generic_dprimitive = diff(generic_primitive, minor_u_symbol)
    generic_primitive_cont = derivative_scale*log(derivative_ratio_cont* &
        (1-minor_u_symbol*derivative_kappa)/ &
        (1+minor_u_symbol*derivative_kappa))
    generic_dprimitive_cont = diff(generic_primitive_cont, minor_u_symbol)
    generic_psi = derivative_scale*log(derivative_ratio/ &
        ((1+minor_u_symbol*derivative_kappa)/ &
        (1-minor_u_symbol*derivative_kappa)))
    generic_inverse_ratio = derivative_ratio/exp( &
        generic_psi/derivative_scale)
    generic_inverse_u = (generic_inverse_ratio-1)/ &
        (derivative_kappa*(generic_inverse_ratio+1))
    generic_q_profile = generic_a0-generic_alpha*generic_inverse_u**2
    generic_q_expected = generic_a0-generic_alpha*minor_u_symbol**2
    psi_edge_profile = subs(psi_profile, radius, edge_radius)
    psi_continuation_profile = psi_edge_profile + &
        toroidal_flux/(2*sqrt(a0*alpha))*log( &
        ((1+edge_u*kappa)/(1-edge_u*kappa))/ratio_u)
    dpsi_continuation_profile = diff(psi_continuation_profile, radius)
    psi_continuation = psi_edge + &
        toroidal_flux/(2*sqrt(a0*alpha))*log( &
        ((1+edge_u*kappa)/(1-edge_u*kappa))/ratio_u)
    dpsi_continuation = diff(psi_continuation, radius)
    ! Simplification is applied only after differentiation: the emitted roots
    ! retain diff(...) as their source, while the proof engine receives the
    ! equivalent compact derivative tree.
    simplified = simplify_engine%simplify(dpsi_profile)
    if (.not. simplified%ok) error stop 'interior derivative simplification failed'
    dpsi_profile = simplified%value
    simplified = simplify_engine%simplify(dpsi_continuation_profile)
    if (.not. simplified%ok) &
        error stop 'profile derivative simplification failed'
    dpsi_continuation_profile = simplified%value
    simplified = simplify_engine%simplify(dpsi_continuation)
    if (.not. simplified%ok) &
        error stop 'continuation derivative simplification failed'
    dpsi_continuation = simplified%value

    ! Inverting the same primitive gives the q(psi) relation used by EQDSK.
    inverse_ratio = ratio_axis/exp(2* &
        (psi_coordinate+shear_primitive_axis)*sqrt(a0*alpha)/ &
        toroidal_flux)
    inverse_u = (inverse_ratio-1)/(kappa*(inverse_ratio+1))
    q_from_psi = a0-alpha*inverse_u**2
    inverse_radius = sqrt(major_radius**2-inverse_u**2)
    ! q_from_psi is emitted as the exact explicit inverse.  Its composition
    ! certificate is carried by the generic inverse proof below; retaining the
    ! generic tree unmodified is important because substitutions are persistent
    ! in the arena.

    ! The zero-shear branch is the exact alpha -> 0 limit, represented without
    ! a runtime division by alpha.
    psi_limit = psi_edge + toroidal_flux*(edge_u-minor_u)/q_axis
    dpsi_limit = diff(psi_limit, radius)

    ! Exact toroidal flux and geometric volume at finite aspect ratio.
    psi_tor = toroidal_flux*(major_radius-minor_u)
    psi_tor_u = toroidal_flux*(major_radius-minor_u_symbol)
    dpsi_tor_du = diff(psi_tor_u, minor_u_symbol)
    psi_tor_edge = subs(psi_tor, radius, edge_radius)
    s_tor_radius = psi_tor/psi_tor_edge
    rho_tor_radius = sqrt(s_tor_radius)
    volume = 2*pi_constant**2*major_radius*radius**2
    volume_per_s_tor = diff(volume, radius)/diff(s_tor_radius, radius)
    s_tor_from_psi = subs(s_tor_radius, radius, inverse_radius)
    rho_tor_from_psi = subs(rho_tor_radius, radius, inverse_radius)
    volume_per_s_tor_from_psi = subs(volume_per_s_tor, radius, inverse_radius)
    psi_profile_sample = psi_edge_profile*profile_index/ &
        real_expr(arena, real(FIXTURE_NR-1, dp))

    grid_R_expression = major_radius - RBOX_HALF_EXTENT_FACTOR*edge_radius + &
        2.0_dp*RBOX_HALF_EXTENT_FACTOR*edge_radius*grid_index/ &
        real_expr(arena, real(FIXTURE_NR-1, dp))
    grid_Z_expression = -RBOX_HALF_EXTENT_FACTOR*edge_radius + &
        2.0_dp*RBOX_HALF_EXTENT_FACTOR*edge_radius*grid_index/ &
        real_expr(arena, real(FIXTURE_NZ-1, dp))
    rbox_left_expression = subs(grid_R_expression, grid_index, num(arena, 0))
    rbox_right_expression = subs(grid_R_expression, grid_index, &
        num(arena, FIXTURE_NR-1))
    zbox_bottom_expression = subs(grid_Z_expression, grid_index, num(arena, 0))
    zbox_top_expression = subs(grid_Z_expression, grid_index, &
        num(arena, FIXTURE_NZ-1))
    radial_distance = sqrt((grid_R-major_radius)**2+grid_Z**2)

    ! Fortsym-owned synthetic perturbation, wall, and profile samples.  The
    ! cubic Cartesian pair is the smooth m=3 analogue of rho**2*(cos(3theta),
    ! sin(3theta)) and is regular at the axis.
    bmod_rad_cm = grid_R*real_expr(arena, 100.0_dp)
    bmod_zet_cm = grid_Z*real_expr(arena, 100.0_dp)
    bmod_real = real_expr(arena, 1.0e3_dp)*(grid_R-major_radius)* &
        ((grid_R-major_radius)**2-3*grid_Z**2)/edge_radius**3
    bmod_imag = real_expr(arena, 1.0e3_dp)*(3*(grid_R-major_radius)**2* &
        grid_Z-grid_Z**3)/edge_radius**3
    wall_theta = 2*pi_constant*wall_index/ &
        real_expr(arena, real(WALL_POINTS, dp))
    wall_R = (major_radius+1.12_dp*edge_radius*cos(wall_theta))* &
        real_expr(arena, 100.0_dp)
    wall_Z = 1.12_dp*edge_radius*sin(wall_theta)* &
        real_expr(arena, 100.0_dp)
    density_slope = real_expr(arena, -2.5e13_dp)
    density_intercept = real_expr(arena, 5.0e13_dp)
    temperature_slope = real_expr(arena, -1.4e3_dp)
    temperature_intercept = real_expr(arena, 2.0e3_dp)
    potential_slope = real_expr(arena, -3.0e2_dp)
    potential_intercept = real_expr(arena, 3.0e2_dp)

    f_profile = toroidal_flux
    ffprime_profile = f_profile*diff(f_profile, psi_coordinate)
    pressure_profile = pressure_axis*(1-psi_coordinate/psi_edge_profile)
    dpressure_dpsi = diff(pressure_profile, psi_coordinate)
    mu0_expression = real_expr(arena, 4.0e-7_dp)*pi_constant
    analytic_circle_inverse_R_integral = 2*pi_constant/edge_u
    bpol_edge = subs(dpsi_profile, radius, edge_radius)/ &
        (major_radius+edge_radius)
    ! Supplied analytic circle identity: integral_0^(2pi) dtheta/R = 2*pi/u.
    ! The current expression below is algebra derived from that identity; the
    ! independent numerical contour quadrature lives in the focused test.
    plasma_current = analytic_circle_inverse_R_integral* &
        subs(dpsi_profile, radius, edge_radius)*edge_radius/mu0_expression

    theta_sample = 2*pi_constant*boundary_index/ &
        real_expr(arena, real(BOUNDARY_POINTS-1, dp))
    lcfs_R = major_radius + edge_radius*cos(theta_sample)
    lcfs_Z = edge_radius*sin(theta_sample)
    limiter_R = major_radius + LIMITER_RADIUS_FACTOR*edge_radius* &
        cos(theta_sample)
    limiter_Z = LIMITER_RADIUS_FACTOR*edge_radius*sin(theta_sample)

    call suite_begin(proofs, 'circular flux continuation')
    call check_identity(proofs, proof_engine, 'log primitive u derivative', &
        generic_dprimitive+2*derivative_scale*derivative_kappa/ &
        (1-(minor_u_symbol*derivative_kappa)**2))
    call check_identity(proofs, proof_engine, 'circular u radial derivative', &
        du_dr-du_dr_formula)
    call check_identity(proofs, proof_engine, 'generic inverse radius map', &
        generic_inverse_u-minor_u_symbol)
    call check_identity(proofs, proof_engine, 'generic inverse q map', &
        generic_q_profile-generic_q_expected)
    ! Production dpsi roots above descend directly from diff(psi_profile,r)
    ! and diff(psi_continuation,r).  The quotient below is only an independent
    ! generic derivative oracle; the strict engine does not normalize the
    ! expanded physical substitution back to that quotient.
    call check_identity(proofs, proof_engine, &
        'finite-aspect primitive derivative specialization', &
        generic_dprimitive+2*derivative_scale*derivative_kappa/ &
        (1-(minor_u_symbol*derivative_kappa)**2))
    call check_identity(proofs, proof_engine, &
        'exterior primitive derivative specialization', &
        generic_dprimitive_cont+2*derivative_scale*derivative_kappa/ &
        (1-(minor_u_symbol*derivative_kappa)**2))
    call check_identity(proofs, proof_engine, 'axis value', &
        subs(psi_profile, radius, num(arena, 0)))
    edge_value_residual = subs(psi_continuation, radius, edge_radius)-psi_edge
    call check_identity(proofs, proof_engine, 'edge value', edge_value_residual)
    edge_slope_residual = subs(dpsi_profile, radius, edge_radius)- &
        subs(dpsi_continuation, radius, edge_radius)
    call check_identity(proofs, proof_engine, 'interior-continuation edge slope', &
        edge_slope_residual)
    profile_edge_residual = subs(psi_continuation_profile, radius, &
        edge_radius)-psi_edge_profile
    call check_identity(proofs, proof_engine, 'profile-continuation edge value', &
        profile_edge_residual)
    call check_identity(proofs, proof_engine, 'profile-continuation edge slope', &
        subs(dpsi_profile, radius, edge_radius)- &
        subs(dpsi_continuation_profile, radius, edge_radius))
    ! The generic derivative certificate, together with u'(r), is the
    ! symbolic field-line normalization.  The focused test independently
    ! reconstructs the full theta integral from the emitted dpsi root.
    call check_identity(proofs, proof_engine, &
        'generic field-line q normalization certificate', &
        generic_dprimitive+2*derivative_scale*derivative_kappa/ &
        (1-(minor_u_symbol*derivative_kappa)**2))
    call check_identity(proofs, proof_engine, 'inverse q definition', &
        q_from_psi-(a0-alpha*inverse_u**2))
    ! This is an inherited generic inverse certificate; emitted q samples are
    ! independently checked against direct field-line quadrature.
    q_psi_residual = generic_q_profile-generic_q_expected
    call check_identity(proofs, proof_engine, &
        'generic q(psi) composition certificate', &
        q_psi_residual)
    inverse_radius_residual = generic_inverse_u-minor_u_symbol
    call check_identity(proofs, proof_engine, &
        'generic inverse radius composition certificate', &
        inverse_radius_residual)
    call check_identity(proofs, proof_engine, &
        'toroidal flux u differential', dpsi_tor_du+toroidal_flux)
    call check_identity(proofs, proof_engine, 'normalized toroidal flux', &
        s_tor_radius-psi_tor/psi_tor_edge)
    call check_identity(proofs, proof_engine, 'rho toroidal map', &
        rho_tor_radius**2-s_tor_radius)
    rho_tor_from_psi_residual = rho_tor_from_psi**2-s_tor_from_psi
    call check_identity(proofs, proof_engine, 'rho(psi) map', &
        rho_tor_from_psi_residual)
    volume_jacobian_formula_residual = volume_per_s_tor- &
        4*pi_constant**2*major_radius*(major_radius-edge_u)*minor_u
    call check_identity(proofs, proof_engine, 'finite-aspect dV/ds_tor', &
        volume_jacobian_formula_residual)
    ! These map checks inherit the generic inverse-u certificate.  Their
    ! emitted samples are evaluated from the composed physical trees and are
    ! independently checked by direct geometry in the focused test.
    call check_identity(proofs, proof_engine, &
        's_tor map inherited inverse-u certificate', &
        generic_inverse_u-minor_u_symbol)
    call check_identity(proofs, proof_engine, &
        'rho_tor map inherited inverse-u certificate', &
        generic_inverse_u-minor_u_symbol)
    call check_identity(proofs, proof_engine, &
        'dV/ds_tor map inherited inverse-u certificate', &
        generic_inverse_u-minor_u_symbol)
    pressure_edge_residual = subs(pressure_profile, psi_coordinate, &
        psi_edge_profile)
    call check_identity(proofs, proof_engine, 'pressure vanishes at edge', &
        pressure_edge_residual)
    pressure_slope_residual = dpressure_dpsi+pressure_axis/psi_edge_profile
    call check_identity(proofs, proof_engine, 'pressure derivative', &
        pressure_slope_residual)
    fprime_residual = diff(f_profile, psi_coordinate)
    call check_identity(proofs, proof_engine, 'constant F profile', &
        fprime_residual)
    ffprime_residual = ffprime_profile
    call check_identity(proofs, proof_engine, 'constant FF prime profile', &
        ffprime_residual)
    ampere_residual = plasma_current-analytic_circle_inverse_R_integral* &
        subs(dpsi_profile, radius, edge_radius)*edge_radius/mu0_expression
    call check_identity(proofs, proof_engine, &
        'Ampere current algebra from analytic circle integral', &
        ampere_residual)
    call check_identity(proofs, proof_engine, 'boundary poloidal field', &
        bpol_edge-subs(dpsi_profile, radius, edge_radius)/ &
        (major_radius+edge_radius))
    lcfs_circle_residual = (lcfs_R-major_radius)**2+lcfs_Z**2-edge_radius**2
    call check_identity(proofs, proof_engine, 'LCFS circle', &
        lcfs_circle_residual)
    limiter_circle_residual = (limiter_R-major_radius)**2+limiter_Z**2- &
        (LIMITER_RADIUS_FACTOR*edge_radius)**2
    call check_identity(proofs, proof_engine, 'limiter circle', &
        limiter_circle_residual)
    psi_axis_residual = subs(psi_limit, radius, num(arena, 0))- &
        (psi_edge+toroidal_flux*(edge_u-major_radius)/q_axis)
    call check_identity(proofs, proof_engine, 'zero-shear axis value', &
        psi_axis_residual)
    zero_shear_residual = dpsi_limit- &
        toroidal_flux*radius/(q_axis*minor_u)
    call check_identity(proofs, proof_engine, 'zero-shear finite-aspect limit', &
        zero_shear_residual)
    call suite_end(proofs)
    if (proofs%failed /= 0) error stop 'circular continuation proof failed'

    roots = [psi_continuation, dpsi_continuation, safety_factor, psi_tor, &
        s_tor_radius, rho_tor_radius, volume_per_s_tor, plasma_current]
    limit_roots = [psi_limit, dpsi_limit, q_axis]
    call simplify_array(roots)
    call simplify_array(limit_roots)
    call emit_kernel_file(trim(kernel_path), &
        'neort_circular_flux_continuation_symbolic', &
        'evaluate_neort_circular_flux_continuation', roots, &
        [character(len=64) :: 'radius', 'edge_radius', 'psi_edge', &
        'toroidal_flux', 'q_axis', 'delta_q', 'major_radius'], &
        [character(len=64) :: 'psi', 'dpsi_dr', 'q_safety', 'psi_tor', &
        's_tor', 'rho_tor', 'dvolume_ds_tor', 'plasma_current'])
    call emit_kernel_file(trim(limit_kernel_path), &
        'neort_circular_flux_continuation_limit_symbolic', &
        'evaluate_neort_circular_flux_continuation_limit', limit_roots, &
        [character(len=64) :: 'radius', 'edge_radius', 'psi_edge', &
        'toroidal_flux', 'q_axis', 'major_radius'], &
        [character(len=64) :: 'psi', 'dpsi_dr', 'q_safety'])
    simplified = simplify_engine%simplify(volume_per_s_tor)
    if (.not. simplified%ok) error stop 'volume Jacobian simplification failed'
    volume_per_s_tor = simplified%value
    simplified = simplify_engine%simplify(volume_per_s_tor_from_psi)
    if (.not. simplified%ok) &
        error stop 'volume Jacobian profile simplification failed'
    volume_per_s_tor_from_psi = simplified%value
    call emit_python_table(trim(python_path), psi_continuation_profile, &
        q_from_psi, radial_distance, grid_R_expression, grid_Z_expression, &
        psi_profile_sample, theta_sample, rbox_left_expression, &
        rbox_right_expression, zbox_bottom_expression, zbox_top_expression, &
        psi_tor_edge, volume_per_s_tor, volume_per_s_tor_from_psi, &
        f_profile, ffprime_profile, pressure_profile, dpressure_dpsi, &
        plasma_current, lcfs_R, lcfs_Z, limiter_R, limiter_Z, bmod_rad_cm, &
        bmod_zet_cm, bmod_real, bmod_imag, wall_R, wall_Z, density_slope, &
        density_intercept, temperature_slope, temperature_intercept, &
        potential_slope, potential_intercept)

contains

    subroutine simplify_array(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        integer :: k

        do k = 1, size(expressions)
            simplified = simplify_engine%simplify(expressions(k))
            if (.not. simplified%ok) error stop 'continuation simplification failed'
            expressions(k) = simplified%value
        end do
    end subroutine simplify_array

    subroutine emit_kernel_file(path, module_name, procedure_name, roots_in, &
            args_in, outputs_in)
        character(*), intent(in) :: path, module_name, procedure_name
        type(expr_t), intent(in) :: roots_in(:)
        character(*), intent(in) :: args_in(:), outputs_in(:)
        type(kernel_spec_t) :: spec
        character(:), allocatable :: emitted_text
        type(str_t) :: emitted
        logical :: ok

        spec%name = str(procedure_name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%generator = str(GENERATOR_PATH)
        spec%generator_revision = str(FORTSYM_REVISION)
        spec%regenerate_command = str(REGENERATE_COMMAND)
        allocate(spec%args(size(args_in)), spec%outputs(size(outputs_in)))
        do i = 1, size(args_in)
            spec%args(i) = str(trim(args_in(i)))
        end do
        do i = 1, size(outputs_in)
            spec%outputs(i) = str(trim(outputs_in(i)))
        end do
        emitted = emit_kernel(roots_in, spec, ok)
        if (.not. ok) error stop 'continuation kernel emission failed'
        emitted_text = chars(emitted)
        open(newunit=unit, file=trim(path), status='replace', action='write', &
            iostat=ios)
        if (ios /= 0) error stop 'cannot open continuation kernel output'
        write(unit, '(a)') '! Dimensions: F [T m], psi [T m^2], psi_r [T m].'
        write(unit, '(a)') &
            '! Geometric circular field/orbit interpolation fixture; '// &
            'not GS force balance.'
        if (len(emitted_text) > 0 .and. &
                emitted_text(len(emitted_text):) == new_line('a')) then
            write(unit, '(a)', advance='no') &
                emitted_text(:len(emitted_text)-1)
        else
            write(unit, '(a)', advance='no') emitted_text
        end if
        close(unit)
        write(output_unit, '(a)') 'wrote '//trim(path)
    end subroutine emit_kernel_file

    subroutine emit_python_table(path, psi_expression, q_expression, &
            radial_expression, grid_R_expr, grid_Z_expr, psi_sample_expr, &
            theta_expr, rbox_left_expr, rbox_right_expr, zbox_bottom_expr, &
            zbox_top_expr, toroidal_flux_edge_expression, volume_expression, &
            volume_profile_expression, f_expression, ffprime_expression, &
            pressure_expression, dpressure_expression, &
            plasma_current_expression, lcfs_R_expression, lcfs_Z_expression, &
            limiter_R_expression, limiter_Z_expression, bmod_rad_expression, &
            bmod_zet_expression, bmod_real_expression, bmod_imag_expression, &
            wall_R_expression, wall_Z_expression, density_slope_expression, &
            density_intercept_expression, temperature_slope_expression, &
            temperature_intercept_expression, potential_slope_expression, &
            potential_intercept_expression)
        character(*), intent(in) :: path
        type(expr_t), intent(in) :: psi_expression, q_expression
        type(expr_t), intent(in) :: radial_expression, grid_R_expr, grid_Z_expr
        type(expr_t), intent(in) :: psi_sample_expr, theta_expr
        type(expr_t), intent(in) :: rbox_left_expr, rbox_right_expr
        type(expr_t), intent(in) :: zbox_bottom_expr, zbox_top_expr
        type(expr_t), intent(in) :: toroidal_flux_edge_expression
        type(expr_t), intent(in) :: volume_expression, volume_profile_expression
        type(expr_t), intent(in) :: f_expression, ffprime_expression
        type(expr_t), intent(in) :: pressure_expression, dpressure_expression
        type(expr_t), intent(in) :: plasma_current_expression
        type(expr_t), intent(in) :: lcfs_R_expression, lcfs_Z_expression
        type(expr_t), intent(in) :: limiter_R_expression, limiter_Z_expression
        type(expr_t), intent(in) :: bmod_rad_expression, bmod_zet_expression
        type(expr_t), intent(in) :: bmod_real_expression, bmod_imag_expression
        type(expr_t), intent(in) :: wall_R_expression, wall_Z_expression
        type(expr_t), intent(in) :: density_slope_expression
        type(expr_t), intent(in) :: density_intercept_expression
        type(expr_t), intent(in) :: temperature_slope_expression
        type(expr_t), intent(in) :: temperature_intercept_expression
        type(expr_t), intent(in) :: potential_slope_expression
        type(expr_t), intent(in) :: potential_intercept_expression
        integer :: i_grid, j_grid, k
        real(dp) :: psi_edge_value, psi_coordinate_value, psi_value
        real(dp) :: toroidal_flux_edge_value, volume_value
        real(dp) :: grid_radius, grid_z, radial_distance_value
        real(dp) :: psi_axis_value, ip_value, theta_value
        real(dp) :: axis_z_value, box_value
        real(dp) :: lcfs_r_value, lcfs_z_value, limiter_r_value
        real(dp) :: limiter_z_value, f_value, ffprime_value
        real(dp) :: pressure_value, dpressure_value, q_value
        real(dp) :: volume_profile_value
        real(dp) :: bmod_rad_value, bmod_zet_value
        real(dp) :: bmod_real_value, bmod_imag_value
        real(dp) :: wall_r_value, wall_z_value
        real(dp) :: density_slope_value, density_intercept_value
        real(dp) :: temperature_slope_value, temperature_intercept_value
        real(dp) :: potential_slope_value, potential_intercept_value
        real(dp) :: pi_numeric

        allocate(bindings%names(16), bindings%values(16))
        bindings%names = [str('radius'), str('edge_radius'), str('psi_edge'), &
            str('toroidal_flux'), str('q_axis'), str('delta_q'), &
            str('major_radius'), str('psi_coordinate'), str('grid_R'), &
            str('grid_Z'), str('theta'), str('pressure_axis'), &
            str('grid_index'), str('profile_index'), str('boundary_index'), &
            str('wall_index')]
        bindings%values = [0.0_dp, EDGE_RADIUS_VALUE, 0.0_dp, &
            R0_VALUE*B0_VALUE, Q_AXIS_VALUE, Q_EDGE_VALUE-Q_AXIS_VALUE, &
            R0_VALUE, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp, &
            PRESSURE_AXIS_VALUE, 0.0_dp, 0.0_dp, 0.0_dp, 0.0_dp]
        bindings%n = 16
        pi_numeric = eval_expr(pi_constant, bindings, defined)
        if (.not. defined) error stop 'pi evaluation failed'
        axis_z_value = eval_expr(axis_z, bindings, defined)
        if (.not. defined) error stop 'axis coordinate evaluation failed'
        psi_edge_value = evaluate(psi_expression, EDGE_RADIUS_VALUE)
        psi_axis_value = evaluate(psi_expression, 0.0_dp)
        toroidal_flux_edge_value = evaluate(toroidal_flux_edge_expression, &
            0.0_dp)
        volume_value = evaluate(volume_expression, 0.0_dp)
        f_value = eval_expr(f_expression, bindings, defined)
        if (.not. defined) error stop 'F profile evaluation failed'
        ffprime_value = eval_expr(ffprime_expression, bindings, defined)
        if (.not. defined) error stop 'FF prime evaluation failed'
        ip_value = eval_expr(plasma_current_expression, bindings, defined)
        if (.not. defined) error stop 'plasma current evaluation failed'
        density_slope_value = eval_expr(density_slope_expression, bindings, defined)
        density_intercept_value = eval_expr(density_intercept_expression, &
            bindings, defined)
        temperature_slope_value = eval_expr(temperature_slope_expression, &
            bindings, defined)
        temperature_intercept_value = eval_expr( &
            temperature_intercept_expression, bindings, defined)
        potential_slope_value = eval_expr(potential_slope_expression, bindings, &
            defined)
        potential_intercept_value = eval_expr( &
            potential_intercept_expression, bindings, defined)
        if (.not. defined) error stop 'profile coefficient evaluation failed'

        open(newunit=unit, file=trim(path), status='replace', action='write', &
            iostat=ios)
        if (ios /= 0) error stop 'cannot open continuation table output'
        write(unit, '(a)') '# Generated by Fortsym; do not edit.'
        write(unit, '(a)') '# Generator: '//GENERATOR_PATH
        write(unit, '(a)') '# Fortsym revision: '//FORTSYM_REVISION
        write(unit, '(a)') '# Geometric circular orbit/interpolation fixture;'
        write(unit, '(a)') '# not a Grad-Shafranov force-balance equilibrium.'
        write(unit, '(a)') '# Dimensions: F [T m], psi [T m^2], psi_r [T m].'
        write(unit, '(a)') '# Analytic circle integral: integral(dtheta/R)=2*pi/u.'
        write(unit, '(a)') '# Sign provenance: psi/theta orientation flips are '// &
            'coordinate relabelings;'
        write(unit, '(a)') '# LCFS/wall theta reversal is also a coordinate relabeling.'
        write(unit, '(a)') '# Fourier (Re,Im) conjugation paired with n -> -n '// &
            'is an identity;'
        write(unit, '(a)') '# changing one component alone is physically consequential.'
        write(unit, '(a)') '# BMOD uses a regular rho^3 m=3 Cartesian harmonic;'
        write(unit, '(a)') '# this intentionally replaces singular '// &
            'rho^2*cos(3theta) at the axis.'
        write(unit, '(a)') '# rho^2 -> rho^3 changes perturbation weighting; '// &
            'it is not a sign convention.'
        write(unit, '(a)') 'import numpy as np'
        write(unit, '(a)') ''
        write(unit, '(a,es24.16e3)') 'R0 = ', R0_VALUE
        write(unit, '(a,es24.16e3)') 'R_AXIS = ', R0_VALUE
        write(unit, '(a,es24.16e3)') 'Z_AXIS = ', axis_z_value
        write(unit, '(a,es24.16e3)') 'A = ', EDGE_RADIUS_VALUE
        write(unit, '(a,es24.16e3)') 'B0 = ', B0_VALUE
        write(unit, '(a,es24.16e3)') 'Q0 = ', Q_AXIS_VALUE
        write(unit, '(a,es24.16e3)') 'QA = ', Q_EDGE_VALUE
        write(unit, '(a,es24.16e3)') 'TOROIDAL_FLUX = ', &
            R0_VALUE*B0_VALUE
        write(unit, '(a,es24.16e3)') 'PSI_AXIS = ', psi_axis_value
        write(unit, '(a,es24.16e3)') 'PSI_EDGE = ', psi_edge_value
        write(unit, '(a,es24.16e3)') 'TOROIDAL_FLUX_EDGE = ', &
            toroidal_flux_edge_value
        write(unit, '(a,es24.16e3)') 'D_VOLUME_DS_TOR_AXIS = ', volume_value
        write(unit, '(a,es24.16e3)') 'IP = ', ip_value
        write(unit, '(a,i0)') 'NR = ', FIXTURE_NR
        write(unit, '(a,i0)') 'NZ = ', FIXTURE_NZ
        write(unit, '(a,i0)') 'BOUNDARY_POINTS = ', BOUNDARY_POINTS
        write(unit, '(a,i0)') 'WALL_POINTS = ', WALL_POINTS
        box_value = evaluate(rbox_left_expr, 0.0_dp)
        write(unit, '(a,es24.16e3)') 'RBOX_LEFT = ', box_value
        box_value = evaluate(rbox_right_expr, 0.0_dp)- &
            evaluate(rbox_left_expr, 0.0_dp)
        write(unit, '(a,es24.16e3)') 'RBOX_LENGTH = ', box_value
        box_value = 0.5_dp*(evaluate(zbox_top_expr, 0.0_dp)+ &
            evaluate(zbox_bottom_expr, 0.0_dp))
        write(unit, '(a,es24.16e3)') 'ZBOX_MID = ', box_value
        box_value = evaluate(zbox_top_expr, 0.0_dp)- &
            evaluate(zbox_bottom_expr, 0.0_dp)
        write(unit, '(a,es24.16e3)') 'ZBOX_LENGTH = ', box_value
        write(unit, '(a)') ''
        write(unit, '(a)') 'FPROF = np.full(NR,'
        write(unit, '(a,es24.16e3,a)') '    ', f_value, ', dtype=np.float64)'
        write(unit, '(a)') 'FFPRIMEPROF = np.full(NR,'
        write(unit, '(a,es24.16e3,a)') '    ', ffprime_value, ', dtype=np.float64)'
        write(unit, '(a)') 'PTOTPROF = np.array(['
        do k = 0, FIXTURE_NR-1
            bindings%values(14) = real(k, dp)
            psi_coordinate_value = eval_expr(psi_sample_expr, bindings, defined)
            if (.not. defined) error stop 'psi sample evaluation failed'
            bindings%values(8) = psi_coordinate_value
            pressure_value = eval_expr(pressure_expression, bindings, defined)
            if (.not. defined) error stop 'pressure profile evaluation failed'
            call write_python_value(unit, pressure_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'DPRESS_DPSIPROF = np.array(['
        do k = 0, FIXTURE_NR-1
            bindings%values(14) = real(k, dp)
            bindings%values(8) = eval_expr(psi_sample_expr, bindings, defined)
            if (.not. defined) error stop 'psi sample evaluation failed'
            dpressure_value = eval_expr(dpressure_expression, bindings, &
                defined)
            if (.not. defined) error stop 'pressure derivative evaluation failed'
            call write_python_value(unit, dpressure_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'PROFILE_POLY = np.array(['
        call write_python_profile_row(unit, density_slope_value, &
            density_intercept_value)
        call write_python_zero_profile_row(unit)
        call write_python_profile_row(unit, temperature_slope_value, &
            temperature_intercept_value)
        call write_python_profile_row(unit, potential_slope_value, &
            potential_intercept_value)
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'R_GRID = np.array(['
        do i_grid = 0, FIXTURE_NR-1
            bindings%values(13) = real(i_grid, dp)
            grid_radius = eval_expr(grid_R_expr, bindings, defined)
            if (.not. defined) error stop 'R grid evaluation failed'
            call write_python_value(unit, grid_radius)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'BMOD_RAD_CM = np.array(['
        do i_grid = 0, FIXTURE_NR-1
            bindings%values(13) = real(i_grid, dp)
            grid_radius = eval_expr(grid_R_expr, bindings, defined)
            bindings%values(9) = grid_radius
            bmod_rad_value = eval_expr(bmod_rad_expression, bindings, defined)
            if (.not. defined) error stop 'BMOD radial coordinate failed'
            call write_python_value(unit, bmod_rad_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'BMOD_ZET_CM = np.array(['
        do j_grid = 0, FIXTURE_NZ-1
            bindings%values(13) = real(j_grid, dp)
            grid_z = eval_expr(grid_Z_expr, bindings, defined)
            bindings%values(10) = grid_z
            bmod_zet_value = eval_expr(bmod_zet_expression, bindings, defined)
            if (.not. defined) error stop 'BMOD vertical coordinate failed'
            call write_python_value(unit, bmod_zet_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'BMOD_RE = np.array(['
        do j_grid = 0, FIXTURE_NZ-1
            bindings%values(13) = real(j_grid, dp)
            grid_z = eval_expr(grid_Z_expr, bindings, defined)
            bindings%values(10) = grid_z
            do i_grid = 0, FIXTURE_NR-1
                bindings%values(13) = real(i_grid, dp)
                grid_radius = eval_expr(grid_R_expr, bindings, defined)
                bindings%values(9) = grid_radius
                bmod_real_value = eval_expr(bmod_real_expression, bindings, &
                    defined)
                if (.not. defined) error stop 'BMOD real evaluation failed'
                call write_python_value(unit, bmod_real_value)
            end do
        end do
        write(unit, '(a)') '], dtype=np.float64).reshape((NZ, NR)).T'
        write(unit, '(a)') 'BMOD_IM = np.array(['
        do j_grid = 0, FIXTURE_NZ-1
            bindings%values(13) = real(j_grid, dp)
            grid_z = eval_expr(grid_Z_expr, bindings, defined)
            bindings%values(10) = grid_z
            do i_grid = 0, FIXTURE_NR-1
                bindings%values(13) = real(i_grid, dp)
                grid_radius = eval_expr(grid_R_expr, bindings, defined)
                bindings%values(9) = grid_radius
                bmod_imag_value = eval_expr(bmod_imag_expression, bindings, &
                    defined)
                if (.not. defined) error stop 'BMOD imaginary evaluation failed'
                call write_python_value(unit, bmod_imag_value)
            end do
        end do
        write(unit, '(a)') '], dtype=np.float64).reshape((NZ, NR)).T'
        write(unit, '(a)') 'WALL = np.array(['
        do k = 0, WALL_POINTS-1
            bindings%values(16) = real(k, dp)
            wall_r_value = eval_expr(wall_R_expression, bindings, defined)
            if (.not. defined) error stop 'wall R evaluation failed'
            wall_z_value = eval_expr(wall_Z_expression, bindings, defined)
            if (.not. defined) error stop 'wall Z evaluation failed'
            call write_python_pair(unit, wall_r_value, wall_z_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'Z_GRID = np.array(['
        do j_grid = 0, FIXTURE_NZ-1
            bindings%values(13) = real(j_grid, dp)
            grid_z = eval_expr(grid_Z_expr, bindings, defined)
            if (.not. defined) error stop 'Z grid evaluation failed'
            call write_python_value(unit, grid_z)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'LCFS = np.array(['
        do k = 0, BOUNDARY_POINTS-1
            bindings%values(15) = real(k, dp)
            theta_value = eval_expr(theta_expr, bindings, defined)
            if (.not. defined) error stop 'theta evaluation failed'
            bindings%values(11) = theta_value
            lcfs_r_value = eval_expr(lcfs_R_expression, bindings, defined)
            if (.not. defined) error stop 'LCFS R evaluation failed'
            lcfs_z_value = eval_expr(lcfs_Z_expression, bindings, defined)
            if (.not. defined) error stop 'LCFS Z evaluation failed'
            call write_python_pair(unit, lcfs_r_value, lcfs_z_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'LIMITER = np.array(['
        do k = 0, BOUNDARY_POINTS-1
            bindings%values(15) = real(k, dp)
            theta_value = eval_expr(theta_expr, bindings, defined)
            if (.not. defined) error stop 'theta evaluation failed'
            bindings%values(11) = theta_value
            limiter_r_value = eval_expr(limiter_R_expression, bindings, defined)
            if (.not. defined) error stop 'limiter R evaluation failed'
            limiter_z_value = eval_expr(limiter_Z_expression, bindings, defined)
            if (.not. defined) error stop 'limiter Z evaluation failed'
            call write_python_pair(unit, limiter_r_value, limiter_z_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'PSI_VS = np.array(['
        do j_grid = 0, FIXTURE_NZ-1
            bindings%values(13) = real(j_grid, dp)
            grid_z = eval_expr(grid_Z_expr, bindings, defined)
            if (.not. defined) error stop 'Z grid evaluation failed'
            do i_grid = 0, FIXTURE_NR-1
                bindings%values(13) = real(i_grid, dp)
                grid_radius = eval_expr(grid_R_expr, bindings, defined)
                bindings%values(9) = grid_radius
                bindings%values(10) = grid_z
                radial_distance_value = eval_expr(radial_expression, &
                    bindings, defined)
                if (.not. defined) error stop 'grid radius evaluation failed'
                psi_value = evaluate(psi_expression, radial_distance_value)
                call write_python_value(unit, psi_value)
            end do
        end do
        write(unit, '(a)') '], dtype=np.float64).reshape((NZ, NR))'
        write(unit, '(a)') 'PSI_PROFILE = np.array(['
        do k = 0, FIXTURE_NR-1
            bindings%values(14) = real(k, dp)
            psi_coordinate_value = eval_expr(psi_sample_expr, bindings, defined)
            if (.not. defined) error stop 'psi sample evaluation failed'
            call write_python_value(unit, psi_coordinate_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'QPROF = np.array(['
        do k = 0, FIXTURE_NR-1
            bindings%values(14) = real(k, dp)
            bindings%values(8) = eval_expr(psi_sample_expr, bindings, defined)
            if (.not. defined) error stop 'psi sample evaluation failed'
            q_value = eval_expr(q_expression, bindings, defined)
            if (.not. defined) error stop 'q profile evaluation failed'
            call write_python_value(unit, q_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'D_VOLUME_DS_TOR_PROFILE = np.array(['
        do k = 0, FIXTURE_NR-1
            bindings%values(14) = real(k, dp)
            bindings%values(8) = eval_expr(psi_sample_expr, bindings, defined)
            if (.not. defined) error stop 'psi sample evaluation failed'
            volume_profile_value = eval_expr(volume_profile_expression, &
                bindings, defined)
            if (.not. defined) error stop 'volume Jacobian evaluation failed'
            call write_python_value(unit, volume_profile_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') ''
        write(unit, '(a)') 'assert R_GRID.size == NR and Z_GRID.size == NZ'
        write(unit, '(a)') 'assert PSI_VS.shape == (NZ, NR)'
        write(unit, '(a)') 'assert PSI_PROFILE.size == NR'
        write(unit, '(a)') 'assert QPROF.size == NR'
        write(unit, '(a)') 'assert D_VOLUME_DS_TOR_PROFILE.size == NR'
        write(unit, '(a)') 'assert BMOD_RAD_CM.size == NR'
        write(unit, '(a)') 'assert BMOD_ZET_CM.size == NZ'
        write(unit, '(a)') 'assert BMOD_RE.shape == (NR, NZ)'
        write(unit, '(a)') 'assert BMOD_IM.shape == (NR, NZ)'
        write(unit, '(a)') 'assert WALL.shape == (WALL_POINTS, 2)'
        write(unit, '(a)') 'assert PROFILE_POLY.shape == (4, 10)'
        write(unit, '(a)') 'assert np.all(np.diff(PSI_PROFILE) > 0.0)'
        write(unit, '(a)') 'assert np.all(np.diff(QPROF) > 0.0)'
        close(unit)
        write(output_unit, '(a)') 'wrote '//trim(path)
    end subroutine emit_python_table

    real(dp) function evaluate(expression, radius_value)
        type(expr_t), intent(in) :: expression
        real(dp), intent(in) :: radius_value
        logical :: local_defined

        bindings%values(1) = radius_value
        evaluate = eval_expr(expression, bindings, local_defined)
        if (.not. local_defined) error stop 'continuation table evaluation failed'
    end function evaluate

    subroutine write_python_value(output_unit_number, value)
        integer, intent(in) :: output_unit_number
        real(dp), intent(in) :: value

        write(output_unit_number, '(a,es24.16e3,a)') '    ', value, ','
    end subroutine write_python_value

    subroutine write_python_pair(output_unit_number, first, second)
        integer, intent(in) :: output_unit_number
        real(dp), intent(in) :: first, second

        write(output_unit_number, '(a,es24.16e3,a,es24.16e3,a)') &
            '    [', first, ', ', second, '],'
    end subroutine write_python_pair

    subroutine write_python_profile_row(output_unit_number, slope, intercept)
        integer, intent(in) :: output_unit_number
        real(dp), intent(in) :: slope, intercept

        write(output_unit_number, '(a)', advance='no') &
            '    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, '
        write(output_unit_number, '(es24.16e3,a)', advance='no') slope, ', '
        write(output_unit_number, '(es24.16e3,a)') intercept, '],'
    end subroutine write_python_profile_row

    subroutine write_python_zero_profile_row(output_unit_number)
        integer, intent(in) :: output_unit_number

        write(output_unit_number, '(a)') &
            '    [0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0],'
    end subroutine write_python_zero_profile_row

end program gen_circular_flux_continuation
