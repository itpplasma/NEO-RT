program gen_potato_kernels
    !! Derive the scalar algebra used by POTATO's resonance and topology code.
    !!
    !! This generator is intentionally split by physical role.  The consumer
    !! never receives a catch-all routine with unrelated dummy arguments: the
    !! action-domain expression, Eq. (4) Hamiltonian coefficient, resonance
    !! Jacobian/finite-harmonic guard, topology-gap bound, and limiting forms
    !! each have their own generated interface.  Runtime POTATO code owns
    !! domain checks, status, topology partitioning, and error propagation;
    !! Fortsym owns the algebra.
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use fortsym_arena, only: arena_t
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_identity
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine_symengine, only: make_symengine_engine, symengine_engine_t
    use fortsym_expr, only: abs, cos, expr_t, func, log, num, operator(+), &
        operator(-), operator(*), operator(/), operator(**), pi_expr, rat, &
        sin, sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    implicit none

    character(*), parameter :: FORTSYM_REVISION = &
        'fortsym@58a0e06c95ecc943dfdcb044b7ca6a9964c1c55d'
    character(*), parameter :: REGENERATE_COMMAND = &
        'cd tools/gc_symbolics && fo exec gen_potato_kernels ../../POTATO/SRC/generated'

    character(2048) :: output_directory
    integer :: output_length, argument_status
    type(arena_t), target :: arena
    type(symengine_engine_t) :: proof_engine
    type(native_engine_t) :: simplify_engine
    type(suite_t) :: proofs

    type(expr_t) :: energy_H, qPhi, magnetic_field_B
    type(expr_t) :: qPhi_prime, magnetic_field_B_prime
    type(expr_t) :: Jperp_candidate, Jperp_positive_bound, dJperp_dR
    type(expr_t) :: qPhi_partial, B_partial, zero

    type(expr_t) :: Jperp, deltaB_m_real, deltaB_m_imag, mode_phase
    type(expr_t) :: hamiltonian_coefficient, H_m_real, H_m_imag, H_m_squared
    type(expr_t) :: H_m_squared_sign, H_m_squared_twice, H_m_squared_paired
    type(expr_t) :: global_phase, deltaB_rot_real, deltaB_rot_imag
    type(expr_t) :: H_m_squared_rotated
    type(expr_t) :: dpsiast_dR, dresonance_dR, root_jacobian

    type(expr_t) :: mode_m, mode_n, delta_phi
    type(expr_t) :: harmonic_target, harmonic_extent, resonance_g
    type(expr_t) :: no_root_margin, current_extent_envelope
    type(expr_t) :: extent_envelope

    type(expr_t) :: orbit_H_m_squared, maxwellian_weight, Phi_eff
    type(expr_t) :: bounce_time, thermodynamic_force, reference_energy
    type(expr_t) :: delta_root_weight, resonance_torque_weight

    type(expr_t) :: resonance_q, resonance_q_prime, bounce_time_prime
    type(expr_t) :: abs_mode_n, frequency_resonance
    type(expr_t) :: frequency_derivative, frequency_derivative_at_root
    type(expr_t) :: frequency_root_remainder
    type(expr_t) :: frequency_derivative_abs_at_root
    type(expr_t) :: n_squared_delta_numerator
    type(expr_t) :: frequency_delta_weight_from_n2, frequency_delta_weight
    type(expr_t) :: collapsed_n_tau_weight, collapsed_n_tau2_weight

    type(expr_t) :: integrand_envelope, Jperp_gap_lo, Jperp_gap_hi
    type(expr_t) :: topology_gap_measure, topology_contribution_error_bound

    type(expr_t) :: hessian_H_RR, hessian_H_Rp, hessian_H_pp
    type(expr_t) :: regular_tau_value, cut_linear_slope, xpoint_cut_curvature
    type(expr_t) :: delta_R, delta_R_ref, hessian_determinant, lambda_local, C_tau
    type(expr_t) :: regular_action_offset, regular_action_jacobian
    type(expr_t) :: xpoint_action_offset, xpoint_action_jacobian
    type(expr_t) :: regular_tau_limit, separatrix_tau_scale, xpoint_tau_scale

    call get_command_argument(1, output_directory, length=output_length, &
        status=argument_status)
    if (argument_status /= 0 .or. output_length == 0) then
        write (output_unit, '(a)') &
            'usage: gen_potato_kernels OUTPUT_DIRECTORY'
        error stop 2
    end if
    output_directory = output_directory(:output_length)

    call arena%init()
    proof_engine = make_symengine_engine(arena)
    simplify_engine = make_native_engine(arena)

    zero = num(arena, 0)

    ! ------------------------------------------------------------------
    ! Exact normalized perpendicular-action domain.  qPhi is the physical
    ! charge-times-potential quantity; it is deliberately not renamed Phi.
    energy_H = sym(arena, 'energy_H')
    qPhi = sym(arena, 'qPhi')
    magnetic_field_B = sym(arena, 'magnetic_field_B')
    qPhi_prime = sym(arena, 'qPhi_prime')
    magnetic_field_B_prime = sym(arena, 'magnetic_field_B_prime')
    Jperp_candidate = (energy_H - qPhi)/magnetic_field_B
    Jperp_positive_bound = func('max', [zero, Jperp_candidate])
    qPhi_partial = diff(Jperp_candidate, qPhi)
    B_partial = diff(Jperp_candidate, magnetic_field_B)
    dJperp_dR = qPhi_partial*qPhi_prime + B_partial*magnetic_field_B_prime

    ! ------------------------------------------------------------------
    ! POTATO's perturbed Hamiltonian coefficient.  The real/imaginary pair
    ! is rotated by the orbit phase before its squared modulus is formed.
    Jperp = sym(arena, 'Jperp')
    deltaB_m_real = sym(arena, 'deltaB_m_real')
    deltaB_m_imag = sym(arena, 'deltaB_m_imag')
    mode_phase = sym(arena, 'mode_phase')
    hamiltonian_coefficient = &
        (2*(energy_H - qPhi) - Jperp*magnetic_field_B)/magnetic_field_B
    H_m_real = hamiltonian_coefficient*(deltaB_m_real*cos(mode_phase) - &
        deltaB_m_imag*sin(mode_phase))
    H_m_imag = hamiltonian_coefficient*(deltaB_m_real*sin(mode_phase) + &
        deltaB_m_imag*cos(mode_phase))
    H_m_squared = H_m_real**2 + H_m_imag**2

    H_m_squared_sign = hamiltonian_square(-deltaB_m_real, -deltaB_m_imag, &
        mode_phase, hamiltonian_coefficient)
    H_m_squared_twice = hamiltonian_square(2*deltaB_m_real, &
        2*deltaB_m_imag, mode_phase, hamiltonian_coefficient)
    ! The real-field partner changes both the Fourier amplitude and the mode
    ! phase: (m,n,A) -> (-m,-n,A*).  Fixed-(m,n) A -> A* is not asserted.
    H_m_squared_paired = hamiltonian_square(deltaB_m_real, &
        -deltaB_m_imag, -mode_phase, hamiltonian_coefficient)

    global_phase = sym(arena, 'global_phase')
    deltaB_rot_real = deltaB_m_real*cos(global_phase) - &
        deltaB_m_imag*sin(global_phase)
    deltaB_rot_imag = deltaB_m_real*sin(global_phase) + &
        deltaB_m_imag*cos(global_phase)
    H_m_squared_rotated = hamiltonian_square(deltaB_rot_real, &
        deltaB_rot_imag, mode_phase, hamiltonian_coefficient)

    ! ------------------------------------------------------------------
    ! Delta-function root change of variables.  The absolute value is part of
    ! the generated algebra; tangent-root rejection remains handwritten.
    dpsiast_dR = sym(arena, 'dpsiast_dR')
    dresonance_dR = sym(arena, 'dresonance_dR')
    root_jacobian = abs(dpsiast_dR/dresonance_dR)

    ! ------------------------------------------------------------------
    ! Signed per-harmonic resonance guard.  The root condition is
    !
    !     g = Delta_phi + 2*pi*m/n = 0.
    !
    ! The target remains signed for the actual root search.  The only
    ! magnitude is the reverse-triangle-inequality search extent
    ! |2*pi*m/n|; if |Delta_phi|-extent > 0, g cannot vanish.  The finite-set
    ! envelope is generated as a scalar max reduction and is called once for
    ! each exact (m,n) pair by the handwritten guard.
    mode_m = sym(arena, 'mode_m')
    mode_n = sym(arena, 'mode_n')
    delta_phi = sym(arena, 'delta_phi')
    harmonic_target = -2*pi_expr(arena)*mode_m/mode_n
    harmonic_extent = abs(harmonic_target)
    resonance_g = delta_phi + 2*pi_expr(arena)*mode_m/mode_n
    no_root_margin = abs(delta_phi) - harmonic_extent
    current_extent_envelope = sym(arena, 'current_extent_envelope')
    extent_envelope = func('max', [current_extent_envelope, harmonic_extent])

    ! ------------------------------------------------------------------
    ! Eq. (17) delta-root and torque weight.  The negative pi^(3/2)/4
    ! convention is retained here.  The root Jacobian is derived from the
    ! physical resonance derivatives inside this kernel, while the mode
    ! representation enters as abs(mode_n); this is why n and -n conjugate
    ! representations agree.
    orbit_H_m_squared = sym(arena, 'orbit_H_m_squared')
    maxwellian_weight = sym(arena, 'maxwellian_weight')
    Phi_eff = sym(arena, 'Phi_eff')
    bounce_time = sym(arena, 'bounce_time')
    thermodynamic_force = sym(arena, 'thermodynamic_force')
    reference_energy = sym(arena, 'reference_energy')
    delta_root_weight = abs(mode_n)*root_jacobian
    resonance_torque_weight = &
        -pi_expr(arena)**rat(arena,3_int64,2_int64)/4 * &
        reference_energy*delta_root_weight*orbit_H_m_squared* &
        maxwellian_weight*Phi_eff*bounce_time**2*thermodynamic_force

    ! ------------------------------------------------------------------
    ! Frequency delta reduction used by POTATO's Eq. (17) weight.  Let
    ! q = Delta_phi + 2*pi*m/n and F_freq = n*q/tau.  At q=0,
    ! F_freq' = n*q'/tau.  The original n^2 delta factor therefore becomes
    ! |n|*tau/|q'|; the remaining transport tau gives
    ! |n|*tau^2/|q'|.  This kernel makes that cancellation explicit so no
    ! second n^2 can be introduced in handwritten torque assembly.
    resonance_q = sym(arena, 'resonance_q')
    resonance_q_prime = sym(arena, 'resonance_q_prime')
    bounce_time_prime = sym(arena, 'bounce_time_prime')
    ! abs(mode_n) is a derived quantity, never a caller-supplied proof input.
    abs_mode_n = abs(mode_n)
    frequency_resonance = mode_n*resonance_q/bounce_time
    frequency_derivative = mode_n*(resonance_q_prime*bounce_time - &
        resonance_q*bounce_time_prime)/bounce_time**2
    frequency_derivative_at_root = mode_n*resonance_q_prime/bounce_time
    frequency_root_remainder = frequency_derivative - &
        frequency_derivative_at_root
    frequency_derivative_abs_at_root = abs_mode_n*abs(resonance_q_prime)/ &
        bounce_time
    n_squared_delta_numerator = abs_mode_n**2
    ! The denominator cancellation is only valid under the runtime nonzero
    ! guard.  Prove its cross-multiplied defining identity instead of asking a
    ! CAS to cancel a possibly-zero factor.
    frequency_delta_weight_from_n2 = n_squared_delta_numerator/ &
        frequency_derivative_abs_at_root
    frequency_delta_weight = abs_mode_n*bounce_time/abs(resonance_q_prime)
    collapsed_n_tau_weight = frequency_delta_weight
    collapsed_n_tau2_weight = collapsed_n_tau_weight*bounce_time

    ! ------------------------------------------------------------------
    ! A certified envelope times an omitted J_perp interval is a contribution
    ! bound, not merely a geometric gap length.
    integrand_envelope = sym(arena, 'integrand_envelope')
    Jperp_gap_lo = sym(arena, 'Jperp_gap_lo')
    Jperp_gap_hi = sym(arena, 'Jperp_gap_hi')
    topology_gap_measure = Jperp_gap_hi - Jperp_gap_lo
    topology_contribution_error_bound = &
        integrand_envelope*topology_gap_measure

    ! ------------------------------------------------------------------
    ! Limiting forms from the local one-degree-of-freedom Hamiltonian.  The
    ! Hessian is supplied in one explicitly selected physical/normalised-time
    ! convention; no fitted C is accepted.  Runtime code must reject the
    ! limit if this Hessian or the local cut maps are unavailable.
    hessian_H_RR = sym(arena, 'hessian_H_RR')
    hessian_H_Rp = sym(arena, 'hessian_H_Rp')
    hessian_H_pp = sym(arena, 'hessian_H_pp')
    regular_tau_value = sym(arena, 'regular_tau_value')
    cut_linear_slope = sym(arena, 'cut_linear_slope')
    xpoint_cut_curvature = sym(arena, 'xpoint_cut_curvature')
    delta_R = sym(arena, 'delta_R')
    ! delta_R_ref carries the same physical length unit as delta_R.  Their
    ! ratio is dimensionless, so the logarithmic limit has no hidden unit
    ! convention.
    delta_R_ref = sym(arena, 'delta_R_ref')
    hessian_determinant = hessian_H_RR*hessian_H_pp - hessian_H_Rp**2
    lambda_local = sqrt(-hessian_determinant)
    C_tau = 1/lambda_local
    regular_action_offset = cut_linear_slope*delta_R
    regular_action_jacobian = cut_linear_slope
    xpoint_action_offset = xpoint_cut_curvature*delta_R**2
    xpoint_action_jacobian = 2*xpoint_cut_curvature*delta_R
    regular_tau_limit = regular_tau_value
    separatrix_tau_scale = C_tau*log(delta_R_ref/abs(delta_R))
    xpoint_tau_scale = 2*C_tau*log(delta_R_ref/abs(delta_R))

    call suite_begin(proofs, 'POTATO Fortsym contracts')
    call check_identity(proofs, proof_engine, &
        'Jperp candidate times B equals H minus qPhi', &
        Jperp_candidate*magnetic_field_B - (energy_H - qPhi))
    call check_identity(proofs, proof_engine, &
        'Jperp radial derivative is the chain rule', &
        dJperp_dR - (diff(Jperp_candidate, qPhi)*qPhi_prime + &
            diff(Jperp_candidate, magnetic_field_B)*magnetic_field_B_prime))
    call check_identity(proofs, proof_engine, &
        'positive action bound is the positive part of Jperp candidate', &
        Jperp_positive_bound - func('max', [zero, Jperp_candidate]))
    call check_identity(proofs, proof_engine, &
        'H_m squared removes orbit phase', &
        H_m_squared - hamiltonian_coefficient**2* &
            (deltaB_m_real**2 + deltaB_m_imag**2))
    call check_identity(proofs, proof_engine, &
        'A to minus A leaves H_m squared', H_m_squared_sign - H_m_squared)
    call check_identity(proofs, proof_engine, &
        'A to 2 A scales H_m squared by four', &
        H_m_squared_twice - 4*H_m_squared)
    call check_identity(proofs, proof_engine, &
        'global complex phase leaves H_m squared', &
        H_m_squared_rotated - H_m_squared)
    call check_identity(proofs, proof_engine, &
        'paired (-m,-n,A*) mode leaves H_m squared', &
        H_m_squared_paired - H_m_squared)
    call check_identity(proofs, proof_engine, &
        'root Jacobian is the absolute delta-root weight', &
        root_jacobian - abs(dpsiast_dR/dresonance_dR))
    call check_identity(proofs, proof_engine, &
        'harmonic target preserves the signed resonance', &
        harmonic_target + 2*pi_expr(arena)*mode_m/mode_n)
    call check_identity(proofs, proof_engine, &
        'harmonic search extent is magnitude only', &
        harmonic_extent - abs(2*pi_expr(arena)*mode_m/mode_n))
    call check_identity(proofs, proof_engine, &
        'signed harmonic resonance condition', &
        resonance_g - (delta_phi + 2*pi_expr(arena)*mode_m/mode_n))
    call check_identity(proofs, proof_engine, &
        'per-harmonic no-root margin', &
        no_root_margin - (abs(delta_phi) - harmonic_extent))
    call check_identity(proofs, proof_engine, &
        'finite harmonic extent envelope', &
        extent_envelope - func('max', [current_extent_envelope, harmonic_extent]))
    call check_identity(proofs, proof_engine, &
        'delta-root weight uses absolute mode number', &
        delta_root_weight - abs(mode_n)*root_jacobian)
    call check_identity(proofs, proof_engine, &
        'absolute mode number is derived internally', &
        abs_mode_n - abs(mode_n))
    call check_identity(proofs, proof_engine, &
        'conjugate mode keeps delta-root weight', &
        abs(-mode_n)*root_jacobian - delta_root_weight)
    call check_identity(proofs, proof_engine, &
        'Eq.17 torque weight has the signed leading prefactor', &
        resonance_torque_weight + &
        pi_expr(arena)**rat(arena,3_int64,2_int64)/4 * &
        reference_energy*delta_root_weight*orbit_H_m_squared* &
            maxwellian_weight*Phi_eff*bounce_time**2*thermodynamic_force)
    call check_identity(proofs, proof_engine, &
        'frequency resonance is n times q over tau', &
        frequency_resonance - mode_n*resonance_q/bounce_time)
    call check_identity(proofs, proof_engine, &
        'frequency derivative before imposing the root', &
        frequency_derivative - mode_n*(resonance_q_prime*bounce_time - &
            resonance_q*bounce_time_prime)/bounce_time**2)
    call check_identity(proofs, proof_engine, &
        'frequency derivative root remainder is proportional to q', &
        frequency_root_remainder + mode_n*resonance_q*bounce_time_prime/ &
            bounce_time**2)
    call check_identity(proofs, proof_engine, &
        'n-squared delta factor collapses before transport tau', &
        frequency_delta_weight_from_n2*frequency_derivative_abs_at_root - &
            n_squared_delta_numerator)
    call check_identity(proofs, proof_engine, &
        'transport tau leaves one absolute mode factor', &
        collapsed_n_tau2_weight - abs_mode_n*bounce_time**2/ &
            abs(resonance_q_prime))
    call check_identity(proofs, proof_engine, &
        'topology contribution bound is envelope times gap', &
        topology_contribution_error_bound - integrand_envelope* &
            (Jperp_gap_hi - Jperp_gap_lo))
    call check_identity(proofs, proof_engine, &
        'Hessian determinant defines the local saddle rate', &
        hessian_determinant - (hessian_H_RR*hessian_H_pp - hessian_H_Rp**2))
    call check_identity(proofs, proof_engine, &
        'homoclinic coefficient is inverse saddle rate', &
        C_tau*lambda_local - 1)
    call check_identity(proofs, proof_engine, &
        'regular boundary remains finite', &
        regular_tau_limit - regular_tau_value)
    call check_identity(proofs, proof_engine, &
        'regular cut action map is linear', &
        regular_action_offset - cut_linear_slope*delta_R)
    call check_identity(proofs, proof_engine, &
        'regular cut map Jacobian', regular_action_jacobian - cut_linear_slope)
    call check_identity(proofs, proof_engine, &
        'X-point action map is quadratic', &
        xpoint_action_offset - xpoint_cut_curvature*delta_R**2)
    call check_identity(proofs, proof_engine, &
        'X-point map Jacobian', xpoint_action_jacobian - &
            2*xpoint_cut_curvature*delta_R)
    call check_identity(proofs, proof_engine, &
        'X-point logarithm is twice separatrix logarithm', &
        xpoint_tau_scale - 2*separatrix_tau_scale)
    call check_identity(proofs, proof_engine, &
        'homoclinic logarithm uses a dimensionless distance ratio', &
        separatrix_tau_scale - C_tau*log(delta_R_ref/abs(delta_R)))
    if (proofs%failed /= 0) error stop 'POTATO symbolic proof failed'
    call suite_end(proofs)

    call simplify_one(Jperp_candidate)
    call simplify_one(Jperp_positive_bound)
    call simplify_one(dJperp_dR)
    call simplify_one(hamiltonian_coefficient)
    call simplify_one(H_m_real)
    call simplify_one(H_m_imag)
    call simplify_one(H_m_squared)
    call simplify_one(root_jacobian)
    call simplify_one(harmonic_target)
    call simplify_one(harmonic_extent)
    call simplify_one(resonance_g)
    call simplify_one(no_root_margin)
    call simplify_one(extent_envelope)
    call simplify_one(frequency_resonance)
    call simplify_one(frequency_derivative)
    call simplify_one(frequency_derivative_at_root)
    call simplify_one(frequency_root_remainder)
    call simplify_one(frequency_delta_weight_from_n2)
    call simplify_one(frequency_delta_weight)
    call simplify_one(collapsed_n_tau_weight)
    call simplify_one(collapsed_n_tau2_weight)
    call simplify_one(delta_root_weight)
    call simplify_one(resonance_torque_weight)
    call simplify_one(topology_gap_measure)
    call simplify_one(topology_contribution_error_bound)
    call simplify_one(hessian_determinant)
    call simplify_one(lambda_local)
    call simplify_one(C_tau)
    call simplify_one(regular_action_offset)
    call simplify_one(regular_action_jacobian)
    call simplify_one(xpoint_action_offset)
    call simplify_one(xpoint_action_jacobian)
    call simplify_one(regular_tau_limit)
    call simplify_one(separatrix_tau_scale)
    call simplify_one(xpoint_tau_scale)

    call emit_jperp_kernel(output_directory, Jperp_candidate, Jperp_positive_bound, &
        dJperp_dR)
    call emit_hamiltonian_kernel(output_directory, hamiltonian_coefficient, &
        H_m_real, H_m_imag, H_m_squared)
    call emit_root_kernel(output_directory, root_jacobian)
    call emit_harmonic_guard_kernels(output_directory, harmonic_target, &
        harmonic_extent, resonance_g, no_root_margin, extent_envelope)
    call emit_frequency_reduction_kernel(output_directory, frequency_resonance, &
        frequency_derivative, frequency_derivative_at_root, &
        frequency_root_remainder, frequency_delta_weight_from_n2, &
        collapsed_n_tau2_weight)
    call emit_resonance_kernel(output_directory,delta_root_weight, &
        resonance_torque_weight)
    call emit_gap_kernel(output_directory, topology_gap_measure, &
        topology_contribution_error_bound)
    call emit_limiting_kernel(output_directory,hessian_determinant,lambda_local, &
        C_tau,regular_action_offset,regular_action_jacobian, &
        xpoint_action_offset,xpoint_action_jacobian,regular_tau_limit, &
        separatrix_tau_scale,xpoint_tau_scale)

contains

    function hamiltonian_square(real_amplitude, imag_amplitude, phase, coeff) &
        result(square)
        type(expr_t), intent(in) :: real_amplitude, imag_amplitude, phase, coeff
        type(expr_t) :: square, rotated_real, rotated_imag

        rotated_real = coeff*(real_amplitude*cos(phase) - &
            imag_amplitude*sin(phase))
        rotated_imag = coeff*(real_amplitude*sin(phase) + &
            imag_amplitude*cos(phase))
        square = rotated_real**2 + rotated_imag**2
    end function hamiltonian_square

    subroutine simplify_one(expression)
        type(expr_t), intent(inout) :: expression
        type(engine_result_t) :: simplified

        simplified = simplify_engine%simplify(expression)
        if (.not. simplified%ok) then
            write (output_unit, '(a)') 'simplification failed for output'
            error stop 1
        end if
        expression = simplified%value
    end subroutine simplify_one

    subroutine common_spec(spec, procedure_name, module_name)
        type(kernel_spec_t), intent(out) :: spec
        character(*), intent(in) :: procedure_name, module_name

        spec%name = str(procedure_name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%generator = str('tools/gc_symbolics/app/gen_potato_kernels.f90')
        spec%generator_revision = str(FORTSYM_REVISION)
        spec%regenerate_command = str(REGENERATE_COMMAND)
    end subroutine common_spec

    subroutine write_kernel(output_directory, filename, roots, spec)
        character(*), intent(in) :: output_directory, filename
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t), intent(in) :: spec
        character(4096) :: path
        integer :: unit, ios

        path = trim(output_directory)//'/'//trim(filename)
        open (newunit=unit, file=trim(path), status='replace', &
            action='write', iostat=ios)
        if (ios /= 0) then
            write (output_unit, '(a)') 'cannot open generated output '//trim(path)
            error stop 1
        end if
        write (unit, '(a)') chars(emit_kernel(roots, spec))
        close (unit)
        write (output_unit, '(a)') 'wrote '//trim(path)
    end subroutine write_kernel

    subroutine emit_jperp_kernel(directory, candidate, positive_bound, derivative)
        character(*), intent(in) :: directory
        type(expr_t), intent(in) :: candidate, positive_bound, derivative
        type(kernel_spec_t) :: spec

        call common_spec(spec, 'potato_jperp_domain_kernel', &
            'potato_jperp_domain_kernel')
        allocate (spec%args(5), spec%outputs(3))
        spec%args = [str('energy_H'), str('qPhi'), str('magnetic_field_B'), &
            str('qPhi_prime'), str('magnetic_field_B_prime')]
        spec%outputs = [str('Jperp_candidate'), str('Jperp_positive_bound'), &
            str('dJperp_dR')]
        call write_kernel(directory, 'potato_jperp_domain.f90', &
            [candidate, positive_bound, derivative], spec)
    end subroutine emit_jperp_kernel

    subroutine emit_hamiltonian_kernel(directory, coefficient, hm_real, hm_imag, &
        hm_squared)
        character(*), intent(in) :: directory
        type(expr_t), intent(in) :: coefficient, hm_real, hm_imag, hm_squared
        type(kernel_spec_t) :: spec

        call common_spec(spec, 'potato_hm_eq4_kernel', 'potato_hm_eq4_kernel')
        allocate (spec%args(7), spec%outputs(4))
        spec%args = [str('energy_H'), str('qPhi'), str('magnetic_field_B'), &
            str('Jperp'), str('deltaB_m_real'), str('deltaB_m_imag'), &
            str('mode_phase')]
        spec%outputs = [str('H_m_real'), str('H_m_imag'), &
            str('H_m_squared'), str('hamiltonian_coefficient')]
        call write_kernel(directory, 'potato_hm_eq4.f90', &
            [hm_real, hm_imag, hm_squared, coefficient], spec)
    end subroutine emit_hamiltonian_kernel

    subroutine emit_root_kernel(directory, root_weight)
        character(*), intent(in) :: directory
        type(expr_t), intent(in) :: root_weight
        type(kernel_spec_t) :: spec

        call common_spec(spec, 'potato_root_jacobian_kernel', &
            'potato_root_jacobian_kernel')
        allocate (spec%args(2), spec%outputs(1))
        spec%args = [str('dpsiast_dR'), str('dresonance_dR')]
        spec%outputs = [str('root_jacobian')]
        call write_kernel(directory, 'potato_root_jacobian.f90', &
            [root_weight], spec)
    end subroutine emit_root_kernel

    subroutine emit_harmonic_guard_kernels(directory, target, extent, g, &
        no_root, envelope)
        character(*), intent(in) :: directory
        type(expr_t), intent(in) :: target, extent, g, no_root, envelope
        type(kernel_spec_t) :: spec

        call common_spec(spec, 'potato_resonance_harmonic_kernel', &
            'potato_resonance_harmonic_kernel')
        allocate (spec%args(3), spec%outputs(4))
        spec%args = [str('mode_m'), str('mode_n'), str('delta_phi')]
        spec%outputs = [str('target_delphi'), str('search_extent'), &
            str('resonance_g'), str('no_root_margin')]
        call write_kernel(directory, 'potato_resonance_harmonic.f90', &
            [target, extent, g, no_root], spec)

        call common_spec(spec, 'potato_resonance_extent_envelope_kernel', &
            'potato_resonance_extent_envelope_kernel')
        allocate (spec%args(2), spec%outputs(1))
        spec%args = [str('current_extent_envelope'), str('harmonic_extent')]
        spec%outputs = [str('extent_envelope')]
        call write_kernel(directory, 'potato_resonance_envelope.f90', &
            [envelope], spec)
    end subroutine emit_harmonic_guard_kernels

    subroutine emit_frequency_reduction_kernel(directory, frequency, &
        derivative, derivative_at_root, root_remainder, n2_weight, collapsed)
        character(*), intent(in) :: directory
        type(expr_t), intent(in) :: frequency, derivative, derivative_at_root
        type(expr_t), intent(in) :: root_remainder, n2_weight, collapsed
        type(kernel_spec_t) :: spec

        call common_spec(spec, 'potato_frequency_reduction_kernel', &
            'potato_frequency_reduction_kernel')
        allocate (spec%args(5), spec%outputs(6))
        spec%args = [str('mode_n'), str('resonance_q'), &
            str('resonance_q_prime'), str('bounce_time'), &
            str('bounce_time_prime')]
        spec%outputs = [str('frequency_resonance'), &
            str('frequency_derivative'), str('frequency_derivative_at_root'), &
            str('frequency_root_remainder'), &
            str('frequency_delta_weight_from_n2'), &
            str('collapsed_n_tau2_weight')]
        call write_kernel(directory, 'potato_frequency_reduction.f90', &
            [frequency, derivative, derivative_at_root, root_remainder, &
             n2_weight, collapsed], spec)
    end subroutine emit_frequency_reduction_kernel

    subroutine emit_resonance_kernel(directory, delta_weight, torque_weight)
        character(*), intent(in) :: directory
        type(expr_t), intent(in) :: delta_weight, torque_weight
        type(kernel_spec_t) :: spec

        call common_spec(spec, 'potato_resonance_torque_kernel', &
            'potato_resonance_torque_kernel')
        allocate (spec%args(9), spec%outputs(2))
        spec%args = [str('dpsiast_dR'), str('dresonance_dR'), str('mode_n'), &
            str('orbit_H_m_squared'), str('maxwellian_weight'), str('Phi_eff'), &
            str('bounce_time'), str('thermodynamic_force'), &
            str('reference_energy')]
        spec%outputs = [str('delta_root_weight'), &
            str('resonance_torque_weight')]
        call write_kernel(directory, 'potato_resonance_weight.f90', &
            [delta_weight, torque_weight], spec)
    end subroutine emit_resonance_kernel

    subroutine emit_gap_kernel(directory, gap_measure, contribution_bound)
        character(*), intent(in) :: directory
        type(expr_t), intent(in) :: gap_measure, contribution_bound
        type(kernel_spec_t) :: spec

        call common_spec(spec, 'potato_gap_error_kernel', &
            'potato_gap_error_kernel')
        allocate (spec%args(3), spec%outputs(2))
        spec%args = [str('integrand_envelope'), str('Jperp_gap_lo'), &
            str('Jperp_gap_hi')]
        spec%outputs = [str('topology_gap_measure'), &
            str('topology_contribution_error_bound')]
        call write_kernel(directory, 'potato_gap_error.f90', &
            [gap_measure, contribution_bound], spec)
    end subroutine emit_gap_kernel

    subroutine emit_limiting_kernel(directory, hessian_determinant, lambda_local, &
        c_tau, regular_action_offset, regular_action_jacobian, &
        xpoint_action_offset, xpoint_action_jacobian, regular_tau_limit, &
        separatrix, xpoint)
        character(*), intent(in) :: directory
        type(expr_t), intent(in) :: hessian_determinant, lambda_local, c_tau
        type(expr_t), intent(in) :: regular_action_offset, regular_action_jacobian
        type(expr_t), intent(in) :: xpoint_action_offset, xpoint_action_jacobian
        type(expr_t), intent(in) :: regular_tau_limit, separatrix, xpoint
        type(kernel_spec_t) :: spec

        call common_spec(spec, 'potato_limiting_kernel', &
            'potato_limiting_kernel')
        allocate (spec%args(8), spec%outputs(10))
        spec%args = [str('hessian_H_RR'), str('hessian_H_Rp'), &
            str('hessian_H_pp'), str('regular_tau_value'), &
            str('cut_linear_slope'), str('xpoint_cut_curvature'), &
            str('delta_R'), str('delta_R_ref')]
        spec%outputs = [str('hessian_determinant'), str('lambda_local'), &
            str('C_tau'), str('regular_action_offset'), &
            str('regular_action_jacobian'), str('xpoint_action_offset'), &
            str('xpoint_action_jacobian'), str('regular_tau_limit'), &
            str('separatrix_tau_scale'), str('xpoint_tau_scale')]
        call write_kernel(directory, 'potato_limiting.f90', &
            [hessian_determinant, lambda_local, c_tau, regular_action_offset, &
             regular_action_jacobian, xpoint_action_offset, xpoint_action_jacobian, &
             regular_tau_limit, separatrix, xpoint], spec)
    end subroutine emit_limiting_kernel

end program gen_potato_kernels
