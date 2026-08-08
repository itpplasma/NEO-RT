program gen_gc_homoclinic_limit
    !! Derive local homoclinic logarithmic coefficients and finite parts.
    !!
    !! POTATO/TEX/equilmaxw.tex distinguishes a usual separatrix crossing
    !! from an X-point by one versus two saddle legs.  Here epsilon is the
    !! positive dimensionless hyperbolic approach coordinate.  A k-leg
    !! quantity has singular part -k*C*log(epsilon), with k=1 or 2.
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use fortsym_arena, only: arena_t
    use fortsym_check, only: check_identity, suite_begin, suite_end, suite_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine_symengine, only: make_symengine_engine, &
        symengine_engine_t
    use fortsym_expr, only: abs, expr_t, log, operator(+), operator(-), &
        operator(*), operator(/), operator(**), rat, sqrt, sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    implicit none

    character(*), parameter :: FORTSYM_REVISION = &
        'fortsym@545788453a204d58705f735b519c3863c2f734c8'
    character(*), parameter :: GENERATOR_PATH = &
        'tools/gc_symbolics/app/gen_gc_homoclinic_limit.f90'
    character(*), parameter :: REGENERATE_COMMAND = &
        'cd tools/gc_symbolics && fo exec gen_gc_homoclinic_limit '// &
        '../../src/generated'

    type(arena_t), target :: arena
    type(symengine_engine_t) :: proof_engine
    type(native_engine_t) :: simplify_engine
    type(suite_t) :: proofs
    type(engine_result_t) :: simplified
    type(expr_t) :: coordinate_q, momentum_p
    type(expr_t) :: h_qq_input, h_qp_input, h_pp_input
    type(expr_t) :: normal_hamiltonian, h_q, h_p
    type(expr_t) :: h_qq, h_qp, h_pp, q_dot, p_dot
    type(expr_t) :: flow_qq, flow_qp, flow_pq, flow_pp
    type(expr_t) :: hessian_determinant, flow_determinant
    type(expr_t) :: saddle_discriminant, saddle_rate
    type(expr_t) :: toroidal_rate_saddle, c_tau, c_phi
    type(expr_t) :: one_leg_tau_coefficient, two_leg_tau_coefficient
    type(expr_t) :: one_leg_phi_coefficient, two_leg_phi_coefficient
    type(expr_t) :: c_tau_section_reverse, c_phi_section_reverse
    type(expr_t) :: c_tau_toroidal_reverse, c_phi_toroidal_reverse
    type(expr_t) :: absolute_phi_reversal_defect
    type(expr_t) :: saddle_rate_floor, epsilon, epsilon_floor
    type(expr_t) :: nondegeneracy_margin, section_reversed_margin
    type(expr_t) :: epsilon_margin, section_reversed_epsilon_margin
    type(expr_t) :: saddle_rate_floor_margin, epsilon_floor_margin
    type(expr_t) :: tau_one_epsilon, tau_two_epsilon
    type(expr_t) :: phi_one_epsilon, phi_two_epsilon, log_epsilon
    type(expr_t) :: tau_finite_one, tau_finite_two
    type(expr_t) :: phi_finite_one, phi_finite_two
    type(expr_t) :: phi_finite_one_toroidal_reverse
    type(expr_t) :: phi_finite_two_toroidal_reverse
    type(expr_t) :: coefficient_roots(14), diagnostic_roots(8)
    type(expr_t) :: finite_part_roots(6)
    character(2048) :: output_directory
    integer :: argument_status, output_length

    call get_command_argument(1, output_directory, length=output_length, &
        status=argument_status)
    if (argument_status /= 0 .or. output_length == 0) then
        write (output_unit, '(a)') &
            'usage: gen_gc_homoclinic_limit OUTPUT_DIRECTORY'
        error stop 2
    end if
    output_directory = output_directory(:output_length)

    call arena%init()
    proof_engine = make_symengine_engine(arena)
    simplify_engine = make_native_engine(arena)

    coordinate_q = sym(arena, 'coordinate_q')
    momentum_p = sym(arena, 'momentum_p')
    h_qq_input = sym(arena, 'h_qq')
    h_qp_input = sym(arena, 'h_qp')
    h_pp_input = sym(arena, 'h_pp')
    toroidal_rate_saddle = sym(arena, 'toroidal_rate_saddle')
    saddle_rate_floor = sym(arena, 'saddle_rate_floor')
    epsilon = sym(arena, 'epsilon')
    epsilon_floor = sym(arena, 'epsilon_floor')
    tau_one_epsilon = sym(arena, 'tau_one_epsilon')
    tau_two_epsilon = sym(arena, 'tau_two_epsilon')
    phi_one_epsilon = sym(arena, 'phi_one_epsilon')
    phi_two_epsilon = sym(arena, 'phi_two_epsilon')

    ! All local derivatives originate from the quadratic Hamiltonian tree.
    normal_hamiltonian = rat(arena, 1_int64, 2_int64)* &
        h_qq_input*coordinate_q**2 + &
        h_qp_input*coordinate_q*momentum_p + &
        rat(arena, 1_int64, 2_int64)*h_pp_input*momentum_p**2
    h_q = diff(normal_hamiltonian, coordinate_q)
    h_p = diff(normal_hamiltonian, momentum_p)
    h_qq = diff(h_q, coordinate_q)
    h_qp = diff(h_q, momentum_p)
    h_pp = diff(h_p, momentum_p)
    q_dot = h_p
    p_dot = -h_q
    flow_qq = diff(q_dot, coordinate_q)
    flow_qp = diff(q_dot, momentum_p)
    flow_pq = diff(p_dot, coordinate_q)
    flow_pp = diff(p_dot, momentum_p)

    hessian_determinant = h_qq*h_pp-h_qp**2
    flow_determinant = flow_qq*flow_pp-flow_qp*flow_pq
    saddle_discriminant = -hessian_determinant
    saddle_rate = sqrt(saddle_discriminant)
    c_tau = 1/saddle_rate
    c_phi = toroidal_rate_saddle*c_tau
    one_leg_tau_coefficient = c_tau
    two_leg_tau_coefficient = 2*c_tau
    one_leg_phi_coefficient = c_phi
    two_leg_phi_coefficient = 2*c_phi

    ! C -> -C changes section orientation only.  It does not transform the
    ! Hamiltonian normal form or the signed toroidal observable.
    c_tau_section_reverse = c_tau
    c_phi_section_reverse = c_phi

    ! phi -> -phi reverses every signed toroidal rate, but not elapsed time.
    c_tau_toroidal_reverse = c_tau
    c_phi_toroidal_reverse = (-toroidal_rate_saddle)*c_tau
    absolute_phi_reversal_defect = &
        abs(c_phi_toroidal_reverse)+abs(c_phi)

    ! These diagnostics contain no square root, division, or logarithm.
    nondegeneracy_margin = &
        saddle_discriminant-saddle_rate_floor**2
    section_reversed_margin = nondegeneracy_margin
    epsilon_margin = epsilon-epsilon_floor
    section_reversed_epsilon_margin = epsilon_margin
    saddle_rate_floor_margin = saddle_rate_floor
    epsilon_floor_margin = epsilon_floor

    log_epsilon = log(epsilon)
    tau_finite_one = tau_one_epsilon+c_tau*log_epsilon
    tau_finite_two = tau_two_epsilon+2*c_tau*log_epsilon
    phi_finite_one = phi_one_epsilon+c_phi*log_epsilon
    phi_finite_two = phi_two_epsilon+2*c_phi*log_epsilon
    phi_finite_one_toroidal_reverse = &
        -phi_one_epsilon-c_phi*log_epsilon
    phi_finite_two_toroidal_reverse = &
        -phi_two_epsilon-2*c_phi*log_epsilon

    call suite_begin(proofs, 'homoclinic limit Fortsym contracts')
    call check_identity(proofs, proof_engine, 'quadratic H qq derivative', &
        h_qq-h_qq_input)
    call check_identity(proofs, proof_engine, 'quadratic H qp derivative', &
        h_qp-h_qp_input)
    call check_identity(proofs, proof_engine, 'quadratic H pp derivative', &
        h_pp-h_pp_input)
    call check_identity(proofs, proof_engine, 'Hamiltonian flow is traceless', &
        flow_qq+flow_pp)
    call check_identity(proofs, proof_engine, 'flow and Hessian determinants', &
        flow_determinant-hessian_determinant)
    call check_identity(proofs, proof_engine, 'saddle discriminant', &
        saddle_discriminant+flow_determinant)
    call check_identity(proofs, proof_engine, 'saddle rate squared', &
        saddle_rate**2-saddle_discriminant)
    call check_identity(proofs, proof_engine, 'positive-time coefficient', &
        c_tau*saddle_rate-1)
    call check_identity(proofs, proof_engine, 'signed toroidal coefficient', &
        c_phi*saddle_rate-toroidal_rate_saddle)
    call check_identity(proofs, proof_engine, 'two saddle time legs', &
        two_leg_tau_coefficient-2*one_leg_tau_coefficient)
    call check_identity(proofs, proof_engine, 'two saddle toroidal legs', &
        two_leg_phi_coefficient-2*one_leg_phi_coefficient)
    call check_identity(proofs, proof_engine, 'section reversal keeps C tau', &
        c_tau_section_reverse-c_tau)
    call check_identity(proofs, proof_engine, 'section reversal keeps C phi', &
        c_phi_section_reverse-c_phi)
    call check_identity(proofs, proof_engine, 'toroidal reversal keeps C tau', &
        c_tau_toroidal_reverse-c_tau)
    call check_identity(proofs, proof_engine, 'toroidal reversal flips C phi', &
        c_phi_toroidal_reverse+c_phi)
    call check_identity(proofs, proof_engine, &
        'absolute C phi violates signed reversal', &
        absolute_phi_reversal_defect-2*abs(c_phi))
    call check_identity(proofs, proof_engine, 'nondegeneracy margin', &
        nondegeneracy_margin- &
        (saddle_discriminant-saddle_rate_floor**2))
    call check_identity(proofs, proof_engine, &
        'section reversal keeps nondegeneracy margin', &
        section_reversed_margin-nondegeneracy_margin)
    call check_identity(proofs, proof_engine, 'epsilon margin', &
        epsilon_margin-(epsilon-epsilon_floor))
    call check_identity(proofs, proof_engine, &
        'section reversal keeps epsilon margin', &
        section_reversed_epsilon_margin-epsilon_margin)
    call check_identity(proofs, proof_engine, 'one-leg tau finite part', &
        tau_finite_one-(tau_one_epsilon+c_tau*log(epsilon)))
    call check_identity(proofs, proof_engine, 'two-leg tau finite part', &
        tau_finite_two-(tau_two_epsilon+2*c_tau*log(epsilon)))
    call check_identity(proofs, proof_engine, 'one-leg phi finite part', &
        phi_finite_one-(phi_one_epsilon+c_phi*log(epsilon)))
    call check_identity(proofs, proof_engine, 'two-leg phi finite part', &
        phi_finite_two-(phi_two_epsilon+2*c_phi*log(epsilon)))
    call check_identity(proofs, proof_engine, 'one-leg phi finite reversal', &
        phi_finite_one_toroidal_reverse+phi_finite_one)
    call check_identity(proofs, proof_engine, 'two-leg phi finite reversal', &
        phi_finite_two_toroidal_reverse+phi_finite_two)
    call suite_end(proofs)
    if (proofs%failed /= 0) error stop 'homoclinic limit proof failed'

    coefficient_roots = [hessian_determinant, saddle_discriminant, &
        saddle_rate, c_tau, c_phi, one_leg_tau_coefficient, &
        two_leg_tau_coefficient, one_leg_phi_coefficient, &
        two_leg_phi_coefficient, c_tau_section_reverse, &
        c_phi_section_reverse, c_tau_toroidal_reverse, &
        c_phi_toroidal_reverse, absolute_phi_reversal_defect]
    diagnostic_roots = [hessian_determinant, saddle_discriminant, &
        nondegeneracy_margin, section_reversed_margin, epsilon_margin, &
        section_reversed_epsilon_margin, saddle_rate_floor_margin, &
        epsilon_floor_margin]
    finite_part_roots = [tau_finite_one, tau_finite_two, phi_finite_one, &
        phi_finite_two, phi_finite_one_toroidal_reverse, &
        phi_finite_two_toroidal_reverse]
    call simplify_array(coefficient_roots)
    call simplify_array(diagnostic_roots)
    call simplify_array(finite_part_roots)

    call emit_kernel_file(trim(output_directory)// &
        '/neort_gc_homoclinic_coefficients_symbolic.f90', &
        'neort_gc_homoclinic_coefficients_symbolic', &
        'evaluate_neort_gc_homoclinic_coefficients', coefficient_roots, &
        [character(len=64) :: 'h_qq', 'h_qp', 'h_pp', &
        'toroidal_rate_saddle'], [character(len=64) :: &
        'hessian_determinant', 'saddle_discriminant', 'saddle_rate', &
        'c_tau', 'c_phi', 'one_leg_tau_coefficient', &
        'two_leg_tau_coefficient', 'one_leg_phi_coefficient', &
        'two_leg_phi_coefficient', 'c_tau_section_reverse', &
        'c_phi_section_reverse', 'c_tau_toroidal_reverse', &
        'c_phi_toroidal_reverse', 'absolute_phi_reversal_defect'], &
        'Require positive diagnostic margins before this kernel.')
    call emit_kernel_file(trim(output_directory)// &
        '/neort_gc_homoclinic_diagnostic_symbolic.f90', &
        'neort_gc_homoclinic_diagnostic_symbolic', &
        'evaluate_neort_gc_homoclinic_diagnostic', diagnostic_roots, &
        [character(len=64) :: 'h_qq', 'h_qp', 'h_pp', &
        'saddle_rate_floor', 'epsilon', 'epsilon_floor'], &
        [character(len=64) :: 'hessian_determinant', &
        'saddle_discriminant', 'nondegeneracy_margin', &
        'section_reversed_margin', 'epsilon_margin', &
        'section_reversed_epsilon_margin', 'saddle_rate_floor_margin', &
        'epsilon_floor_margin'], &
        'Division-free diagnostics; require every margin to be positive.')
    call emit_kernel_file(trim(output_directory)// &
        '/neort_gc_homoclinic_finite_part_symbolic.f90', &
        'neort_gc_homoclinic_finite_part_symbolic', &
        'evaluate_neort_gc_homoclinic_finite_part', finite_part_roots, &
        [character(len=64) :: 'tau_one_epsilon', 'tau_two_epsilon', &
        'phi_one_epsilon', 'phi_two_epsilon', 'epsilon', 'c_tau', 'c_phi'], &
        [character(len=64) :: 'tau_finite_one', 'tau_finite_two', &
        'phi_finite_one', 'phi_finite_two', &
        'phi_finite_one_toroidal_reverse', &
        'phi_finite_two_toroidal_reverse'], &
        'Require positive epsilon diagnostic before evaluating log(epsilon).')

contains

    subroutine simplify_array(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        integer :: i

        do i = 1, size(expressions)
            simplified = simplify_engine%simplify(expressions(i))
            if (.not. simplified%ok) then
                error stop 'homoclinic limit simplification failed'
            end if
            expressions(i) = simplified%value
        end do
    end subroutine simplify_array

    subroutine emit_kernel_file(path, module_name, procedure_name, roots, &
            argument_names, output_names, contract_comment)
        character(*), intent(in) :: path, module_name, procedure_name
        character(*), intent(in) :: argument_names(:), output_names(:)
        character(*), intent(in) :: contract_comment
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t) :: spec
        type(str_t) :: emitted
        character(:), allocatable :: emitted_text
        integer :: i, ios, unit
        logical :: ok

        if (size(roots) /= size(output_names)) then
            error stop 'homoclinic limit output arity mismatch'
        end if
        spec%name = str(procedure_name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%generator = str(GENERATOR_PATH)
        spec%generator_revision = str(FORTSYM_REVISION)
        spec%regenerate_command = str(REGENERATE_COMMAND)
        allocate (spec%args(size(argument_names)), &
            spec%outputs(size(output_names)))
        do i = 1, size(argument_names)
            spec%args(i) = str(trim(argument_names(i)))
        end do
        do i = 1, size(output_names)
            spec%outputs(i) = str(trim(output_names(i)))
        end do
        emitted = emit_kernel(roots, spec, ok)
        if (.not. ok) error stop 'homoclinic limit kernel emission failed'
        emitted_text = chars(emitted)
        open (newunit=unit, file=trim(path), status='replace', action='write', &
            iostat=ios)
        if (ios /= 0) error stop 'cannot open homoclinic limit output'
        write (unit, '(a)') &
            '! Fortsym-generated homoclinic limit; do not hand-edit.'
        write (unit, '(a)') '! '//trim(contract_comment)
        write (unit, '(a)', advance='no') emitted_text
        close (unit)
        write (output_unit, '(a)') 'wrote '//trim(path)
    end subroutine emit_kernel_file

end program gen_gc_homoclinic_limit
