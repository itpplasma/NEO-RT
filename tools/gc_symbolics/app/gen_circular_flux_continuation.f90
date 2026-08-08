program gen_circular_flux_continuation
    !! Derive the smooth exterior continuation of the circular golden record.
    !!
    !! The circular fixture prescribes the sheared safety factor
    !!
    !!     q(r) = q_axis + delta_q*(r/edge_radius)**2,
    !!
    !! and the equilibrium relation d psi / d r = F*r/q(r).  Fortsym derives
    !! the logarithmic primitive, its first derivative, and the zero-shear
    !! limiting expression.  The sampled Python fixture data are emitted from
    !! those same expression trees; the fixture never carries a second copy of
    !! the continuation algebra.
    use, intrinsic :: iso_fortran_env, only: output_unit, real64
    use fortsym_arena, only: arena_t
    use fortsym_check, only: check_identity, suite_begin, suite_end, suite_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine_symengine, only: make_symengine_engine, symengine_engine_t
    use fortsym_eval, only: binding_t, eval_expr
    use fortsym_expr, only: expr_t, log, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    use fortsym_subs, only: subs
    implicit none

    integer, parameter :: dp = real64
    integer, parameter :: N_INTERIOR = 1024
    integer, parameter :: N_EXTERIOR = 1024
    real(dp), parameter :: R0_VALUE = 1.60_dp
    real(dp), parameter :: EDGE_RADIUS_VALUE = 0.50_dp
    real(dp), parameter :: B0_VALUE = 2.0_dp
    real(dp), parameter :: Q_AXIS_VALUE = 1.5_dp
    real(dp), parameter :: Q_EDGE_VALUE = 4.0_dp
    real(dp), parameter :: RADIUS_EXTENT_FACTOR = 2.4_dp
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
    type(expr_t) :: radius, edge_radius, psi_edge, toroidal_flux
    type(expr_t) :: q_axis, delta_q, q_edge, safety_factor
    type(expr_t) :: psi_continuation, dpsi_continuation
    type(expr_t) :: psi_limit, dpsi_limit, psi_profile, dpsi_profile
    type(expr_t) :: roots(3), limit_roots(3), edge_value_residual
    type(expr_t) :: edge_slope_residual
    type(engine_result_t) :: simplified
    character(2048) :: kernel_path, limit_kernel_path, python_path
    integer :: kernel_length, limit_kernel_length, python_length, argument_status
    integer :: unit, ios, i
    real(dp) :: r_value, psi_value, q_value
    logical :: defined
    type(binding_t) :: bindings

    call get_command_argument(1, kernel_path, length=kernel_length, &
        status=argument_status)
    if (argument_status /= 0 .or. kernel_length == 0) error stop &
        'usage: gen_circular_flux_continuation KERNEL_PATH LIMIT_KERNEL_PATH PYTHON_PATH'
    call get_command_argument(2, limit_kernel_path, length=limit_kernel_length, &
        status=argument_status)
    if (argument_status /= 0 .or. limit_kernel_length == 0) error stop &
        'usage: gen_circular_flux_continuation KERNEL_PATH LIMIT_KERNEL_PATH PYTHON_PATH'
    call get_command_argument(3, python_path, length=python_length, &
        status=argument_status)
    if (argument_status /= 0 .or. python_length == 0) error stop &
        'usage: gen_circular_flux_continuation KERNEL_PATH LIMIT_KERNEL_PATH PYTHON_PATH'
    kernel_path = kernel_path(:kernel_length)
    limit_kernel_path = limit_kernel_path(:limit_kernel_length)
    python_path = python_path(:python_length)

    call arena%init()
    proof_engine = make_symengine_engine(arena)
    simplify_engine = make_native_engine(arena)

    radius = sym(arena, 'radius')
    edge_radius = sym(arena, 'edge_radius')
    psi_edge = sym(arena, 'psi_edge')
    toroidal_flux = sym(arena, 'toroidal_flux')
    q_axis = sym(arena, 'q_axis')
    delta_q = sym(arena, 'delta_q')
    q_edge = q_axis + delta_q
    safety_factor = q_axis + delta_q*(radius/edge_radius)**2

    ! The additive constant is the exact LCFS value.  Thus the continuation
    ! is value-continuous at radius=edge_radius by construction.
    psi_continuation = psi_edge + toroidal_flux*edge_radius**2/(2*delta_q)* &
        log(safety_factor/q_edge)
    dpsi_continuation = diff(psi_continuation, radius)
    psi_profile = toroidal_flux*edge_radius**2/(2*delta_q)* &
        log(safety_factor/q_axis)
    dpsi_profile = diff(psi_profile, radius)

    ! The delta_q=0 branch is derived as the finite zero-shear limit of the
    ! same radial flux law.  It is emitted separately so no runtime division by
    ! delta_q is ever used at the limiting parameter value.
    psi_limit = psi_edge + toroidal_flux*(radius**2-edge_radius**2)/ &
        (2*q_axis)
    dpsi_limit = diff(psi_limit, radius)

    call suite_begin(proofs, 'circular flux continuation')
    call check_identity(proofs, proof_engine, 'continuation derivative', &
        dpsi_continuation - toroidal_flux*radius/safety_factor)
    call check_identity(proofs, proof_engine, 'limit derivative', &
        dpsi_limit - toroidal_flux*radius/q_axis)
    call check_identity(proofs, proof_engine, 'shared interior flux law', &
        dpsi_profile - toroidal_flux*radius/safety_factor)
    edge_value_residual = subs(psi_continuation, radius, edge_radius) - psi_edge
    call check_identity(proofs, proof_engine, 'edge value', edge_value_residual)
    edge_slope_residual = subs(dpsi_profile, radius, edge_radius) - &
        subs(dpsi_continuation, radius, edge_radius)
    call check_identity(proofs, proof_engine, 'interior-continuation edge slope', &
        edge_slope_residual)
    call suite_end(proofs)
    if (proofs%failed /= 0) error stop 'circular continuation proof failed'

    roots = [psi_continuation, dpsi_continuation, safety_factor]
    limit_roots = [psi_limit, dpsi_limit, q_axis]
    call simplify_array(roots)
    call simplify_array(limit_roots)
    call emit_kernel_file(trim(kernel_path), &
        'neort_circular_flux_continuation_symbolic', &
        'evaluate_neort_circular_flux_continuation', roots, &
        [character(len=64) :: 'radius', 'edge_radius', 'psi_edge', &
        'toroidal_flux', 'q_axis', 'delta_q'], &
        [character(len=64) :: 'psi', 'dpsi_dr', 'q_safety'])
    call emit_kernel_file(trim(limit_kernel_path), &
        'neort_circular_flux_continuation_limit_symbolic', &
        'evaluate_neort_circular_flux_continuation_limit', limit_roots, &
        [character(len=64) :: 'radius', 'edge_radius', 'psi_edge', &
        'toroidal_flux', 'q_axis'], &
        [character(len=64) :: 'psi', 'dpsi_dr', 'q_safety'])
    call emit_python_table(trim(python_path), psi_profile, safety_factor, &
        q_axis, delta_q, edge_radius, toroidal_flux)

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
        write(unit, '(a)') emitted_text
        close(unit)
        write(output_unit, '(a)') 'wrote '//trim(path)
    end subroutine emit_kernel_file

    subroutine emit_python_table(path, psi_expression, q_expression, q0, dq, &
            a, flux)
        character(*), intent(in) :: path
        type(expr_t), intent(in) :: psi_expression, q_expression, q0, dq, a, flux
        integer :: count, k
        real(dp) :: psi_edge_value

        count = N_INTERIOR + N_EXTERIOR + 1
        allocate(bindings%names(6), bindings%values(6))
        bindings%names = [str('radius'), str('edge_radius'), str('psi_edge'), &
            str('toroidal_flux'), str('q_axis'), str('delta_q')]
        bindings%values = [0.0_dp, EDGE_RADIUS_VALUE, 0.0_dp, &
            R0_VALUE*B0_VALUE, Q_AXIS_VALUE, Q_EDGE_VALUE-Q_AXIS_VALUE]
        bindings%n = 6
        bindings%values(2) = EDGE_RADIUS_VALUE
        psi_edge_value = evaluate(psi_expression, EDGE_RADIUS_VALUE)
        bindings%values(3) = psi_edge_value

        open(newunit=unit, file=trim(path), status='replace', action='write', &
            iostat=ios)
        if (ios /= 0) error stop 'cannot open continuation table output'
        write(unit, '(a)') '# Generated by Fortsym; do not edit.'
        write(unit, '(a)') '# Generator: '//GENERATOR_PATH
        write(unit, '(a)') '# Fortsym revision: '//FORTSYM_REVISION
        write(unit, '(a)') 'import numpy as np'
        write(unit, '(a)') ''
        write(unit, '(a,es24.16e3)') 'R0 = ', R0_VALUE
        write(unit, '(a,es24.16e3)') 'A = ', EDGE_RADIUS_VALUE
        write(unit, '(a,es24.16e3)') 'B0 = ', B0_VALUE
        write(unit, '(a,es24.16e3)') 'Q0 = ', Q_AXIS_VALUE
        write(unit, '(a,es24.16e3)') 'QA = ', Q_EDGE_VALUE
        write(unit, '(a,i0)') 'EDGE_INDEX = ', N_INTERIOR
        write(unit, '(a,es24.16e3)') 'TOROIDAL_FLUX = ', &
            R0_VALUE*B0_VALUE
        write(unit, '(a,es24.16e3)') 'PSI_EDGE = ', psi_edge_value
        write(unit, '(a)') ''
        write(unit, '(a)') 'RADIUS = np.array(['
        do k = 0, count-1
            r_value = radial_value(k)
            call write_python_value(unit, r_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'PSI = np.array(['
        do k = 0, count-1
            r_value = radial_value(k)
            psi_value = evaluate(psi_expression, r_value)
            call write_python_value(unit, psi_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') 'Q = np.array(['
        do k = 0, count-1
            r_value = radial_value(k)
            q_value = evaluate(q_expression, r_value)
            call write_python_value(unit, q_value)
        end do
        write(unit, '(a)') '], dtype=np.float64)'
        write(unit, '(a)') ''
        write(unit, '(a)') 'assert RADIUS[EDGE_INDEX] == A'
        write(unit, '(a)') 'assert np.all(np.diff(PSI) > 0.0)'
        close(unit)
        write(output_unit, '(a)') 'wrote '//trim(path)
        associate(unused_q0 => q0, unused_dq => dq, unused_a => a, &
                unused_flux => flux)
        end associate
    end subroutine emit_python_table

    real(dp) function radial_value(index)
        integer, intent(in) :: index

        if (index <= N_INTERIOR) then
            radial_value = EDGE_RADIUS_VALUE*real(index, dp)/real(N_INTERIOR, dp)
        else
            radial_value = EDGE_RADIUS_VALUE + &
                (RADIUS_EXTENT_FACTOR-1.0_dp)*EDGE_RADIUS_VALUE* &
                real(index-N_INTERIOR, dp)/real(N_EXTERIOR, dp)
        end if
    end function radial_value

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

end program gen_circular_flux_continuation
