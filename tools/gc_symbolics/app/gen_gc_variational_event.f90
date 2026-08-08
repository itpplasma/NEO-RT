program gen_gc_variational_event
    !! Derive the first variational return-map event kernel.
    !!
    !! The regular kernel is used only after the generated diagnostic margin
    !! is positive.  The diagnostic kernel exposes D and all numerator
    !! quantities without dividing by D, so a rejected event is fail-closed.
    !! Every derivative and orientation identity below is an expression tree
    !! owned by Fortsym; the consumer only applies the returned certificate.
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_arena, only: arena_t
    use fortsym_check, only: check_identity, suite_begin, suite_end, suite_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine_symengine, only: make_symengine_engine, &
        symengine_engine_t
    use fortsym_expr, only: abs, expr_t, operator(+), operator(-), &
        operator(*), operator(/), sym
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str, str_t
    implicit none

    character(*), parameter :: FORTSYM_REVISION = &
        'fortsym@545788453a204d58705f735b519c3863c2f734c8'
    character(*), parameter :: GENERATOR_PATH = &
        'tools/gc_symbolics/app/gen_gc_variational_event.f90'
    character(*), parameter :: REGENERATE_COMMAND = &
        'cd tools/gc_symbolics && fo exec gen_gc_variational_event '// &
        '../../src/generated'

    type(arena_t), target :: arena
    type(symengine_engine_t) :: proof_engine
    type(native_engine_t) :: simplify_engine
    type(suite_t) :: proofs
    type(engine_result_t) :: simplified
    type(expr_t) :: f_x, f_y, c_x, c_y, c_t, c_lambda
    type(expr_t) :: s_x, s_y, y_x, y_y, y_t, y_lambda_explicit
    type(expr_t) :: d, event_numerator, tau_lambda
    type(expr_t) :: y_explicit, y_transport, y_lambda
    type(expr_t) :: c_x_reverse, c_y_reverse, c_t_reverse
    type(expr_t) :: c_lambda_reverse, d_reverse
    type(expr_t) :: event_numerator_reverse, tau_lambda_reverse
    type(expr_t) :: y_lambda_reverse, abs_d, d_squared
    type(expr_t) :: transversality_floor, transversality_margin
    type(expr_t) :: transversality_margin_reverse
    type(expr_t) :: tau_orientation_residual
    type(expr_t) :: y_orientation_residual
    type(expr_t) :: regular_roots(14), diagnostic_roots(10)
    character(2048) :: output_directory
    integer :: argument_status, output_length

    call get_command_argument(1, output_directory, length=output_length, &
        status=argument_status)
    if (argument_status /= 0 .or. output_length == 0) then
        write (output_unit, '(a)') &
            'usage: gen_gc_variational_event OUTPUT_DIRECTORY'
        error stop 2
    end if
    output_directory = output_directory(:output_length)

    call arena%init()
    proof_engine = make_symengine_engine(arena)
    simplify_engine = make_native_engine(arena)

    f_x = sym(arena, 'f_x')
    f_y = sym(arena, 'f_y')
    c_x = sym(arena, 'c_x')
    c_y = sym(arena, 'c_y')
    c_t = sym(arena, 'c_t')
    c_lambda = sym(arena, 'c_lambda')
    s_x = sym(arena, 's_x')
    s_y = sym(arena, 's_y')
    y_x = sym(arena, 'y_x')
    y_y = sym(arena, 'y_y')
    y_t = sym(arena, 'y_t')
    y_lambda_explicit = sym(arena, 'y_lambda_explicit')
    transversality_floor = sym(arena, 'transversality_floor')

    ! D = C_z f + C_t, and N is the fixed-time event sensitivity.
    d = c_x*f_x + c_y*f_y + c_t
    event_numerator = c_x*s_x + c_y*s_y + c_lambda
    tau_lambda = -event_numerator/d

    y_explicit = y_lambda_explicit
    y_transport = y_x*f_x + y_y*f_y + y_t
    y_lambda = y_x*s_x + y_y*s_y + y_explicit + y_transport*tau_lambda

    ! Reverse only the section orientation: C -> -C.  The observable and
    ! trajectory are unchanged, so the event-time and Y derivatives survive.
    c_x_reverse = -c_x
    c_y_reverse = -c_y
    c_t_reverse = -c_t
    c_lambda_reverse = -c_lambda
    d_reverse = c_x_reverse*f_x + c_y_reverse*f_y + c_t_reverse
    event_numerator_reverse = c_x_reverse*s_x + c_y_reverse*s_y + &
        c_lambda_reverse
    tau_lambda_reverse = -event_numerator_reverse/d_reverse
    y_lambda_reverse = y_x*s_x + y_y*s_y + y_explicit + &
        y_transport*tau_lambda_reverse
    tau_orientation_residual = tau_lambda_reverse - tau_lambda
    y_orientation_residual = y_lambda_reverse - y_lambda
    abs_d = abs(d)
    d_squared = d*d
    transversality_margin = d_squared- &
        transversality_floor*transversality_floor
    transversality_margin_reverse = &
        d_reverse*d_reverse-transversality_floor*transversality_floor

    call suite_begin(proofs, 'variational event Fortsym contracts')
    call check_identity(proofs, proof_engine, 'event transversality D', &
        d-(c_x*f_x+c_y*f_y+c_t))
    call check_identity(proofs, proof_engine, 'event sensitivity numerator', &
        event_numerator-(c_x*s_x+c_y*s_y+c_lambda))
    call check_identity(proofs, proof_engine, 'implicit event-time derivative', &
        tau_lambda*d+event_numerator)
    call check_identity(proofs, proof_engine, 'explicit terminal derivative', &
        y_explicit-y_lambda_explicit)
    call check_identity(proofs, proof_engine, 'terminal transport derivative', &
        y_transport-(y_x*f_x+y_y*f_y+y_t))
    call check_identity(proofs, proof_engine, 'event-corrected terminal derivative', &
        y_lambda-(y_x*s_x+y_y*s_y+y_lambda_explicit+ &
        (y_x*f_x+y_y*f_y+y_t)*tau_lambda))
    call check_identity(proofs, proof_engine, 'reversed D changes sign', &
        d_reverse+d)
    call check_identity(proofs, proof_engine, 'reversed numerator changes sign', &
        event_numerator_reverse+event_numerator)
    call check_identity(proofs, proof_engine, 'section reversal preserves tau', &
        tau_orientation_residual)
    call check_identity(proofs, proof_engine, 'section reversal preserves Y', &
        y_orientation_residual)
    call check_identity(proofs, proof_engine, 'absolute transversality output', &
        abs_d-abs(d))
    call check_identity(proofs, proof_engine, 'squared transversality output', &
        d_squared-d*d)
    call check_identity(proofs, proof_engine, 'transversality margin', &
        transversality_margin-(d*d- &
        transversality_floor*transversality_floor))
    call check_identity(proofs, proof_engine, &
        'section reversal preserves transversality margin', &
        transversality_margin_reverse-transversality_margin)
    call suite_end(proofs)
    if (proofs%failed /= 0) error stop 'variational event proof failed'

    regular_roots = [d, event_numerator, tau_lambda, y_explicit, y_transport, &
        y_lambda, d_reverse, event_numerator_reverse, tau_lambda_reverse, &
        y_lambda_reverse, tau_orientation_residual, y_orientation_residual, &
        abs_d, d_squared]
    diagnostic_roots = [d, event_numerator, y_explicit, y_transport, &
        d_reverse, event_numerator_reverse, abs_d, d_squared, &
        transversality_margin, transversality_margin_reverse]
    call simplify_array(regular_roots)
    call simplify_array(diagnostic_roots)

    call emit_kernel_file(trim(output_directory)// &
        '/neort_gc_variational_event_symbolic.f90', &
        'neort_gc_variational_event_symbolic', &
        'evaluate_neort_gc_variational_event', regular_roots, .false.)
    call emit_kernel_file(trim(output_directory)// &
        '/neort_gc_variational_event_diagnostic_symbolic.f90', &
        'neort_gc_variational_event_diagnostic_symbolic', &
        'evaluate_neort_gc_variational_event_diagnostic', diagnostic_roots, &
        .true.)

contains

    subroutine simplify_array(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        integer :: i

        do i = 1, size(expressions)
            simplified = simplify_engine%simplify(expressions(i))
            if (.not. simplified%ok) then
                error stop 'variational event simplification failed'
            end if
            expressions(i) = simplified%value
        end do
    end subroutine simplify_array

    subroutine emit_kernel_file(path, module_name, procedure_name, roots, &
            diagnostic)
        character(*), intent(in) :: path, module_name, procedure_name
        type(expr_t), intent(in) :: roots(:)
        logical, intent(in) :: diagnostic
        character(len=64), parameter :: regular_args(12) = [character(len=64) :: &
            'f_x', 'f_y', 'c_x', 'c_y', 'c_t', 'c_lambda', 's_x', 's_y', &
            'y_x', 'y_y', 'y_t', 'y_lambda_explicit']
        character(len=64), parameter :: diagnostic_args(13) = &
            [character(len=64) :: 'f_x', 'f_y', 'c_x', 'c_y', 'c_t', &
            'c_lambda', 's_x', 's_y', 'y_x', 'y_y', 'y_t', &
            'y_lambda_explicit', 'transversality_floor']
        character(len=64) :: outputs(14)
        type(kernel_spec_t) :: spec
        type(str_t) :: emitted
        character(:), allocatable :: emitted_text
        integer :: i, ios, unit
        logical :: ok

        if (diagnostic) then
            outputs(1:10) = [character(len=64) :: 'd', 'event_numerator', &
                'y_explicit', 'y_transport', 'd_reverse', &
                'event_numerator_reverse', 'abs_d', 'd_squared', &
                'transversality_margin', 'transversality_margin_reverse']
        else
            outputs = [character(len=64) :: 'd', 'event_numerator', &
                'tau_lambda', 'y_explicit', 'y_transport', 'y_lambda', &
                'd_reverse', 'event_numerator_reverse', 'tau_lambda_reverse', &
                'y_lambda_reverse', 'tau_orientation_residual', &
                'y_orientation_residual', 'abs_d', 'd_squared']
        end if
        spec%name = str(procedure_name)
        spec%module_name = str(module_name)
        spec%mode = KERNEL_SUBROUTINE
        spec%pure_procedure = .true.
        spec%generator = str(GENERATOR_PATH)
        spec%generator_revision = str(FORTSYM_REVISION)
        spec%regenerate_command = str(REGENERATE_COMMAND)
        if (diagnostic) then
            allocate (spec%args(size(diagnostic_args)))
            do i = 1, size(diagnostic_args)
                spec%args(i) = str(diagnostic_args(i))
            end do
        else
            allocate (spec%args(size(regular_args)))
            do i = 1, size(regular_args)
                spec%args(i) = str(regular_args(i))
            end do
        end if
        allocate (spec%outputs(size(roots)))
        do i = 1, size(roots)
            spec%outputs(i) = str(outputs(i))
        end do
        emitted = emit_kernel(roots, spec, ok)
        if (.not. ok) error stop 'variational event kernel emission failed'
        emitted_text = chars(emitted)
        open (newunit=unit, file=trim(path), status='replace', action='write', &
            iostat=ios)
        if (ios /= 0) error stop 'cannot open variational event output'
        write (unit, '(a)') &
            '! Fortsym-generated; do not hand-edit this kernel.'
        if (diagnostic) then
            write (unit, '(a)') &
                '! This diagnostic kernel never divides by D.'
        else
            write (unit, '(a)') &
                '! Require diagnostic transversality_margin > 0 first.'
        end if
        write (unit, '(a)', advance='no') emitted_text
        close (unit)
        write (output_unit, '(a)') 'wrote '//trim(path)
    end subroutine emit_kernel_file

end program gen_gc_variational_event
