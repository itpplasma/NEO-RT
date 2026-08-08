program gen_gc_class_coordinate_maps
    !! Derive the bounded class-coordinate quadrature maps used by POTATO's
    !! Poincare cut.  This slice is separate from physical homoclinic time
    !! limits and from topology certificates.
    !!
    !! Boundary types 1 and 2 are regular and turning endpoints.  Types 3 and
    !! 4 have the same map shape but different logarithmic truncations: the
    !! regular separatrix uses one half and the singular X-point one quarter.
    !! The generated interface is deliberately a flat value/derivative table;
    !! topology and class assembly remain runtime responsibilities.
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use fortsym_arena, only: arena_t
    use fortsym_check, only: check_identity, suite_begin, suite_end, suite_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine_symengine, only: make_symengine_engine, &
        symengine_engine_t
    use fortsym_expr, only: cos, expr_t, log, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**), pi_expr, rat, sin, sym, tanh
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    use fortsym_trigrewrite, only: trig_expand
    implicit none

    character(*), parameter :: FORTSYM_REVISION = &
        'fortsym@545788453a204d58705f735b519c3863c2f734c8'
    character(*), parameter :: REGENERATE_COMMAND = &
        'cd tools/gc_symbolics && fo exec gen_gc_class_coordinate_maps '// &
        'OUTPUT_PATH'

    character(2048) :: output_path
    integer :: output_length, argument_status
    type(arena_t), target :: arena
    type(symengine_engine_t) :: proof_engine
    type(native_engine_t) :: simplify_engine
    type(suite_t) :: proofs
    type(expr_t) :: x, relmargin, widthclass, log_ratio
    type(expr_t) :: xi(4, 4), dxi(4, 4), xbeg(4, 4), xend(4, 4)
    type(expr_t) :: section(4, 4), section_derivative(4, 4)
    type(expr_t), allocatable :: roots(:)
    type(kernel_spec_t) :: spec
    integer :: left_type, right_type, cursor

    call get_command_argument(1, output_path, length=output_length, &
        status=argument_status)
    if (argument_status /= 0 .or. output_length == 0) then
        write (output_unit, '(a)') &
            'usage: gen_gc_class_coordinate_maps OUTPUT_PATH'
        error stop 2
    end if
    output_path = output_path(:output_length)

    call arena%init()
    proof_engine = make_symengine_engine(arena)
    simplify_engine = make_native_engine(arena)
    x = sym(arena, 'class_coordinate')
    relmargin = sym(arena, 'relative_margin')
    widthclass = sym(arena, 'relative_class_width')
    log_ratio = log(relmargin/widthclass)

    do left_type = 1, 4
        do right_type = 1, 4
            xi(left_type, right_type) = map_expression(arena, left_type, &
                right_type, x)
            dxi(left_type, right_type) = diff(xi(left_type, right_type), x)
            xbeg(left_type, right_type) = left_bound(arena, left_type, &
                relmargin, log_ratio)
            xend(left_type, right_type) = right_bound(arena, left_type, &
                right_type, relmargin, log_ratio)
            section(left_type, right_type) = sym(arena, 'section_left') + &
                (sym(arena, 'section_right')-sym(arena, 'section_left')) &
                *xi(left_type, right_type)
            section_derivative(left_type, right_type) = &
                diff(section(left_type, right_type), x)
        end do
    end do

    call prove_contracts(arena, proof_engine, proofs, x, relmargin, &
        log_ratio, xi, dxi, xbeg, xend, section, &
        section_derivative)
    call simplify_all(simplify_engine, xi)
    call simplify_all(simplify_engine, dxi)
    call simplify_all(simplify_engine, xbeg)
    call simplify_all(simplify_engine, xend)
    call simplify_all(simplify_engine, section)
    call simplify_all(simplify_engine, section_derivative)

    allocate (roots(4*4*6))
    cursor = 0
    call append_table(roots, cursor, xi)
    call append_table(roots, cursor, dxi)
    call append_table(roots, cursor, xbeg)
    call append_table(roots, cursor, xend)
    call append_table(roots, cursor, section)
    call append_table(roots, cursor, section_derivative)
    call emit_output(trim(output_path), roots)

contains

    function map_expression(arena, left_type, right_type, argument) &
            result(value)
        type(arena_t), target, intent(inout) :: arena
        integer, intent(in) :: left_type, right_type
        type(expr_t), intent(in) :: argument
        type(expr_t) :: value
        type(expr_t) :: half

        half = rat(arena, 1_int64, 2_int64)
        select case (left_type)
        case (1)
            select case (right_type)
            case (1); value = argument
            case (2); value = 1- (1-argument)**2
            case default; value = tanh(argument)
            end select
        case (2)
            select case (right_type)
            case (1); value = argument**2
            case (2); value = half*(1-cos(pi_expr(arena)*argument))
            case default; value = tanh(argument)**2
            end select
        case (3, 4)
            select case (right_type)
            case (1); value = 1+tanh(argument)
            case (2); value = 1-tanh(argument)**2
            case default; value = half*(1+tanh(argument))
            end select
        end select
    end function map_expression

    function left_bound(arena, boundary_type, margin, log_ratio) result(value)
        type(arena_t), target, intent(inout) :: arena
        integer, intent(in) :: boundary_type
        type(expr_t), intent(in) :: margin, log_ratio
        type(expr_t) :: value

        select case (boundary_type)
        case (1); value = num(arena, 0)
        case (2); value = margin
        case (3); value = rat(arena, 1_int64, 2_int64)*log_ratio
        case (4); value = rat(arena, 1_int64, 4_int64)*log_ratio
        end select
    end function left_bound

    function right_bound(arena, left_type, boundary_type, margin, log_ratio) &
            result(value)
        type(arena_t), target, intent(inout) :: arena
        integer, intent(in) :: left_type, boundary_type
        type(expr_t), intent(in) :: margin, log_ratio
        type(expr_t) :: value

        if (boundary_type >= 3) then
            select case (boundary_type)
            case (3); value = -rat(arena, 1_int64, 2_int64)*log_ratio
            case (4); value = -rat(arena, 1_int64, 4_int64)*log_ratio
            end select
        else if (left_type >= 3) then
            select case (boundary_type)
            case (1); value = num(arena, 0)
            case (2); value = -margin
            end select
        else
            select case (boundary_type)
            case (1); value = num(arena, 1)
            case (2); value = 1-margin
            end select
        end if
    end function right_bound

    subroutine prove_contracts(arena, engine, checks, x, margin, &
            log_ratio, xi, dxi, xbeg, xend, section, section_derivative)
        type(arena_t), target, intent(inout) :: arena
        type(symengine_engine_t) :: engine
        type(suite_t), intent(out) :: checks
        type(expr_t), intent(in) :: x, margin, log_ratio
        type(expr_t), intent(in) :: xi(4, 4), dxi(4, 4), xbeg(4, 4), xend(4, 4)
        type(expr_t), intent(in) :: section(4, 4), section_derivative(4, 4)
        integer :: left_type, right_type
        type(expr_t) :: reverse_x, reversed_map, reversed_section

        call suite_begin(checks, 'POTATO class-coordinate Fortsym contracts')
        call check_bound_contracts(checks, engine, arena, margin, &
            log_ratio, xbeg, xend)
        do left_type = 1, 4
            do right_type = 1, 4
                call check_identity(checks, engine, 'map derivative', &
                    dxi(left_type, right_type)-diff(xi(left_type, right_type), x))
                reverse_x = reverse_argument(arena, left_type, right_type, x)
                reversed_map = map_expression(arena, right_type, left_type, &
                    reverse_x)
                reversed_section = sym(arena, 'section_right') + &
                    (sym(arena, 'section_left')-sym(arena, 'section_right')) &
                    *reversed_map
                call check_orientation(checks, engine, left_type, right_type, &
                    xi(left_type, right_type)+reversed_map-1, &
                    'section orientation map')
                call check_orientation(checks, engine, left_type, right_type, &
                    section(left_type, right_type)-reversed_section, &
                    'section orientation')
                call check_orientation(checks, engine, left_type, right_type, &
                    section_derivative(left_type, right_type)- &
                    diff(reversed_section, x), 'section orientation derivative')
            end do
        end do
        if (checks%failed /= 0) error stop 'class-coordinate proof failed'
        call suite_end(checks)
    end subroutine prove_contracts

    subroutine check_orientation(checks, engine, left_type, right_type, &
            residual, label)
        type(suite_t), intent(inout) :: checks
        type(symengine_engine_t) :: engine
        integer, intent(in) :: left_type, right_type
        type(expr_t), intent(in) :: residual
        character(*), intent(in) :: label
        type(expr_t) :: expanded, rewritten
        type(engine_result_t) :: expansion
        logical :: rewrite_ok
        character(:), allocatable :: rewrite_message

        if (left_type == 2 .and. right_type == 2) then
            expansion = simplify_engine%expand(residual)
            if (.not. expansion%ok) error stop 'class-map trig expansion failed'
            expanded = expansion%value
            call trig_expand(expanded, rewritten, rewrite_ok, rewrite_message)
            if (.not. rewrite_ok) error stop 'class-map trig rewrite failed'
            call check_identity(checks, engine, label, rewritten)
        else
            call check_identity(checks, engine, label, residual)
        end if
    end subroutine check_orientation

    subroutine check_bound_contracts(checks, engine, arena, margin, &
            log_ratio, xbeg, xend)
        type(suite_t), intent(inout) :: checks
        type(symengine_engine_t) :: engine
        type(arena_t), target, intent(inout) :: arena
        type(expr_t), intent(in) :: margin, log_ratio
        type(expr_t), intent(in) :: xbeg(4, 4), xend(4, 4)
        integer :: left_type, right_type
        type(expr_t) :: left_expected, right_expected

        do left_type = 1, 4
            do right_type = 1, 4
                left_expected = left_bound(arena, left_type, margin, log_ratio)
                right_expected = right_bound(arena, left_type, right_type, &
                    margin, log_ratio)
                call check_identity(checks, engine, 'left bound', &
                    xbeg(left_type, right_type)-left_expected)
                call check_identity(checks, engine, 'right bound', &
                    xend(left_type, right_type)-right_expected)
            end do
        end do
    end subroutine check_bound_contracts

    function reverse_argument(arena, left_type, right_type, argument) &
            result(reversed)
        type(arena_t), target, intent(inout) :: arena
        integer, intent(in) :: left_type, right_type
        type(expr_t), intent(in) :: argument
        type(expr_t) :: reversed

        if (left_type <= 2 .and. right_type <= 2) then
            reversed = 1-argument
        else
            reversed = -argument
        end if
    end function reverse_argument

    subroutine simplify_all(engine, expressions)
        type(native_engine_t) :: engine
        type(expr_t), intent(inout) :: expressions(:, :)
        type(engine_result_t) :: result
        integer :: i, j

        do i = 1, size(expressions, 1)
            do j = 1, size(expressions, 2)
                result = engine%simplify(expressions(i, j))
                if (.not. result%ok) error stop 'class-map simplification failed'
                expressions(i, j) = result%value
            end do
        end do
    end subroutine simplify_all

    subroutine append_table(roots, cursor, table)
        type(expr_t), intent(inout) :: roots(:)
        integer, intent(inout) :: cursor
        type(expr_t), intent(in) :: table(:, :)
        integer :: i, j

        do i = 1, size(table, 1)
            do j = 1, size(table, 2)
                cursor = cursor+1
                roots(cursor) = table(i, j)
            end do
        end do
    end subroutine append_table

    subroutine emit_output(path, roots)
        character(*), intent(in) :: path
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t) :: kernel
        integer :: unit, ios, i, j, group, cursor
        character(32) :: name

        kernel%name = str( &
            'evaluate_neort_gc_class_coordinate_quadrature_maps')
        kernel%module_name = str( &
            'neort_gc_class_coordinate_quadrature_maps_generated_mod')
        kernel%mode = KERNEL_SUBROUTINE
        kernel%pure_procedure = .true.
        kernel%generator = str( &
            'tools/gc_symbolics/app/gen_gc_class_coordinate_maps.f90')
        kernel%generator_revision = str(FORTSYM_REVISION)
        kernel%regenerate_command = str(REGENERATE_COMMAND)
        allocate (kernel%args(5), kernel%outputs(size(roots)))
        kernel%args = [str('class_coordinate'), str('relative_margin'), &
            str('relative_class_width'), str('section_left'), &
            str('section_right')]
        cursor = 0
        do group = 1, 6
            do i = 1, 4
                do j = 1, 4
                    cursor = cursor+1
                    write (name, '(a,i1,i1)') trim(table_prefix(group)), i, j
                    kernel%outputs(cursor) = str(trim(name))
                end do
            end do
        end do
        open (newunit=unit, file=trim(path), status='replace', action='write', &
            iostat=ios)
        if (ios /= 0) error stop 'cannot open class-map generated output'
        write (unit, '(a)') chars(emit_kernel(roots, kernel))
        close (unit)
        write (output_unit, '(a)') 'wrote '//trim(path)
    end subroutine emit_output

    character(8) function table_prefix(table_number)
        integer, intent(in) :: table_number

        select case (table_number)
        case (1); table_prefix = 'xi'
        case (2); table_prefix = 'dxi'
        case (3); table_prefix = 'xbeg'
        case (4); table_prefix = 'xend'
        case (5); table_prefix = 'section'
        case (6); table_prefix = 'dsection'
        end select
    end function table_prefix

end program gen_gc_class_coordinate_maps
