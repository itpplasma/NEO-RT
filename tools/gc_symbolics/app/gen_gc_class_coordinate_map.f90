program gen_gc_class_coordinate_map
    !! Generate the bounded POTATO class-coordinate map and its derivative.
    !!
    !! The four boundary families are the pinned POTATO conventions:
    !! 1 regular rho boundary, 2 inner/turning boundary, 3 regular
    !! separatrix, and 4 X-point.  The runtime only validates inputs and
    !! dispatches this generated table.
    use, intrinsic :: iso_fortran_env, only: int64, output_unit
    use fortsym_arena, only: arena_t
    use fortsym_check, only: check_identity, suite_begin, suite_end, suite_t
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: make_native_engine, native_engine_t
    use fortsym_engine_symengine, only: make_symengine_engine, symengine_engine_t
    use fortsym_expr, only: acos, cos, expr_t, log, num, operator(+), &
        operator(-), operator(*), operator(/), operator(**), pi_expr, rat, &
        sqrt, sym, tanh, atanh
    use fortsym_kernel, only: emit_kernel, kernel_spec_t, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    implicit none

    character(*), parameter :: FORTSYM_REVISION = &
        'fortsym@545788453a204d58705f735b519c3863c2f734c8'
    character(*), parameter :: REGENERATE_COMMAND = &
        'cd tools/gc_symbolics && fo exec gen_gc_class_coordinate_map '// &
        'OUTPUT_FORTRAN OUTPUT_INVENTORY'

    character(2048) :: output_fortran, output_inventory
    integer :: output_length, argument_status
    type(arena_t), target :: arena
    type(symengine_engine_t) :: proof_engine
    type(native_engine_t) :: simplify_engine
    type(suite_t) :: proofs
    type(expr_t) :: x, relmargin, widthclass, log_ratio
    type(expr_t) :: xi(4, 4), dxi(4, 4), xbeg(4, 4), xend(4, 4)
    type(expr_t) :: xpoint_left_limit(4, 4), xpoint_right_limit(4, 4)
    type(expr_t), allocatable :: roots(:)
    integer :: left_type, right_type, cursor

    call get_command_argument(1, output_fortran, length=output_length, &
        status=argument_status)
    if (argument_status /= 0 .or. output_length == 0) then
        write (output_unit, '(a)') &
            'usage: gen_gc_class_coordinate_map OUTPUT_FORTRAN OUTPUT_INVENTORY'
        error stop 2
    end if
    output_fortran = output_fortran(:output_length)
    call get_command_argument(2, output_inventory, length=output_length, &
        status=argument_status)
    if (argument_status /= 0 .or. output_length == 0) then
        write (output_unit, '(a)') &
            'usage: gen_gc_class_coordinate_map OUTPUT_FORTRAN OUTPUT_INVENTORY'
        error stop 2
    end if
    output_inventory = output_inventory(:output_length)

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
                right_type, relmargin, log_ratio)
            xend(left_type, right_type) = right_bound(arena, left_type, &
                right_type, relmargin, log_ratio)
            xpoint_left_limit(left_type, right_type) = num(arena, 0)
            xpoint_right_limit(left_type, right_type) = num(arena, 1)
        end do
    end do

    call prove_contracts(arena, proof_engine, proofs, x, relmargin, &
        widthclass, log_ratio, xi, dxi, xbeg, xend, xpoint_left_limit, &
        xpoint_right_limit)
    call simplify_all(simplify_engine, xi)
    call simplify_all(simplify_engine, dxi)
    call simplify_all(simplify_engine, xbeg)
    call simplify_all(simplify_engine, xend)
    call simplify_all(simplify_engine, xpoint_left_limit)
    call simplify_all(simplify_engine, xpoint_right_limit)

    allocate (roots(6*4*4))
    cursor = 0
    call append_table(roots, cursor, xi)
    call append_table(roots, cursor, dxi)
    call append_table(roots, cursor, xbeg)
    call append_table(roots, cursor, xend)
    call append_table(roots, cursor, xpoint_left_limit)
    call append_table(roots, cursor, xpoint_right_limit)
    call emit_fortran(trim(output_fortran), roots)
    call emit_inventory(trim(output_inventory))

contains

    function map_expression(arena, left_type, right_type, argument) result(value)
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
            case (2); value = 1-(1-argument)**2
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

    function left_bound(arena, left_type, right_type, margin, log_ratio) result(value)
        type(arena_t), target, intent(inout) :: arena
        integer, intent(in) :: left_type, right_type
        type(expr_t), intent(in) :: margin, log_ratio
        type(expr_t) :: value

        select case (left_type)
        case (1); value = num(arena, 0)
        case (2)
            select case (right_type)
            case (1); value = sqrt(margin)
            case (2); value = acos(1-2*margin)/pi_expr(arena)
            case default; value = atanh(sqrt(margin))
            end select
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
        else if (boundary_type == 2) then
            select case (left_type)
            case (1); value = 1-sqrt(margin)
            case (2); value = 1-acos(1-2*margin)/pi_expr(arena)
            case default; value = -atanh(sqrt(margin))
            end select
        else if (left_type >= 3) then
            value = num(arena, 0)
        else
            select case (boundary_type)
            case (1); value = num(arena, 1)
            case default; value = num(arena, 0)
            end select
        end if
    end function right_bound

    subroutine prove_contracts(arena, engine, checks, x, margin, width, &
            log_ratio, xi, dxi, xbeg, xend, xpoint_left_limit, &
            xpoint_right_limit)
        type(arena_t), target, intent(inout) :: arena
        type(symengine_engine_t) :: engine
        type(suite_t), intent(out) :: checks
        type(expr_t), intent(in) :: x, margin, width, log_ratio
        type(expr_t), intent(in) :: xi(4, 4), dxi(4, 4), xbeg(4, 4), xend(4, 4)
        type(expr_t), intent(in) :: xpoint_left_limit(4, 4), &
            xpoint_right_limit(4, 4)
        integer :: left_type, right_type
        type(expr_t) :: expected_left, expected_right

        call suite_begin(checks, 'POTATO class-coordinate Fortsym contracts')
        do left_type = 1, 4
            do right_type = 1, 4
                call check_identity(checks, engine, 'map derivative', &
                    dxi(left_type, right_type)-diff(xi(left_type, right_type), x))
                expected_left = left_bound(arena, left_type, right_type, &
                    margin, log_ratio)
                expected_right = right_bound(arena, left_type, right_type, &
                    margin, log_ratio)
                call check_identity(checks, engine, 'left bound', &
                    xbeg(left_type, right_type)-expected_left)
                call check_identity(checks, engine, 'right bound', &
                    xend(left_type, right_type)-expected_right)
                call check_identity(checks, engine, 'left X-point limit', &
                    xpoint_left_limit(left_type, right_type)-num(arena, 0))
                call check_identity(checks, engine, 'right X-point limit', &
                    xpoint_right_limit(left_type, right_type)-num(arena, 1))
            end do
        end do
        if (checks%failed /= 0) error stop 'class-coordinate proof failed'
        call suite_end(checks)
    end subroutine prove_contracts

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

    subroutine emit_fortran(path, roots)
        character(*), intent(in) :: path
        type(expr_t), intent(in) :: roots(:)
        type(kernel_spec_t) :: kernel
        integer :: unit, ios, i, j, group, cursor
        character(32) :: name

        kernel%name = str('evaluate_neort_gc_class_coordinate_map')
        kernel%module_name = str('neort_gc_class_coordinate_map_generated')
        kernel%mode = KERNEL_SUBROUTINE
        kernel%pure_procedure = .true.
        kernel%generator = str( &
            'tools/gc_symbolics/app/gen_gc_class_coordinate_map.f90')
        kernel%generator_revision = str(FORTSYM_REVISION)
        kernel%regenerate_command = str(REGENERATE_COMMAND)
        allocate (kernel%args(3), kernel%outputs(size(roots)))
        kernel%args = [str('class_coordinate'), str('relative_margin'), &
            str('relative_class_width')]
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
    end subroutine emit_fortran

    subroutine emit_inventory(path)
        character(*), intent(in) :: path
        integer :: unit, ios, left_type, right_type

        open (newunit=unit, file=trim(path), status='replace', action='write', &
            iostat=ios)
        if (ios /= 0) error stop 'cannot open class-map inventory'
        write (unit, '(a)') '# neort_gc_class_coordinate_map_inventory_v2'
        write (unit, '(a)') 'generator=tools/gc_symbolics/app/'// &
            'gen_gc_class_coordinate_map.f90'
        write (unit, '(a)') 'fortsym_revision='//FORTSYM_REVISION
        write (unit, '(a)') 'kernel=neort_gc_class_coordinate_map_generated:'// &
            'evaluate_neort_gc_class_coordinate_map'
        write (unit, '(a)') 'boundary_convention=1:rho;2:inner;3:separatrix;'// &
            '4:xpoint'
        write (unit, '(a)') 'orientation=ifuntype=10*left+right;34 and 43 '// &
            'use the same increasing 1/2*(1+tanh(x)) map'
        write (unit, '(a)') 'log_truncation=type3:0.5*log(margin/width);'// &
            'type4:0.25*log(margin/width)'
        write (unit, '(a)') 'inner_margin=quadratic maps use exact inverse '// &
            'sqrt/acos/atanh charts'
        write (unit, '(a)') 'xpoint_limits=generated left=0;right=1 coordinate limits'
        write (unit, '(a)') 'outputs=xi,dxi_dx,xbeg,xend,xpoint_limits for all 16 pairs'
        do left_type = 1, 4
            do right_type = 1, 4
                write (unit, '(a,i1,a,i1,a,i2)') 'pair=', left_type, ':', &
                    right_type, ':ifuntype=', 10*left_type+right_type
            end do
        end do
        close (unit)
        write (output_unit, '(a)') 'wrote '//trim(path)
    end subroutine emit_inventory

    character(8) function table_prefix(table_number)
        integer, intent(in) :: table_number

        select case (table_number)
        case (1); table_prefix = 'xi'
        case (2); table_prefix = 'dxi'
        case (3); table_prefix = 'xbeg'
        case (4); table_prefix = 'xend'
        case (5); table_prefix = 'xpoint_left'
        case (6); table_prefix = 'xpoint_right'
        end select
    end function table_prefix

end program gen_gc_class_coordinate_map
