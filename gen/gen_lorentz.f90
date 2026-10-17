program gen_lorentz
    use iso_fortran_env, only: output_unit, real64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(*), operator(/), &
                            operator(+), operator(-)
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: node, span, nud, eta, tau, inverse_b
    type(expr_t) :: shape_roots(4), diffusion_root(1)
    type(kernel_ir_t) :: ir
    type(kernel_emit_spec_t) :: spec
    type(str_t) :: source
    character(len=2048) :: output, revision
    character(:), allocatable :: message
    integer :: unit
    logical :: ok

    call get_command_argument(1, output)
    call get_command_argument(2, revision)
    if (len_trim(output) == 0 .or. len_trim(revision) == 0) &
        error stop 'usage: gen_lorentz OUTPUT FORTSYM_REVISION'
    call arena%init()
    node = sym(arena, 'node'); span = sym(arena, 'span')
    shape_roots(1) = (1.0_real64 - node)/2.0_real64
    shape_roots(2) = (1.0_real64 + node)/2.0_real64
    shape_roots(3) = -1.0_real64/span
    shape_roots(4) = 1.0_real64/span
    call lower_kernel_ir(shape_roots, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('neort_lorentz_shape')
    spec%pure_procedure = .true.
    spec%args = [str('node'), str('span')]
    spec%outputs = [str('p1'), str('p2'), str('g1'), str('g2')]
    spec%generator = str('gen/gen_lorentz.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = &
        str('gen_lorentz src/generated/lorentz_leaf.inc FORTSYM_REVISION')
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    open (newunit=unit, file=trim(output), status='replace')
    call write_source(unit, source)
    write (unit, '(a)') ''

    nud = sym(arena, 'nud'); eta = sym(arena, 'eta')
    tau = sym(arena, 'tau'); inverse_b = sym(arena, 'inverse_b')
    diffusion_root(1) = 2.0_real64*nud*eta*tau*(inverse_b - eta)
    call lower_kernel_ir(diffusion_root, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('neort_lorentz_diffusion')
    spec%args = [str('nud'), str('eta'), str('tau'), str('inverse_b')]
    spec%outputs = [str('a')]
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    call write_source(unit, source)
    close (unit)
    write (output_unit, '(a)') 'Generated '//trim(output)
contains
    subroutine write_source(output_unit, source)
        integer, intent(in) :: output_unit
        type(str_t), intent(in) :: source
        character(:), allocatable :: emitted
        integer :: n

        emitted = chars(source)
        do
            n = len(emitted)
            if (n == 0) exit
            if (emitted(n:n) /= achar(10) .and. emitted(n:n) /= achar(13)) exit
            emitted = emitted(:n - 1)
        end do
        write (output_unit, '(a)') emitted
    end subroutine write_source
end program gen_lorentz
