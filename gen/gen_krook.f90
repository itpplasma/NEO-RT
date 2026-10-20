program gen_krook
    use iso_fortran_env, only: output_unit, real64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, sin, cos, operator(*), operator(/), &
                            operator(+), operator(-), operator(**)
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: r, nu, scale, pi, base, speed, harmonic, geometry, roots(3)
    type(expr_t) :: angle, node, span, abs_span, angle_roots(2)
    type(expr_t) :: da, db, secant
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
        error stop 'usage: gen_krook OUTPUT FORTSYM_REVISION'
    call arena%init()
    r = sym(arena, 'r'); nu = sym(arena, 'nu'); scale = sym(arena, 'scale')
    pi = sym(arena, 'pi'); base = sym(arena, 'base'); speed = sym(arena, 'speed')
    harmonic = sym(arena, 'harmonic'); geometry = sym(arena, 'geometry')
    roots(1) = nu/(pi*scale*(r**2 + nu**2))
    roots(2) = base*geometry/speed**3
    roots(3) = roots(2)*(1.0_real64 + harmonic**2/4.0_real64)
    call lower_kernel_ir(roots, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('krook_leaf')
    spec%args = [str('r'), str('nu'), str('scale'), str('pi'), str('base'), &
                 str('speed'), str('harmonic'), str('geometry')]
    spec%outputs = [str('kernel'), str('energy_rate'), str('harmonic_rate')]
    spec%generator = str('gen/gen_krook.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = str('gen_krook src/generated/krook_leaf.inc FORTSYM_REVISION')
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    open (newunit=unit, file=trim(output), status='replace')
    write (unit, '(a)') chars(source)
    node = sym(arena, 'node'); span = sym(arena, 'span')
    abs_span = sym(arena, 'abs_span')
    angle = span*sin(pi*node/2.0_real64)**2
    angle_roots(1) = scale*(r*cos(angle) + nu*sin(angle))/ &
                     (nu*cos(angle) - r*sin(angle))
    angle_roots(2) = abs_span*sin(pi*node)/2.0_real64
    call lower_kernel_ir(angle_roots, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('krook_angle_leaf')
    spec%args = [str('r'), str('nu'), str('scale'), str('pi'), str('node'), &
                 str('span'), str('abs_span')]
    spec%outputs = [str('target'), str('jacobian')]
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    write (unit, '(a)') chars(source)
    da = sym(arena, 'da'); db = sym(arena, 'db'); secant = sym(arena, 'secant')
    angle_roots(1) = 3.0_real64*da + 3.0_real64*db - 6.0_real64*secant
    angle_roots(2) = 6.0_real64*secant - 4.0_real64*da - 2.0_real64*db
    call lower_kernel_ir(angle_roots, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('krook_hermite_leaf')
    spec%args = [str('da'), str('db'), str('secant')]
    spec%outputs = [str('quadratic'), str('linear')]
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    write (unit, '(a)') chars(source)
    close (unit)
    write (output_unit, '(a)') 'Generated '//trim(output)
end program gen_krook
