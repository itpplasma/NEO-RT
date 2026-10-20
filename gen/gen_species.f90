program gen_species
    use iso_fortran_env, only: output_unit
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(*), operator(/), operator(+), operator(-)
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: base, gradient, delta, radius, speed, dspeed, omega, domega
    type(expr_t) :: affine(1), rotation(2)
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
        error stop 'usage: gen_species OUTPUT FORTSYM_REVISION'
    call arena%init()
    base = sym(arena, 'base'); gradient = sym(arena, 'gradient'); delta = sym(arena, 'delta')
    affine(1) = base + gradient*delta
    call lower_kernel_ir(affine, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('species_affine_leaf')
    spec%args = [str('base'), str('gradient'), str('delta')]
    spec%outputs = [str('value')]
    spec%generator = str('gen/gen_species.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = str('gen_species src/generated/species_leaf.inc FORTSYM_REVISION')
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    open(newunit=unit, file=trim(output), status='replace')
    write(unit, '(a)') chars(source)
    radius = sym(arena, 'radius'); speed = sym(arena, 'speed')
    dspeed = sym(arena, 'dspeed'); omega = sym(arena, 'omega'); domega = sym(arena, 'domega')
    rotation(1) = radius*omega/speed
    rotation(2) = radius*(domega/speed - omega*dspeed/(speed*speed))
    call lower_kernel_ir(rotation, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('species_rotation_leaf')
    spec%args = [str('radius'), str('speed'), str('dspeed'), str('omega'), str('domega')]
    spec%outputs = [str('mach'), str('dmach')]
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    write(unit, '(a)', advance='no') chars(source)
    close(unit)
    write(output_unit, '(a)') 'Generated '//trim(output)
end program gen_species
