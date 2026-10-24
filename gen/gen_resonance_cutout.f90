program gen_resonance_cutout
    use iso_fortran_env, only: output_unit, real64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, sqrt, operator(*), operator(/), &
                            operator(+), operator(-)
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: eta, tp, dt, q, inv_b, gap, ratio, m, nu, average
    type(expr_t) :: diffusion, curvature, layer, dlo, dhi, h
    type(expr_t) :: roots(4)
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
        error stop 'usage: gen_resonance_cutout OUTPUT FORTSYM_REVISION'
    call arena%init()
    eta = sym(arena, 'eta'); tp = sym(arena, 'etatp'); dt = sym(arena, 'etadt')
    q = sym(arena, 'q'); inv_b = sym(arena, 'inv_b0')
    gap = sym(arena, 'gap'); ratio = sym(arena, 'ratio'); m = sym(arena, 'm')
    nu = sym(arena, 'nu_d'); average = sym(arena, 'average')
    diffusion = sym(arena, 'diffusion'); curvature = sym(arena, 'curvature')
    layer = sym(arena, 'layer')
    dlo = sym(arena, 'dlo'); dhi = sym(arena, 'dhi'); h = sym(arena, 'h')
    spec%generator = str('gen/gen_resonance_cutout.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = &
        str('gen_resonance_cutout src/generated/resonance_cutout_leaf.inc '// &
            'FORTSYM_REVISION')
    spec%pure_procedure = .true.
    open (newunit=unit, file=trim(output), status='replace')

    ! The large-aspect-ratio orbit parameters; branch guards stay in the caller.
    roots(1) = 2.0_real64 / (1.0_real64 / tp + 1.0_real64 / dt)
    roots(2) = eta * (1.0_real64 / tp - 1.0_real64 / dt)
    roots(3) = (dt - eta) / (dt - tp)
    roots(4) = 1.0_real64 - eta / dt
    call emit('cutout_pitch_parameters', [str('eta'), str('etatp'), str('etadt')], &
              [str('inv_b0'), str('q'), str('trapped_m'), str('gap')], 4)
    roots(1) = q / gap
    call emit('cutout_passing_modulus', [str('q'), str('gap')], [str('k2')], 1)
    roots(1) = q * inv_b * (ratio - (1.0_real64 - m))
    call emit('cutout_trapped_average', &
              [str('q'), str('inv_b0'), str('ratio'), str('m')], [str('average')], 1)
    roots(1) = gap * inv_b * ratio
    call emit('cutout_passing_average', &
              [str('gap'), str('inv_b0'), str('ratio')], [str('average')], 1)
    roots(1) = 2.0_real64 * nu * eta * average
    call emit('cutout_pitch_diffusion', &
              [str('nu_d'), str('eta'), str('average')], [str('diffusion')], 1)
    ! Characteristic fold width normalized by the linear-layer-fit criterion.
    roots(1) = sqrt(sqrt(diffusion / curvature))
    call emit('cutout_fold_width', &
              [str('diffusion'), str('curvature')], [str('layer')], 1)
    roots(1) = 1.0_real64 / (curvature * layer)
    call emit('cutout_boundary_fill', &
              [str('curvature'), str('layer')], [str('weight')], 1)
    roots(1) = (dhi - dlo) / (2.0_real64 * h)
    call emit('cutout_curvature', &
              [str('dlo'), str('dhi'), str('h')], [str('curvature')], 1)
    close (unit)
    write (output_unit, '(a)') 'Generated '//trim(output)
contains
    subroutine emit(name, args, outputs, count)
        character(len=*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        integer, intent(in) :: count
        call lower_kernel_ir(roots(:count), ir, ok, message)
        if (.not. ok) error stop message
        spec%name = str(name)
        spec%args = args
        spec%outputs = outputs
        source = emit_fortran_kernel_ir(ir, spec, ok, message)
        if (.not. ok) error stop message
        write (unit, '(a)') chars(source)
    end subroutine emit
end program gen_resonance_cutout
