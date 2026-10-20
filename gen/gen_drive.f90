program gen_drive
    ! Drive weights of the first-order guiding-centre Hamiltonian per unit energy,
    !   H/E = lagrangian*dB_L/B + curvature*kappa.xi_perp,  etab = eta*B = mu*B/E,
    ! i.e. H = mu dB_L + m v_par^2 kappa.xi_perp for deltaA = xi_perp x B0 in the
    ! unperturbed canonical variables (doc/hamiltonian_drive.md).
    use iso_fortran_env, only: output_unit
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(*), operator(-)
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: etab, weights(2)
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
        error stop 'usage: gen_drive OUTPUT FORTSYM_REVISION'
    call arena%init()
    etab = sym(arena, 'etab')
    weights(1) = etab                      ! mu B / E           (mirror force on dB_L)
    weights(2) = 2*(1 - etab)              ! m v_par^2 / E      (curvature)
    call lower_kernel_ir(weights, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('drive_weight_leaf')
    spec%args = [str('etab')]
    spec%outputs = [str('lagrangian'), str('curvature')]
    spec%generator = str('gen/gen_drive.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = str('gen_drive src/generated/drive_leaf.inc FORTSYM_REVISION')
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    open(newunit=unit, file=trim(output), status='replace')
    write(unit, '(a)', advance='no') chars(source)
    close(unit)
    write(output_unit, '(a)') 'Generated '//trim(output)
end program gen_drive
