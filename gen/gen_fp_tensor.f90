program gen_fp_tensor
    use iso_fortran_env, only: output_unit, real64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, exp, operator(*), operator(/), &
                            operator(+), operator(-), operator(**)
    use fortsym_diff, only: diff
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: node_u, node_eta, span_u, span_eta, lo_u, hi_u
    type(expr_t) :: lo_eta, hi_eta, shapes(12), coefficients(3)
    type(expr_t) :: u, tau_u1, nud, duu, eta, inverse_b, weight
    type(expr_t) :: hreal_u1, himag_u1, r_u0, r_u2, phase(3)
    type(kernel_ir_t) :: ir
    type(kernel_emit_spec_t) :: spec
    type(str_t) :: source
    character(len=2048) :: output, revision
    character(:), allocatable :: message
    integer :: unit, i
    logical :: ok

    call get_command_argument(1, output)
    call get_command_argument(2, revision)
    if (len_trim(output) == 0 .or. len_trim(revision) == 0) &
        error stop 'usage: gen_fp_tensor OUTPUT FORTSYM_REVISION'
    call arena%init()
    node_u = sym(arena, 'node_u'); node_eta = sym(arena, 'node_eta')
    span_u = sym(arena, 'span_u'); span_eta = sym(arena, 'span_eta')
    lo_u = (1.0_real64 - node_u)/2.0_real64
    hi_u = (1.0_real64 + node_u)/2.0_real64
    lo_eta = (1.0_real64 - node_eta)/2.0_real64
    hi_eta = (1.0_real64 + node_eta)/2.0_real64
    ! Corners: (u-,eta-), (u+,eta-), (u+,eta+), (u-,eta+).
    shapes(1) = lo_u*lo_eta
    shapes(2) = hi_u*lo_eta
    shapes(3) = hi_u*hi_eta
    shapes(4) = lo_u*hi_eta
    do i = 1, 4
        shapes(4 + i) = 2.0_real64*diff(shapes(i), node_u)/span_u
        shapes(8 + i) = 2.0_real64*diff(shapes(i), node_eta)/span_eta
    end do
    call lower_kernel_ir(shapes, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('neort_fp_q1_shape_scalar')
    spec%pure_procedure = .true.
    spec%args = [str('node_u'), str('node_eta'), str('span_u'), str('span_eta')]
    spec%outputs = [str('p1'), str('p2'), str('p3'), str('p4'), &
                    str('du1'), str('du2'), str('du3'), str('du4'), &
                    str('de1'), str('de2'), str('de3'), str('de4')]
    spec%generator = str('gen/gen_fp_tensor.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = &
        str('gen_fp_tensor src/generated/fp_tensor_leaf.inc FORTSYM_REVISION')
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    open (newunit=unit, file=trim(output), status='replace')
    call write_source(unit, source)
    write (unit, '(a)') ''

    u = sym(arena, 'u'); tau_u1 = sym(arena, 'tau_u1')
    nud = sym(arena, 'nud'); duu = sym(arena, 'duu')
    eta = sym(arena, 'eta'); inverse_b = sym(arena, 'inverse_b')
    weight = u**2*exp(-u**2)*tau_u1
    coefficients(1) = weight
    coefficients(2) = 2.0_real64*eta*weight*nud*(inverse_b - eta)
    coefficients(3) = weight*duu
    call lower_kernel_ir(coefficients, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('neort_fp_tensor_coeff')
    spec%args = [str('u'), str('tau_u1'), str('nud'), str('duu'), &
                 str('eta'), str('inverse_b')]
    spec%outputs = [str('w'), str('p'), str('e')]
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    call write_source(unit, source)
    write (unit, '(a)') ''

    hreal_u1 = sym(arena, 'hreal_u1'); himag_u1 = sym(arena, 'himag_u1')
    r_u0 = sym(arena, 'r_u0'); r_u2 = sym(arena, 'r_u2')
    phase(1) = u**2*hreal_u1
    phase(2) = u**2*himag_u1
    phase(3) = r_u0 + u**2*r_u2
    call lower_kernel_ir(phase, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('neort_fp_phase_scale')
    spec%args = [str('u'), str('hreal_u1'), str('himag_u1'), &
                 str('r_u0'), str('r_u2')]
    spec%outputs = [str('hreal'), str('himag'), str('r')]
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
end program gen_fp_tensor
