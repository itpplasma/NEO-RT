program gen_conserving_basis
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
    type(expr_t) :: vx, vy, vz, vth, b, gu, ge, ht, db
    type(expr_t) :: velocity(3), speed, xi, eta_gradient(3), outputs(10)
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
        error stop 'usage: gen_conserving_basis OUTPUT FORTSYM_REVISION'
    call arena%init()
    vx = sym(arena, 'vx'); vy = sym(arena, 'vy'); vz = sym(arena, 'vz')
    vth = sym(arena, 'vth'); b = sym(arena, 'b')
    gu = sym(arena, 'gu'); ge = sym(arena, 'ge')
    ht = sym(arena, 'ht'); db = sym(arena, 'db')
    velocity = [vx, vy, vz]
    speed = (vx**2 + vy**2 + vz**2)**0.5_real64
    xi = vz/speed
    do i = 1, 3
        eta_gradient(i) = diff((1.0_real64 - xi**2)/b, velocity(i))
        outputs(i) = velocity(i)/speed*gu/vth + eta_gradient(i)*ge
    end do
    outputs(4) = vz/vth
    outputs(5) = (vx**2 + vy**2 + vz**2)/vth**2 - 1.5_real64
    do i = 1, 3
        outputs(5 + i) = diff(outputs(5), velocity(i))
    end do
    outputs(9) = diff(outputs(4), vz)
    outputs(10) = -50.0_real64*(vx**2 + vy**2)*ht*db/vth
    call lower_kernel_ir(outputs, ir, ok, message)
    if (.not. ok) error stop message
    spec%name = str('conserving_velocity_leaf')
    spec%pure_procedure = .true.
    spec%args = [str('vx'), str('vy'), str('vz'), str('vth'), str('b'), &
                 str('gu'), str('ge'), str('ht'), str('db')]
    spec%outputs = [str('gx'), str('gy'), str('gz'), str('p'), str('e'), &
                    str('ex'), str('ey'), str('ez'), str('pz'), str('dpdt')]
    spec%generator = str('gen/gen_conserving_basis.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = &
        str('gen_conserving_basis src/generated/conserving_basis_leaf.inc '// &
            'FORTSYM_REVISION')
    source = emit_fortran_kernel_ir(ir, spec, ok, message)
    if (.not. ok) error stop message
    open (newunit=unit, file=trim(output), status='replace')
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
end program gen_conserving_basis
