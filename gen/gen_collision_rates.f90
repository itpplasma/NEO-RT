program gen_collision_rates
    use iso_fortran_env, only: output_unit, real64, int64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, rat, erf, exp, operator(*), operator(/), &
                            operator(+), operator(-), operator(**)
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none

    type(arena_t), target :: arena
    type(expr_t) :: x, x2, inverse_x, sqrt_pi, pi, density, mass, qa, qb
    type(expr_t) :: charge, epsilon0, loglam, polynomial, root(1)
    type(kernel_ir_t) :: ir
    type(kernel_emit_spec_t) :: spec
    type(str_t) :: source
    character(len=2048) :: output, revision
    character(:), allocatable :: message
    integer(int64) :: factorial
    integer :: unit, k
    logical :: ok

    call get_command_argument(1, output)
    call get_command_argument(2, revision)
    if (len_trim(output) == 0 .or. len_trim(revision) == 0) &
        error stop 'usage: gen_collision_rates OUTPUT FORTSYM_REVISION'
    call arena%init()
    x = sym(arena, 'x'); x2 = sym(arena, 'x2')
    inverse_x = sym(arena, 'inverse_x'); sqrt_pi = sym(arena, 'sqrt_pi')
    pi = sym(arena, 'pi'); density = sym(arena, 'density')
    mass = sym(arena, 'mass'); qa = sym(arena, 'qa'); qb = sym(arena, 'qb')
    charge = sym(arena, 'charge'); epsilon0 = sym(arena, 'epsilon0')
    loglam = sym(arena, 'loglam')
    spec%generator = str('gen/gen_collision_rates.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = &
        str('gen_collision_rates src/generated/collision_rate_leaf.inc FORTSYM_REVISION')
    open (newunit=unit, file=trim(output), status='replace')

    ! B(x)/x = 4/sqrt(pi) sum_k (-x^2)^k/[k!(2k+1)(2k+3)].
    ! This follows by integrating (1-z^2) exp[-(xz)^2] on [0,1].
    factorial = 1_int64
    do k = 1, 12
        factorial = factorial*k
    end do
    polynomial = rat(arena, 1_int64, factorial*25_int64*27_int64)
    do k = 11, 0, -1
        factorial = factorial/(k + 1)
        polynomial = rat(arena, (-1_int64)**k, &
                         factorial*int((2*k + 1)*(2*k + 3), int64)) + x2*polynomial
    end do
    root(1) = 4.0_real64*polynomial/sqrt_pi
    call emit('collision_small_shape', [str('x2'), str('sqrt_pi')], 'b_over_x')

    root(1) = (erf(x)*(1.0_real64 - inverse_x**2/2.0_real64) + &
               exp(-x**2)*inverse_x/sqrt_pi)*inverse_x
    call emit('collision_middle_shape', &
              [str('x'), str('inverse_x'), str('sqrt_pi')], 'b_over_x')
    root(1) = inverse_x*(1.0_real64 - inverse_x**2/2.0_real64)
    call emit('collision_large_shape', [str('inverse_x')], 'b_over_x')

    ! G(x)/x = 2/sqrt(pi) integral_0^1 z^2 exp[-(xz)^2] dz.
    ! Its degree-12 series has relative error below 8e-19 on x <= 1/2.
    factorial = 1_int64
    do k = 1, 12
        factorial = factorial*k
    end do
    polynomial = rat(arena, 1_int64, factorial*27_int64)
    do k = 11, 0, -1
        factorial = factorial/(k + 1)
        polynomial = rat(arena, (-1_int64)**k, &
                         factorial*int(2*k + 3, int64)) + x2*polynomial
    end do
    root(1) = 2.0_real64*polynomial/sqrt_pi
    call emit('collision_small_parallel_shape', [str('x2'), str('sqrt_pi')], 'g_over_x')
    root(1) = (erf(x) - 2.0_real64*x*exp(-x**2)/sqrt_pi)*inverse_x**3/2.0_real64
    call emit('collision_middle_parallel_shape', &
              [str('x'), str('inverse_x'), str('sqrt_pi')], 'g_over_x')
    root(1) = inverse_x**3/2.0_real64
    call emit('collision_large_parallel_shape', [str('inverse_x')], 'g_over_x')

    root(1) = 4.0_real64*pi*density*loglam* &
              (qa*qb*charge**2/(4.0_real64*pi*epsilon0*mass))**2
    call emit('collision_prefactor', [str('pi'), str('density'), str('loglam'), &
              str('qa'), str('qb'), str('charge'), str('epsilon0'), str('mass')], &
              'gamma_ab', last=.true.)
    close (unit)
    write (output_unit, '(a)') 'Generated '//trim(output)
contains
    subroutine emit(name, args, result_name, last)
        character(len=*), intent(in) :: name, result_name
        type(str_t), intent(in) :: args(:)
        logical, optional, intent(in) :: last
        character(:), allocatable :: text

        call lower_kernel_ir(root, ir, ok, message)
        if (.not. ok) error stop message
        spec%name = str(name)
        spec%args = args
        spec%outputs = [str(result_name)]
        source = emit_fortran_kernel_ir(ir, spec, ok, message)
        if (.not. ok) error stop message
        if (present(last)) then
            if (last) then
                text = chars(source)
                write (unit, '(a)') text(:len(text) - 1)
                return
            end if
        end if
        write (unit, '(a)') chars(source)
    end subroutine emit
end program gen_collision_rates
