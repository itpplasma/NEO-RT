program gen_collision_rates
    use iso_fortran_env, only: output_unit, real64, int64
    use fortsym_expr, only: expr_t, sym, rat, erf, exp, operator(*), operator(/), &
                            operator(+), operator(-), operator(**)
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t, words
    implicit none

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
    x = sym('x'); x2 = sym('x2')
    inverse_x = sym('inverse_x'); sqrt_pi = sym('sqrt_pi')
    pi = sym('pi'); density = sym('density')
    mass = sym('mass'); qa = sym('qa'); qb = sym('qb')
    charge = sym('charge'); epsilon0 = sym('epsilon0')
    loglam = sym('loglam')
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
    polynomial = rat(1_int64, factorial*25_int64*27_int64)
    do k = 11, 0, -1
        factorial = factorial/(k + 1)
        polynomial = rat((-1_int64)**k, &
                         factorial*int((2*k + 1)*(2*k + 3), int64)) + x2*polynomial
    end do
    root(1) = 4.0_real64*polynomial/sqrt_pi
    call emit('collision_small_shape', words('x2 sqrt_pi'), 'b_over_x')

    root(1) = (erf(x)*(1.0_real64 - inverse_x**2/2.0_real64) + &
               exp(-x**2)*inverse_x/sqrt_pi)*inverse_x
    call emit('collision_middle_shape', &
        words('x inverse_x sqrt_pi'), 'b_over_x')
    root(1) = inverse_x*(1.0_real64 - inverse_x**2/2.0_real64)
    call emit('collision_large_shape', [str('inverse_x')], 'b_over_x')
    root(1) = 4.0_real64*pi*density*loglam* &
              (qa*qb*charge**2/(4.0_real64*pi*epsilon0*mass))**2
    call emit('collision_prefactor', &
        words('pi density loglam qa qb charge epsilon0 mass'), 'gamma_ab')
    close (unit)
    write (output_unit, '(a)') 'Generated '//trim(output)
contains
    subroutine emit(name, args, result_name)
        character(len=*), intent(in) :: name, result_name
        type(str_t), intent(in) :: args(:)

        call lower_kernel_ir(root, ir, ok, message)
        if (.not. ok) error stop message
        spec%name = str(name)
        spec%args = args
        spec%outputs = [str(result_name)]
        source = emit_fortran_kernel_ir(ir, spec, ok, message)
        if (.not. ok) error stop message
        write (unit, '(a)') chars(source)
    end subroutine emit
end program gen_collision_rates
