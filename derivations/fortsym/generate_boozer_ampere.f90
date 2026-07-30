program generate_boozer_ampere
    ! Derive the Fourier-amplitude form of Ampere's law in coordinates
    ! (s, phi, theta).  The physical real field is
    !
    !   Re[B_hat(s) exp(i*(n*phi + m*theta))].
    !
    ! fortsym differentiates that field through its chart curl operator.  The
    ! result is sampled at zero phase to obtain the real and imaginary parts of
    ! the complex current amplitude, then emitted as a reviewable Fortran
    ! kernel.  "scale" in the generated leaf is 1/(mu0*J), where J is the
    ! signed coordinate Jacobian.
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_arena, only: arena_t
    use fortsym_chart, only: DIM, chart_t, chart_create, curl
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, num, sym, operator(+), operator(-), operator(*)
    use fortsym_kernel, only: KERNEL_SUBROUTINE, emit_kernel, kernel_spec_t
    use fortsym_string, only: chars, str
    use fortsym_subs, only: subs_many
    implicit none

    character(*), parameter :: DEFAULT_OUTPUT = &
        "../../src/generated/boozer_ampere_kernel.f90"
    character(*), parameter :: FORTSYM_REVISION = &
        "fortsym@0cb719ca71c5cb840cd9a83abb8d05c23bc6bbd7"
    type(arena_t), target :: arena
    type(chart_t) :: cartesian
    type(symengine_engine_t) :: engine
    type(expr_t) :: u(DIM), x(DIM), phase, scale, m, n
    type(expr_t) :: bs_re, bs_im, bp_re, bp_im, bt_re, bt_im
    type(expr_t) :: dbp_re, dbp_im, dbt_re, dbt_im
    type(expr_t) :: b_real(DIM), b_quadrature(DIM)
    type(expr_t) :: current_real(DIM), current_imag(DIM), roots(6)
    type(expr_t) :: zero
    type(kernel_spec_t) :: spec
    character(1024) :: output
    integer :: k, unit, status

    output = DEFAULT_OUTPUT
    if (command_argument_count() > 0) call get_command_argument(1, output)
    call arena%init()
    engine = make_symengine_engine(arena)
    u = [sym(arena, "s"), sym(arena, "phi"), sym(arena, "theta")]
    x = u
    cartesian = chart_create(arena, u, x)
    zero = num(arena, 0)

    scale = sym(arena, "scale")
    m = sym(arena, "m")
    n = sym(arena, "n")
    bs_re = sym(arena, "bs_re")
    bs_im = sym(arena, "bs_im")
    bp_re = sym(arena, "bphi_re")
    bp_im = sym(arena, "bphi_im")
    bt_re = sym(arena, "btheta_re")
    bt_im = sym(arena, "btheta_im")
    dbp_re = sym(arena, "dbphi_ds_re")
    dbp_im = sym(arena, "dbphi_ds_im")
    dbt_re = sym(arena, "dbtheta_ds_re")
    dbt_im = sym(arena, "dbtheta_ds_im")
    phase = n*u(2) + m*u(3)

    call make_phase_fields(.false., b_real)
    call make_phase_fields(.true., b_quadrature)
    current_real = curl(cartesian, b_real)
    current_imag = curl(cartesian, b_quadrature)
    call evaluate_at_origin(current_real)
    call evaluate_at_origin(current_imag)

    roots = [ &
        current_real(1), current_imag(1), &
        current_real(2), current_imag(2), &
        current_real(3), current_imag(3)]
    do k = 1, size(roots)
        roots(k) = scale*roots(k)
    end do
    call simplify_all(roots)

    spec%name = str("boozer_ampere_kernel")
    spec%module_name = str("generated_boozer_ampere_kernel")
    spec%mode = KERNEL_SUBROUTINE
    spec%pure_procedure = .true.
    spec%temp_prefix = str("t")
    spec%generator = str("derivations/fortsym/generate_boozer_ampere.f90")
    spec%generator_revision = str(FORTSYM_REVISION)
    spec%regenerate_command = str( &
        "cd derivations/fortsym && fpm run --profile release")
    allocate (spec%args(13), spec%outputs(6))
    spec%args = [ &
        str("scale"), str("m"), str("n"), str("bs_re"), str("bs_im"), &
        str("bphi_re"), str("bphi_im"), str("btheta_re"), str("btheta_im"), &
        str("dbphi_ds_re"), str("dbphi_ds_im"), &
        str("dbtheta_ds_re"), str("dbtheta_ds_im")]
    spec%outputs = [ &
        str("js_re"), str("js_im"), str("jphi_re"), str("jphi_im"), &
        str("jtheta_re"), str("jtheta_im")]

    open (newunit=unit, file=trim(output), status="replace", action="write", iostat=status)
    if (status /= 0) error stop "cannot open generated Ampere kernel"
    write (unit, "(a)") chars(emit_kernel(roots, spec))
    close (unit)
    write (output_unit, "(a)") "generated "//trim(output)

contains

    subroutine make_phase_fields(quadrature, field)
        use fortsym_expr, only: sin, cos
        logical, intent(in) :: quadrature
        type(expr_t), intent(out) :: field(DIM)
        type(expr_t) :: br(DIM), bi(DIM)
        integer :: k

        br = [bs_re, bp_re + dbp_re*u(1), bt_re + dbt_re*u(1)]
        bi = [bs_im, bp_im + dbp_im*u(1), bt_im + dbt_im*u(1)]
        do k = 1, DIM
            if (quadrature) then
                field(k) = bi(k)*cos(phase) + br(k)*sin(phase)
            else
                field(k) = br(k)*cos(phase) - bi(k)*sin(phase)
            end if
        end do
    end subroutine make_phase_fields

    subroutine evaluate_at_origin(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        integer :: k

        do k = 1, size(expressions)
            expressions(k) = subs_many(expressions(k), u, [zero, zero, zero])
        end do
    end subroutine evaluate_at_origin

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        type(engine_result_t) :: result
        integer :: k

        do k = 1, size(expressions)
            result = engine%simplify(expressions(k))
            if (.not. result%ok) error stop "fortsym failed to simplify Ampere kernel"
            expressions(k) = result%value
        end do
    end subroutine simplify_all

end program generate_boozer_ampere
