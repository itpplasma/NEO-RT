program gen_neort_separatrix_frequency
    ! Derive the singular-orbit fit and emit its value/derivative kernel.
    !
    ! The generator belongs to NEO-RT.  fortsym supplies the expression DAG,
    ! native differentiation/simplification, and Fortran emitter; no symbolic
    ! result is transcribed into the consumer by hand.
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym
    implicit none

    character(*), parameter :: DEFAULT_OUTPUT = &
        "src/generated/neort_separatrix_frequency_symbolic.f90"
    type(expr_t) :: z, tau_a, tau_b, tau_c, tau_d
    type(expr_t) :: num_a, num_b, num_c, num_d
    type(expr_t) :: log_distance, tau, numerator, drift
    type(expr_t) :: roots(6)
    type(kernel_spec_t) :: spec
    type(engine_result_t) :: result
    character(1024) :: output
    character(:), allocatable :: source
    logical :: ok
    integer :: unit, ios

    output = DEFAULT_OUTPUT
    call get_command_argument(1, output)
    if (len_trim(output) == 0) output = DEFAULT_OUTPUT

    call reset()
    z = "z"
    tau_a = "tau_a"
    tau_b = "tau_b"
    tau_c = "tau_c"
    tau_d = "tau_d"
    num_a = "num_a"
    num_b = "num_b"
    num_c = "num_c"
    num_d = "num_d"

    ! Both orbit integrals share one endpoint form.  Around a smooth
    ! non-degenerate magnetic maximum the near-saddle passage contributes
    ! log(1/x) and every remaining piece of the period is analytic in x, so
    !
    !   I(x) = A log(1/x) + B + x [C log(1/x) + D] + O(x**2 log(1/x)).
    !
    ! tau is the bounce time; the numerator is the same integral weighted by
    ! the drift integrand.  Their leading coefficients differ only by the
    ! value of that integrand at the barrier top, so the bounce-averaged
    ! drift is the ratio of two functions of this form and its separatrix
    ! limit num_a/tau_a is finite.  The ratio is kept intact rather than
    ! re-expanded: the 1/log series it generates falls only by
    ! tau_b/(tau_a log(1/x)) per order, which is not small at the orbit
    ! distances the integrator resolves, so a truncated 1/log model cannot be
    ! extrapolated from there to the boundary.
    log_distance = log(1/z)
    tau = tau_a*log_distance + tau_b + z*(tau_c*log_distance + tau_d)
    numerator = num_a*log_distance + num_b + z*(num_c*log_distance + num_d)
    drift = numerator/tau

    roots(1) = tau
    roots(2) = simplified_derivative(tau, z)
    roots(3) = simplified(2*pi_expr(default_arena())/tau)
    roots(4) = simplified_derivative(roots(3), z)
    roots(5) = simplified(drift)
    roots(6) = simplified_derivative(drift, z)

    spec%name = str("neort_separatrix_frequency_kernel")
    spec%mode = KERNEL_SUBROUTINE
    spec%module_name = str("neort_separatrix_frequency_symbolic")
    spec%generator = str("gen_neort_separatrix_frequency")
    spec%regenerate_command = str(&
        "cmake --build build --target gen_neort_separatrix_frequency && "// &
        "build/gen_neort_separatrix_frequency")
    spec%pure_procedure = .true.
    spec%temp_prefix = str("sep_t")
    allocate (spec%args(9), spec%outputs(6))
    spec%args = [str("z"), str("tau_a"), str("tau_b"), str("tau_c"), &
        str("tau_d"), str("num_a"), str("num_b"), str("num_c"), &
        str("num_d")]
    spec%outputs = [str("tau"), str("dtau_dz"), str("omega"), &
        str("domega_dz"), str("drift"), str("ddrift_dz")]

    source = chars(emit_kernel(roots, spec, ok))
    if (.not. ok) then
        write (output_unit, '(a)') "gen_neort_separatrix_frequency: emission failed"
        error stop 1
    end if
    if (len(source) > 0) then
        if (source(len(source):len(source)) == new_line("a")) then
            source = source(:len(source) - 1)
        end if
    end if

    open (newunit=unit, file=trim(output), status="replace", action="write", &
        iostat=ios)
    if (ios /= 0) then
        write (output_unit, '(a)') "gen_neort_separatrix_frequency: cannot write "//trim(output)
        error stop 1
    end if
    write (unit, '(a)', advance='no') source
    close (unit)
    write (output_unit, '(a)') "gen_neort_separatrix_frequency: wrote "//trim(output)

contains

    function simplified(expression) result(value)
        type(expr_t), intent(in) :: expression
        type(expr_t) :: value

        result = simplify(expression)
        if (.not. result%ok) then
            write (output_unit, '(a)') "gen_neort_separatrix_frequency: simplify failed"
            error stop 1
        end if
        value = result%value
    end function simplified

    function simplified_derivative(expression, variable) result(value)
        type(expr_t), intent(in) :: expression, variable
        type(expr_t) :: value
        type(engine_result_t) :: derivative

        derivative = diff(expression, variable)
        if (.not. derivative%ok) then
            write (output_unit, '(a)') "gen_neort_separatrix_frequency: differentiation failed"
            error stop 1
        end if
        value = simplified(derivative%value)
    end function simplified_derivative

end program gen_neort_separatrix_frequency
