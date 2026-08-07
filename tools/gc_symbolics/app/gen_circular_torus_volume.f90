program gen_circular_torus_volume
    !! Derive dV/ds_tor for concentric circular surfaces at finite aspect ratio.
    !! The generated kernel is the independent metric oracle used by the direct
    !! EQDSK regression; no large-aspect approximation enters the result.
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_arena, only: arena_t
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_zero
    use fortsym_diff, only: diff
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, num, pi_expr, operator(+), operator(-), &
        operator(*), operator(/), operator(**), cos, sqrt
    use fortsym_integrate, only: integrate
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str
    use fortsym_subs, only: subs

    implicit none

    character(2048) :: output_path
    character(:), allocatable :: generated
    character(:), allocatable :: why
    integer :: output_length, argument_status, unit, ios
    logical :: integrated
    type(arena_t), target :: arena
    type(symengine_engine_t) :: proof_engine
    type(native_engine_t) :: simplify_engine
    type(engine_result_t) :: simplified
    type(suite_t) :: checks
    type(expr_t) :: s_tor, radius, theta, major_radius, minor_radius
    type(expr_t) :: pi, cylindrical_radius, sqrtg, theta_primitive
    type(expr_t) :: volume_rate, edge_gap, s_from_radius, ds_dr
    type(expr_t) :: dVds, dVds_on_radius
    type(kernel_spec_t) :: spec

    call get_command_argument(1, output_path, length=output_length, &
        status=argument_status)
    if (argument_status /= 0 .or. output_length == 0) then
        write(output_unit, '(a)') 'usage: gen_circular_torus_volume OUTPUT_PATH'
        error stop 2
    end if
    output_path = output_path(:output_length)

    call arena%init()
    proof_engine = make_symengine_engine(arena)
    simplify_engine = make_native_engine(arena)

    s_tor = sym(arena, 's_tor')
    radius = sym(arena, 'radius')
    theta = sym(arena, 'theta')
    major_radius = sym(arena, 'major_radius')
    minor_radius = sym(arena, 'minor_radius')
    pi = pi_expr(arena)

    cylindrical_radius = major_radius + radius*cos(theta)
    sqrtg = cylindrical_radius*radius
    theta_primitive = integrate(arena, sqrtg, theta, integrated, why)
    if (.not. integrated) then
        write(output_unit, '(a)') 'Jacobian integration failed: '//why
        error stop 1
    end if
    volume_rate = 2*pi*(subs(theta_primitive, theta, pi) &
        - subs(theta_primitive, theta, -pi))

    ! For B_phi=F/R, the toroidal flux inside radius r is proportional to
    ! R0-sqrt(R0**2-r**2).  Normalizing at r=a defines s_tor exactly.
    edge_gap = major_radius - sqrt(major_radius**2 - minor_radius**2)
    s_from_radius = (major_radius - sqrt(major_radius**2 - radius**2))/edge_gap
    ds_dr = diff(s_from_radius, radius)

    ! V(r)=2*pi**2*R0*r**2, while
    ! sqrt(R0**2-r**2)=R0-s_tor*edge_gap.
    dVds = 4*pi**2*major_radius*edge_gap &
        *(major_radius - s_tor*edge_gap)
    dVds_on_radius = subs(dVds, s_tor, s_from_radius)

    call suite_begin(checks, 'finite-aspect circular torus volume metric')
    call check_zero(checks, proof_engine, 'integrated cylindrical Jacobian', &
        volume_rate - 4*pi**2*major_radius*radius)
    call check_zero(checks, proof_engine, 'toroidal-flux radial derivative', &
        ds_dr - radius/(edge_gap*sqrt(major_radius**2 - radius**2)))
    call check_zero(checks, proof_engine, 'volume/flux chain rule', &
        dVds_on_radius*ds_dr - volume_rate)
    call check_zero(checks, proof_engine, 'finite magnetic-axis limit', &
        subs(dVds, s_tor, num(arena, 0)) &
        - 4*pi**2*major_radius**2*edge_gap)
    if (checks%failed /= 0) error stop 'circular torus volume proof failed'
    call suite_end(checks)

    simplified = simplify_engine%simplify(dVds)
    if (.not. simplified%ok) error stop 'volume-metric simplification failed'
    dVds = simplified%value
    spec%name = str('evaluate_circular_torus_dvds')
    spec%module_name = str('circular_torus_volume_kernel')
    spec%mode = KERNEL_SUBROUTINE
    spec%pure_procedure = .true.
    spec%generator = str('tools/gc_symbolics/app/gen_circular_torus_volume.f90')
    spec%generator_revision = str( &
        'fortsym@77b031204c76fa88872ddface3af6ac3a25fbb00')
    spec%regenerate_command = str( &
        'cd tools/gc_symbolics && fo exec gen_circular_torus_volume '// &
        '../../test/generated_circular_torus_volume.f90')
    allocate(spec%args(3), spec%outputs(1))
    spec%args = [str('s_tor'), str('major_radius'), str('minor_radius')]
    spec%outputs = [str('dVds')]

    open(newunit=unit, file=trim(output_path), status='replace', &
        action='write', iostat=ios)
    if (ios /= 0) error stop 'cannot open generated volume-metric output'
    generated = chars(emit_kernel([dVds], spec))
    if (generated(len(generated):len(generated)) == new_line('a')) then
        generated = generated(:len(generated) - 1)
    end if
    write(unit, '(a)', advance='no') generated
    close(unit)
    write(output_unit, '(a)') 'wrote '//trim(output_path)
end program gen_circular_torus_volume
