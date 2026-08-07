program gen_manufactured_torus
    !! Derive the manufactured circular-torus field used by the guiding-center
    !! behavioral tests.  fortsym obtains the metric, Jacobian, grad(log B),
    !! and curl(h) from the Cartesian chart and proves the identities that make
    !! the field a valid axisymmetric oracle before emitting Fortran.
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_arena, only: arena_t
    use fortsym_chart, only: DIM, chart_t, chart_create, metric_covariant, &
        jacobian, divergence, curl
    use fortsym_check, only: suite_t, suite_begin, suite_end, check_zero
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_native, only: native_engine_t, make_native_engine
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, sym, num, operator(+), operator(-), &
        operator(*), operator(/), operator(**), sin, cos, sqrt
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    use fortsym_string, only: chars, str

    implicit none

    character(2048) :: output_path
    integer :: output_length, argument_status, unit, ios, i, j
    type(arena_t), target :: arena
    type(symengine_engine_t) :: proof_engine
    type(native_engine_t) :: simplify_engine
    type(suite_t) :: checks
    type(chart_t) :: chart
    type(expr_t) :: u(DIM), xyz(DIM), metric(DIM, DIM)
    type(expr_t) :: radius, phi, theta, major_radius, toroidal_flux
    type(expr_t) :: q_axis, shear_shape, cylindrical_radius, poloidal_flux_rate
    type(expr_t) :: b_con(DIM), b_cov(DIM), h_con(DIM), h_cov(DIM)
    type(expr_t) :: grad_log_b(DIM), curl_h(DIM), chart_curl(DIM), grad_psi(DIM)
    type(expr_t) :: bmod, b_squared, sqrtg, psi, q_safety, magnetic_shear
    type(expr_t) :: field_norm
    type(expr_t) :: roots(20), norm_h
    type(kernel_spec_t) :: spec

    call get_command_argument(1, output_path, length=output_length, &
        status=argument_status)
    if (argument_status /= 0 .or. output_length == 0) then
        write(output_unit, '(a)') &
            'usage: gen_manufactured_torus OUTPUT_PATH'
        error stop 2
    end if
    output_path = output_path(:output_length)

    call arena%init()
    proof_engine = make_symengine_engine(arena)
    simplify_engine = make_native_engine(arena)

    radius = sym(arena, 'radius')
    phi = sym(arena, 'phi')
    theta = sym(arena, 'theta')
    major_radius = sym(arena, 'major_radius')
    toroidal_flux = sym(arena, 'toroidal_flux')
    q_axis = sym(arena, 'q_axis')
    shear_shape = sym(arena, 'shear_shape')
    u = [radius, phi, theta]

    cylindrical_radius = major_radius + radius*cos(theta)
    xyz(1) = cylindrical_radius*cos(phi)
    xyz(2) = cylindrical_radius*sin(phi)
    xyz(3) = radius*sin(theta)
    chart = chart_create(arena, u, xyz)
    metric = metric_covariant(chart)
    sqrtg = cylindrical_radius*radius

    ! J B^theta=d(psi)/dr is a flux function, so div(B)=0 exactly.
    poloidal_flux_rate = toroidal_flux*radius &
        *(1 + shear_shape*(radius/major_radius)**2)/(q_axis*major_radius)
    b_con(1) = num(arena, 0)
    b_con(2) = toroidal_flux/cylindrical_radius**2
    b_con(3) = poloidal_flux_rate/(cylindrical_radius*radius)
    do i = 1, DIM
        b_cov(i) = metric(i, 1)*b_con(1)
        do j = 2, DIM
            b_cov(i) = b_cov(i) + metric(i, j)*b_con(j)
        end do
    end do
    b_squared = b_cov(1)*b_con(1) + b_cov(2)*b_con(2) &
        + b_cov(3)*b_con(3)
    field_norm = sqrt(toroidal_flux**2 + poloidal_flux_rate**2)
    bmod = field_norm/cylindrical_radius
    h_con(1) = num(arena, 0)
    h_con(2) = toroidal_flux/(cylindrical_radius*field_norm)
    h_con(3) = poloidal_flux_rate/(radius*field_norm)
    h_cov(1) = num(arena, 0)
    h_cov(2) = toroidal_flux*cylindrical_radius/field_norm
    h_cov(3) = radius*poloidal_flux_rate/field_norm
    do i = 1, DIM
        grad_log_b(i) = diff(bmod, u(i))/bmod
    end do
    curl_h(1) = -diff(h_cov(2), theta)/sqrtg
    curl_h(2) = -diff(h_cov(3), radius)/sqrtg
    curl_h(3) = diff(h_cov(2), radius)/sqrtg
    chart_curl = curl(chart, h_cov)

    psi = toroidal_flux/(q_axis*major_radius) &
        *(radius**2/2 + shear_shape*radius**4/(4*major_radius**2))
    do i = 1, DIM
        grad_psi(i) = diff(psi, u(i))
    end do
    q_safety = toroidal_flux*radius &
        /(poloidal_flux_rate*sqrt(major_radius**2 - radius**2))
    magnetic_shear = radius*diff(q_safety, radius)/q_safety

    call suite_begin(checks, 'manufactured torus geometry')
    call check_zero(checks, proof_engine, 'signed chart Jacobian', &
        jacobian(chart) - sqrtg)
    call check_zero(checks, proof_engine, 'field is divergence free', &
        divergence(chart, b_con))
    call check_zero(checks, proof_engine, 'field magnitude from metric', &
        b_squared - bmod**2)
    norm_h = h_cov(1)*h_con(1) + h_cov(2)*h_con(2) &
        + h_cov(3)*h_con(3)
    call check_zero(checks, proof_engine, 'unit field normalization', norm_h - 1)
    call check_zero(checks, proof_engine, 'canonical flux generates poloidal field', &
        grad_psi(1) - sqrtg*b_con(3))
    do i = 1, DIM
        call check_zero(checks, proof_engine, 'manual curl matches chart curl', &
            curl_h(i) - chart_curl(i))
    end do
    if (checks%failed /= 0) error stop 'manufactured torus proof failed'
    call suite_end(checks)

    roots(1) = bmod
    roots(2) = sqrtg
    roots(3:5) = grad_log_b
    roots(6:8) = h_cov
    roots(9:11) = h_con
    roots(12:14) = curl_h
    roots(15) = psi
    roots(16:18) = grad_psi
    roots(19) = q_safety
    roots(20) = magnetic_shear
    call simplify_all(roots)

    spec%name = str('evaluate_manufactured_torus')
    spec%module_name = str('manufactured_torus_kernel')
    spec%mode = KERNEL_SUBROUTINE
    spec%pure_procedure = .true.
    spec%generator = str('tools/gc_symbolics/app/gen_manufactured_torus.f90')
    spec%generator_revision = str( &
        'fortsym@a2b9bb353816b03445282fc71d0584b34a787549')
    spec%regenerate_command = str( &
        'cd tools/gc_symbolics && fo exec gen_manufactured_torus ../../test/generated_manufactured_torus.f90')
    allocate(spec%args(7), spec%outputs(20))
    spec%args = [str('radius'), str('phi'), str('theta'), &
        str('major_radius'), str('toroidal_flux'), str('q_axis'), &
        str('shear_shape')]
    spec%outputs = [str('bmod'), str('sqrtg'), &
        str('grad_log_b_1'), str('grad_log_b_2'), str('grad_log_b_3'), &
        str('h_cov_1'), str('h_cov_2'), str('h_cov_3'), &
        str('h_con_1'), str('h_con_2'), str('h_con_3'), &
        str('curl_h_1'), str('curl_h_2'), str('curl_h_3'), &
        str('psi'), str('grad_psi_1'), str('grad_psi_2'), &
        str('grad_psi_3'), str('q_safety'), str('magnetic_shear')]

    open(newunit=unit, file=trim(output_path), status='replace', &
        action='write', iostat=ios)
    if (ios /= 0) error stop 'cannot open generated torus output'
    write(unit, '(a)') chars(emit_kernel(roots, spec))
    close(unit)
    write(output_unit, '(a)') 'wrote '//trim(output_path)

contains

    subroutine simplify_all(expressions)
        type(expr_t), intent(inout) :: expressions(:)
        type(engine_result_t) :: simplified
        integer :: k

        do k = 1, size(expressions)
            simplified = simplify_engine%simplify(expressions(k))
            if (.not. simplified%ok) then
                write(output_unit, '(a,i0)') 'simplification failed for output ', k
                error stop 1
            end if
            expressions(k) = simplified%value
        end do
    end subroutine simplify_all

end program gen_manufactured_torus
