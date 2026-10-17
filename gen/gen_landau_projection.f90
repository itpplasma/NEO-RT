program gen_landau_projection
    use iso_fortran_env, only: output_unit, real64
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, sqrt, sin, cos, operator(*), &
        operator(/), operator(+), operator(-), operator(**)
    use fortsym_kernel_emit, only: kernel_emit_spec_t, emit_fortran_kernel_ir
    use fortsym_kernel_ir, only: kernel_ir_t, lower_kernel_ir
    use fortsym_string, only: chars, str, str_t
    implicit none
    type(arena_t), target :: arena
    type(kernel_emit_spec_t) :: spec
    character(len=2048) :: output, revision
    integer :: unit

    call get_command_argument(1, output)
    call get_command_argument(2, revision)
    if (len_trim(output) == 0 .or. len_trim(revision) == 0) &
        error stop 'usage: gen_landau_projection OUTPUT FORTSYM_REVISION'
    call arena%init()
    spec%generator = str('gen/gen_landau_projection.f90')
    spec%generator_revision = str(trim(revision))
    spec%regenerate_command = &
        str('gen_landau_projection src/generated/landau_projection_leaf.inc REVISION')
    spec%pure_procedure = .true.
    open (newunit=unit, file=trim(output), status='replace')
    call generate_scales()
    call generate_velocities()
    call generate_transverse()
    close (unit)
    write (output_unit, '(a)') 'Generated '//trim(output)
contains
    subroutine generate_scales()
        type(expr_t) :: ma, mb, energy, za, zb, e, eps, loglam, na, nb, pi
        type(expr_t) :: root(5), relative, coupling

        ma = sym(arena, 'ma'); mb = sym(arena, 'mb')
        energy = sym(arena, 'energy'); pi = sym(arena, 'pi')
        za = sym(arena, 'za'); zb = sym(arena, 'zb')
        e = sym(arena, 'e'); eps = sym(arena, 'eps')
        loglam = sym(arena, 'loglam'); na = sym(arena, 'na'); nb = sym(arena, 'nb')
        relative = sqrt(2.0_real64*energy*(ma + mb)/(ma*mb))
        coupling = 2.0_real64*pi*loglam*(za*zb*e**2/(4.0_real64*pi*eps))**2
        root = [sqrt(2.0_real64*energy/(ma + mb)), relative, &
            ma/(ma + mb), mb/(ma + mb), coupling*na*nb/(2.0_real64*pi**2*relative)]
        call emit('landau_pair_scales', [str('ma'), str('mb'), str('energy'), &
            str('za'), str('zb'), str('e'), str('eps'), str('loglam'), &
            str('na'), str('nb'), str('pi')], &
            [str('vcm'), str('vrel'), str('fa'), str('fb'), str('prefactor')], root)
    end subroutine generate_scales

    subroutine generate_velocities()
        type(expr_t) :: vcm, vrel, fa, fb, mu_cm, mu_rel, phi
        type(expr_t) :: center(3), direction(3), root(9)
        integer :: k

        vcm = sym(arena, 'vcm'); vrel = sym(arena, 'vrel')
        fa = sym(arena, 'fa'); fb = sym(arena, 'fb')
        mu_cm = sym(arena, 'mu_cm'); mu_rel = sym(arena, 'mu_rel')
        phi = sym(arena, 'phi')
        center = [vcm*sqrt(1.0_real64 - mu_cm**2), 0.0_real64*vcm, vcm*mu_cm]
        direction = [sqrt(1.0_real64 - mu_rel**2)*cos(phi), &
            sqrt(1.0_real64 - mu_rel**2)*sin(phi), mu_rel]
        do k = 1, 3
            root(k) = center(k) + fb*vrel*direction(k)
            root(k + 3) = center(k) - fa*vrel*direction(k)
            root(k + 6) = direction(k)
        end do
        call emit('landau_pair_velocities', [str('vcm'), str('vrel'), str('fa'), &
            str('fb'), str('mu_cm'), str('mu_rel'), str('phi')], &
            [str('va1'), str('va2'), str('va3'), str('vb1'), str('vb2'), str('vb3'), &
            str('n1'), str('n2'), str('n3')], root)
    end subroutine generate_velocities

    subroutine generate_transverse()
        type(expr_t) :: n(3), g(3), root(3)

        n = [sym(arena, 'n1'), sym(arena, 'n2'), sym(arena, 'n3')]
        g = [sym(arena, 'g1'), sym(arena, 'g2'), sym(arena, 'g3')]
        ! |n cross g|^2 = conj(g).(I-n*n^T).g for a unit real n.
        ! The Gram construction keeps the finite quadrature matrix positive.
        root = [n(2)*g(3) - n(3)*g(2), n(3)*g(1) - n(1)*g(3), &
                n(1)*g(2) - n(2)*g(1)]
        call emit('landau_transverse', [str('n1'), str('n2'), str('n3'), &
            str('g1'), str('g2'), str('g3')], [str('p1'), str('p2'), str('p3')], root)
    end subroutine generate_transverse

    subroutine emit(name, args, outputs, roots)
        character(len=*), intent(in) :: name
        type(str_t), intent(in) :: args(:), outputs(:)
        type(expr_t), intent(in) :: roots(:)
        type(kernel_ir_t) :: ir
        type(str_t) :: source
        character(:), allocatable :: message, text
        integer :: n
        logical :: ok

        call lower_kernel_ir(roots, ir, ok, message)
        if (.not. ok) error stop message
        spec%name = str(name); spec%args = args; spec%outputs = outputs
        source = emit_fortran_kernel_ir(ir, spec, ok, message)
        if (.not. ok) error stop message
        text = chars(source)
        do
            n = len(text)
            if (n == 0) exit
            if (text(n:n) /= achar(10)) exit
            text = text(:n - 1)
        end do
        write (unit, '(a)') text
    end subroutine emit
end program gen_landau_projection
