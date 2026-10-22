program gen_native_harmonic
    use fortsym_string, only: str, chars
    use fortsym_arena, only: arena_t
    use fortsym_expr, only: expr_t, sym, operator(+), operator(*), cos, sin
    use fortsym_diff, only: diff
    use fortsym_kernel, only: kernel_spec_t, emit_kernel, KERNEL_SUBROUTINE
    implicit none
    type(arena_t), target :: arena
    type(expr_t) :: x, theta, m, a(4), b(4), pa, pb, f, roots(3)
    type(kernel_spec_t) :: spec
    character(2) :: name
    integer :: k, unit

    call arena%init()
    x=sym(arena,'x'); theta=sym(arena,'theta'); m=sym(arena,'m')
    do k=1,4
        write(name,'(a,i1)') 'a',k
        a(k)=sym(arena,name)
        write(name,'(a,i1)') 'b',k
        b(k)=sym(arena,name)
    end do
    pa=((a(1)*x+a(2))*x+a(3))*x+a(4)
    pb=((b(1)*x+b(2))*x+b(3))*x+b(4)
    f=pa*cos(m*theta)+pb*sin(m*theta)
    roots=[f,diff(f,x),diff(f,theta)]
    spec%name=str('native_harmonic')
    spec%mode=KERNEL_SUBROUTINE
    spec%generator=str('POTATO/gen/gen_native_harmonic.f90')
    allocate(spec%args(11),spec%outputs(3))
    spec%args(1)=str('x'); spec%args(2)=str('theta'); spec%args(3)=str('m')
    do k=1,4
        write(name,'(a,i1)') 'a',k
        spec%args(3+k)=str(name)
        write(name,'(a,i1)') 'b',k
        spec%args(7+k)=str(name)
    end do
    spec%outputs(1)=str('v'); spec%outputs(2)=str('dr'); spec%outputs(3)=str('dt')
    open(newunit=unit,file='POTATO/SRC/generated/native_harmonic.inc',status='replace')
    write(unit,'(a)') chars(emit_kernel(roots,spec))
    close(unit)
end program gen_native_harmonic
