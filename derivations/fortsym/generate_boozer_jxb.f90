program generate_boozer_jxb
    ! Generate the metric raise and phase-averaged toroidal Lorentz-torque
    ! density for coordinates ordered (s, phi, theta).
    !
    ! Let K^i = Jcal J^i and C^i = Jcal B^i.  The volume-weighted covariant
    ! toroidal force is
    !
    !   Jcal (J cross B)_phi = K^theta C^s - K^s C^theta.
    !
    ! The input magnetic perturbation is covariant.  Since det(g)=Jcal**2,
    !
    !   C^i = Jcal g^ij B_j = cof(g)_ij B_j / Jcal.
    !
    ! The generated real kernel performs that raise and the 1/2 Re phasor
    ! contraction.  Angular integration is deliberately left to NEO-RT.
    use, intrinsic :: iso_fortran_env, only: output_unit
    use fortsym_arena, only: arena_t
    use fortsym_engine, only: engine_result_t
    use fortsym_engine_symengine, only: symengine_engine_t, make_symengine_engine
    use fortsym_expr, only: expr_t, num, sym, &
        operator(+), operator(-), operator(*), operator(/)
    use fortsym_kernel, only: KERNEL_SUBROUTINE, emit_kernel, kernel_spec_t
    use fortsym_string, only: chars, str
    implicit none

    character(*), parameter :: DEFAULT_OUTPUT = &
        "../../src/generated/boozer_jxb_kernel.f90"
    character(*), parameter :: FORTSYM_REVISION = &
        "fortsym@0cb719ca71c5cb840cd9a83abb8d05c23bc6bbd7"
    type(arena_t), target :: arena
    type(symengine_engine_t) :: engine
    type(expr_t) :: inv_jacobian
    type(expr_t) :: gss, gsp, gst, gpp, gpt, gtt
    type(expr_t) :: bs_re, bs_im, bp_re, bp_im, bt_re, bt_im
    type(expr_t) :: ks_re, ks_im, kt_re, kt_im
    type(expr_t) :: cofactor_ss, cofactor_sp, cofactor_st
    type(expr_t) :: cofactor_pt, cofactor_tt
    type(expr_t) :: cs_re, cs_im, ct_re, ct_im, torque
    type(expr_t) :: roots(5)
    type(kernel_spec_t) :: spec
    character(:), allocatable :: generated
    character(1024) :: output
    integer :: k, status, unit

    output = DEFAULT_OUTPUT
    if (command_argument_count() > 0) call get_command_argument(1, output)
    call arena%init()
    engine = make_symengine_engine(arena)

    inv_jacobian = sym(arena, "inv_jacobian")
    gss = sym(arena, "g_ss")
    gsp = sym(arena, "g_sphi")
    gst = sym(arena, "g_stheta")
    gpp = sym(arena, "g_phiphi")
    gpt = sym(arena, "g_phitheta")
    gtt = sym(arena, "g_thetatheta")
    bs_re = sym(arena, "bs_re")
    bs_im = sym(arena, "bs_im")
    bp_re = sym(arena, "bphi_re")
    bp_im = sym(arena, "bphi_im")
    bt_re = sym(arena, "btheta_re")
    bt_im = sym(arena, "btheta_im")
    ks_re = sym(arena, "ks_re")
    ks_im = sym(arena, "ks_im")
    kt_re = sym(arena, "ktheta_re")
    kt_im = sym(arena, "ktheta_im")

    cofactor_ss = gpp*gtt - gpt*gpt
    cofactor_sp = gst*gpt - gsp*gtt
    cofactor_st = gsp*gpt - gst*gpp
    cofactor_pt = gsp*gst - gss*gpt
    cofactor_tt = gss*gpp - gsp*gsp

    cs_re = inv_jacobian*( &
        cofactor_ss*bs_re + cofactor_sp*bp_re + cofactor_st*bt_re)
    cs_im = inv_jacobian*( &
        cofactor_ss*bs_im + cofactor_sp*bp_im + cofactor_st*bt_im)
    ct_re = inv_jacobian*( &
        cofactor_st*bs_re + cofactor_pt*bp_re + cofactor_tt*bt_re)
    ct_im = inv_jacobian*( &
        cofactor_st*bs_im + cofactor_pt*bp_im + cofactor_tt*bt_im)
    torque = ( &
        kt_re*cs_re + kt_im*cs_im - ks_re*ct_re - ks_im*ct_im)/ &
        num(arena, 2)

    roots = [cs_re, cs_im, ct_re, ct_im, torque]
    do k = 1, size(roots)
        call simplify(roots(k))
    end do

    spec%name = str("boozer_jxb_kernel")
    spec%module_name = str("generated_boozer_jxb_kernel")
    spec%mode = KERNEL_SUBROUTINE
    spec%pure_procedure = .true.
    spec%temp_prefix = str("t")
    spec%generator = str("derivations/fortsym/generate_boozer_jxb.f90")
    spec%generator_revision = str(FORTSYM_REVISION)
    spec%regenerate_command = str( &
        "cd derivations/fortsym && fpm run --profile release generate_boozer_jxb")
    allocate (spec%args(17), spec%outputs(5))
    spec%args = [ &
        str("inv_jacobian"), str("g_ss"), str("g_sphi"), str("g_stheta"), &
        str("g_phiphi"), str("g_phitheta"), str("g_thetatheta"), &
        str("bs_re"), str("bs_im"), str("bphi_re"), str("bphi_im"), &
        str("btheta_re"), str("btheta_im"), str("ks_re"), str("ks_im"), &
        str("ktheta_re"), str("ktheta_im")]
    spec%outputs = [ &
        str("cs_re"), str("cs_im"), str("ctheta_re"), str("ctheta_im"), &
        str("torque_density")]

    open (newunit=unit, file=trim(output), status="replace", action="write", iostat=status)
    if (status /= 0) error stop "cannot open generated Boozer jxb kernel"
    generated = chars(emit_kernel(roots, spec))
    do while (len(generated) > 0)
        if (generated(len(generated):) /= new_line("a")) exit
        generated = generated(:len(generated) - 1)
    end do
    write (unit, "(a)") generated
    close (unit)
    write (output_unit, "(a)") "generated "//trim(output)

contains

    subroutine simplify(expression)
        type(expr_t), intent(inout) :: expression
        type(engine_result_t) :: result

        result = engine%simplify(expression)
        if (.not. result%ok) error stop "fortsym failed to simplify Boozer jxb kernel"
        expression = result%value
    end subroutine simplify

end program generate_boozer_jxb
