program test_gc_cylindrical_eqdsk
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: inp_swi, read_boozer_file, R0, a
    use neort_gc_cylindrical_model, only: GC_CYL_SUCCESS, &
        gc_cylindrical_field_sample_t
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        configure_eqdsk_cylindrical_field, eqdsk_cylindrical_field_t
    use util_for_test, only: pass_test

    implicit none

    character(len=1024) :: eqdsk_file
    type(eqdsk_cylindrical_field_t) :: field
    type(gc_cylindrical_field_sample_t) :: sample
    integer :: status

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    if (len_trim(eqdsk_file) == 0) error stop 'EQDSK_FILE is required'
    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))
    call configure_eqdsk_cylindrical_field(1.0_dp, field, status)
    if (status /= GC_CYL_SUCCESS) error stop 'direct EQDSK field setup failed'
    call field%evaluate([R0, 0.0_dp, 0.0_dp], sample, status)
    if (status /= GC_CYL_SUCCESS) error stop 'direct EQDSK axis sample failed'
    if (.not. all(ieee_is_finite([sample%b, sample%grad_b, &
        sample%curl_bhat, sample%psi, sample%grad_psi]))) then
        error stop 'direct EQDSK cylindrical sample is non-finite'
    end if
    if (sample%bmod <= 0.0_dp) error stop 'direct EQDSK axis field is zero'
    ! A full guiding-centre orbit is not clipped at the LCFS.  The point is
    ! outside the nominal R0 +/- a plasma radius but inside the EQDSK box.
    call field%evaluate([R0 + 1.1_dp*a, 0.0_dp, 0.0_dp], &
        sample, status)
    if (status /= GC_CYL_SUCCESS) then
        error stop 'direct EQDSK finite-width sample was clipped at the LCFS'
    end if
    call field%evaluate([field%domain_R_max + 1.0e-6_dp, 0.0_dp, 0.0_dp], &
        sample, status)
    if (status == GC_CYL_SUCCESS) then
        error stop 'direct EQDSK box boundary was not enforced'
    end if
    call pass_test
end program test_gc_cylindrical_eqdsk
