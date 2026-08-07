program test_gc_cylindrical_eqdsk
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: inp_swi, read_boozer_file, R0
    use neort_gc_cylindrical_model, only: GC_CYL_SUCCESS, &
        gc_cylindrical_field_sample_t
    use neort_gc_eqdsk_cylindrical_adapter, only: eqdsk_cylindrical_field_t
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
    field%field_scale = 1.0_dp
    call field%evaluate([R0, 0.0_dp, 0.0_dp], sample, status)
    if (status /= GC_CYL_SUCCESS) error stop 'direct EQDSK axis sample failed'
    if (.not. all(ieee_is_finite([sample%b, sample%grad_b, &
        sample%curl_bhat, sample%psi, sample%grad_psi]))) then
        error stop 'direct EQDSK cylindrical sample is non-finite'
    end if
    if (sample%bmod <= 0.0_dp) error stop 'direct EQDSK axis field is zero'
    call pass_test
end program test_gc_cylindrical_eqdsk
