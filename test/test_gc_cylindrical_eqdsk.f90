program test_gc_cylindrical_eqdsk
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: inp_swi, read_boozer_file, R0
    use field_eq_mod, only: hrad, hzet, rad, zet
    use neort_gc_cylindrical_model, only: GC_CYL_SUCCESS, &
        gc_cylindrical_field_sample_t
    use neort_gc_eqdsk_cylindrical_adapter, only: eqdsk_cylindrical_field_t
    use neort_gc_eqdsk_cut_jet, only: EQDSK_CUT_JET_SUCCESS, &
        eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
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
    call test_cut_jet_adapter(field)
    call pass_test

contains

    subroutine test_cut_jet_adapter(eqdsk_field)
        type(eqdsk_cylindrical_field_t), intent(inout) :: eqdsk_field
        type(eqdsk_cut_jet_t) :: cut, reversed
        type(gc_cylindrical_field_sample_t) :: local_sample
        real(dp) :: position(3), velocity(3), plus_position(3), minus_position(3)
        real(dp) :: direct_value, plus_value, minus_value
        real(dp) :: finite_d_R, finite_d_Z, finite_rate, step
        integer :: i_R, i_Z, local_status

        if (size(rad) < 4 .or. size(zet) < 4) then
            error stop 'EQDSK cut-jet fixture grid is too small'
        end if
        i_R = max(2, min(size(rad)-2, size(rad)/2))
        i_Z = max(2, min(size(zet)-2, 3*size(zet)/5))
        position = [0.5_dp*(rad(i_R)+rad(i_R+1)), &
            0.5_dp*(zet(i_Z)+zet(i_Z+1)), 0.37_dp]
        velocity = [0.13_dp, -0.17_dp, 0.09_dp]
        eqdsk_field%field_scale = 1.7_dp

        call eqdsk_field%evaluate(position, local_sample, local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            error stop 'scaled direct EQDSK sample failed'
        end if
        call evaluate_eqdsk_cut_jet(position, eqdsk_field%field_scale, 1, &
            velocity, cut, local_status)
        if (local_status /= EQDSK_CUT_JET_SUCCESS) then
            error stop 'production EQDSK cut-jet adapter failed'
        end if
        call direct_cut_value(eqdsk_field, position, direct_value)
        call require_close(cut%cut_value, direct_value, 2.0e-12_dp, &
            'EQDSK cut value versus direct Eq. 13 contraction')
        call require_close(eqdsk_field%field_scale*cut%psi_jet(1), &
            local_sample%psi, 2.0e-12_dp, 'EQDSK cut psi cell selection')
        call require_close(eqdsk_field%field_scale*cut%f_jet(1), &
            position(1)*local_sample%b(2), 2.0e-12_dp, &
            'EQDSK cut F-profile selection')

        step = 1.0e-4_dp*min(hrad, hzet)
        plus_position = position
        minus_position = position
        plus_position(1) = position(1)+step
        minus_position(1) = position(1)-step
        call direct_cut_value(eqdsk_field, plus_position, plus_value)
        call direct_cut_value(eqdsk_field, minus_position, minus_value)
        finite_d_R = (plus_value-minus_value)/(2.0_dp*step)
        plus_position = position
        minus_position = position
        plus_position(2) = position(2)+step
        minus_position(2) = position(2)-step
        call direct_cut_value(eqdsk_field, plus_position, plus_value)
        call direct_cut_value(eqdsk_field, minus_position, minus_value)
        finite_d_Z = (plus_value-minus_value)/(2.0_dp*step)
        finite_rate = finite_d_R*velocity(1)+finite_d_Z*velocity(2)
        call require_close(cut%d_cut_d_R, finite_d_R, 2.0e-7_dp, &
            'EQDSK cut R derivative')
        call require_close(cut%d_cut_d_Z, finite_d_Z, 2.0e-7_dp, &
            'EQDSK cut Z derivative')
        call require_close(cut%d_cut_d_arc_phi, 0.0_dp, 0.0_dp, &
            'axisymmetric EQDSK cut phi derivative')
        call require_close(cut%cut_rate, finite_rate, 2.0e-7_dp, &
            'EQDSK cut directional rate')

        call evaluate_eqdsk_cut_jet(position, eqdsk_field%field_scale, -1, &
            velocity, reversed, local_status)
        if (local_status /= EQDSK_CUT_JET_SUCCESS) then
            error stop 'reversed EQDSK cut-jet adapter failed'
        end if
        call require_close(reversed%cut_value, -cut%cut_value, 2.0e-12_dp, &
            'EQDSK cut orientation reversal')
        call require_close(reversed%cut_rate, -cut%cut_rate, 2.0e-12_dp, &
            'EQDSK cut-rate orientation reversal')
        call require_close(reversed%absolute_cut_rate, cut%absolute_cut_rate, &
            2.0e-12_dp, 'EQDSK absolute cut-rate orientation invariance')
    end subroutine test_cut_jet_adapter

    subroutine direct_cut_value(eqdsk_field, position, value)
        type(eqdsk_cylindrical_field_t), intent(in) :: eqdsk_field
        real(dp), intent(in) :: position(3)
        real(dp), intent(out) :: value
        type(gc_cylindrical_field_sample_t) :: local_sample
        integer :: local_status

        call eqdsk_field%evaluate(position, local_sample, local_status)
        if (local_status /= GC_CYL_SUCCESS) then
            error stop 'direct Eq. 13 field sample failed'
        end if
        value = (local_sample%grad_b(3)*local_sample%grad_psi(1) - &
            local_sample%grad_b(1)*local_sample%grad_psi(3))/position(1)
    end subroutine direct_cut_value

    subroutine require_close(actual, expected, relative_tolerance, label)
        real(dp), intent(in) :: actual, expected, relative_tolerance
        character(*), intent(in) :: label

        if (abs(actual-expected) > relative_tolerance* &
                max(1.0_dp, abs(expected))) then
            write (*, '(a,2(1x,es24.16))') trim(label), actual, expected
            error stop 'direct EQDSK cut-jet oracle failed'
        end if
    end subroutine require_close
end program test_gc_cylindrical_eqdsk
