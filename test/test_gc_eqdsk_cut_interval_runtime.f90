program test_gc_eqdsk_cut_interval_runtime
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use field_eq_mod, only: nrad, nzet, psi_sep, rad, zet
    use neort_gc_eqdsk_cut_interval, only: &
        EQDSK_CUT_INTERVAL_CELL_MISMATCH, EQDSK_CUT_INTERVAL_SUCCESS, &
        eqdsk_cut_interval_result_t, evaluate_eqdsk_cut_interval_box
    use neort_gc_eqdsk_cut_jet, only: EQDSK_CUT_JET_SUCCESS, &
        eqdsk_cut_jet_t, evaluate_eqdsk_cut_jet
    use neort_eqdsk_cut_r_flux_chart_symbolic, only: &
        evaluate_neort_eqdsk_cut_r_flux_chart
    use neort_gc_eqdsk_cylindrical_adapter, only: &
        eqdsk_cylindrical_field_t, initialize_eqdsk_cylindrical_field
    implicit none

    type(eqdsk_cylindrical_field_t) :: field
    type(eqdsk_cut_interval_result_t) :: box
    type(eqdsk_cut_jet_t) :: point
    character(len=1024) :: path
    real(dp) :: R_value, Z_value, position(3), denominator
    real(dp) :: point_dZ_dR, point_dpsihat_dR
    integer :: status, point_status, cell_R, cell_Z, zero_Z, i, j

    call get_environment_variable('EQDSK_FILE', path)
    call require(len_trim(path) > 0, 'EQDSK_FILE is required')
    call initialize_eqdsk_cylindrical_field(trim(path), 1.0_dp, field, status)
    call require(status == 0, 'failed to initialize circular EQDSK')
    call require(nrad >= 3 .and. nzet >= 3, 'EQDSK grid is too small')

    zero_Z = 0
    do cell_Z = 1, nzet
        if (zet(cell_Z) == 0.0_dp) then
            zero_Z = cell_Z
            exit
        end if
    end do
    call require(zero_Z > 0, 'circular EQDSK has no exact midplane row')
    cell_R = max(1, min(nrad-1, nrad/2))
    cell_Z = max(1, min(nzet-1, zero_Z))
    call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, rad(cell_R), &
        rad(cell_R+1), 0.0_dp, 0.0_dp, box, status)
    call require(status == EQDSK_CUT_INTERVAL_SUCCESS, &
        'midplane interval evaluation failed')
    call require(box%numerator%lo <= 0.0_dp .and. box%numerator%hi >= 0.0_dp, &
        'exact symmetric midplane was excluded from Eq.13')
    call require(box%denominator_positive_certified, &
        'positive Eq.13 denominator was not certified')

    call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, rad(cell_R), &
        rad(cell_R+1), zet(cell_Z), zet(cell_Z+1), box, status)
    call require(status == EQDSK_CUT_INTERVAL_SUCCESS, &
        'nondegenerate interval evaluation failed')
    do i = 0, 2
        R_value = rad(cell_R)+real(i,dp)*(rad(cell_R+1)-rad(cell_R))/2.0_dp
        do j = 0, 2
            Z_value = zet(cell_Z)+real(j,dp)*(zet(cell_Z+1)-zet(cell_Z))/2.0_dp
            position = [R_value, Z_value, 0.0_dp]
            call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, &
                [0.0_dp,0.0_dp,0.0_dp], point, point_status)
            call require(point_status == EQDSK_CUT_JET_SUCCESS, &
                'sampled scalar cut evaluation failed')
            call require(encloses(box%numerator, point%cut_numerator), &
                'sampled numerator is outside tightened interval result')
        end do
    end do

    ! Regression box retained from the first global-atlas failure.  Its scalar
    ! midpoint is an independent oracle for interval tightening diagnostics.
    cell_R = 12
    cell_Z = 26
    call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, &
        1.0997585571929329e2_dp, 1.0997586161389307e2_dp, &
        -1.5955200195312504e1_dp, -1.5954589843750004e1_dp, box, status)
    call require(status == EQDSK_CUT_INTERVAL_SUCCESS, &
        'atlas regression-box interval evaluation failed')
    R_value = 0.5_dp*(1.0997585571929329e2_dp+1.0997586161389307e2_dp)
    Z_value = 0.5_dp*(-1.5955200195312504e1_dp-1.5954589843750004e1_dp)
    position = [R_value, Z_value, 0.0_dp]
    call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, &
        [0.0_dp,0.0_dp,0.0_dp], point, point_status)
    call require(point_status == EQDSK_CUT_JET_SUCCESS, &
        'atlas regression-box scalar evaluation failed')
    call require(encloses(box%numerator, point%cut_numerator), &
        'atlas regression-box midpoint escaped its enclosure')
    write (*, '(a,2(1x,es24.16))') 'regression N interval', &
        box%numerator%lo, box%numerator%hi
    write (*, '(a,3(1x,es24.16))') 'regression psi_hat interval/scalar', &
        box%psi_hat%lo, box%psi_hat%hi, point%psi_jet(1)/psi_sep
    write (*, '(a,3(1x,es24.16))') 'regression scalar N/NR/NZ', &
        point%cut_numerator, point%d_cut_numerator_d_R, &
        point%d_cut_numerator_d_Z
    write (*, '(a,2(1x,es24.16))') 'regression NR interval', &
        box%numerator_R%lo, box%numerator_R%hi
    write (*, '(a,2(1x,es24.16))') 'regression NZ interval', &
        box%numerator_Z%lo, box%numerator_Z%hi
    do j = 0, 8
        Z_value = -1.5955200195312504e1_dp+real(j,dp) &
            *(-1.5953979492187504e1_dp+1.5955200195312504e1_dp)/8.0_dp
        position = [R_value, Z_value, 0.0_dp]
        call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, &
            [0.0_dp,0.0_dp,0.0_dp], point, point_status)
        call require(point_status == EQDSK_CUT_JET_SUCCESS, &
            'atlas regression scalar branch evaluation failed')
        write (*, '(a,2(1x,es24.16))') 'regression scalar Z/N', &
            Z_value, point%cut_numerator
    end do

    cell_R = max(1, min(nrad-1, nrad/2))
    cell_Z = max(1, min(nzet-1, zero_Z))
    Z_value = 0.5_dp*(zet(cell_Z)+zet(cell_Z+1))
    R_value = 0.5_dp*(rad(cell_R)+rad(cell_R+1))
    call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, R_value, R_value, &
        Z_value, Z_value, box, status)
    call require(status == EQDSK_CUT_INTERVAL_SUCCESS, &
        'point interval evaluation failed')
    position = [R_value, Z_value, 0.0_dp]
    call evaluate_eqdsk_cut_jet(position, 1.0_dp, 1, [0.0_dp,0.0_dp,0.0_dp], &
        point, point_status)
    call require(point_status == EQDSK_CUT_JET_SUCCESS, &
        'independent scalar cut evaluation failed')
    denominator = point%f_jet(1)**2+point%psi_jet(2)**2+point%psi_jet(3)**2
    call require(encloses(box%numerator, point%cut_numerator), &
        'point numerator is outside interval kernel result')
    call require(encloses(box%numerator_R, point%d_cut_numerator_d_R), &
        'point radial numerator derivative is outside interval result')
    call require(encloses(box%numerator_Z, point%d_cut_numerator_d_Z), &
        'point vertical numerator derivative is outside interval result')
    call require(encloses(box%positive_denominator, denominator), &
        'point positive denominator is outside interval result')
    call require(box%r_chart_certified, &
        'regular point did not produce an interval R chart')
    call evaluate_neort_eqdsk_cut_r_flux_chart( &
        point%d_cut_numerator_d_R, point%d_cut_numerator_d_Z, &
        point%psi_jet(2), point%psi_jet(3), &
        psi_sep, point_dZ_dR, point_dpsihat_dR)
    call require(encloses(box%dZ_dR, point_dZ_dR), &
        'point cut slope is outside interval chart result')
    call require(encloses(box%dpsihat_dR, point_dpsihat_dR), &
        'point flux derivative is outside interval chart result')

    call evaluate_eqdsk_cut_interval_box(cell_R, cell_Z, rad(cell_R), &
        rad(cell_R+2), zet(cell_Z), zet(cell_Z+1), box, status)
    call require(status == EQDSK_CUT_INTERVAL_CELL_MISMATCH, &
        'cross-cell interval was accepted by one-cell evaluator')

    write (*, '(a)') 'test_gc_eqdsk_cut_interval_runtime OK'

contains

    logical function encloses(interval, value)
        use neort_gc_outward_interval, only: gc_outward_interval_t
        type(gc_outward_interval_t), intent(in) :: interval
        real(dp), intent(in) :: value

        encloses = value >= interval%lo .and. value <= interval%hi
    end function encloses

    subroutine require(condition, message)
        logical, intent(in) :: condition
        character(*), intent(in) :: message

        if (.not. condition) error stop message
    end subroutine require

end program test_gc_eqdsk_cut_interval_runtime
