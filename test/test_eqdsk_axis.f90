program test_eqdsk_axis
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, &
        init_magfie_at_s, do_magfie, q
    use do_magfie_pert_mod, only: inp_swi_pert, read_boozer_pert_file, &
        init_magfie_pert_at_s, do_magfie_pert_amp, set_mph
    use util, only: pi
    use util_for_test, only: pass_test, fail_test

    implicit none

    integer, parameter :: ntheta = 4096
    character(len=1024) :: eqdsk_file, perturbation_file
    real(dp) :: x(3), bmod, sqrtg, bder(3), hcov(3), hcon(3), hcurl(3)
    real(dp) :: theta, dtheta, q_integral
    complex(dp) :: bamp
    integer :: k

    call get_environment_variable('EQDSK_FILE', eqdsk_file)
    call get_environment_variable('RZ_PERT_FILE', perturbation_file)
    if (len_trim(eqdsk_file) == 0 .or. len_trim(perturbation_file) == 0) then
        call fail_test
    end if

    inp_swi = 11
    call read_boozer_file(trim(eqdsk_file))
    call set_s(0.25_dp)
    call init_magfie_at_s()

    dtheta = 2.0_dp*pi/real(ntheta, dp)
    q_integral = 0.0_dp
    do k = 0, ntheta - 1
        theta = -pi + (real(k, dp) + 0.5_dp)*dtheta
        x = [0.25_dp, 0.0_dp, theta]
        call do_magfie(x, bmod, sqrtg, bder, hcov, hcon, hcurl)
        if (.not. all(ieee_is_finite([bmod, sqrtg, bder, hcov, hcon, hcurl]))) then
            call fail_test
        end if
        if (bmod <= 0.0_dp .or. abs(sqrtg) <= tiny(sqrtg) .or. &
            abs(hcon(3)) <= tiny(hcon(3))) then
            call fail_test
        end if
        q_integral = q_integral + hcon(2)/hcon(3)*dtheta/(2.0_dp*pi)
    end do
    if (q*q_integral <= 0.0_dp .or. abs(q_integral/q) < 0.5_dp .or. &
        abs(q_integral/q) > 2.0_dp) then
        write(*,*) 'GEQDSK q and integrated field-line pitch are inconsistent:', q, q_integral
        call fail_test
    end if

    inp_swi_pert = 11
    call set_mph(-3)
    call read_boozer_pert_file(trim(perturbation_file))
    call init_magfie_pert_at_s()
    x = [0.25_dp, 0.0_dp, 0.0_dp]
    call do_magfie_pert_amp(x, bamp)
    if (.not. ieee_is_finite(real(bamp)) .or. &
        .not. ieee_is_finite(aimag(bamp))) call fail_test

    call pass_test
end program test_eqdsk_axis
