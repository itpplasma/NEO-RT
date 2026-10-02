program test_boozer_angle_map
    ! Boozer angles on the direct GEQDSK chart (inp_swi=11), checked against
    ! the field line of the direct field: phi_B - q*theta_B is constant along
    ! a field line, with phi integrated as the local pitch over theta_geo.
    ! The map is built from the Boozer file's own Fourier series of R, Z and
    ! nu; the field line comes from libneo's cylindrical field, so neither
    ! side reproduces the other.  The oracle must also reject a flipped sign
    ! of phi_B - phi and theta_B = theta_geo (both checked below).
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, bfac, read_boozer_file
    use do_magfie_pert_mod, only: inp_swi_pert, pert_angle_map, read_boozer_pert_file
    use neort_eqdsk_field, only: eqdsk_local_pitch
    use neort_boozer_angle_map, only: boozer_angles
    use util, only: pi

    implicit none

    integer, parameter :: nth = 2048
    real(dp), parameter :: s_test(3) = [0.15_dp, 0.4_dp, 0.75_dp]
    character(len=1024) :: geqdsk, boozer, pert
    real(dp) :: residual, swing, worst, worst_wrong_sign, worst_no_map, winding
    character(len=32) :: orientation
    integer :: is

    call get_environment_variable("EQDSK_SOLOVEV_FILE", geqdsk)
    call get_environment_variable("BOOZER_SOLOVEV_FILE", boozer)
    call get_environment_variable("PERT_SOLOVEV_FILE", pert)
    call get_environment_variable("BOOZER_WINDING", orientation)
    winding = 1.0_dp
    if (trim(orientation) == "-1") winding = -1.0_dp
    if (len_trim(pert) == 0) error stop "EQDSK/BOOZER/PERT_SOLOVEV_FILE must be set"
    inp_swi = 11
    inp_swi_pert = 9
    bfac = 1.0_dp
    pert_angle_map = boozer
    call read_boozer_file(trim(geqdsk))
    call read_boozer_pert_file(trim(pert))

    worst = 0.0_dp
    worst_wrong_sign = 0.0_dp
    worst_no_map = 0.0_dp
    do is = 1, size(s_test)
        call label_residual(s_test(is), 1.0_dp, .true., residual, swing)
        print '(a,f5.2,a,es10.3,a,es10.3)', "s =", s_test(is), &
            "  residual [rad]:", residual, "  max|phi_B - phi|:", swing
        worst = max(worst, residual)
        call label_residual(s_test(is), -1.0_dp, .true., residual, swing)
        worst_wrong_sign = max(worst_wrong_sign, residual)
        call label_residual(s_test(is), 1.0_dp, .false., residual, swing)
        worst_no_map = max(worst_no_map, residual)
    end do
    print '(a,es10.3)', "residual with phi_B - phi sign flipped:", worst_wrong_sign
    print '(a,es10.3)', "residual with theta_B = theta_geo:     ", worst_no_map
    if (worst > 2.0e-4_dp) error stop "field-line label not constant"
    if (worst_wrong_sign < 50.0_dp*2.0e-4_dp) error stop "oracle blind to phi sign"
    if (worst_no_map < 50.0_dp*2.0e-4_dp) error stop "oracle blind to theta_B"
    print *, "PASS test_boozer_angle_map"

contains

    subroutine label_residual(s, dphi_sign, use_map, residual, swing)
        ! max |phi + dphi - q_line*theta_B - const| over one poloidal turn,
        ! q_line being the field-line q of the direct field on this surface.
        real(dp), intent(in) :: s, dphi_sign
        logical, intent(in) :: use_map
        real(dp), intent(out) :: residual, swing
        real(dp) :: theta(0:nth), phi(0:nth), pitch(0:nth), label(0:nth)
        real(dp) :: s_b, theta_b, dphi, q_line
        integer :: i

        do i = 0, nth
            theta(i) = 2.0_dp*pi*i/nth
            pitch(i) = eqdsk_local_pitch(s, theta(i))
        end do
        phi(0) = 0.0_dp
        do i = 1, nth
            phi(i) = phi(i - 1) + 0.5_dp*(pitch(i) + pitch(i - 1))*(theta(i) - theta(i - 1))
        end do
        q_line = phi(nth)/(2.0_dp*pi)
        swing = 0.0_dp
        do i = 0, nth
            call boozer_angles(s, theta(i), s_b, theta_b, dphi)
            if (.not. use_map) then
                theta_b = theta(i)
                dphi = 0.0_dp
            end if
            label(i) = phi(i) + dphi_sign*dphi - winding*q_line*theta_b
            swing = max(swing, abs(dphi))
        end do
        residual = maxval(abs(label - label(0)))
    end subroutine label_residual

end program test_boozer_angle_map
