module diag_pitch_action
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use neort_lib, only: neort_init, neort_prepare_splines, neort_setup_at_s
    use neort, only: set_to_passing_region, set_to_trapped_region, vmax_over_vth
    use neort_profiles, only: vth, Om_tE, A1, A2
    use neort_freq, only: Om_th, Om_ph, Om_tB
    use neort_transport, only: timestep_transport, Tphi_int, transport_omth => Omth, &
        transport_domthdv => dOmthdv, transport_domthdeta => dOmthdeta
    use neort_orbit, only: bounce_fast, bounce_fast_toleranced, nvar, noshear
    use driftorbit, only: mth, mph, mi, sign_vpar, nonlin, supban, nopassing, &
        comptorque, magdrift, magdrift_passing, etatp, etadt
    use do_magfie_mod, only: s, q, iota, psi_pr, sign_theta
    implicit none
    private
    public :: run_pitch_action_diag, read_pitch_points, resonance_coefficients
    public :: pitch_point_t

    type :: pitch_point_t
        real(dp) :: surface, eta, ux
        integer :: branch, harmonic
    end type pitch_point_t

contains

    pure function resonance_coefficients(harmonic, toroidal, passing, rotational, &
            electric, unit_theta, unit_drift) result(coeff)
        integer, intent(in) :: harmonic, toroidal
        logical, intent(in) :: passing
        real(dp), intent(in) :: rotational, electric, unit_theta, unit_drift
        real(dp) :: coeff(3), transit

        transit = real(harmonic, dp)
        if (passing) then
            if (rotational == 0.0_dp) error stop "zero passing rotational transform"
            transit = transit + real(toroidal, dp)/rotational
        end if
        coeff(1) = real(toroidal, dp)*unit_drift
        coeff(2) = transit*unit_theta
        coeff(3) = real(toroidal, dp)*electric
    end function resonance_coefficients

    subroutine parse_point(line, point)
        character(*), intent(in) :: line
        type(pitch_point_t), intent(out) :: point
        character(len=256) :: extra
        integer :: ios

        read (line, *, iostat=ios) point%surface, point%branch, point%harmonic, &
            point%eta, point%ux
        if (ios /= 0) error stop "pitch point requires s branch mth eta ux"
        read (line, *, iostat=ios) point%surface, point%branch, point%harmonic, &
            point%eta, point%ux, extra
        if (ios == 0) error stop "extra pitch-point column"
        if (.not. ieee_is_finite(point%surface)) error stop "nonfinite surface"
        if (.not. ieee_is_finite(point%eta)) error stop "nonfinite pitch"
        if (.not. ieee_is_finite(point%ux)) error stop "nonfinite speed"
        if (point%surface <= 0.0_dp .or. point%surface >= 1.0_dp) &
            error stop "pitch surface must lie inside (0,1)"
        if (point%eta <= 0.0_dp .or. point%ux <= 0.0_dp) &
            error stop "pitch and normalized speed must be positive"
        if (point%branch < 1 .or. point%branch > 3) &
            error stop "branch must be 1=co, 2=counter, or 3=trapped"
    end subroutine parse_point

    subroutine read_pitch_points(path, points)
        character(*), intent(in) :: path
        type(pitch_point_t), allocatable, intent(out) :: points(:)
        character(len=1024) :: line
        type(pitch_point_t) :: point
        integer :: unit, ios, count, k

        open (newunit=unit, file=path, status="old", action="read")
        count = 0
        do
            read (unit, '(A)', iostat=ios) line
            if (ios < 0) exit
            if (ios /= 0) error stop "cannot read pitch point"
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == "#") cycle
            call parse_point(line, point)
            count = count + 1
        end do
        if (count == 0) error stop "empty pitch-point file"
        rewind (unit)
        allocate (points(count))
        k = 0
        do
            read (unit, '(A)', iostat=ios) line
            if (ios < 0) exit
            if (ios /= 0) error stop "cannot reread pitch point"
            line = adjustl(line)
            if (len_trim(line) == 0) cycle
            if (line(1:1) == "#") cycle
            k = k + 1
            call parse_point(line, points(k))
        end do
        close (unit)
    end subroutine read_pitch_points

    subroutine run_pitch_action_diag(runname, point_file, tight)
        character(*), intent(in) :: runname, point_file
        logical, intent(in), optional :: tight
        type(pitch_point_t), allocatable :: points(:)
        real(dp) :: current_surface
        integer :: unit, k
        logical :: tight_mode
        character(len=256) :: output_name

        tight_mode = .false.
        if (present(tight)) tight_mode = tight
        call read_pitch_points(point_file, points)
        call neort_init(trim(runname)//".in", "in_file", "in_file_pert")
        if (nonlin) error stop "pitch action requires nonlin=false"
        if (supban) error stop "pitch action requires supban=false"
        if (.not. comptorque) error stop "pitch action requires comptorque=true"
        call neort_prepare_splines("plasma.in", "profile.in")
        if (tight_mode) then
            output_name = trim(runname)//"_pitch_action_tight.dat"
        else
            output_name = trim(runname)//"_pitch_action.dat"
        end if
        open (newunit=unit, file=trim(output_name), &
            status="replace", action="write")
        write (unit, '(A)') "# schema neort-pitch-action-v1"
        write (unit, '(A)') "# branch: 1=passing_co 2=passing_ctr 3=trapped"
        write (unit, '(A)') "# eta is inverse gauss; ux=v/vth; frequencies in s^-1"
        write (unit, '(A,L1,A,I0,A,L1)') "# magdrift=", magdrift, &
            " magdrift_passing=", magdrift_passing, " noshear=", noshear
        write (unit, '(A)') "# nonlin=false supban=false spline_init_sign=+1"
        if (tight_mode) write (unit, '(A)') &
            "# integration_rtol=1e-12 integration_atol=1e-14 diagnostic_only"
        write (unit, '(A)') "# columns: point branch mth mph istate "// &
            "s_tor eta ux vth sign_vpar eta_min eta_max umin umax "// &
            "q iota psi_pr sign_theta A1 A2 OmE Omth Omph dg_du "// &
            "dg_deta a b c g_direct g_quadratic dg_du_quadratic "// &
            "taub bounce_re bounce_im Hmn2 Tphi_int etatp etadt"
        current_surface = -1.0_dp
        do k = 1, size(points)
            if (points(k)%surface /= current_surface) then
                sign_vpar = 1.0_dp
                call neort_setup_at_s(points(k)%surface)
                if (.not. ieee_is_finite(q*iota)) error stop "nonfinite q*iota"
                if (abs(q*iota - 1.0_dp) > 1.0e-10_dp) &
                    error stop "pitch action requires native q*iota=1"
                current_surface = points(k)%surface
            end if
            call write_pitch_point(unit, k, points(k), tight_mode)
        end do
        close (unit)
    end subroutine run_pitch_action_diag

    subroutine write_pitch_point(unit, index, point, tight)
        integer, intent(in) :: unit, index
        type(pitch_point_t), intent(in) :: point
        logical, intent(in) :: tight
        real(dp) :: eta_min, eta_max, unit_theta, unit_drift, d1, d2
        real(dp) :: omph, domphdv, domphdeta, v, coeff(3), g, dgdu, dgdeta
        real(dp) :: taub, bounceavg(nvar), hmn2, weight, values(33)
        integer :: istate
        logical :: passing

        passing = point%branch /= 3
        if (passing) then
            if (nopassing) error stop "passing point forbidden by nopassing"
            call set_to_passing_region(eta_min, eta_max)
        else
            call set_to_trapped_region(eta_min, eta_max)
        end if
        sign_vpar = 1.0_dp
        if (point%branch == 2) sign_vpar = -1.0_dp
        mth = point%harmonic
        if (point%eta < eta_min .or. point%eta > eta_max) &
            error stop "point outside native pitch support"
        if (point%ux < 1.0e-6_dp .or. point%ux > vmax_over_vth) &
            error stop "point outside native velocity support"
        call Om_th(vth, point%eta, unit_theta, d1, d2)
        unit_drift = 0.0_dp
        if (magdrift) call Om_tB(vth, point%eta, unit_drift, d1, d2)
        coeff = resonance_coefficients(mth, mph, passing, iota, Om_tE, &
            unit_theta, unit_drift)
        v = point%ux*vth
        call Om_th(v, point%eta, transport_omth, transport_domthdv, &
            transport_domthdeta)
        call Om_ph(v, point%eta, omph, domphdv, domphdeta)
        if (transport_omth == 0.0_dp) error stop "zero orbit frequency"
        g = real(mth, dp)*transport_omth + real(mph, dp)*omph
        dgdu = vth*(real(mth, dp)*transport_domthdv + &
            real(mph, dp)*domphdv)
        dgdeta = real(mth, dp)*transport_domthdeta + real(mph, dp)*domphdeta
        taub = 2.0_dp*acos(-1.0_dp)/abs(transport_omth)
        if (tight) then
            call bounce_fast_toleranced(v, point%eta, taub, bounceavg, &
                timestep_transport, istate, 1.0e-12_dp, 1.0e-14_dp)
        else
            call bounce_fast(v, point%eta, taub, bounceavg, timestep_transport, istate)
        end if
        if (istate /= 2) error stop "pitch action bounce solver failed"
        hmn2 = (bounceavg(3)**2 + bounceavg(4)**2)*(mi*v**2/2.0_dp)**2
        weight = Tphi_int(point%ux, taub, hmn2)
        values(1) = s
        values(2) = point%eta
        values(3) = point%ux
        values(4) = vth
        values(5) = sign_vpar
        values(6) = eta_min
        values(7) = eta_max
        values(8) = 1.0e-6_dp
        values(9) = vmax_over_vth
        values(10) = q
        values(11) = iota
        values(12) = psi_pr
        values(13) = sign_theta
        values(14) = A1
        values(15) = A2
        values(16) = Om_tE
        values(17) = transport_omth
        values(18) = omph
        values(19) = dgdu
        values(20) = dgdeta
        values(21:23) = coeff
        values(24) = g
        values(25) = (coeff(1)*point%ux + coeff(2))*point%ux + coeff(3)
        values(26) = 2.0_dp*coeff(1)*point%ux + coeff(2)
        values(27) = taub
        values(28) = bounceavg(3)
        values(29) = bounceavg(4)
        values(30) = hmn2
        values(31) = weight
        values(32) = etatp
        values(33) = etadt
        if (.not. all(ieee_is_finite(values))) error stop "nonfinite pitch action"
        write (unit, '(5(I0,1X),*(ES24.16,1X))') index, point%branch, mth, mph, &
            istate, values
    end subroutine write_pitch_point

end module diag_pitch_action
