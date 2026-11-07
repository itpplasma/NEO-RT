program test_trapped_precession_sign
    ! Usage: test_trapped_precession_sign.x <orientation> <iota_sign> <field_sign>
    !
    ! Oracle from cylindrical geometry, independent of NEO-RT's chart algebra.
    ! Large-aspect-ratio circular tokamak with R = R0 + r cos(theta),
    ! Z = orientation*r sin(theta), |B| = B0 (1 - (r/R0) cos(theta)), toroidal
    ! field field_sign*B0 R0/R along the right-handed cylindrical azimuth and
    ! chart rotational transform iota_sign*|iota|. A deeply trapped ion
    ! (v_par = 0) at theta = 0 drifts vertically with
    !   v_Z = field_sign * c mu B0/(qi B R0),  mu = m v**2/(2 B),
    ! and grad(theta) = orientation*e_Z/r there, so the field-line label
    ! alpha = phi - q theta precesses at
    !   Omega = -q * orientation * v_Z / r.
    ! NEO-RT's bounce-averaged Om_tB must have this sign and magnitude.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: s, R0, init_magfie_at_s, set_s
    use do_magfie_pert_mod, only: init_mph_from_shared
    use neort_lib, only: neort_init
    use neort, only: set_to_trapped_region
    use neort_freq, only: Om_tB, init_canon_freq_trapped_spline
    use neort_magfie, only: init_flux_surface_average
    use neort_profiles, only: init_profiles
    use driftorbit, only: vth, etadt, etatp, etamin, etamax, sign_vpar
    use util, only: c, mi, qi, pi

    implicit none

    real(dp), parameter :: field_t = 2.0_dp
    real(dp), parameter :: major_radius_m = 50.0_dp
    real(dp), parameter :: minor_radius_m = 0.5_dp
    real(dp), parameter :: iota_abs = 0.5_dp
    real(dp), parameter :: evaluation_s = 0.5_dp
    real(dp), parameter :: deep_fraction = 0.999_dp
    real(dp), parameter :: tolerance = 0.03_dp
    real(dp), parameter :: current_to_covar = 0.2_dp
    integer, parameter :: nsurf = 17

    character(len=16) :: arg
    integer :: orientation, iota_sign, field_sign
    real(dp) :: OmtB, d1, d2, eta, v, radius_cm, b_inner, expected, q_chart

    call get_command_argument(1, arg); read (arg, *) orientation
    call get_command_argument(2, arg); read (arg, *) iota_sign
    call get_command_argument(3, arg); read (arg, *) field_sign

    call write_chart("in_file", orientation, iota_sign, field_sign)
    call write_config("precession.in")

    call neort_init("precession.in", "in_file")
    call set_s(evaluation_s)
    call init_magfie_at_s()
    call init_mph_from_shared()
    call init_profiles(R0)
    call init_flux_surface_average(s)
    call init_canon_freq_trapped_spline()
    sign_vpar = 1
    call set_to_trapped_region(etamin, etamax)

    v = vth
    eta = etatp + deep_fraction*(etadt - etatp)
    call Om_tB(v, eta, OmtB, d1, d2)

    radius_cm = 1.0e2_dp*minor_radius_m*sqrt(evaluation_s)
    b_inner = 1.0e4_dp*field_t*(1.0_dp - radius_cm/(1.0e2_dp*major_radius_m))
    q_chart = 1.0_dp/(iota_sign*iota_abs)
    expected = -q_chart*orientation*field_sign*c*mi*v**2*1.0e4_dp*field_t &
        /(2.0_dp*qi*b_inner**2*radius_cm*1.0e2_dp*major_radius_m)

    print '(a,3i3,a,es13.5,a,es13.5)', "orientation, iota sign, field sign", &
        orientation, iota_sign, field_sign, ": Om_tB =", OmtB, "  oracle =", expected
    if (sign(1.0_dp, OmtB) /= sign(1.0_dp, expected)) then
        error stop "trapped precession has the wrong toroidal sign"
    end if
    if (abs(OmtB - expected) > tolerance*abs(expected)) then
        error stop "trapped precession magnitude disagrees with the oracle"
    end if

contains

    subroutine write_chart(path, orientation, iota_sign, field_sign)
        ! Coordinate phi is the cylindrical azimuth, so (s, theta, phi) has
        ! handedness -orientation, and the stored flux follows NEO-RT's
        ! convention flux = -handedness*field_sign*pi*a**2*B0.
        character(len=*), intent(in) :: path
        integer, intent(in) :: orientation, iota_sign, field_sign

        integer :: u, k, handedness
        real(dp) :: s_k, r_m, iota_chart, Bph_si, Bth_si, Jpol, Itor, dvds, flux

        handedness = -orientation
        iota_chart = iota_sign*iota_abs
        flux = -handedness*field_sign*pi*minor_radius_m**2*field_t
        Bph_si = field_sign*field_t*major_radius_m

        open (newunit=u, file=path, status="replace", action="write")
        write (u, "(a)") "CC manufactured large-aspect-ratio circular chart"
        write (u, "(a)") "CC"
        write (u, "(a)") "CC"
        write (u, "(a)") "CC"
        write (u, "(a)") "CC m0b n0b nsurf nper flux a R"
        write (u, *) 1, 0, nsurf, 1, flux, minor_radius_m, major_radius_m
        do k = 1, nsurf
            s_k = real(k, dp)/real(nsurf + 1, dp)
            r_m = minor_radius_m*sqrt(s_k)
            Bth_si = field_sign*field_t*r_m**2*iota_chart/major_radius_m
            Jpol = Bph_si*1.0e6_dp/(current_to_covar*handedness)
            Itor = Bth_si*1.0e6_dp/(current_to_covar*handedness)
            dvds = 4.0_dp*pi**2*handedness*major_radius_m*minor_radius_m**2/2.0_dp
            write (u, "(a)") "CC s iota Jpol Itor pprime sqrtg00"
            write (u, "(a)") "CC"
            write (u, *) s_k, iota_chart, Jpol, Itor, 0.0_dp, dvds
            write (u, "(a)") "CC m n rmnc rmns zmnc zmns vmnc vmns bmnc bmns"
            write (u, *) 0, 0, major_radius_m, 0.0_dp, 0.0_dp, 0.0_dp, &
                0.0_dp, 0.0_dp, field_t, 0.0_dp
            write (u, *) 1, 0, r_m, 0.0_dp, 0.0_dp, orientation*r_m, &
                0.0_dp, 0.0_dp, -field_t*r_m/major_radius_m, 0.0_dp
        end do
        close (u)
    end subroutine write_chart

    subroutine write_config(path)
        character(len=*), intent(in) :: path
        integer :: u

        open (newunit=u, file=path, status="replace", action="write")
        write (u, "(a)") "&params"
        write (u, "(a,es24.16)") "s = ", evaluation_s
        write (u, "(a)") "M_t = 0.0d0, qs = 1.0d0, ms = 2.014d0, vth = 1.0d8"
        write (u, "(a)") "epsmn = 1.0d-3, m0 = 0, mph = 1, magdrift = .true."
        write (u, "(a)") "nopassing = .true., noshear = .false., pertfile = .false."
        write (u, "(a)") "nonlin = .false., comptorque = .false., bfac = 1.0d0"
        write (u, "(a)") "efac = 1.0d0, inp_swi = 9, vsteps = 64, log_level = -1"
        write (u, "(a)") "/"
        close (u)
    end subroutine write_config

end program test_trapped_precession_sign
