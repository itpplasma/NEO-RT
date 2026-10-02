program diagnose_line_drive
    !! Small S3 evidence tables on examples/base, using the ideal_helical field.
    !! Run via fo exec diagnose_line_drive.x. The output directory must exist.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: do_magfie_init, R0, s, q
    use driftorbit, only: efac, bfac, mth, m0, epsmn, etatp, etadt, &
                          sign_vpar, sign_vpar_htheta
    use neort, only: init
    use neort_config, only: read_and_set_config
    use neort_orbit, only: noshear, vpar, th0, evaluate_bfield_local
    use neort_freq, only: Om_th
    use neort_profiles, only: init_profiles, read_and_init_plasma_input, &
                              read_and_init_profile_input, vth
    use neort_line_drive, only: set_line_drive_options, line_bounce, &
                                line_harmonics_t, exact_orbit_period, &
                                gauge_amp, collapse_eta, line_point, gauge_fix_t, &
                                local_field_t, local_field, NCOMP
    use do_magfie_pert_mod, only: mph
    use fortnum_ode_vode, only: vode_state_t, vode_init, vode_integrate_to
    use fortnum_status, only: fortnum_status_t, FORTNUM_OK
    use util, only: pi, imun, mi

    implicit none

    integer, parameter :: nsample = 256, npitch = 18
    character(len=1024) :: output_dir
    integer :: status, unit, iclass, j, mb
    real(dp) :: fraction, eta, period, omega, om_path, om_te, max_error(2), sampled_eta
    type(line_harmonics_t) :: harmonics

    call get_environment_variable('LINE_DRIVE_DIAGNOSTICS_DIR', output_dir, &
                                  status=status)
    if (status /= 0) output_dir = '.'
    if (len_trim(output_dir) == 0) output_dir = '.'
    call read_and_set_config('driftorbit.in')
    noshear = .false.
    m0 = 2
    epsmn = 1.0e-3_dp
    gauge_amp = 0
    collapse_eta = .true.
    call set_line_drive_options('line', 'ideal_helical', .false., 8, noshear)
    call do_magfie_init('in_file')
    call init_profiles(R0)
    call read_and_init_plasma_input('plasma.in', s)
    call read_and_init_profile_input('profile.in', s, R0, efac, bfac)
    call init
    sign_vpar = 1.0_dp
    max_error = 0

    open (newunit=unit, file=trim(output_dir)//'/harmonics.csv', status='replace')
    write (unit, '(a)') 'orbit,pitch_fraction,eta_per_gauss,m_b,n,'// &
        'period_s,omega_rad_per_s,om_te_rad_per_s,line_re,line_im,'// &
        'boozer_re,boozer_im,line_abs2,boozer_abs2,relative_error'
    do iclass = 1, 2
        mb = -4
        if (iclass == 2) mb = 0
        do j = 1, npitch
            fraction = 0.05_dp + 0.4_dp*real(j - 1, dp)/real(npitch - 1, dp)
            if (iclass == 2) fraction = fraction + 0.5_dp
            call resonant_harmonic(fraction, mb, eta, period, omega, &
                                   om_path, om_te, harmonics)
            max_error(iclass) = max(max_error(iclass), relative_error(harmonics))
            write (unit, '(a,",",es24.16,",",es24.16,",",i0,",",i0,'// &
                   '10(",",es24.16))') orbit_name(iclass), fraction, eta, mb, mph, &
                period, omega, om_te, real(harmonics%h_line), &
                aimag(harmonics%h_line), real(harmonics%h_boozer), &
                aimag(harmonics%h_boozer), abs(harmonics%h_line)**2, &
                abs(harmonics%h_boozer)**2, relative_error(harmonics)
        end do
    end do
    close (unit)

    call sample_orbit(0.4_dp, -4, 'passing')
    call sample_orbit(0.6_dp, 0, 'trapped')
    open (newunit=unit, file=trim(output_dir)//'/metadata.txt', status='replace')
    write (unit, '(a)') 'S3 ideal_helical, examples/base, no gauge shift'
    write (unit, '(a)') 'All harmonics and pointwise drives are divided '// &
        'by E=mi*vth**2/2.'
    write (unit, '(a)') 'Harmonic scan uses a different resonant E x B '// &
        'rotation at each pitch.'
    write (unit, '(a)') 'Pitch fraction: eta/(2*etatp) passing; '// &
        '0.5+(eta-etatp)/(2*(etadt-etatp)) trapped.'
    write (unit, '(a)') 'Scan excludes eta=0, the trapped-passing boundary and etadt.'
    write (unit, '(a)') 'Pointwise integrands include the same Fourier '// &
        'phase as line_bounce.'
    write (unit, '(a,es24.16)') 's = ', s
    write (unit, '(a,es24.16)') 'vth [cm/s] = ', vth
    write (unit, '(a,es24.16)') 'epsmn = ', epsmn
    write (unit, '(a,i0)') 'm0 = ', m0
    write (unit, '(a,i0)') 'n = ', mph
    write (unit, '(a,es24.16)') 'passing max relative harmonic error = ', max_error(1)
    write (unit, '(a,es24.16)') 'trapped max relative harmonic error = ', max_error(2)
    close (unit)
    print '(a)', 'Line-drive diagnostics written to '//trim(output_dir)
    print '(a,2es12.3)', 'Max relative harmonic errors (passing, trapped): ', max_error

contains

    function orbit_name(iclass) result(name)
        integer, intent(in) :: iclass
        character(len=:), allocatable :: name

        name = 'passing'
        if (iclass == 2) name = 'trapped'
    end function orbit_name

    real(dp) function relative_error(res)
        type(line_harmonics_t), intent(in) :: res

        if (abs(res%h_boozer) <= tiny(1.0_dp)) then
            error stop 'diagnostics: vacuous Boozer harmonic'
        end if
        relative_error = abs(res%h_line - res%h_boozer)/abs(res%h_boozer)
    end function relative_error

    subroutine resonant_harmonic(frac, mb, eta, period, omega, om_path, om_te, res)
        real(dp), intent(in) :: frac
        integer, intent(in) :: mb
        real(dp), intent(out) :: eta, period, omega, om_path, om_te
        type(line_harmonics_t), intent(out) :: res
        real(dp) :: derivative_v, derivative_eta, mdotom
        integer :: istate

        eta = 2.0_dp*frac*etatp
        if (frac > 0.5_dp) eta = etatp + 2.0_dp*(frac - 0.5_dp)*(etadt - etatp)
        mth = mb
        call Om_th(vth, eta, omega, derivative_v, derivative_eta)
        call exact_orbit_period(vth, eta, 2.0_dp*pi/abs(omega), period, omega)
        mdotom = mb*omega
        if (eta <= etatp) mdotom = (mb + q*mph)*omega
        om_path = -mdotom/mph
        call line_bounce(vth, eta, period, omega, om_path, 0.0_dp, res, istate)
        if (istate /= 2) error stop 'diagnostics: drift integration failed'
        om_te = om_path - res%om_drift
        call line_bounce(vth, eta, period, omega, om_path, om_te, res, istate)
        if (istate /= 2) error stop 'diagnostics: harmonic integration failed'
    end subroutine resonant_harmonic

    subroutine sample_orbit(frac, mb, name)
        real(dp), intent(in) :: frac
        integer, intent(in) :: mb
        character(len=*), intent(in) :: name
        type(vode_state_t) :: state
        type(fortnum_status_t) :: solver_status
        type(gauge_fix_t) :: gf
        type(local_field_t) :: lf
        type(line_harmonics_t) :: res
        real(dp), allocatable :: y(:)
        real(dp) :: initial(2), eta, period, omega, om_path, om_te, t, energy
        real(dp) :: bmod, htheta, drift_rate, phase
        complex(dp) :: values(NCOMP), fourier_phase
        integer :: unit, k

        call resonant_harmonic(frac, mb, eta, period, omega, om_path, om_te, res)
        sampled_eta = eta
        call evaluate_bfield_local(bmod, htheta)
        sign_vpar_htheta = sign(1.0_dp, htheta)*sign_vpar
        initial = [th0, sign_vpar_htheta*vpar(vth, eta, bmod)]
        y = initial
        energy = 0.5_dp*mi*vth**2
        ! ideal_helical has dA_parallel=0, so its gauge-fix coefficients vanish.
        gf = gauge_fix_t()
        call vode_init(state, 2, 0.0_dp, initial)
        open (newunit=unit, file=trim(output_dir)//'/integrand_'//name//'.csv', &
              status='replace')
        write (unit, '(a)') 'time_over_period,time_s,theta_rad,vpar_cm_per_s,'// &
            'bmod_gauss,line_re,line_im,boozer_re,boozer_im'
        do k = 0, nsample
            t = period*real(k, dp)/real(nsample, dp)
            if (k > 0) then
                call vode_integrate_to(rhs, state, t, 1.0e-11_dp, &
                                       [1.0e-13_dp, 1.0e-11_dp*vth], y, solver_status)
                if (solver_status%code /= FORTNUM_OK) then
                    error stop 'diagnostics: sampled orbit integration failed'
                end if
            end if
            call line_point(vth, eta, y(1), y(2), om_path, om_te, &
                            gf, lf, values, drift_rate)
            phase = mph*q*y(1) - mb*omega*t
            if (eta <= etatp) phase = mph*q*y(1) - (mb + q*mph)*omega*t
            fourier_phase = exp(imun*phase)/energy
            values = values*fourier_phase
            write (unit, '(es24.16,8(",",es24.16))') &
                t/period, t, y(1), y(2), lf%bmod, real(values(1)), &
                aimag(values(1)), real(values(2)), aimag(values(2))
        end do
        close (unit)
    end subroutine sample_orbit

    subroutine rhs(t, y, dydt, context)
        real(dp), intent(in) :: t, y(:)
        real(dp), intent(out) :: dydt(:)
        class(*), intent(in), optional :: context
        type(local_field_t) :: lf
        associate (unused_t => t, unused_context => context)
        end associate
        call local_field(s, y(1), lf)
        dydt(1) = y(2)*lf%hth
        dydt(2) = -0.5_dp*vth**2*sampled_eta*lf%hth*lf%dB_dth
    end subroutine rhs

end program diagnose_line_drive
