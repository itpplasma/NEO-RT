program probe_frequency_comparison
    !! Matched-equilibrium diagnostic for Boozer thin and native real-space
    !! finite-width frequencies.  POTATO writes its own canonical-frequency
    !! table from the same EQDSK; plotting and root finding stay outside the
    !! production solver.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, read_boozer_file, set_s, init_magfie_at_s, R0
    use driftorbit, only: FREQUENCY_MODEL_LEGACY, frequency_model, magdrift, &
        magdrift_passing, sign_vpar, etatp, etadt
    use neort_freq, only: init_canon_freq_trapped_spline, &
        init_canon_freq_passing_spline, Om_th, Om_ph
    use neort_gc_frequency_provider, only: GC_FREQUENCY_SUCCESS, &
        gc_frequency_context_t, gc_full_orbit_frequency_result_t, &
        initialize_gc_frequency_context, evaluate_gc_full_orbit_frequency
    use neort_gc_orbit_integrator, only: GC_ORBIT_PASSING, GC_ORBIT_TRAPPED
    use neort_magfie, only: init_flux_surface_average
    use neort_orbit, only: bounce_time, th0
    use neort_profiles, only: vth, Om_tE
    use util, only: qe, mu, qi, mi

    implicit none

    integer, parameter :: npoints = 180
    real(dp), parameter :: reference_speed = 6.9199e7_dp ! 5 keV D [cm/s]
    real(dp), parameter :: mass = 2.014_dp*mu
    character(len=1024) :: chartmap, eqdsk, output, wall_file, argument
    type(gc_frequency_context_t) :: context
    type(gc_full_orbit_frequency_result_t) :: full
    real(dp) :: boozer_surface, direct_surface, eta_ratio, eta
    real(dp) :: omega_b, omega_phi, d1, d2, period_estimate, speed, energy_kev
    integer :: unit, class_index, k, status, orbit_class
    character(len=8) :: class_name

    call get_command_argument(1, chartmap)
    call get_command_argument(2, eqdsk)
    call get_command_argument(3, argument)
    read(argument, *) boozer_surface
    call get_command_argument(4, argument)
    read(argument, *) direct_surface
    call get_command_argument(5, output)
    call get_command_argument(6, wall_file)
    if (len_trim(chartmap) == 0 .or. len_trim(eqdsk) == 0 &
        .or. len_trim(output) == 0 .or. len_trim(wall_file) == 0) then
        error stop 'six arguments are required: chartmap eqdsk boozer_s direct_s output wall_file'
    end if
    energy_kev = 5.0_dp
    call get_command_argument(7, argument)
    if (len_trim(argument) > 0) read(argument, *) energy_kev
    if (energy_kev <= 0.0_dp) error stop 'energy_kev must be positive'
    speed = reference_speed*sqrt(energy_kev/5.0_dp)

    qi = qe
    mi = mass
    vth = speed
    Om_tE = 0.0_dp
    sign_vpar = 1.0_dp
    magdrift = .true.
    magdrift_passing = 1
    frequency_model = FREQUENCY_MODEL_LEGACY
    open(newunit=unit, file=trim(output), status='replace', action='write')
    write(unit, '(a)') '# model class eta_over_etatp omega_b omega_phi status orbit_status'
    write(unit, '(a,es20.12)') '# energy_kev ', energy_kev

    inp_swi = 10
    call read_boozer_file(trim(chartmap))
    call set_s(boozer_surface)
    call init_magfie_at_s()
    call init_flux_surface_average(boozer_surface)
    call init_canon_freq_trapped_spline()
    call init_canon_freq_passing_spline()
    write(unit, '(a,es20.12)') '# boozer_surface_s_tor ', boozer_surface
    write(unit, '(a,es20.12)') '# boozer_eta_tp ', etatp
    do class_index = 1, 2
        class_name = merge('passing ', 'trapped ', class_index == 1)
        do k = 1, npoints
            eta_ratio = class_eta_ratio(class_index, k)
            eta = eta_ratio*etatp
            call Om_th(speed, eta, omega_b, d1, d2)
            call Om_ph(speed, eta, omega_phi, d1, d2)
            write(unit, '(a,1x,a,1x,3(es20.12,1x),2(i4,1x))') &
                'boozer_thin', trim(class_name), eta_ratio, omega_b, &
                omega_phi, 0, 0
        end do
    end do

    inp_swi = 11
    call read_boozer_file(trim(eqdsk))
    call set_s(direct_surface)
    call init_magfie_at_s()
    call init_flux_surface_average(direct_surface)
    call initialize_gc_frequency_context(direct_surface, th0, 1.0_dp, 0.0_dp, &
        mass, qe, speed, context, status, selected_frequency_model=2, &
        wall_file=trim(wall_file), wall_units='cm')
    if (status /= GC_FREQUENCY_SUCCESS) error stop 'full-orbit context failed'
    write(unit, '(a,es20.12)') '# direct_surface_s_tor ', direct_surface
    write(unit, '(a,es20.12)') '# direct_eta_tp ', etatp
    write(unit, '(a,es20.12)') '# direct_q_fieldline ', context%q_fieldline
    do class_index = 1, 2
        class_name = merge('passing ', 'trapped ', class_index == 1)
        orbit_class = merge(GC_ORBIT_PASSING, GC_ORBIT_TRAPPED, class_index == 1)
        do k = 1, npoints
            eta_ratio = class_eta_ratio(class_index, k)
            eta = eta_ratio*etatp
            ! The direct backend must not inherit the legacy chart's passing
            ! period estimate.  Its physical q is computed from the same
            ! EQDSK field that drives the return map.  Retain the legacy
            ! estimate only as a conservative trapped-orbit timeout, where
            ! the near-separatrix bounce period is singular in the thin
            ! chart and a short transit estimate would truncate the return.
            period_estimate = 12.0_dp*abs(context%q_fieldline)*R0/speed
            if (orbit_class == GC_ORBIT_TRAPPED) then
                period_estimate = max(period_estimate, bounce_time(speed, eta))
            end if
            call evaluate_gc_full_orbit_frequency(context, eta, 1, orbit_class, &
                abs(period_estimate), full, status)
            write(unit, '(a,1x,a,1x,3(es20.12,1x),2(i4,1x))') &
                'neort_full', trim(class_name), eta_ratio, full%omega_b, &
                full%omega_phi, status, full%orbit_status
        end do
    end do
    close(unit)

contains

    real(dp) function class_eta_ratio(class_value, point)
        integer, intent(in) :: class_value, point
        real(dp) :: fraction, trapped_max

        fraction = real(point - 1, dp)/real(npoints - 1, dp)
        if (class_value == 1) then
            class_eta_ratio = 0.01_dp + 0.989_dp*fraction
        else
            trapped_max = etadt/etatp
            class_eta_ratio = 1.001_dp &
                + (trapped_max - 1.002_dp)*fraction
        end if
    end function class_eta_ratio

end program probe_frequency_comparison
