program test_eqdsk_boozer_transport
    ! The direct chart's complete perturbation transport chain is compared
    ! with the independent Boozer Fourier backend on the same physical orbits.
    ! The latter needs neither phi_B-phi, the y(7) toroidal integral, nor the
    ! direct-chart phase. Several positive and negative bounce harmonics
    ! prevent a conjugate phase from giving a vacuous |H_mn|^2 comparison.
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use do_magfie_mod, only: inp_swi, bfac, read_boozer_file, set_s, &
        init_magfie_at_s, psi_pr
    use do_magfie_pert_mod, only: inp_swi_pert, pert_angle_map, &
        read_boozer_pert_file, init_magfie_pert_at_s, mph
    use driftorbit, only: etatp, etadt, mth, sign_vpar, sign_vpar_htheta, &
        pertfile, pertfile_scale, nonlin
    use neort_magfie, only: init_flux_surface_average
    use neort_orbit, only: bounce_time, bounce_fast, nvar, noshear
    use neort_transport, only: timestep_transport, Omth
    use util, only: pi, qi, mi, qe, mu

    implicit none

    integer, parameter :: nk = 2, nclass = 2, nharm = 5
    integer, parameter :: harmonics(nharm) = [-2, -1, 0, 1, 2]
    real(dp), parameter :: surfaces(2) = [0.35_dp, 0.75_dp]
    real(dp), parameter :: kappa(nk) = [0.35_dp, 0.7_dp]
    real(dp), parameter :: v = 1.0e8_dp
    ! The unmodified trapped turn locator has a ~1e-2 period error; do not
    ! require its fix here. Passing trajectories do not have turning points.
    real(dp), parameter :: power_tolerance(nclass) = [3.0e-2_dp, 2.0e-3_dp]
    character(len=1024) :: geqdsk, boozer, perturbation
    real(dp) :: flux_ratio, eta(nk, nclass), power(nharm, nk, nclass, 2)
    integer :: isurf, failures

    call read_fixture_paths()
    call initialize_physics()
    inp_swi = 11
    call read_boozer_file(trim(geqdsk))
    flux_ratio = abs(psi_pr)
    inp_swi = 9
    call read_boozer_file(trim(boozer))
    flux_ratio = flux_ratio/abs(psi_pr)

    failures = 0
    print '(a)', "HARMONIC_CSV,s,class,kappa,eta,mth,H2_Boozer_erg2,H2_direct_erg2"
    do isurf = 1, size(surfaces)
        call initialize_backend(9, trim(boozer), surfaces(isurf)*flux_ratio)
        call choose_physical_pitches(eta)
        call measure_harmonic_power(eta, power(:, :, :, 1))
        call initialize_backend(11, trim(geqdsk), surfaces(isurf))
        call measure_harmonic_power(eta, power(:, :, :, 2))
        call write_harmonic_csv(surfaces(isurf), eta, power)
        call compare_harmonic_power(surfaces(isurf), power, failures)
    end do
    if (failures > 0) error stop "EQDSK and Boozer perturbation transport disagree"
    print *, "PASS test_eqdsk_boozer_transport"

contains

    subroutine read_fixture_paths()
        call get_environment_variable("EQDSK_SOLOVEV_FILE", geqdsk)
        call get_environment_variable("BOOZER_SOLOVEV_FILE", boozer)
        call get_environment_variable("PERT_SOLOVEV_FILE", perturbation)
        if (len_trim(geqdsk) == 0) error stop "EQDSK_SOLOVEV_FILE must be set"
        if (len_trim(boozer) == 0) error stop "BOOZER_SOLOVEV_FILE must be set"
        if (len_trim(perturbation) == 0) error stop "PERT_SOLOVEV_FILE must be set"
    end subroutine read_fixture_paths

    subroutine initialize_physics()
        qi = qe
        mi = 2.0_dp*mu
        bfac = 1.0_dp
        inp_swi_pert = 9
        pertfile = .true.
        pertfile_scale = 1.0_dp
        nonlin = .false.
        noshear = .false.
        sign_vpar = 1.0_dp
    end subroutine initialize_physics

    subroutine initialize_backend(switch, path, s_surface)
        integer, intent(in) :: switch
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: s_surface

        inp_swi = switch
        pert_angle_map = ""
        if (switch == 11) pert_angle_map = boozer
        call read_boozer_file(path)
        call set_s(s_surface)
        call init_magfie_at_s()
        call init_flux_surface_average(s_surface)
        call read_boozer_pert_file(trim(perturbation))
        call init_magfie_pert_at_s()
        if (mph /= -2) error stop "Solov'ev perturbation must have toroidal mode -2"
    end subroutine initialize_backend

    subroutine choose_physical_pitches(pitches)
        real(dp), intent(out) :: pitches(nk, nclass)

        ! Eta is a physical invariant, independent of the flux chart. Use
        ! exactly these values in both backends, rather than matching a
        ! normalized pitch derived from their slightly different extrema.
        pitches(:, 1) = etatp + kappa*(etadt - etatp)
        pitches(:, 2) = kappa*etatp
    end subroutine choose_physical_pitches

    subroutine measure_harmonic_power(pitches, hmn2)
        real(dp), intent(in) :: pitches(nk, nclass)
        real(dp), intent(out) :: hmn2(nharm, nk, nclass)
        real(dp) :: taub, bounceavg(nvar), energy_squared
        integer :: iclass, k, ih, istate

        energy_squared = (mi*v**2/2.0_dp)**2
        do iclass = 1, nclass
            do k = 1, nk
                taub = bounce_time(v, pitches(k, iclass))
                if (.not. ieee_is_finite(taub)) error stop "nonfinite bounce period"
                if (taub <= 0.0_dp) error stop "nonpositive bounce period"
                Omth = sign_vpar_htheta*2.0_dp*pi/taub
                do ih = 1, nharm
                    mth = harmonics(ih)
                    call bounce_fast(v, pitches(k, iclass), taub, bounceavg, &
                        timestep_transport, istate)
                    if (istate /= 2) error stop "transport orbit integration failed"
                    hmn2(ih, k, iclass) = &
                        sum(bounceavg(3:4)**2)*energy_squared
                end do
            end do
        end do
        if (.not. all(ieee_is_finite(hmn2))) error stop "nonfinite harmonic power"
    end subroutine measure_harmonic_power

    subroutine write_harmonic_csv(s_surface, pitches, hmn2)
        real(dp), intent(in) :: s_surface, pitches(nk, nclass)
        real(dp), intent(in) :: hmn2(nharm, nk, nclass, 2)
        integer :: iclass, k, ih

        do iclass = 1, nclass
            do k = 1, nk
                do ih = 1, nharm
                    print '(a,f6.3,",",i1,",",f6.3,",",es17.9,",",i3,' &
                        //'2(",",es17.9))', "HARMONIC_CSV,", s_surface, iclass, &
                        kappa(k), pitches(k, iclass), harmonics(ih), &
                        hmn2(ih, k, iclass, :)
                end do
            end do
        end do
    end subroutine write_harmonic_csv

    subroutine compare_harmonic_power(s_surface, hmn2, nfail)
        real(dp), intent(in) :: s_surface, hmn2(nharm, nk, nclass, 2)
        integer, intent(inout) :: nfail
        real(dp) :: reference_scale, error
        integer :: iclass, k

        do iclass = 1, nclass
            do k = 1, nk
                reference_scale = maxval(hmn2(:, k, iclass, 1))
                if (reference_scale < 1.0e-30_dp) then
                    error stop "Boozer harmonic-power oracle is vacuous"
                end if
                ! Normalize by peak power on this orbit: near-zero Fourier
                ! coefficients must not dominate the comparison.
                error = maxval(abs(hmn2(:, k, iclass, 2) &
                    - hmn2(:, k, iclass, 1)))/reference_scale
                print '(a,f5.2,a,i1,a,f5.2,a,es10.3,a,es10.3)', &
                    "s=", s_surface, " class=", iclass, " kappa=", kappa(k), &
                    " max|H|^2 [erg^2]=", reference_scale, " relative error=", error
                if (error > power_tolerance(iclass)) nfail = nfail + 1
            end do
        end do
    end subroutine compare_harmonic_power

end program test_eqdsk_boozer_transport
