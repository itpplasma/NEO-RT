program test_eqdsk_boozer_backend
    ! The same Solov'ev field read two ways: directly from GEQDSK (inp_swi=11)
    ! and as the Boozer file that libneo's efit_to_boozer.x makes of it
    ! (inp_swi=9).  Surface quantities, the bounce/transit frequency and the
    ! bounce-averaged magnetic precession Om_tB must agree on trapped and
    ! co-passing orbits; the frequencies are also checked against the exact
    ! field.  The Boozer path is an
    ! independent implementation (Fourier series in Boozer angles, closed-form
    ! precession), so this is a behavioural check of the direct-chart drift
    ! projection including the radial field-line-label term (issue #85).
    !
    ! efit_to_boozer.x places its last surface slightly inside the GEQDSK
    ! boundary, so its s is normalized to a smaller edge flux.  Surfaces are
    ! therefore matched by absolute toroidal flux, s_B = s*psi_pr/psi_pr_B.
    ! Counter-passing orbits are left out: bounce() only detects a turn with
    ! increasing theta, in either backend.
    use iso_fortran_env, only: dp => real64
    use do_magfie_mod, only: inp_swi, bfac, read_boozer_file, set_s, &
        init_magfie_at_s, q, psi_pr, sign_theta
    use driftorbit, only: B0, etatp, etadt, dVds
    use neort_magfie, only: init_flux_surface_average
    use neort_orbit, only: bounce, nvar, noshear
    use util, only: pi, qi, mi, qe, mu

    implicit none

    integer, parameter :: nk = 7, nclass = 2
    real(dp), parameter :: s0 = 0.35_dp, v = 1.0e8_dp
    real(dp), parameter :: kappa(nk) = [0.05_dp, 0.2_dp, 0.35_dp, 0.5_dp, &
        0.65_dp, 0.8_dp, 0.95_dp]
    ! Bounce (trapped) and transit (co-passing) frequencies in rad/s of the
    ! exact Solov'ev field at s0, v, from an independent integration of the
    ! same orbit equations on the analytic surface (scipy DOP853, rtol 1e-12,
    ! geometry by bisection and spectral differentiation, no libneo).
    real(dp), parameter :: om_exact(nk, nclass) = reshape([ &
        7.843362e4_dp, 1.071284e5_dp, 1.250950e5_dp, 1.396712e5_dp, &
        1.527002e5_dp, 1.653318e5_dp, 1.786544e5_dp, &
        4.316009e5_dp, 4.028513e5_dp, 3.716217e5_dp, 3.370387e5_dp, &
        2.974866e5_dp, 2.491611e5_dp, 1.750738e5_dp], [nk, nclass])
    ! bounce() resolves the trapped turn only to ~1e-2 in either backend
    ! (Boozer 9.3e-3, direct 3.1e-3 here); transit times have no turning
    ! point and are resolved to 1e-6.
    real(dp), parameter :: tol_om(nclass) = [2.0e-2_dp, 1.0e-5_dp]
    character(len=1024) :: geqdsk, boozer
    real(dp) :: om_b(nk, nclass, 2), om_t(nk, nclass, 2), surf(5, 2)
    real(dp) :: om_t_noshear(nk, nclass), dummy_b(nk, nclass), dummy_s(5)
    real(dp) :: om_t_again(nk, nclass)
    real(dp) :: err_s, shear_size, psi_pr_direct, s_boozer
    integer :: iclass, ib, nfail

    call get_environment_variable("EQDSK_SOLOVEV_FILE", geqdsk)
    call get_environment_variable("BOOZER_SOLOVEV_FILE", boozer)
    if (len_trim(geqdsk) == 0 .or. len_trim(boozer) == 0) then
        error stop "EQDSK_SOLOVEV_FILE and BOOZER_SOLOVEV_FILE must be set"
    end if
    qi = qe
    mi = 2.0_dp*mu
    bfac = 1.0_dp

    inp_swi = 11
    call read_boozer_file(trim(geqdsk))
    psi_pr_direct = psi_pr
    ! Initialize the direct input first: the Fourier buffers have never been
    ! allocated when the following Boozer scan begins.
    call set_s(s0)
    call init_magfie_at_s()
    inp_swi = 9
    call read_boozer_file(trim(boozer))
    s_boozer = s0*abs(psi_pr_direct/psi_pr)
    print '(a,f10.6)', "Boozer-file s of the matched surface: ", s_boozer

    call scan_backend(9, trim(boozer), s_boozer, om_b(:, :, 1), om_t(:, :, 1), &
        surf(:, 1))
    noshear = .true.
    call scan_backend(9, trim(boozer), s_boozer, dummy_b, om_t_noshear, dummy_s)
    noshear = .false.
    call scan_backend(11, trim(geqdsk), s0, om_b(:, :, 2), om_t(:, :, 2), &
        surf(:, 2))
    ! Switching back from the direct field must reallocate the Fourier work
    ! buffers and reproduce the first Boozer scan exactly.
    call scan_backend(9, trim(boozer), s_boozer, dummy_b, om_t_again, dummy_s)

    call print_table(om_b, om_t)
    print '(a,5es14.6)', "surface booz:  ", surf(:, 1)
    print '(a,5es14.6)', "surface direct:", surf(:, 2)
    shear_size = maxval(abs(om_t(:, :, 1) - om_t_noshear))/maxval(abs(om_t(:, :, 1)))
    print '(a,es10.3)', "Boozer shear term / max|Om_tB|: ", shear_size

    nfail = 0
    err_s = maxval(abs(surf(:, 2)/surf(:, 1) - 1.0_dp))
    call report("q, sign(sign_theta*psi_pr), B0, Bmax, dV/dpsi", err_s, &
        2.0e-3_dp, nfail)
    do iclass = 1, nclass
        do ib = 1, 2
            call report(trim(merge("trapped Om_b  ", "passing Om_b  ", iclass == 1)) &
                //merge(" Boozer", " direct", ib == 1), &
                maxval(abs(om_b(:, iclass, ib)/om_exact(:, iclass) - 1.0_dp)), &
                tol_om(iclass), nfail)
        end do
    end do
    call report("trapped Om_tB direct-Boozer / max|Om_tB|", &
        maxval(abs(om_t(:, 1, 2) - om_t(:, 1, 1)))/maxval(abs(om_t(:, 1, 1))), &
        1.0e-2_dp, nfail)
    call report("Boozer rescan after direct field, max |diff|", &
        maxval(abs(om_t_again - om_t(:, :, 1))), 0.0_dp, nfail)
    call report("passing Om_tB direct-Boozer / max|Om_tB|", &
        maxval(abs(om_t(:, 2, 2) - om_t(:, 2, 1)))/maxval(abs(om_t(:, 2, 1))), &
        1.0e-4_dp, nfail)
    if (nfail > 0) error stop "test_eqdsk_boozer_backend failed"
    print *, "PASS test_eqdsk_boozer_backend"

contains

    subroutine scan_backend(switch, path, s_surf, omega_b, omega_t, surface)
        integer, intent(in) :: switch
        character(len=*), intent(in) :: path
        real(dp), intent(in) :: s_surf
        real(dp), intent(out) :: omega_b(nk, nclass), omega_t(nk, nclass)
        real(dp), intent(out) :: surface(5)
        real(dp) :: eta, taub, bounceavg(nvar)
        integer :: k, iclass

        inp_swi = switch
        call read_boozer_file(path)
        call set_s(s_surf)
        call init_magfie_at_s()
        call init_flux_surface_average(s_surf)
        ! dV/d(psi_tor) rather than dV/ds, whose normalization differs.
        surface = [q, sign(1.0_dp, sign_theta*psi_pr), B0, 1.0_dp/etatp, &
            dVds/abs(psi_pr)]
        do iclass = 1, nclass
            do k = 1, nk
                ! Trapped and co-passing.
                if (iclass == 1) then
                    eta = etatp + kappa(k)*(etadt - etatp)
                else
                    eta = kappa(k)*etatp
                end if
                call bounce(v, eta, taub, bounceavg)
                omega_b(k, iclass) = 2.0_dp*pi/taub
                omega_t(k, iclass) = bounceavg(3)*v**2
            end do
        end do
    end subroutine scan_backend

    subroutine print_table(omega_b, omega_t)
        real(dp), intent(in) :: omega_b(nk, nclass, 2), omega_t(nk, nclass, 2)
        integer :: k, iclass

        print '(a)', "class kappa Om_b: exact Boozer direct, Om_tB: Boozer direct"
        do iclass = 1, nclass
            do k = 1, nk
                print '(i3,f7.2,5es13.5)', iclass, kappa(k), om_exact(k, iclass), &
                    omega_b(k, iclass, :), omega_t(k, iclass, :)
            end do
        end do
    end subroutine print_table

    subroutine report(label, value, tol, nfail)
        character(len=*), intent(in) :: label
        real(dp), intent(in) :: value, tol
        integer, intent(inout) :: nfail

        if (value <= tol) then
            print '(a,": ",es10.3," <= ",es9.2)', label, value, tol
        else
            print '("FAIL ",a,": ",es10.3," > ",es9.2)', label, value, tol
            nfail = nfail + 1
        end if
    end subroutine report

end program test_eqdsk_boozer_backend
