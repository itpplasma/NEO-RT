module boozer_perturbation_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use netcdf, only: nf90_open, nf90_nowrite, nf90_get_att, nf90_global, &
        nf90_close, nf90_noerr
    use libneo_coordinates_base, only: coordinate_system_t, PSI_TOR_NORM, BOOZER
    use libneo_coordinates_chartmap, only: make_chartmap_coordinate_system, &
        chartmap_coordinate_system_t
    use do_magfie_pert_mod, only: read_boozer_pert_file, init_magfie_pert_at_s, &
        do_magfie_pert_amp, inp_swi_pert, mph_shared
    use fortnum_multiroot, only: multiroot_hybrid
    use fortnum_status, only: fortnum_status_t
    implicit none
    private
    public :: use_boozer_perturbation, boozer_geometry_file, &
        boozer_perturbation_file, boozer_perturbation

    logical :: use_boozer_perturbation = .false.
    character(len=256) :: boozer_geometry_file = 'boozer_geometry.nc'
    character(len=256) :: boozer_perturbation_file = 'in_file_pert'
    class(coordinate_system_t), allocatable :: geometry
    logical :: initialized = .false.
    real(dp) :: support_s(2), zero_disks(6)
    logical :: thread_initialized = .false.
    real(dp) :: previous_u(3) = 0.0_dp
    logical :: have_previous_u = .false.
    !$omp threadprivate(thread_initialized, previous_u, have_previous_u)

    type :: inverse_target_t
        real(dp) :: x(3)
    end type inverse_target_t

contains

    subroutine initialize_boozer_perturbation(n_tor)
        integer, intent(in) :: n_tor
        logical :: ready
        integer :: ncid, status

        !$omp atomic read
        ready = initialized
        if (ready) return
        !$omp critical (potato_boozer_init)
        if (.not. initialized) then
            status = nf90_open(trim(boozer_geometry_file), nf90_nowrite, ncid)
            if (status /= nf90_noerr) error stop 'POTATO Boozer geometry open failed'
            status = nf90_get_att(ncid, nf90_global, 'support_s', support_s)
            if (status /= nf90_noerr) error stop 'Boozer geometry needs support_s'
            status = nf90_get_att(ncid, nf90_global, 'zero_disks', zero_disks)
            if (status /= nf90_noerr) error stop 'Boozer geometry needs zero_disks'
            status = nf90_close(ncid)
            if (status /= nf90_noerr) error stop 'POTATO Boozer geometry close failed'
            if (.not. all(ieee_is_finite(support_s)) .or. &
                .not. all(ieee_is_finite(zero_disks))) &
                error stop 'POTATO Boozer geometry support must be finite'
            if (support_s(1) >= support_s(2) .or. &
                minval(zero_disks([3, 6])) < 0.0_dp) &
                error stop 'POTATO Boozer geometry support is invalid'
            call make_chartmap_coordinate_system(geometry, trim(boozer_geometry_file))
            select type (g => geometry)
            type is (chartmap_coordinate_system_t)
                if (g%rho_convention /= PSI_TOR_NORM) &
                    error stop 'POTATO Boozer geometry requires rho_convention=s'
                if (g%zeta_convention /= BOOZER) &
                    error stop 'POTATO geometry requires zeta_convention=boozer'
            end select
            inp_swi_pert = 9
            call read_boozer_pert_file(trim(boozer_perturbation_file))
            if (mph_shared /= n_tor) error stop 'POTATO/BC toroidal modes differ'
            !$omp atomic write
            initialized = .true.
        end if
        !$omp end critical (potato_boozer_init)
    end subroutine initialize_boozer_perturbation

    subroutine boozer_perturbation(r, z, n_tor, amplitude)
        real(dp), intent(in) :: r, z
        integer, intent(in) :: n_tor
        complex(dp), intent(out) :: amplitude
        real(dp) :: u(3), x(3), guess(3), roundtrip(3)
        integer :: ierr
        type(inverse_target_t) :: target
        type(fortnum_status_t) :: status

        if (.not. all(ieee_is_finite([r, z]))) &
            error stop 'POTATO Boozer inverse: nonfinite R,Z'
        call initialize_boozer_perturbation(n_tor)
        amplitude = cmplx(0.0_dp, 0.0_dp, dp)
        ! Certified disks from native BC Fourier contours avoid inversion in the
        ! uncharted core and far exterior. There the opt-in field is zero.
        if (sum(([r, z] - zero_disks(1:2))**2) < zero_disks(3)) return
        if (sum(([r, z] - zero_disks(4:5))**2) > zero_disks(6)) return
        if (.not. thread_initialized) then
            call init_magfie_pert_at_s()
            thread_initialized = .true.
        end if
        ! libneo returns (s, theta_B, phi_B) at cylindrical phi_geom=0.
        ! This also avoids identifying poloidal flux with normalized toroidal flux.
        if (have_previous_u) then
            guess = previous_u
            ! Request the full interior solve, including theta, at edge points.
            ! The libneo edge shortcut holds the previous angular ray fixed.
            guess(1) = min(guess(1), 0.89_dp)
            select type (g => geometry)
            type is (chartmap_coordinate_system_t)
                call g%invert_cart([r, 0.0_dp, z], guess, u, ierr)
            end select
        else
            call geometry%from_cyl([r, 0.0_dp, z], u, ierr)
        end if
        roundtrip = huge(1.0_dp)
        if (ierr == 0) call geometry%evaluate_cart(u, roundtrip)
        ! The warm solver's LOCATED status allows a fraction of a radial cell.
        ! Near the axis that can exceed our field-coordinate accuracy contract.
        ! Retry from independent seeds rather than accept a displaced BC surface.
        if (ierr /= 0 .or. &
            maxval(abs(roundtrip - [r, 0.0_dp, z])) > 1.0e-8_dp) then
            call geometry%from_cyl([r, 0.0_dp, z], u, ierr)
            if (ierr == 0) call geometry%evaluate_cart(u, roundtrip)
        end if
        if (ierr /= 0) then
            print *, 'POTATO Boozer inverse failed at R,Z/status:', r, z, ierr
            error stop 'POTATO Boozer coordinate inverse failed'
        end if
        if (maxval(abs(roundtrip - [r, 0.0_dp, z])) > 1.0e-8_dp) then
            ! The cold inverse also permits a step-size exit at chart accuracy.
            ! Refine that located root to POTATO's stricter Cartesian tolerance.
            ! Keep the final round-trip check independent of solver status.
            target%x = [r, 0.0_dp, z]
            guess = u
            call multiroot_hybrid(inverse_residual, 3, guess, u, status, &
                xtol=epsilon(1.0_dp), ftol=1.0e-10_dp, max_iter=8, ctx=target)
            call geometry%evaluate_cart(u, roundtrip)
        end if
        if (.not. all(ieee_is_finite(roundtrip))) &
            error stop 'POTATO Boozer inverse returned nonfinite coordinates'
        if (maxval(abs(roundtrip - [r, 0.0_dp, z])) > 1.0e-8_dp) &
            error stop 'POTATO Boozer inverse residual exceeds 1e-8 cm'
        previous_u = u
        have_previous_u = .true.
        if (u(1) < support_s(1) .or. u(1) > support_s(2)) return
        x = [u(1), 0.0_dp, u(2)]
        call do_magfie_pert_amp(x, amplitude)
        amplitude = amplitude*exp(cmplx(0.0_dp, n_tor*u(3), dp))
    end subroutine boozer_perturbation

    subroutine inverse_residual(u, residual, jacobian, ctx)
        real(dp), intent(in) :: u(:)
        real(dp), intent(out) :: residual(:), jacobian(:, :)
        class(*), intent(in), optional :: ctx

        if (.not. present(ctx)) error stop 'missing POTATO inverse target'
        select type (ctx)
        type is (inverse_target_t)
            call geometry%evaluate_cart(u, residual)
            residual = residual - ctx%x
            call geometry%covariant_basis(u, jacobian)
        class default
            error stop 'invalid POTATO inverse target'
        end select
    end subroutine inverse_residual
end module boozer_perturbation_mod
