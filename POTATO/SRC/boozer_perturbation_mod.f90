module boozer_perturbation_mod
    use, intrinsic :: iso_fortran_env, only: dp => real64
    use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
    use libneo_coordinates_base, only: coordinate_system_t, PSI_TOR_NORM, BOOZER
    use libneo_coordinates_chartmap, only: make_chartmap_coordinate_system, &
        chartmap_coordinate_system_t
    use do_magfie_pert_mod, only: read_boozer_pert_file, init_magfie_pert_at_s, &
        do_magfie_pert_amp, inp_swi_pert, mph_shared
    implicit none
    private
    public :: use_boozer_perturbation, boozer_geometry_file, &
        boozer_perturbation_file, boozer_perturbation

    logical :: use_boozer_perturbation = .false.
    character(len=256) :: boozer_geometry_file = 'boozer_geometry.nc'
    character(len=256) :: boozer_perturbation_file = 'in_file_pert'
    class(coordinate_system_t), allocatable :: geometry
    logical :: initialized = .false.
    logical :: thread_initialized = .false.
    real(dp) :: previous_u(3) = 0.0_dp
    logical :: have_previous_u = .false.
    !$omp threadprivate(thread_initialized, previous_u, have_previous_u)

contains

    subroutine initialize_boozer_perturbation(n_tor)
        integer, intent(in) :: n_tor
        logical :: ready

        !$omp atomic read
        ready = initialized
        if (ready) return
        !$omp critical (potato_boozer_init)
        if (.not. initialized) then
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

        if (.not. all(ieee_is_finite([r, z]))) &
            error stop 'POTATO Boozer inverse: nonfinite R,Z'
        call initialize_boozer_perturbation(n_tor)
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
        if (ierr /= 0) then
            print *, 'POTATO Boozer inverse failed at R,Z/status:', r, z, ierr
            error stop 'POTATO Boozer coordinate inverse failed'
        end if
        call geometry%evaluate_cart(u, roundtrip)
        if (maxval(abs(roundtrip - [r, 0.0_dp, z])) > 1.0e-8_dp) &
            error stop 'POTATO Boozer inverse residual exceeds 1e-8 cm'
        previous_u = u
        have_previous_u = .true.
        x = [u(1), 0.0_dp, u(2)]
        call do_magfie_pert_amp(x, amplitude)
        amplitude = amplitude*exp(cmplx(0.0_dp, n_tor*u(3), dp))
    end subroutine boozer_perturbation
end module boozer_perturbation_mod
