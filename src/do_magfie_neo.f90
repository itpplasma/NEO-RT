! The optional provider shares NEO-RT's public field interface. NEO-2 owns
! Fourier interpolation and field equations; this module owns unit conversion,
! lifecycle and serialization of the provider's single retained global state.
module do_magfie_mod
    use iso_fortran_env, only: dp => real64
    use util, only: pi
    use neo_field_provider, only: neo_field_load, neo_field_prepare, &
        field_minor_radius_m, field_s_min, field_s_max
    use magfie_mod, only: magfie
    use neo_magfie, only: boozer_curr_tor_hat, boozer_curr_pol_hat, &
        boozer_curr_tor_hat_s, boozer_curr_pol_hat_s, boozer_psi_pr_hat, &
        boozer_iota, boozer_iota_s
    use neo_exchange, only: rt0
    use partpa_mod, only: bmod0
    implicit none
    real(dp), parameter :: sign_theta = -1.0_dp
    real(dp), parameter, private :: tesla_to_gauss = 1.0e4_dp
    real(dp), parameter, private :: metre_to_cm = 100.0_dp
    real(dp), parameter, private :: reference_bmod_t = 1.0e-4_dp
    private :: copy_field_metadata
    logical, parameter :: hcurl_toroidal_is_available = .false.
    real(dp) :: s = 0.0_dp, psi_pr = 0.0_dp, Bthcov = 0.0_dp, Bphcov = 0.0_dp
    real(dp) :: dBthcovds = 0.0_dp, dBphcovds = 0.0_dp
    real(dp) :: q = 0.0_dp, dqds = 0.0_dp, iota = 0.0_dp
    real(dp) :: R0 = 0.0_dp, a = 0.0_dp, eps = 0.0_dp, bfac = 1.0_dp
    integer :: inp_swi = 0
    real(dp), private :: prepared_s = -huge(1.0_dp)
    logical, private :: field_loaded = .false.
    !$omp threadprivate(s, psi_pr, Bthcov, Bphcov, dBthcovds, dBphcovds)
    !$omp threadprivate(q, dqds, iota, eps)
contains
    subroutine set_s(value)
        real(dp), intent(in) :: value
        s = value
    end subroutine set_s

    subroutine magfie_thread_init()
        ! NEO-2's retained state is serialized below; no per-thread allocation.
    end subroutine magfie_thread_init

    subroutine read_boozer_file(path)
        use neo_input, only: ixn
        character(*), intent(in) :: path
        !$omp critical(neort_neo_field)
        call neo_field_load(path, inp_swi)
        if (any(ixn /= 0)) error stop "NEO-RT requires an axisymmetric background"
        bmod0 = reference_bmod_t ! provider hats -> Gauss, as in the historic adapter
        a = metre_to_cm*field_minor_radius_m
        R0 = metre_to_cm*rt0
        prepared_s = -huge(1.0_dp)
        field_loaded = .true.
        !$omp end critical(neort_neo_field)
    end subroutine read_boozer_file

    logical function s_outside_data(value)
        use, intrinsic :: ieee_arithmetic, only: ieee_is_finite
        real(dp), intent(in) :: value
        s_outside_data = .true.
        if (.not. ieee_is_finite(value)) return
        s_outside_data = value < field_s_min .or. value > field_s_max
    end function s_outside_data

    subroutine init_magfie_at_s()
        real(dp) :: bmin, bmax, jacobian, x(3), bder(3), hcovar(3), hctrvr(3), hcurl(3)
        if (s_outside_data(s)) error stop "s outside Boozer file radial range"
        x = [s, 0.0_dp, 0.0_dp]
        call do_magfie(x, bmin, jacobian, bder, hcovar, hctrvr, hcurl)
        x(3) = pi
        call do_magfie(x, bmax, jacobian, bder, hcovar, hctrvr, hcurl)
        eps = (bmax - bmin)/(bmax + bmin)
    end subroutine init_magfie_at_s

    subroutine do_magfie_init(path)
        character(*), intent(in) :: path
        call read_boozer_file(path)
        call init_magfie_at_s()
    end subroutine do_magfie_init

    subroutine do_magfie(x, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: bmod, sqrtg
        real(dp), intent(out) :: bder(size(x)), hcovar(size(x))
        real(dp), intent(out) :: hctrvr(size(x)), hcurl(size(x))
        real(dp) :: position(3)
        position = x(1:3)
        position(1) = max(field_s_min, min(field_s_max, position(1)))
        !$omp critical(neort_neo_field)
        if (.not. field_loaded) error stop "NEO-2 field evaluation before read"
        if (position(1) /= prepared_s) then
            call neo_field_prepare(position(1))
            prepared_s = position(1)
        end if
        call magfie(position, bmod, sqrtg, bder, hcovar, hctrvr, hcurl)
        bmod = tesla_to_gauss*bmod*bfac
        call copy_field_metadata()
        ! Keep the same explicit radial-field approximation as standalone.
        hcovar(1) = 0.0_dp
        hcurl(2) = 0.0_dp
        !$omp end critical(neort_neo_field)
    end subroutine do_magfie
    subroutine copy_field_metadata()
        psi_pr = boozer_psi_pr_hat*bfac
        Bthcov = boozer_curr_tor_hat*bfac
        Bphcov = boozer_curr_pol_hat*bfac
        dBthcovds = boozer_curr_tor_hat_s*bfac
        dBphcovds = boozer_curr_pol_hat_s*bfac
        iota = boozer_iota
        q = 1.0_dp/iota
        dqds = -boozer_iota_s/iota**2
    end subroutine copy_field_metadata
end module do_magfie_mod

module do_magfie_pert_mod
    use iso_fortran_env, only: dp => real64
    use neo_magfie_perturbation, only: neo_load_pert, neo_magfie_pert_amplitude
    use ntv_mod, only: m_phi
    use do_magfie_mod, only: s, bfac, inp_swi
    implicit none
    integer :: inp_swi_pert = -1, mph = 0
    integer, private :: mph_shared = 0
    logical, private :: perturbation_loaded = .false.
    !$omp threadprivate(mph)
contains
    subroutine set_mph(value)
        integer, intent(in) :: value
        mph_shared = value
        mph = value
    end subroutine set_mph

    subroutine init_mph_from_shared()
        mph = mph_shared
    end subroutine init_mph_from_shared

    subroutine read_boozer_pert_file(path)
        use neo_control, only: lab_swi
        use partpa_mod, only: bmod0
        character(*), intent(in) :: path
        integer :: input_switch
        input_switch = inp_swi_pert
        if (input_switch < 0) input_switch = inp_swi
        !$omp critical(neort_neo_field)
        ! Match the background provider's chart convention even when the
        ! perturbation is loaded first (the legacy reader otherwise STOPs).
        lab_swi = 10
        bmod0 = 1.0e-4_dp ! same Gauss normalization as the background loader
        call neo_load_pert(path, input_switch)
        mph_shared = m_phi
        perturbation_loaded = .true.
        !$omp end critical(neort_neo_field)
        call init_mph_from_shared()
    end subroutine read_boozer_pert_file

    subroutine init_magfie_pert_at_s()
        real(dp) :: x(3)
        complex(dp) :: amplitude
        call init_mph_from_shared()
        x = [s, 0.0_dp, 0.0_dp]
        call do_magfie_pert_amp(x, amplitude)
    end subroutine init_magfie_pert_at_s

    subroutine do_magfie_pert_init(path)
        character(*), intent(in) :: path
        call read_boozer_pert_file(path)
        call init_magfie_pert_at_s()
    end subroutine do_magfie_pert_init

    subroutine do_magfie_pert_amp(x, amplitude)
        real(dp), intent(in) :: x(:)
        complex(dp), intent(out) :: amplitude
        integer :: input_switch
        input_switch = inp_swi_pert
        if (input_switch < 0) input_switch = inp_swi
        !$omp critical(neort_neo_field)
        if (.not. perturbation_loaded) error stop "NEO-2 perturbation before read"
        call neo_magfie_pert_amplitude(x, amplitude)
        !$omp end critical(neort_neo_field)
        ! Type 8 is NEO-RT's cosine envelope. Type 9 retains the signed-m
        ! complex Fourier coefficient. NEO-2 uses bmod0=1e-4, hence Gauss.
        if (input_switch == 8) amplitude = cmplx(real(amplitude, dp), 0.0_dp, dp)
        amplitude = amplitude*bfac
    end subroutine do_magfie_pert_amp

    subroutine do_magfie_pert(x, bmod)
        use util, only: imun
        real(dp), intent(in) :: x(:)
        real(dp), intent(out) :: bmod
        complex(dp) :: amplitude
        call do_magfie_pert_amp(x, amplitude)
        bmod = real(amplitude*exp(imun*mph*x(2)), dp)
    end subroutine do_magfie_pert
end module do_magfie_pert_mod
