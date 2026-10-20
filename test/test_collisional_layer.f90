program test_collisional_layer
    ! Unit tests for the collisional boundary-layer factor:
    !   (1) complete elliptic integrals against known values,
    !   (2) independently checked interpolation, rejection and asymptotic limits,
    !   (3) physics: an interior resonance root approaches the collisionless
    !       (delta) result as collisionality -> 0 and is suppressed as
    !       the narrow-domain asymptote at high collisionality.
    use iso_fortran_env, only: dp => real64
    use neort_config, only: read_and_set_config
    use neort_collisional_layer, only: layer_ratio, layer_delta_eta, &
        comp_ellint_k, comp_ellint_e, deflection_frequency, collisional_layer

    implicit none

    character(len=16) :: mode
    character(len=64) :: argument
    real(dp) :: a, b
    integer :: ia, ib

    real(dp), parameter :: PI_HALF = 1.5707963267948966_dp

    call get_command_argument(1, mode)
    if (mode == 'config-reset') then
        collisional_layer = .true.
        call read_and_set_config('driftorbit_test.in')
        if (collisional_layer) error stop 'absent collisional_layer retained previous true value'
        stop
    end if
    if (mode == 'ratio') then
        call get_command_argument(2, argument)
        read(argument, *) a
        call get_command_argument(3, argument)
        read(argument, *) b
        print *, layer_ratio(a, b)
        stop
    end if
    if (mode == 'edges') then
        do ia = -16, 23
            a = 2.0_dp**((real(ia, dp) + 0.5_dp)/4.0_dp)
            print *, a, 0.0625_dp, layer_ratio(a, 0.0625_dp)
            print *, 0.0625_dp, a, layer_ratio(0.0625_dp, a)
            print *, a, 64.0_dp, layer_ratio(a, 64.0_dp)
            print *, 64.0_dp, a, layer_ratio(64.0_dp, a)
            print *, 0.0_dp, a, layer_ratio(0.0_dp, a)
            print *, a, 0.0_dp, layer_ratio(a, 0.0_dp)
            do ib = 1, 3
                b = real(ib, dp)/64.0_dp
                print *, a, b, layer_ratio(a, b)
                print *, b, a, layer_ratio(b, a)
            end do
        end do
        stop
    end if
    if (mode == 'table' .or. mode == 'interior') then
        do ia = -16, 24
            do ib = -16, 24
                a = 2.0_dp**(real(ia, dp)/4.0_dp)
                b = 2.0_dp**(real(ib, dp)/4.0_dp)
                if (mode == 'interior') then
                    if (ia == 24 .or. ib == 24) cycle
                    a = a*2.0_dp**0.125_dp
                    b = b*2.0_dp**0.125_dp
                end if
                print *, a, b, layer_ratio(a, b)
            end do
        end do
        if (mode == 'table') then
            do ia = -16, 24
                a = 2.0_dp**(real(ia, dp)/4.0_dp)
                print *, 0.0_dp, a, layer_ratio(0.0_dp, a)
                print *, a, 0.0_dp, layer_ratio(a, 0.0_dp)
            end do
            print *, 0.0_dp, 0.0_dp, layer_ratio(0.0_dp, 0.0_dp)
        else
            do ia = -16, 23
                a = 2.0_dp**((real(ia, dp) + 0.5_dp)/4.0_dp)
                print *, 0.03125_dp, a, layer_ratio(0.03125_dp, a)
                print *, a, 0.03125_dp, layer_ratio(a, 0.03125_dp)
            end do
            print *, 0.03125_dp, 0.03125_dp, layer_ratio(0.03125_dp, 0.03125_dp)
        end if
        stop
    end if

    call check('narrow domain independent solver', layer_ratio(0.1_dp, 0.1_dp), &
        0.00084882598725_dp, 1.0e-9_dp)
    call test_background_rates
    call test_mixture_rates
    call test_root_factor
    call test_bounce_diffusion
    call test_elliptic_integrals
    call test_table_interpolation
    call test_table_far_limit
    call test_interior_root_collisionless_limit
    call test_interior_root_suppression

    print *, "test_collisional_layer ... OK"

contains

    subroutine check(name, got, expected, tol)
        character(*), intent(in) :: name
        real(dp), intent(in) :: got, expected, tol

        if (abs(got - expected) > tol) then
            print *, "FAIL ", name, " got=", got, " expected=", expected
            error stop
        end if
    end subroutine check

    subroutine test_background_rates
        use neort_profiles, only: ni1, ni2, Ti1, Ti2, Te, vth, profile_species, &
            am1_global, am2_global, Z1_global, Z2_global
        use util, only: qi, mi, qe, mu
        real(dp), parameter :: speeds(0:2) = &
            [2297054857.3951044_dp, 31062098.147308018_dp, 11342274.559146686_dp]
        ! SI Coulomb prefactor and erf/G independently evaluated with SciPy;
        ! n_e=n_D+6*n_C, T_e=1500, T_D=1000, T_C=800 eV, lnLambda=16.111415468.
        real(dp), parameter :: expected(0:2) = &
            [60067.24793336906_dp, 1305.705509393793_dp, 16991.290847333712_dp]
        integer :: species
        ni1 = 1.0e13_dp
        ni2 = 1.0e12_dp
        Ti1 = 1000.0_dp
        Ti2 = 800.0_dp
        Te = 1500.0_dp
        am1_global = 2.0_dp
        am2_global = 12.0_dp
        Z1_global = 1.0_dp
        Z2_global = 6.0_dp
        do species = 0, 2
            profile_species = species
            vth = speeds(species)
            qi = qe
            mi = 2.0_dp*mu
            if (species == 0) then
                qi = -qe
                mi = 9.109382e-28_dp
            else if (species == 2) then
                qi = 6.0_dp*qe
                mi = 12.0_dp*mu
            end if
            call check("all-background species rate", &
                deflection_frequency(vth), expected(species), 2.0e-7_dp)
        end do
    end subroutine test_background_rates

    subroutine test_mixture_rates
        use neort_species, only: mixture_t, species_state_t, configure_mixture, &
            evaluate_mixture_at_s, selected_species, reset_mixture
        use util, only: mu
        type(mixture_t) :: mixture
        integer :: i
        real(dp), parameter :: speed(3) = &
            [2297054857.3951044_dp, 31062098.147308018_dp, 11342274.559146686_dp]
        real(dp), parameter :: expected(3) = &
            [63380.10567073069_dp, 1377.7184074095176_dp, 17928.402689277795_dp]
        mixture%coulomb_log = 17.0_dp
        allocate(mixture%species(3))
        mixture%species(1) = species_state_t('e', 9.109382e-31_dp, -1.0_dp, &
            1.6e19_dp, 1500.0_dp, 0.0_dp, 0.0_dp)
        mixture%species(2) = species_state_t('D', 2.0_dp*mu*0.001_dp, 1.0_dp, &
            1.0e19_dp, 1000.0_dp, 0.0_dp, 0.0_dp)
        mixture%species(3) = species_state_t('C', 12.0_dp*mu*0.001_dp, 6.0_dp, &
            1.0e18_dp, 800.0_dp, 0.0_dp, 0.0_dp)
        call configure_mixture(mixture)
        call evaluate_mixture_at_s(0.0_dp)
        do i = 1, 3
            selected_species = i
            call check('explicit mixture rate', deflection_frequency(speed(i)), &
                expected(i), 2.0e-6_dp)
        end do
        call reset_mixture()
    end subroutine test_mixture_rates

    subroutine test_root_factor
        use neort_collisional_layer, only: collisional_layer_factor
        real(dp) :: got, width
        ! D_eta = 4*nu*eps*eta^2*<m*cos(chi)^2>, verified below by quadrature.
        ! Choose |Omega'| so delta=.05; root distances are then grid node (2,1).
        width = layer_delta_eta(10.0_dp, 0.95_dp, 0.1_dp, 0.9_dp, 1.05_dp, 1.0_dp)
        got = collisional_layer_factor(10.0_dp, 0.95_dp, 0.1_dp, 0.9_dp, &
            1.05_dp, (width/0.05_dp)**3)
        call check('interior root maps to reflecting/absorbing axes', got, &
            layer_ratio(2.0_dp, 1.0_dp), 1.0e-12_dp)
        got = collisional_layer_factor(0.0_dp, 0.95_dp, 0.1_dp, 0.9_dp, 1.05_dp, 1.0_dp)
        call check('zero rate is exactly collisionless', got, 1.0_dp, 0.0_dp)
        got = collisional_layer_factor(10.0_dp, 0.8_dp, 0.1_dp, 0.9_dp, 1.05_dp, 1.0_dp)
        call check('passing factor is exactly one', got, 1.0_dp, 0.0_dp)
    end subroutine test_root_factor

    subroutine test_bounce_diffusion
        ! Independent bounce integral in chi, sin(theta/2)=sqrt(k2)*sin(chi).
        ! <(1-eta B)/B> = 2 eps eta <k2 cos(chi)^2> to leading order.
        real(dp) :: eta, k2, chi, weight, numerator, denominator, expected
        real(dp) :: got
        integer :: i, j
        integer, parameter :: n = 10000
        real(dp), parameter :: eps = 0.01_dp, etatp = 0.99_dp, etadt = 1.01_dp
        do j = 1, 3
            k2 = 0.2_dp*real(j, dp)
            eta = etadt - k2*(etadt - etatp)
            numerator = 0.0_dp
            denominator = 0.0_dp
            do i = 1, n
                chi = PI_HALF*(real(i, dp) - 0.5_dp)/real(n, dp)
                weight = 1.0_dp/sqrt(1.0_dp - k2*sin(chi)**2)
                numerator = numerator + weight*k2*cos(chi)**2
                denominator = denominator + weight
            end do
            expected = (4.0_dp*eps*eta**2*numerator/denominator)**(1.0_dp/3.0_dp)
            got = layer_delta_eta(1.0_dp, eta, eps, etatp, etadt, 1.0_dp)
            call check("bounce diffusion width", got, expected, 2.0e-8_dp)
        end do
    end subroutine test_bounce_diffusion

    subroutine test_elliptic_integrals
        ! K(0) = E(0) = pi/2; reference values at m = 0.5.
        call check("K(0)", comp_ellint_k(0.0_dp), PI_HALF, 1.0e-7_dp)
        call check("E(0)", comp_ellint_e(0.0_dp), PI_HALF, 1.0e-7_dp)
        call check("K(0.5)", comp_ellint_k(0.5_dp), 1.8540746773_dp, 1.0e-7_dp)
        call check("E(0.5)", comp_ellint_e(0.5_dp), 1.3506438810_dp, 1.0e-7_dp)
    end subroutine test_elliptic_integrals

    subroutine test_table_interpolation
        real(dp) :: r_lo, r_hi, r_mid

        ! Interpolation between two nodes stays between them (monotone edge).
        r_lo = layer_ratio(0.5_dp, 0.5_dp)
        r_hi = layer_ratio(1.0_dp, 0.5_dp)
        r_mid = layer_ratio(sqrt(2.0_dp)/2.0_dp, 0.5_dp)  ! log2 = -0.5
        if (r_mid < min(r_lo, r_hi) .or. r_mid > max(r_lo, r_hi)) then
            print *, "FAIL interpolation not bracketed", r_lo, r_mid, r_hi
            error stop
        end if
    end subroutine test_table_interpolation

    subroutine test_table_far_limit
        ! An interior root many layer widths from both boundaries recovers the
        ! collisionless plateau R -> 1.
        call check("far-far R -> 1", layer_ratio(64.0_dp, 64.0_dp), &
                   1.0_dp, 0.03_dp)
    end subroutine test_table_far_limit

    subroutine test_interior_root_collisionless_limit
        real(dp) :: delta_eta, ratio
        real(dp), parameter :: ETATP = 0.9_dp, ETADT = 1.1_dp, ETA = 1.0_dp
        real(dp), parameter :: EPS = 0.1_dp, ABS_DOM = 1.0e5_dp
        real(dp), parameter :: NU_TINY = 1.0e-20_dp

        ! Vanishing collisionality: layer width -> 0, root sits many widths from
        ! both boundaries, so R must approach the delta result (1) within 2%.
        delta_eta = layer_delta_eta(NU_TINY, ETA, EPS, ETATP, ETADT, ABS_DOM)
        ratio = layer_ratio((ETADT - ETA)/delta_eta, (ETA - ETATP)/delta_eta)
        call check("interior collisionless limit R -> 1", ratio, 1.0_dp, 0.02_dp)
    end subroutine test_interior_root_collisionless_limit

    subroutine test_interior_root_suppression
        real(dp) :: delta_eta, ratio
        real(dp), parameter :: ETATP = 0.9_dp, ETADT = 1.1_dp, ETA = 1.0_dp
        real(dp), parameter :: EPS = 0.1_dp, ABS_DOM = 1.0e5_dp
        real(dp), parameter :: NU_HUGE = 1.0e30_dp

        ! Huge collisionality: the independently validated narrow-domain limit.
        delta_eta = layer_delta_eta(NU_HUGE, ETA, EPS, ETATP, ETADT, ABS_DOM)
        ratio = layer_ratio((ETADT - ETA)/delta_eta, (ETA - ETATP)/delta_eta)
        call check('narrow-domain cubic suppression', ratio, &
            ((ETADT - ETATP)/delta_eta)**3/(3.0_dp*2.0_dp*PI_HALF), 1.0e-30_dp)
        if (ratio >= 0.6_dp) then
            print *, "FAIL suppression: R =", ratio, " not < 0.6"
            error stop
        end if
    end subroutine test_interior_root_suppression

end program test_collisional_layer
