! Collisional boundary-layer factor for the resonance quadrature.
!
! Optional closure that multiplies each resonance-root contribution of the
! collisionless (delta-function) quadrature in compute_transport_integral by a
! factor R in [0, ~1.7] that accounts for the finite collisional width of the
! resonance layer and its proximity to the trapped-domain boundaries.
!
! Derivation and verification: CAS-verified two-boundary Airy layer,
! project evidence repo review/ell0_layer/neort_layer_kernel.wl
! (7/7 CAS checks, 2026-07-20).  The scaled two-boundary layer equation is
!   i z G - G'' = 1 ,  reflecting (Neumann) deeply-trapped boundary at z=-d_dt,
!   absorbing (Dirichlet) trapped-passing seam at z=+d_tp,  R = Re<int G dz>/pi.
! R is tabulated versus the two boundary distances scaled by the layer width
! (R_TABLE below) and interpolated at run time.
!
! The NTVTOK detrapping enhancement nu_d/(2 eps) is NOT inserted by hand; it
! emerges from the eta-space geometry through the bounce-averaged pitch
! diffusion D_eta_ba and the relative layer width (see .wl section 1).
!
! Enable with the namelist logical collisional_layer = .true. (default off is
! bit-identical to the collisionless quadrature).
module neort_collisional_layer
    use iso_fortran_env, only: dp => real64
    use util, only: pi

    implicit none
    private

    public :: collisional_layer_factor
    public :: layer_ratio, layer_delta_eta, deflection_frequency
    public :: comp_ellint_k, comp_ellint_e

    ! Namelist switch (shared read-only configuration, mirrors driftorbit::supban).
    logical, public :: collisional_layer = .false.

    ! Electron mass (g).  util.f90 keeps only ion-scale constants, so the
    ! electron thermal velocity needs its own documented constant here.
    real(dp), parameter :: me = 9.109382e-28_dp

    ! --- Two-boundary layer ratio table R(d_dt/delta, d_tp/delta) ------------
    ! Generated from review/ell0_layer/layer_ratio_table.csv, itself produced by
    ! review/ell0_layer/neort_layer_kernel.wl (LinearSolve of the scaled layer
    ! equation on an 8x8 grid of scaled boundary distances 2^{-1..6}).
    ! First index: d_dt/delta (deeply-trapped, reflecting side).
    ! Second index: d_tp/delta (trapped-passing, absorbing side).
    integer, parameter :: NGRID = 8

    ! Grid nodes are log2-uniform: {0.5, 1, 2, 4, 8, 16, 32, 64}.
    real(dp), parameter :: LOG2_MIN = -1.0_dp, LOG2_MAX = 6.0_dp
    real(dp), parameter :: LOG2_OF_2 = 0.6931471805599453_dp

    real(dp), parameter :: R_TABLE(NGRID, NGRID) = reshape([ &
        1.053948050794e-01_dp, 3.544311527926e-01_dp, 1.063785128436e+00_dp, &
        6.800688934433e-01_dp, 6.574541303283e-01_dp, 6.507516485138e-01_dp, &
        6.494987635793e-01_dp, 6.509365949275e-01_dp, &  ! d_dt/delta = 0.5
        2.879599612852e-01_dp, 6.023514061161e-01_dp, 1.741294514705e+00_dp, &
        1.611866899033e+00_dp, 1.580619620580e+00_dp, 1.574183018069e+00_dp, &
        1.573460458078e+00_dp, 1.575928545181e+00_dp, &  ! d_dt/delta = 1
        2.372962674487e-01_dp, 3.838590065884e-01_dp, 9.947134200297e-01_dp, &
        9.891541412040e-01_dp, 9.547558484416e-01_dp, 9.475536929149e-01_dp, &
        9.453121161076e-01_dp, 9.447590265197e-01_dp, &  ! d_dt/delta = 2
        2.614473790623e-01_dp, 4.405375043125e-01_dp, 1.090062658789e+00_dp, &
        1.043967231877e+00_dp, 1.010773845016e+00_dp, 1.003478736500e+00_dp, &
        1.001049442976e+00_dp, 1.000120471143e+00_dp, &  ! d_dt/delta = 4
        2.626083475719e-01_dp, 4.411038131436e-01_dp, 1.089662794289e+00_dp, &
        1.043854867353e+00_dp, 1.010668737630e+00_dp, 1.003400221074e+00_dp, &
        1.001024156048e+00_dp, 1.000201758308e+00_dp, &  ! d_dt/delta = 8
        2.628079030285e-01_dp, 4.413033120399e-01_dp, 1.089861684066e+00_dp, &
        1.044055035228e+00_dp, 1.010869401259e+00_dp, 1.003602099551e+00_dp, &
        1.001228442690e+00_dp, 1.000410856873e+00_dp, &  ! d_dt/delta = 16
        2.628303243460e-01_dp, 4.413255903606e-01_dp, 1.089881716365e+00_dp, &
        1.044078619773e+00_dp, 1.010892751389e+00_dp, 1.003625707650e+00_dp, &
        1.001252439560e+00_dp, 1.000435560611e+00_dp, &  ! d_dt/delta = 32
        2.628278437524e-01_dp, 4.413225582707e-01_dp, 1.089870071366e+00_dp, &
        1.044080500378e+00_dp, 1.010893556594e+00_dp, 1.003626946623e+00_dp, &
        1.001253953035e+00_dp, 1.000437228250e+00_dp], &  ! d_dt/delta = 64
        [NGRID, NGRID], order=[2, 1])

contains

    function collisional_layer_factor(v, eta, abs_dOmega_deta) result(ratio)
        ! Layer factor R for one resonance root at pitch eta and velocity v.
        ! Passing particles have no trapped-domain boundaries: R = 1 (documented
        ! limitation; the passing/loss coupling is not modelled here).
        use do_magfie_mod, only: eps
        use driftorbit, only: etatp, etadt

        real(dp), intent(in) :: v            ! particle velocity (cm/s)
        real(dp), intent(in) :: eta          ! resonance pitch (1-xi^2)/B at root
        real(dp), intent(in) :: abs_dOmega_deta  ! |dOmega/deta| at root = |eta_res(2)|
        real(dp) :: ratio
        real(dp) :: nu_d, delta_eta, d_tp, d_dt

        ratio = 1.0_dp
        if (eta <= etatp) return  ! passing particle

        nu_d = deflection_frequency(v)
        delta_eta = layer_delta_eta(nu_d, eta, eps, etatp, etadt, abs_dOmega_deta)
        if (delta_eta <= 0.0_dp) return  ! degenerate width -> collisionless limit

        d_tp = eta - etatp
        d_dt = etadt - eta
        ratio = layer_ratio(d_dt/delta_eta, d_tp/delta_eta)
    end function collisional_layer_factor

    pure function layer_delta_eta(nu_d, eta, eps, etatp, etadt, abs_dOmega_deta) &
        result(delta_eta)
        ! Collisional layer width in eta-space,
        !   delta_eta = (D_eta_ba / |dOmega/deta|)^(1/3),
        ! with the large-aspect closed-form bounce-averaged pitch diffusion
        !   D_eta_ba = 2 nu_d eta <(1 - eta B)/B>_b
        !            = 2 nu_d eta (2 eps eta) [E(m)/K(m) - (1 - m)],
        !   m = (eta - etatp)/(etadt - etatp).
        ! The bounce average is evaluated in the large-aspect-ratio limit
        ! (derivation review/ell0_layer/neort_layer_kernel.wl section 1); this is
        ! consistent with the large-aspect two-boundary problem defining R and
        ! avoids adding an extra integrand to the orbit bounce-average machinery.
        real(dp), intent(in) :: nu_d, eta, eps, etatp, etadt, abs_dOmega_deta
        real(dp) :: delta_eta
        real(dp) :: m, shape_factor, D_eta_ba
        real(dp), parameter :: M_MIN = 1.0e-6_dp, M_MAX = 1.0_dp - 1.0e-6_dp
        real(dp), parameter :: ONE_THIRD = 1.0_dp/3.0_dp

        m = (eta - etatp)/(etadt - etatp)
        m = min(max(m, M_MIN), M_MAX)
        shape_factor = comp_ellint_e(m)/comp_ellint_k(m) - (1.0_dp - m)
        D_eta_ba = 2.0_dp*nu_d*eta*(2.0_dp*eps*eta)*shape_factor

        if (D_eta_ba <= 0.0_dp .or. abs_dOmega_deta <= 0.0_dp) then
            delta_eta = 0.0_dp
            return
        end if
        delta_eta = (D_eta_ba/abs_dOmega_deta)**ONE_THIRD
    end function layer_delta_eta

    pure function layer_ratio(d_dt_scaled, d_tp_scaled) result(ratio)
        ! Bilinear interpolation of R_TABLE in log2 of the scaled boundary
        ! distances.  Outside the grid the scaled distance is clamped to the
        ! nearest node, so the far-far corner returns R ~ 1 and near-boundary
        ! roots return the tabulated suppression.
        real(dp), intent(in) :: d_dt_scaled, d_tp_scaled
        real(dp) :: ratio
        real(dp) :: f1, f2
        integer :: i1, i2

        call locate_log2(d_dt_scaled, i1, f1)
        call locate_log2(d_tp_scaled, i2, f2)
        ratio = (1.0_dp - f1)*(1.0_dp - f2)*R_TABLE(i1, i2) &
                + f1*(1.0_dp - f2)*R_TABLE(i1 + 1, i2) &
                + (1.0_dp - f1)*f2*R_TABLE(i1, i2 + 1) &
                + f1*f2*R_TABLE(i1 + 1, i2 + 1)
    end function layer_ratio

    pure subroutine locate_log2(scaled_distance, idx, frac)
        ! Locate a scaled distance on the log2-uniform grid, clamped to
        ! [LOG2_MIN, LOG2_MAX].  Returns lower node index idx (1..NGRID-1) and
        ! interpolation weight frac in [0, 1].
        real(dp), intent(in) :: scaled_distance
        integer, intent(out) :: idx
        real(dp), intent(out) :: frac
        real(dp), parameter :: TINY_DISTANCE = 1.0e-300_dp
        real(dp) :: log2_x, t

        log2_x = log(max(scaled_distance, TINY_DISTANCE))/LOG2_OF_2
        log2_x = min(max(log2_x, LOG2_MIN), LOG2_MAX)
        t = log2_x - LOG2_MIN  ! in [0, NGRID-1]
        idx = min(int(t) + 1, NGRID - 1)
        frac = t - real(idx - 1, dp)
    end subroutine locate_log2

    function deflection_frequency(v) result(nu_d)
        ! Energy-dependent deuterium pitch-angle deflection frequency (1/s) from
        ! both background species (deuterium self-collisions and electrons),
        ! Helander-Sigmar form (Gaussian units):
        !   nu_D(v) = sum_b nu_hat_ab [erf(x_b) - G(x_b)] / x_a^3,
        !   nu_hat_ab = 4 pi n_b (Z_a e)^2 (Z_b e)^2 lnLambda / (m_a^2 v_Ta^3),
        !   x_a = v/v_th (deuteron), x_b = v/v_Tb, G = Chandrasekhar function.
        ! Equivalent to NTVTOK's nu_d = sum_b nu_ab G(x_b)/x^3
        ! (Sun et al., PoP 26, 072504 (2019)).
        use util, only: qe, mi, qi, ev
        use neort_profiles, only: vth, ni1, ni2, Ti1, Te
        use neort_profiles, only: Z1_global, Z2_global

        real(dp), intent(in) :: v
        real(dp) :: nu_d
        real(dp) :: Z1, ne, vthe, ln_lambda, nu_hat_prefac
        real(dp) :: nu_hat_D, nu_hat_e, xa, xe, chan_D, chan_e

        ! Without loaded plasma profiles (e.g. torque-free runs) collisionality
        ! is undefined; treat the layer as inert (nu_d = 0 -> R = 1).
        nu_d = 0.0_dp
        if (ni1 <= 0.0_dp .or. Te <= 0.0_dp .or. vth <= 0.0_dp) return

        Z1 = qi/qe                       ! deuteron charge number
        ne = Z1_global*ni1 + Z2_global*ni2  ! electron density from quasineutrality
        vthe = sqrt(2.0_dp*Te*ev/me)

        ! NRL/Helander electron Coulomb logarithm (n in cm^-3, T in eV); the same
        ! form NEO-RT already uses in collis_nbi.f90 (loacol_nbi: alame).
        ln_lambda = max(1.0_dp, 24.0_dp - log(sqrt(ne)/Te))

        nu_hat_prefac = 4.0_dp*pi*qe**4*ln_lambda/(mi**2*vth**3)
        nu_hat_D = nu_hat_prefac*ni1*Z1**4      ! Z_a^2 Z_b^2, both deuterons
        nu_hat_e = nu_hat_prefac*ne*Z1**2       ! Z_a^2 * Z_e^2, Z_e^2 = 1

        xa = v/vth
        xe = v/vthe
        chan_D = erf(xa) - chandrasekhar(xa)
        chan_e = erf(xe) - chandrasekhar(xe)

        nu_d = (nu_hat_D*chan_D + nu_hat_e*chan_e)/xa**3
    end function deflection_frequency

    pure elemental function chandrasekhar(x) result(g)
        ! Chandrasekhar function G(x) = [erf(x) - x erf'(x)] / (2 x^2),
        ! erf'(x) = (2/sqrt(pi)) exp(-x^2).  Small-x series avoids 0/0.
        real(dp), intent(in) :: x
        real(dp) :: g
        real(dp), parameter :: two_over_sqrt_pi = 1.1283791670955126_dp
        real(dp), parameter :: X_SMALL = 1.0e-4_dp
        real(dp) :: erf_prime

        if (abs(x) < X_SMALL) then
            g = two_over_sqrt_pi*x/3.0_dp  ! G(x) -> 2x/(3 sqrt(pi))
            return
        end if
        erf_prime = two_over_sqrt_pi*exp(-x*x)
        g = (erf(x) - x*erf_prime)/(2.0_dp*x*x)
    end function chandrasekhar

    pure elemental function comp_ellint_k(m) result(k)
        ! Complete elliptic integral K(m), parameter (m = k^2) convention.
        ! Abramowitz & Stegun 17.3.34, |error| < 3e-8 for m in [0, 1).
        real(dp), intent(in) :: m
        real(dp) :: k
        real(dp) :: m1, log_term
        real(dp), parameter :: a0 = 1.38629436112_dp, a1 = 0.09666344259_dp, &
            a2 = 0.03590092383_dp, a3 = 0.03742563713_dp, a4 = 0.01451196212_dp
        real(dp), parameter :: b0 = 0.5_dp, b1 = 0.12498593597_dp, &
            b2 = 0.06880248576_dp, b3 = 0.03328355346_dp, b4 = 0.00441787012_dp

        m1 = 1.0_dp - m
        log_term = -log(m1)
        k = (a0 + m1*(a1 + m1*(a2 + m1*(a3 + m1*a4)))) &
            + log_term*(b0 + m1*(b1 + m1*(b2 + m1*(b3 + m1*b4))))
    end function comp_ellint_k

    pure elemental function comp_ellint_e(m) result(e)
        ! Complete elliptic integral E(m), parameter (m = k^2) convention.
        ! Abramowitz & Stegun 17.3.36, |error| < 2e-8 for m in [0, 1).
        real(dp), intent(in) :: m
        real(dp) :: e
        real(dp) :: m1, log_term
        real(dp), parameter :: a1 = 0.44325141463_dp, a2 = 0.06260601220_dp, &
            a3 = 0.04757383546_dp, a4 = 0.01736506451_dp
        real(dp), parameter :: b1 = 0.24998368310_dp, b2 = 0.09200180037_dp, &
            b3 = 0.04069697526_dp, b4 = 0.00526449639_dp

        m1 = 1.0_dp - m
        log_term = -log(m1)
        e = (1.0_dp + m1*(a1 + m1*(a2 + m1*(a3 + m1*a4)))) &
            + log_term*(m1*(b1 + m1*(b2 + m1*(b3 + m1*b4))))
    end function comp_ellint_e

end module neort_collisional_layer
