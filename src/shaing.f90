! Shaing superbanana-plateau (SBP) analytic model.
!
! Resurrected and modernized from the historical src/shaing.f90 (removed in
! 2025-02-23, commits 54aaf57 and 651acf7).  Implements the bounce-averaged
! toroidal drift frequency of Shaing 2009 PPCF 51 035009 Eq. (8),
!
!     Omph = Om_tE + <Om_tB>,
!     <Om_tB> = -(c*mu*B0)/(qi*chi') * eps' * (2*E(kappa)/K(kappa) - 1),
!
! which supplies the resonant toroidal drift frequency of the ell=0
! superbanana resonance for trapped particles.  See the "Superbanana plateau"
! subsection of doc/driftorbit.lyx for the g(kappa) normalization.
!
! Field and profile inputs are taken through the current module APIs
! (do_magfie_mod, driftorbit, neort_profiles, util); the 2015-era globals are
! not resurrected.
module shaing
    use iso_fortran_env, only: dp => real64
    use util, only: c, qi, mi
    use do_magfie_mod, only: eps, s, a, psi_pr, sign_theta
    use driftorbit, only: B0
    use neort_profiles, only: Om_tE

    implicit none

    private
    public :: kappa2, omph_shaing, comelp

contains

    pure function kappa2(eta, B_ref, epsilon) result(k2)
        ! Trapping parameter kappa^2 (doc/driftorbit.lyx eq:kappa2):
        !     kappa^2 = (1 - eta*B_ref*(1 - epsilon))/(2*eta*B_ref*epsilon)
        ! kappa^2 -> 0 for deeply trapped, kappa^2 -> 1 at the trapped-passing
        ! boundary.
        real(dp), intent(in) :: eta, B_ref, epsilon
        real(dp) :: k2

        k2 = (1.0_dp - eta*B_ref*(1.0_dp - epsilon))/(2.0_dp*eta*B_ref*epsilon)
    end function kappa2

    function omph_shaing(v, eta) result(omph)
        ! Bounce-averaged toroidal drift frequency, Shaing 2009 Eq. (8).
        ! Only meaningful in the trapped region (eta > etatp).
        real(dp), intent(in) :: v, eta
        real(dp) :: omph

        real(dp) :: kappa, k2, Kell, Eell, depsdr, chi_prime, mu, drift

        ! Clamp kappa^2 into the open interval (0, 1) so the complete elliptic
        ! integrals stay well defined at the trapped-region boundaries, where
        ! the analytic B0*(1 +- eps) model and the numeric Bmin/Bmax can differ
        ! by rounding.
        k2 = kappa2(eta, B0, eps)
        if (k2 < 1.0e-12_dp) k2 = 1.0e-12_dp
        if (k2 > 1.0_dp - 1.0e-12_dp) k2 = 1.0_dp - 1.0e-12_dp
        kappa = sqrt(k2)

        call comelp(kappa, Kell, Eell)

        ! Magnetic moment mu = m*v_perp^2/(2*B) = mi*eta*v^2/2 (eta = v_perp^2/(v^2*B)).
        mu = 0.5_dp*mi*eta*v**2

        ! Radial derivative of the inverse aspect ratio for a large-aspect-ratio
        ! circular flux surface: eps = r/R0 with r = a*sqrt(s), so
        ! deps/dr = eps/r = eps/(a*sqrt(s)) = 1/R0.
        depsdr = eps/(a*sqrt(s))

        ! Poloidal-flux derivative chi' = sign_theta*psi_pr, consistent with the
        ! reference drift frequency Om_tBref = c*mi*vth^2/(2*qi*chi') used in
        ! neort::check_magfie.
        chi_prime = sign_theta*psi_pr

        drift = -(c*mu*B0)/(qi*chi_prime)*depsdr*(2.0_dp*Eell/Kell - 1.0_dp)

        omph = Om_tE + drift
    end function omph_shaing

    pure subroutine comelp(hk, ck, ce)
        ! Complete elliptic integrals K(k) and E(k) of modulus hk in [0, 1].
        !
        ! Modernized from Zhang & Jin, "Computation of Special Functions"
        ! (Wiley, 1996), routine COMELP.  Incorporation permitted by the
        ! authors provided the copyright is acknowledged.
        real(dp), intent(in) :: hk
        real(dp), intent(out) :: ck, ce

        real(dp) :: pk, ak, bk, ae, be

        pk = 1.0_dp - hk*hk

        if (hk == 1.0_dp) then
            ck = 1.0e300_dp
            ce = 1.0_dp
        else
            ak = (((0.01451196212_dp*pk + 0.03742563713_dp)*pk &
                   + 0.03590092383_dp)*pk + 0.09666344259_dp)*pk &
                 + 1.38629436112_dp
            bk = (((0.00441787012_dp*pk + 0.03328355346_dp)*pk &
                   + 0.06880248576_dp)*pk + 0.12498593597_dp)*pk &
                 + 0.5_dp
            ck = ak - bk*log(pk)

            ae = (((0.01736506451_dp*pk + 0.04757383546_dp)*pk &
                   + 0.0626060122_dp)*pk + 0.44325141463_dp)*pk &
                 + 1.0_dp
            be = (((0.00526449639_dp*pk + 0.04069697526_dp)*pk &
                   + 0.09200180037_dp)*pk + 0.2499836831_dp)*pk
            ce = ae - be*log(pk)
        end if
    end subroutine comelp

end module shaing
