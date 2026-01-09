
module collis_alp
    use iso_fortran_env, only: dp => real64

    implicit none

    integer, parameter :: nsorts = 3, ns = 10000
    logical :: swcoll = .false.
    real(dp), dimension(nsorts) :: efcolf, velrat, enrat

    !$omp threadprivate (swcoll, efcolf, velrat, enrat)

contains

    subroutine coleff(p, dpp, dhh, fpeff)
    !
    !  Computes local values of dimensionless contravariant components
    !  of collisional diffusion tensor and friction force for nonrelativistic
    !  plasma. Backgound temperature is the same for all sorts.
    !
    !     Input variables:
    !        formal: p      - dimensionless momentum module (p/(sqrt(2)*p_T)
    !        common: efcolf - dmls collision frequencies
    !                velrat - ratio of test species thermal velocity to
    !                         background species thermal velocity
    !     Output variables:
    !        formal: dpp    - dimensionless momentum module diffusion
    !                         coefficient
    !                dhh    - dimensionless pitch angle diffusion coeff.
    !                fpeff  - effective dimensionless drag force (prop. to linear
    !                         deviation in Fokker-Planck eq.)

        implicit none

        integer i
        real(dp), intent(in) :: p
        real(dp), intent(out) :: dpp, dhh, fpeff
        real(dp) :: plim, xbeta, d_p, dh, dpd

        plim = max(p, 1.0e-8_dp)

        dpp = 0.0_dp
        dhh = 0.0_dp
        fpeff = 0.0_dp

        do i = 1, nsorts
            xbeta = p*velrat(i)

            call onseff(xbeta, d_p, dh, dpd)

            dpp = dpp + d_p*efcolf(i)
            dhh = dhh + dh*efcolf(i)
            fpeff = fpeff + (dpd/plim - 2.0_dp*d_p*p*enrat(i))*efcolf(i)
        end do

        dhh = dhh/plim**2

        return
    end subroutine coleff

    subroutine onseff(v, d_p, dh, dpd)
    !  d_p - dimensionless dpp
    !  dh - dhh*p^2     (p - dmls)
    !  dpd - (1/p)(d/dp)p^2*d_p   (p - dmls)

        implicit none

        ! square root of pi
        real(dp), parameter :: sqp = 1.7724538_dp
        ! cons=4./(3.*sqrt(pi))
        real(dp), parameter :: cons = 0.75225278_dp
        real(dp), intent(in) :: v
        real(dp), intent(out) :: d_p, dh, dpd
        real(dp) :: v2, v3, ex, er

        v2 = v**2
        v3 = v2*v
        if (v < 0.01_dp) then
            d_p = cons*(1.0_dp - 0.6_dp*v2)
            dh = cons*(1.0_dp - 0.2_dp*v2)
            dpd = 2.0_dp*cons*(1.0_dp - 1.2_dp*v2)
        elseif (v > 6.0_dp) then
            d_p = 1.0_dp/v3
            dh = (1.0_dp - 0.5_dp/v2)/v
            dpd = -1.0_dp/v3
        else
            ex = exp(-v2)/sqp
            er = erf(v)
            d_p = er/v3 - 2.0_dp*ex/v2
            dh = er*(1.0_dp - 0.5_dp/v2)/v + ex/v2
            dpd = 4.0_dp*ex - d_p
        end if

        return
    end subroutine onseff

    subroutine loacol_nbi(amb, am1, am2, Zb, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, ebeam, &
                          v0, dchichi, slowrate, dchichi_norm, slowrate_norm)
    !
    !   Performs precomputation of the constants for Coulomb collision
    !   operator for alpha-particles colliding with 2 sorts of ions and electrons
    !
    !   Normalisation: test particle (alpha) velocity is normalized by v0, alpha-particle
    !   birth velocity, time is multiplied with v0 and has a meaning of free path of alpha
    !   particle with v0. Respectively, collision frequencies have the meaning of inverse
    !   mean free paths.
    !
    !   Input variables:
    !        formal: am1,am2       - mass numbers of the first and second background ion species
    !                Z1,Z2         - charge numbers of these species
    !                densi1,densi2 - densities of these species, 1/cm**3
    !                tempi1,tempi2,tempe - temperatures of two ion species and electrons, eV
    !                ebeam         - initial beam particle energy, eV
    !   Output variables:
    !        formal: v0            - initial alpha particle velocity, cm/s
    !                dchichi       - pitch angle scattering frequency, $D^{\chi\chi}$, of alpha
    !                                particle with initial velocity, 1/s
    !                slowrate      - slowing down rate, $F^v / v_0$, of alpha particle with
    !                                initial velocity, 1/s
    !                dchichi_norim - normalized pitch angle scattering frequency, 1/cm
    !                slowrate_norm - normalized slowing down rate, 1/cm
    !        module collis_alp:
    !                efcolf - normalized collision frequencies
    !                velrat - ratio of initial alpha particle velocity v0 to the
    !                         specific background particle thermal velocity $v_{t}=\sqrt(2T/m)$
    !                enrat  - ratio of initial alpha particle energy to the background species
    !                         energy

        implicit none

        real(dp) :: amb, am1, am2, Zb, Z1, Z2, densi1, densi2, tempi1, tempi2, tempe, ebeam, dense
        real(dp) :: v0, dchichi, slowrate, dchichi_norm, slowrate_norm, vti1, vti2, vte
        real(dp) :: pi, pmass, emass, e, ev, alame, frecol_base, alami1, alami2

        pi = 3.14159265358979_dp
        pmass = 1.6726e-24_dp
        emass = 9.1094e-28_dp
        e = 4.8032e-10_dp
        ev = 1.6022e-12_dp

        enrat(1) = ebeam/tempi1
        enrat(2) = ebeam/tempi2
        enrat(3) = ebeam/tempe

        v0 = sqrt(2.0_dp*ebeam*ev/(amb*pmass))
        vti1 = sqrt(2.0_dp*tempi1*ev/(pmass*am1))
        vti2 = sqrt(2.0_dp*tempi2*ev/(pmass*am2))
        vte = sqrt(2.0_dp*tempe*ev/emass)

        velrat(1) = v0/vti1
        velrat(2) = v0/vti2
        velrat(3) = v0/vte

        dense = densi1*Z1 + densi2*Z2
        alami1 = 23.0_dp - log(max(epsilon(1.0_dp), &
                                 sqrt(densi1*Z1**2/tempi1)*Zb*Z1*(amb + am1)/(amb*tempi1 + am1*ebeam)))
        alami2 = 23.0_dp - log(max(epsilon(1.0_dp), &
                                 sqrt(densi2*Z2**2/tempi2)*Zb*Z2*(amb + am2)/(amb*tempi2 + am2*ebeam)))
        alame = 24.0_dp - log(sqrt(dense)/tempe)
        frecol_base = 2.0_dp*pi*dense*e**4*Zb**2/((amb*pmass)**2*v0**3) ! usual
        frecol_base = frecol_base/v0                                  ! normalized

        efcolf(1) = frecol_base*Z1**2*alami1*densi1/dense
        efcolf(2) = frecol_base*Z2**2*alami2*densi2/dense
        efcolf(3) = frecol_base*alame

        efcolf = efcolf*velrat
    end subroutine loacol_nbi

end module collis_alp
