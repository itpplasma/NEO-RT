module neort_nonlin
    use iso_fortran_env, only: dp => real64
    use logger, only: debug, trace, get_log_level, LOG_TRACE
    use util, only: c, mi, pi
    use collis_alp, only: coleff
    use neort_freq, only: Om_th, Om_ph, d_Om_ds
    use neort_orbit, only: nvar
    use do_magfie_mod, only: sign_theta
    use driftorbit, only: nonlin, vth, mth, mph, qi, iota, psi_pr

    implicit none

contains

    function nonlinear_attenuation(ux, eta, bounceavg, Omth, dOmthdv, dOmthdeta, Hmn2)
        real(dp), intent(in) :: ux, eta, bounceavg(nvar), Omth, dOmthdv, dOmthdeta, Hmn2
        real(dp) :: nonlinear_attenuation

        real(dp) :: dpp, dhh, fpeff, dres, dnorm, Omph, dOmphdv, dOmphdeta, dOmdv, &
                   dOmdeta, Ompr, dOmphds, dOmthds, dOmdpph, v

        call trace('nonlinear_attenuation')

        if (.not. nonlin) then
            nonlinear_attenuation = 1.0_dp
            return
        end if

        v = ux * vth

        if (nonlin) then
            call Om_ph(v, eta, Omph, dOmphdv, dOmphdeta)
            call d_Om_ds(v, eta, 2.0_dp * pi / Omth, dOmthds, dOmphds)
            dOmdv = mth * dOmthdv + mph * dOmphdv
            dOmdeta = mth * dOmthdeta + mph * dOmphdeta
            dOmdpph = -(qi / c * iota * sign_theta * psi_pr)**(-1) * (mth * dOmthds + mph * dOmphds)
            Ompr = omega_prime(ux, eta, bounceavg, Omth, dOmdv, dOmdeta, dOmdpph)
            call coleff(ux, dpp, dhh, fpeff)
            dhh = vth * dhh
            dpp = vth**3 * dpp
            dres = dpp * (dOmdv / Ompr)**2 + dhh * eta * (bounceavg(5) - eta) * (dOmdeta / Ompr)**2
            dnorm = dres * sqrt(abs(Ompr)) / sqrt(abs(Hmn2))**(3.0_dp / 2.0_dp)
            call attenuation_factor(dnorm, nonlinear_attenuation)
        end if

        call trace('nonlinear_attenuation complete')
    end function nonlinear_attenuation

    pure function omega_prime(ux, eta, bounceavg, Omth, dOmdv, dOmdeta, dOmdpph)
        real(dp), intent(in) :: ux, eta, bounceavg(nvar), Omth, dOmdv, dOmdeta, dOmdpph
        real(dp) :: omega_prime

        real(dp) :: ma, mb, mc, md, me, mf, dvdJ, detadJ

        ma = mi * (ux * vth) * mi * c / qi * eta
        mb = mi * (ux * vth)**2 / 2 * mi * c / qi
        mc = mi / (2.0_dp * Omth) * (ux * vth) * (1.0_dp - eta * bounceavg(6))
        md = mi * (ux * vth)**2 / 2.0_dp * Omth
        me = -dble(mth) / mph
        mf = 1.0_dp / mph

        dvdJ = mb * me / (ma * md * mf - mb * mc * mf)
        detadJ = ma * me / (mb * mc * mf - ma * md * mf)

        omega_prime = dOmdv * dvdJ + dOmdeta * detadJ + mph * dOmdpph
    end function omega_prime

end module neort_nonlin
