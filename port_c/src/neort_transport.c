#include "neort_transport.h"
#include "neort_field.h"
#include "neort_profiles.h"
#include "neort_orbit.h"
#include "neort_freq.h"
#include "neort_resonance.h"
#include "neort_nonlin.h"
#include "neort_driftorbit.h"
#include "neort_util.h"
#include <complex.h>
#include <math.h>

/* Module-shared Omth: set by the Om_th call per resonance, read by
 * timestep_transport (mirrors the neort_transport module variable Omth). */
static double transport_Omth = 0.0;

static double D11int(double ux, double taub, double Hmn2)
{
    return pow(NEORT_PI, 1.5) * do_mph * do_mph * NEORT_C * NEORT_C * field_q * prof_vth /
           (util_qi * util_qi * do_dVds * fabs(field_psi_pr)) * ux * ux * ux *
           exp(-ux * ux) * taub * Hmn2;
}

static double D12int(double ux, double taub, double Hmn2)
{
    return D11int(ux, taub, Hmn2) * ux * ux;
}

static double Tphi_int(double ux, double taub, double Hmn2)
{
    double sgn = (field_psi_pr * field_q * FIELD_SIGN_THETA >= 0) ? 1.0 : -1.0;
    return sgn * pow(NEORT_PI, 1.5) * do_mph * do_mph * prof_ni1 * NEORT_C * prof_vth /
           util_qi * ux * ux * ux * exp(-ux * ux) * taub * Hmn2 *
           (prof_A1 + prof_A2 * ux * ux);
}

static void timestep_transport(double v, double eta, int neq, double t,
                               const double *y, double *ydot)
{
    (void)neq;
    double x[3] = {orbit_get_s(), 0.0, y[0]}, bmod, sqrtg, hder[3], hcovar[3],
           hctrvr[3], hcurl[3];
    do_magfie(x, &bmod, &sqrtg, hder, hcovar, hctrvr, hcurl);
    poloidal_velocity(v, eta, bmod, hctrvr[2], hder[2], y[1], ydot);

    double _Complex epsn;
    if (do_pertfile) {
        do_magfie_pert_amp(x, &epsn);
        epsn = do_epsmn * epsn / bmod;
    } else {
        epsn = do_epsmn * cexp(I * do_m0 * y[0]);
    }

    double _Complex Hn;
    if (eta > do_etatp) {
        double t0 = 0.0;
        Hn = (2.0 - eta * bmod) * epsn *
             cexp(I * (field_q * do_mph * y[0] - do_mth * (t - t0) * transport_Omth));
    } else {
        Hn = (2.0 - eta * bmod) * epsn *
             cexp(I * (field_q * do_mph * y[0] -
                       (do_mth + field_q * do_mph) * t * transport_Omth));
    }
    ydot[2] = creal(Hn);
    ydot[3] = cimag(Hn);
    if (do_nonlin) {
        ydot[4] = 1.0 / bmod;
        ydot[5] = bmod;
    } else {
        ydot[4] = 0.0;
        ydot[5] = 0.0;
    }
    ydot[6] = 0.0; /* y(7) unused downstream; keep derivative defined */
}

void compute_transport_integral(double vmin, double vmax, int vsteps, double D[2],
                                double *T)
{
    double dOmthdv, dOmthdeta;
    double roots[DO_NLEV * 3];
    int nroots;

    D[0] = 0.0;
    D[1] = 0.0;
    *T = 0.0;
    double du = (vmax - vmin) / (vsteps * prof_vth);
    double ux = vmin / prof_vth + du / 2.0;

    for (int ku = 1; ku <= vsteps; ku++) {
        double v = ux * prof_vth;
        driftorbit_coarse(v, do_etamin, do_etamax, roots, DO_NLEV, &nroots);
        for (int kr = 0; kr < nroots; kr++) {
            double eta_res[2];
            driftorbit_root(v, 1.0e-8 * fabs(prof_Om_tE), roots[kr * 3 + 0],
                            roots[kr * 3 + 1], eta_res);
            double eta = eta_res[0];

            Om_th(v, eta, &transport_Omth, &dOmthdv, &dOmthdeta);
            double taub = 2.0 * NEORT_PI / fabs(transport_Omth);
            double bounceavg[ORBIT_NVAR];
            int istate;
            bounce_fast_ext(v, eta, taub, bounceavg, timestep_transport, &istate);

            double Hmn2 = (bounceavg[2] * bounceavg[2] + bounceavg[3] * bounceavg[3]) *
                          pow(util_mi * (ux * prof_vth) * (ux * prof_vth) / 2.0, 2.0);
            double atten = nonlinear_attenuation(ux, eta, bounceavg, transport_Omth,
                                                 dOmthdv, dOmthdeta, Hmn2);
            double dD11 = du * D11int(ux, taub, Hmn2) / fabs(eta_res[1]);
            double dD12 = du * D12int(ux, taub, Hmn2) / fabs(eta_res[1]);
            D[0] += dD11 * atten;
            D[1] += dD12 * atten;
            if (do_comptorque) {
                double dT = du * Tphi_int(ux, taub, Hmn2) / fabs(eta_res[1]);
                *T += dT * atten;
            }
        }
        ux += du;
    }

    double D_plateau = NEORT_PI * pow(prof_vth, 3.0) /
                       (16.0 * field_R0 * field_iota *
                        pow(util_qi * do_B0 / (util_mi * NEORT_C), 2.0));
    double dsdreff = 2.0 / field_a * sqrt(orbit_get_s());
    D[0] = pow(dsdreff, -2.0) * D[0] / D_plateau;
    D[1] = pow(dsdreff, -2.0) * D[1] / D_plateau;
}
