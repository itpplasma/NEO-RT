#include "neort_nonlin.h"
#include "neort_profiles.h"
#include "neort_driftorbit.h"
#include "neort_util.h"
#include <stdio.h>
#include <stdlib.h>

double omega_prime(double ux, double eta, const double bounceavg[ORBIT_NVAR],
                   double Omth, double dOmdv, double dOmdeta, double dOmdpph)
{
    double v = ux * prof_vth, mi = util_mi, qi = util_qi, c = NEORT_C;
    double ma = mi * v * mi * c / qi * eta;
    double mb = mi * v * v / 2.0 * mi * c / qi;
    double mc = mi / (2.0 * Omth) * v * (1.0 - eta * bounceavg[5]);
    double md = mi * v * v / 2.0 * Omth;
    double me = -(double)do_mth / do_mph;
    double mf = 1.0 / do_mph;
    double dvdJ = mb * me / (ma * md * mf - mb * mc * mf);
    double detadJ = ma * me / (mb * mc * mf - ma * md * mf);
    return dOmdv * dvdJ + dOmdeta * detadJ + do_mph * dOmdpph;
}

double nonlinear_attenuation(double ux, double eta, const double bounceavg[ORBIT_NVAR],
                             double Omth, double dOmthdv, double dOmthdeta, double Hmn2)
{
    (void)ux; (void)eta; (void)bounceavg; (void)Omth; (void)dOmthdv;
    (void)dOmthdeta; (void)Hmn2;
    if (!do_nonlin)
        return 1.0;
    /* nonlin=true path needs attenuation_factor (thetafun_inp.dat + polylag_3),
     * which is outside the golden gate (nonlin=false) and not yet ported. */
    fprintf(stderr, "nonlinear_attenuation: nonlin=true path not ported "
                    "(needs attenuation_factor lookup)\n");
    exit(1);
}
