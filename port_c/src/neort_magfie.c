#include "neort_magfie.h"
#include "neort_field.h"
#include "neort_driftorbit.h"
#include "neort_orbit.h"
#include "neort_util.h"
#include <math.h>

void init_flux_surface_average(double s)
{
    const int nth = 1000;
    double dth = 2.0 * NEORT_PI / nth;
    double x[3] = {s, 0.0, 0.0};
    double bmod, sqrtg, hder[3], hcovar[3], hctrvr[3], hcurl[3];

    do_dVds = 0.0;
    do_B0 = 0.0;
    field_eps = 0.0;
    do_Bmin = -1.0;
    do_Bmax = 0.0;

    for (int k = 1; k <= nth; k++) {
        double th = -NEORT_PI + k * 2.0 * NEORT_PI / nth;
        x[2] = th;
        do_magfie(x, &bmod, &sqrtg, hder, hcovar, hctrvr, hcurl);
        do_dVds += fabs(sqrtg) * dth;
        do_B0 += bmod * dth;
        field_eps -= cos(th) * bmod * dth;
        if (do_Bmin < 0 || bmod < do_Bmin) {
            do_Bmin = bmod;
            orbit_th0 = th;
        }
        if (bmod > do_Bmax)
            do_Bmax = bmod;
    }

    do_dVds = 2.0 * NEORT_PI * do_dVds;
    do_B0 = do_B0 / (2.0 * NEORT_PI);
    field_eps = field_eps / (do_B0 * NEORT_PI);

    do_etatp = 1.0 / do_Bmax;
    do_etadt = 1.0 / do_Bmin;
}
