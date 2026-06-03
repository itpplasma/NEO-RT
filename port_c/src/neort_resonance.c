#include "neort_resonance.h"
#include "neort_freq.h"
#include "neort_driftorbit.h"
#include <math.h>

#define RST 3 /* roots row stride, mirrors Fortran roots(nlev,3) */

void driftorbit_coarse(double v, double eta_min, double eta_max, double *roots,
                       int ninterv, int *nroots)
{
    double deta = (eta_max - eta_min) * 1.0 / ninterv;
    double resold = 0.0;
    *nroots = 0;
    for (int k = 0; k <= ninterv; k++) {
        double eta = eta_min + k * deta;
        double Omth, dOmthdv, dOmthdeta, Omph, dOmphdv, dOmphdeta;
        Om_th(v, eta, &Omth, &dOmthdv, &dOmthdeta);
        Om_ph(v, eta, &Omph, &dOmphdv, &dOmphdeta);
        double res = do_mth * Omth + do_mph * Omph;
        if (k > 0) {
            if (((res < 0) ? -1.0 : 1.0) != ((resold < 0) ? -1.0 : 1.0)) {
                roots[(*nroots) * RST + 0] = eta - deta;
                roots[(*nroots) * RST + 1] = eta;
                (*nroots)++;
            }
        }
        resold = res;
    }
}

int driftorbit_nroot(double v, double eta_min, double eta_max)
{
    double roots[DO_NLEV * RST];
    int n;
    driftorbit_coarse(v, eta_min, eta_max, roots, DO_NLEV, &n);
    return n;
}

void driftorbit_root(double v, double tol, double eta_min, double eta_max,
                     double out[2])
{
    int maxit = 100, state = -2;
    double etamin2 = eta_min, etamax2 = eta_max;
    double Omth, dOmthdv, dOmthdeta, Omph, dOmphdv, dOmphdeta;
    double res, resmin, resmax;
    double eta = etamin2;

    Om_ph(v, eta, &Omph, &dOmphdv, &dOmphdeta);
    Om_th(v, eta, &Omth, &dOmthdv, &dOmthdeta);
    res = do_mph * Omph + do_mth * Omth;
    resmin = res;

    eta = etamax2;
    Om_ph(v, eta, &Omph, &dOmphdv, &dOmphdeta);
    Om_th(v, eta, &Omth, &dOmthdv, &dOmthdeta);
    resmax = do_mph * Omph + do_mth * Omth;
    int slope_pos = (resmax - resmin > 0);

    if (driftorbit_nroot(v, etamin2, etamax2) == 0) {
        out[0] = eta;
        out[1] = do_mph * dOmphdeta + do_mth * dOmthdeta;
        return;
    }

    for (int k = 1; k <= maxit; k++) {
        Om_ph(v, eta, &Omph, &dOmphdv, &dOmphdeta);
        Om_th(v, eta, &Omth, &dOmthdv, &dOmthdeta);
        res = do_mph * Omph + do_mth * Omth;
        out[0] = eta;
        if (fabs(res) < tol) {
            state = 1;
            out[1] = do_mph * dOmphdeta + do_mth * dOmthdeta;
            break;
        } else if ((slope_pos && res > 0) || (!slope_pos && res < 0)) {
            etamax2 = eta;
            eta = (eta + etamin2) / 2.0;
        } else {
            etamin2 = eta;
            eta = (eta + etamax2) / 2.0;
        }
    }
    if (state < 0)
        out[1] = do_mph * dOmphdeta + do_mth * dOmthdeta;
}
