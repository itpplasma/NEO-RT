#include "neort_collis.h"
#include <math.h>

double collis_efcolf[COLLIS_NSORTS] = {0};
double collis_velrat[COLLIS_NSORTS] = {0};
double collis_enrat[COLLIS_NSORTS] = {0};

static void onseff(double v, double *d_p, double *dh, double *dpd)
{
    const double sqp = 1.7724538;     /* sqrt(pi) */
    const double cons = 0.75225278;   /* 4/(3 sqrt(pi)) */
    double v2 = v * v, v3 = v2 * v, ex, er;
    if (v < 0.01) {
        *d_p = cons * (1.0 - 0.6 * v2);
        *dh = cons * (1.0 - 0.2 * v2);
        *dpd = 2.0 * cons * (1.0 - 1.2 * v2);
    } else if (v > 6.0) {
        *d_p = 1.0 / v3;
        *dh = (1.0 - 0.5 / v2) / v;
        *dpd = -1.0 / v3;
    } else {
        ex = exp(-v2) / sqp;
        er = erf(v);
        *d_p = er / v3 - 2.0 * ex / v2;
        *dh = er * (1.0 - 0.5 / v2) / v + ex / v2;
        *dpd = 4.0 * ex - *d_p;
    }
}

void coleff(double p, double *dpp, double *dhh, double *fpeff)
{
    double plim = p > 1.0e-8 ? p : 1.0e-8;
    *dpp = 0.0;
    *dhh = 0.0;
    *fpeff = 0.0;
    for (int i = 0; i < COLLIS_NSORTS; i++) {
        double xbeta = p * collis_velrat[i], d_p, dh, dpd;
        onseff(xbeta, &d_p, &dh, &dpd);
        *dpp += d_p * collis_efcolf[i];
        *dhh += dh * collis_efcolf[i];
        *fpeff += (dpd / plim - 2.0 * d_p * p * collis_enrat[i]) * collis_efcolf[i];
    }
    *dhh /= plim * plim;
}

void loacol_nbi(double amb, double am1, double am2, double Zb, double Z1, double Z2,
                double densi1, double densi2, double tempi1, double tempi2,
                double tempe, double ebeam, double *v0, double *dchichi,
                double *slowrate, double *dchichi_norm, double *slowrate_norm)
{
    const double pi = 3.14159265358979;
    const double pmass = 1.6726e-24, emass = 9.1094e-28, e = 4.8032e-10, ev = 1.6022e-12;
    double dense, vti1, vti2, vte, alame, frecol_base, alami1, alami2;
    double eps = 2.220446049250313e-16; /* epsilon(1.0_dp) for real64 */

    collis_enrat[0] = ebeam / tempi1;
    collis_enrat[1] = ebeam / tempi2;
    collis_enrat[2] = ebeam / tempe;

    *v0 = sqrt(2.0 * ebeam * ev / (amb * pmass));
    vti1 = sqrt(2.0 * tempi1 * ev / (pmass * am1));
    vti2 = sqrt(2.0 * tempi2 * ev / (pmass * am2));
    vte = sqrt(2.0 * tempe * ev / emass);

    collis_velrat[0] = *v0 / vti1;
    collis_velrat[1] = *v0 / vti2;
    collis_velrat[2] = *v0 / vte;

    dense = densi1 * Z1 + densi2 * Z2;
    {
        double a1 = sqrt(densi1 * Z1 * Z1 / tempi1) * Zb * Z1 * (amb + am1) /
                    (amb * tempi1 + am1 * ebeam);
        double a2 = sqrt(densi2 * Z2 * Z2 / tempi2) * Zb * Z2 * (amb + am2) /
                    (amb * tempi2 + am2 * ebeam);
        alami1 = 23.0 - log(a1 > eps ? a1 : eps);
        alami2 = 23.0 - log(a2 > eps ? a2 : eps);
    }
    alame = 24.0 - log(sqrt(dense) / tempe);
    frecol_base = 2.0 * pi * dense * e * e * e * e * Zb * Zb / ((amb * pmass) * (amb * pmass) * (*v0) * (*v0) * (*v0));
    frecol_base = frecol_base / *v0;

    collis_efcolf[0] = frecol_base * Z1 * Z1 * alami1 * densi1 / dense;
    collis_efcolf[1] = frecol_base * Z2 * Z2 * alami2 * densi2 / dense;
    collis_efcolf[2] = frecol_base * alame;

    for (int i = 0; i < COLLIS_NSORTS; i++)
        collis_efcolf[i] *= collis_velrat[i];

    /* dchichi/slowrate outputs are computed by the Fortran but unused by callers;
     * keep the signature and set them to 0 (the Fortran leaves them undefined too
     * for the NEO-RT path -- it only reads efcolf/velrat/enrat). */
    *dchichi = 0.0;
    *slowrate = 0.0;
    *dchichi_norm = 0.0;
    *slowrate_norm = 0.0;
}
