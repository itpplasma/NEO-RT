#include "neort_profiles.h"
#include "neort_util.h"
#include "neort_collis.h"
#include "neort_spline.h"
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

double util_qi = 1.0 * NEORT_QE;
double util_mi = 2.014 * NEORT_MU;

double prof_vth = 0, prof_dvthds = 0, prof_M_t = 0, prof_dM_tds = 0;
double prof_Om_tE = 0, prof_dOm_tEds = 0;
double prof_ni1 = 0, prof_ni2 = 0, prof_Ti1 = 0, prof_Ti2 = 0, prof_Te = 0;
double prof_dni1ds = 0, prof_dni2ds = 0, prof_dTi1ds = 0, prof_dTi2ds = 0, prof_dTeds = 0;
double prof_A1 = 0, prof_A2 = 0;

static const double SIGN_THETA = -1.0;

static double *plasma_spl = NULL; /* 5 splines, each (nplasma-1)*5 */
static int plasma_nseg = 0;
static double am1_g, am2_g, Z1_g, Z2_g;

static double *Mt_spl = NULL; /* 1 spline, (nprof-1)*5 */
static int prof_nseg = 0;
static int jstart = 1;

static void skip_line(FILE *f)
{
    int c;
    while ((c = fgetc(f)) != '\n' && c != EOF)
        ;
}

void read_and_init_plasma_input(const char *path, double s)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "cannot open %s\n", path);
        exit(1);
    }
    int nplasma;
    skip_line(f); /* "% N am1 am2 Z1 Z2" */
    if (fscanf(f, "%d %lf %lf %lf %lf", &nplasma, &am1_g, &am2_g, &Z1_g, &Z2_g) != 5) {
        fprintf(stderr, "bad plasma header\n");
        exit(1);
    }
    skip_line(f); /* finish header value line */
    skip_line(f); /* "% s ni_1 ..." */
    double *plasma = (double *)malloc((size_t)nplasma * 6 * sizeof(double));
    for (int k = 0; k < nplasma; k++)
        for (int j = 0; j < 6; j++)
            if (fscanf(f, "%lf", &plasma[k * 6 + j]) != 1) {
                fprintf(stderr, "bad plasma row %d\n", k);
                exit(1);
            }
    fclose(f);

    plasma_nseg = nplasma - 1;
    double *x = (double *)malloc(nplasma * sizeof(double));
    double *y = (double *)malloc(nplasma * sizeof(double));
    for (int i = 0; i < nplasma; i++)
        x[i] = plasma[i * 6 + 0];
    plasma_spl = (double *)malloc((size_t)5 * plasma_nseg * 5 * sizeof(double));
    for (int kk = 0; kk < 5; kk++) {
        for (int i = 0; i < nplasma; i++)
            y[i] = plasma[i * 6 + (kk + 1)];
        spline_coeff(x, y, nplasma, plasma_spl + kk * plasma_nseg * 5);
    }
    free(x);
    free(y);
    free(plasma);

    /* init_plasma_at_s */
    double sv[3];
    double *S0 = plasma_spl + 0 * plasma_nseg * 5;
    spline_val_0(S0, plasma_nseg, s, &jstart, sv);
    prof_ni1 = sv[0];
    prof_dni1ds = sv[1];
    spline_val_0(plasma_spl + 1 * plasma_nseg * 5, plasma_nseg, s, &jstart, sv);
    prof_ni2 = sv[0];
    prof_dni2ds = sv[1];
    spline_val_0(plasma_spl + 2 * plasma_nseg * 5, plasma_nseg, s, &jstart, sv);
    prof_Ti1 = sv[0];
    prof_dTi1ds = sv[1];
    spline_val_0(plasma_spl + 3 * plasma_nseg * 5, plasma_nseg, s, &jstart, sv);
    prof_Ti2 = sv[0];
    prof_dTi2ds = sv[1];
    spline_val_0(plasma_spl + 4 * plasma_nseg * 5, plasma_nseg, s, &jstart, sv);
    prof_Te = sv[0];
    prof_dTeds = sv[1];

    util_qi = Z1_g * NEORT_QE;
    util_mi = am1_g * NEORT_MU;
    prof_vth = sqrt(2.0 * prof_Ti1 * NEORT_EV / util_mi);
    prof_dvthds = 0.5 * sqrt(2.0 * NEORT_EV / (util_mi * prof_Ti1)) * prof_dTi1ds;

    const double pmass = 1.6726e-24;
    double v0 = prof_vth, amb = 2.0, Zb = 1.0;
    double ebeam = amb * pmass * v0 * v0 / (2.0 * NEORT_EV);
    double dchichi, slowrate, dn, sn;
    loacol_nbi(amb, am1_g, am2_g, Zb, Z1_g, Z2_g, prof_ni1, prof_ni2, prof_Ti1,
               prof_Ti2, prof_Te, ebeam, &v0, &dchichi, &slowrate, &dn, &sn);
}

void read_and_init_profile_input(const char *path, double s, double R0, double efac,
                                 double bfac)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "cannot open %s\n", path);
        exit(1);
    }
    /* readdata(path, 2): count rows, read first 2 columns of each. */
    int cap = 64, n = 0;
    double *sx = (double *)malloc(cap * sizeof(double));
    double *mt = (double *)malloc(cap * sizeof(double));
    while (1) {
        double a, b, c;
        int r = fscanf(f, "%lf %lf %lf", &a, &b, &c);
        if (r < 3)
            break;
        if (n >= cap) {
            cap *= 2;
            sx = (double *)realloc(sx, cap * sizeof(double));
            mt = (double *)realloc(mt, cap * sizeof(double));
        }
        sx[n] = a;
        mt[n] = b;
        n++;
    }
    fclose(f);

    prof_nseg = n - 1;
    Mt_spl = (double *)malloc((size_t)prof_nseg * 5 * sizeof(double));
    spline_coeff(sx, mt, n, Mt_spl);
    free(sx);
    free(mt);

    double sv[3];
    spline_val_0(Mt_spl, prof_nseg, s, &jstart, sv);
    prof_M_t = sv[0] * efac / bfac;
    prof_dM_tds = sv[1] * efac / bfac;
    prof_Om_tE = prof_vth * prof_M_t / R0;
    prof_dOm_tEds = prof_vth * prof_dM_tds / R0 + prof_M_t * prof_dvthds / R0;
}

void init_thermodynamic_forces(double psi_pr, double q)
{
    prof_A1 = prof_dni1ds / prof_ni1 -
              util_qi / (prof_Ti1 * NEORT_EV) * SIGN_THETA * psi_pr / (q * NEORT_C) *
                  prof_Om_tE -
              1.5 * prof_dTi1ds / prof_Ti1;
    prof_A2 = prof_dTi1ds / prof_Ti1;
}
