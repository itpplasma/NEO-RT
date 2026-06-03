#include "neort_field.h"
#include "neort_spline.h"
#include <complex.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

static const double PI = 3.14159265358979323846;
static const double SIGN_THETA = -1.0;        /* left-handed */
static const double ITOB = 2.0e-1 * (-1.0);    /* 2e-1 * sign_theta */

double field_bfac = 1.0;
int field_inp_swi = 0;
double field_psi_pr = 0, field_R0 = 0, field_a = 0, field_B00 = 0;
double field_Bthcov = 0, field_Bphcov = 0, field_dBthcovds = 0, field_dBphcovds = 0;
double field_q = 0, field_dqds = 0, field_iota = 0, field_eps = 0, field_B0h = 0;
int field_nflux = 0, field_nmode = 0;

/* ---- module state for the axisymmetric field (do_magfie_mod) ---- */
static int ncol1 = 5, ncol2 = 0;
static double *params0 = NULL;   /* nflux x (ncol1+1) row-major */
static double *modes0 = NULL;    /* nflux x nmode x (ncol2+2) row-major */
static double *spl1 = NULL;      /* ncol1 splines, each (nflux-1)*5 */
static double *spl2 = NULL;      /* ncol2*nmode splines, each (nflux-1)*5 */
static int nseg = 0;             /* nflux-1 */
static int jstart = 1;           /* shared spline index guess (value-independent) */

/* ---- module state for the perturbation field (do_magfie_pert_mod) ---- */
static int p_ncol1 = 5, p_ncol2 = 0, p_nflux = 0, p_nmode = 0, p_nfp = 0, p_nseg = 0;
static double *p_params = NULL, *p_modes = NULL, *p_spl2 = NULL;

#define PAR(i, k) params0[(i) * (ncol1 + 1) + (k)]
#define MOD(i, j, k) modes0[((i) * field_nmode + (j)) * (ncol2 + 2) + (k)]
#define S1(k) (spl1 + (k) * nseg * 5)
#define S2(k, j) (spl2 + (((k) * field_nmode + (j)) * nseg * 5))

#define PPAR(i, k) p_params[(i) * (p_ncol1 + 1) + (k)]
#define PMOD(i, j, k) p_modes[((i) * p_nmode + (j)) * (p_ncol2 + 2) + (k)]
#define PS2(k, j) (p_spl2 + (((k) * p_nmode + (j)) * p_nseg * 5))

static void skip_lines(FILE *f, int n)
{
    int c;
    for (int k = 0; k < n; k++)
        while ((c = fgetc(f)) != '\n' && c != EOF)
            ;
}

/* Generic Boozer reader. Fills caller-provided arrays; returns via out params.
 * Mirrors boozer_read / boozer_read_pert exactly (line skips empirically
 * matched: header skips 5, per-surface skips 2 before params, 1 before modes). */
static void boozer_read(const char *path, int nc2, int *nflux_out, int *nmode_out,
                        int *nfp_out, double *flux_out, double *a_out, double *R0_out,
                        double **params_out, double **modes_out)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "cannot open %s\n", path);
        exit(1);
    }
    int m0b, n0b, nflux, nfp;
    double flux, a, R0;
    skip_lines(f, 5);
    if (fscanf(f, "%d %d %d %d %lf %lf %lf", &m0b, &n0b, &nflux, &nfp, &flux, &a,
               &R0) != 7) {
        fprintf(stderr, "bad Boozer header in %s\n", path);
        exit(1);
    }
    skip_lines(f, 1); /* finish the header value line */
    int nmode = (m0b + 1) * (n0b + 1);
    int pcols = ncol1 + 1, mcols = nc2 + 2;
    double *params = (double *)malloc((size_t)nflux * pcols * sizeof(double));
    double *modes = (double *)malloc((size_t)nflux * nmode * mcols * sizeof(double));

    for (int ks = 0; ks < nflux; ks++) {
        skip_lines(f, 2); /* two per-surface label lines */
        for (int k = 0; k < pcols; k++)
            if (fscanf(f, "%lf", &params[ks * pcols + k]) != 1) {
                fprintf(stderr, "bad params surface %d in %s\n", ks, path);
                exit(1);
            }
        skip_lines(f, 2); /* finish param line + skip mode-label line */
        for (int j = 0; j < nmode; j++) {
            for (int k = 0; k < mcols; k++)
                if (fscanf(f, "%lf", &modes[(ks * nmode + j) * mcols + k]) != 1) {
                    fprintf(stderr, "bad mode (%d,%d) in %s\n", ks, j, path);
                    exit(1);
                }
            skip_lines(f, 1); /* finish this mode line */
        }
    }
    fclose(f);
    *nflux_out = nflux;
    *nmode_out = nmode;
    *nfp_out = nfp;
    *flux_out = flux;
    *a_out = a;
    *R0_out = R0;
    *params_out = params;
    *modes_out = modes;
}

static void build_splines(const double *params, const double *modes, int nflux,
                          int nmode, int nc2, int pcols, int mcols, double **spl1_out,
                          double **spl2_out)
{
    int seg = nflux - 1;
    double *x = (double *)malloc(nflux * sizeof(double));
    double *y = (double *)malloc(nflux * sizeof(double));
    for (int i = 0; i < nflux; i++)
        x[i] = params[i * pcols + 0];

    double *s1 = (double *)malloc((size_t)ncol1 * seg * 5 * sizeof(double));
    for (int k = 0; k < ncol1; k++) {
        for (int i = 0; i < nflux; i++)
            y[i] = params[i * pcols + (k + 1)];
        spline_coeff(x, y, nflux, s1 + k * seg * 5);
    }
    double *s2 = (double *)malloc((size_t)nc2 * nmode * seg * 5 * sizeof(double));
    for (int j = 0; j < nmode; j++)
        for (int k = 0; k < nc2; k++) {
            for (int i = 0; i < nflux; i++)
                y[i] = modes[(i * nmode + j) * mcols + (k + 2)];
            spline_coeff(x, y, nflux, s2 + (k * nmode + j) * seg * 5);
        }
    free(x);
    free(y);
    *spl1_out = s1;
    *spl2_out = s2;
}

void do_magfie_init(const char *path)
{
    double flux, a, R0;
    int nfp_dummy;
    ncol1 = 5;
    ncol2 = (field_inp_swi == 8) ? 4 : 8;
    boozer_read(path, ncol2, &field_nflux, &field_nmode, &nfp_dummy, &flux, &a, &R0,
                &params0, &modes0);
    field_a = 100.0 * a;
    field_R0 = 100.0 * R0;
    field_psi_pr = 1.0e8 * flux / (2.0 * PI) * field_bfac;
    nseg = field_nflux - 1;
    build_splines(params0, modes0, field_nflux, field_nmode, ncol2, ncol1 + 1,
                  ncol2 + 2, &spl1, &spl2);
    /* R0 from first harmonic rmnc, B00 from modes0(1,1,6) (col index 5) */
    field_R0 = MOD(0, 0, 2) * 100.0;
    field_B00 = 1.0e4 * MOD(0, 0, 5) * field_bfac;
}

static void fast_sin_cos(double m0, double dm, double x, int n, const double *mvals,
                         double *sinterm, double *costerm)
{
    (void)mvals;
    double _Complex ff = cexp(I * m0 * x);
    double _Complex rot = cexp(I * dm * x);
    for (int j = 0; j < n; j++) {
        costerm[j] = creal(ff);
        sinterm[j] = cimag(ff);
        ff *= rot;
    }
}

void do_magfie(const double x[3], double *bmod, double *sqrtg, double bder[3],
               double hcovar[3], double hctrvr[3], double hcurl[3])
{
    double sv[3], x1, bf = field_bfac;
    x1 = PAR(0, 0) > x[0] ? PAR(0, 0) : x[0];
    if (PAR(field_nflux - 1, 0) < x1)
        x1 = PAR(field_nflux - 1, 0);

    spline_val_0(S1(2), nseg, x1, &jstart, sv);
    field_Bthcov = ITOB * sv[0] * bf;
    field_dBthcovds = ITOB * sv[1] * bf;
    spline_val_0(S1(1), nseg, x1, &jstart, sv);
    field_Bphcov = ITOB * sv[0] * bf;
    field_dBphcovds = ITOB * sv[1] * bf;
    spline_val_0(S1(0), nseg, x1, &jstart, sv);
    field_iota = sv[0];
    field_q = 1.0 / field_iota;
    field_dqds = -sv[1] / (field_iota * field_iota);

    int nm = field_nmode;
    double *cost = (double *)malloc(nm * sizeof(double));
    double *sint = (double *)malloc(nm * sizeof(double));
    double m0 = MOD(0, 0, 0), m1 = MOD(0, 1, 0);
    fast_sin_cos(m0, m1 - m0, x[2], nm, NULL, sint, cost);

    double bm = 0.0, db = 0.0, bth = 0.0;
    if (field_inp_swi == 8) {
        for (int j = 0; j < nm; j++) {
            spline_val_0(S2(3, j), nseg, x1, &jstart, sv);
            double bmnc = 1.0e4 * sv[0] * bf, dbmnc = 1.0e4 * sv[1] * bf;
            if (j == 0)
                field_B0h = bmnc;
            bm += bmnc * cost[j];
            db += dbmnc * cost[j];
            bth += -MOD(0, j, 0) * bmnc * sint[j];
        }
    } else { /* inp_swi == 9 */
        for (int j = 0; j < nm; j++) {
            spline_val_0(S2(6, j), nseg, x1, &jstart, sv);
            double bmnc = 1.0e4 * sv[0] * bf;
            double dbmnc = 1.0e4 * sv[1] * bf;
            spline_val_0(S2(7, j), nseg, x1, &jstart, sv);
            double bmns = 1.0e4 * sv[0] * bf;
            double dbmns = 1.0e4 * sv[1] * bf;
            if (j == 0)
                field_B0h = bmnc;
            double mj = MOD(0, j, 0);
            bm += bmnc * cost[j] + bmns * sint[j];
            db += dbmnc * cost[j] + dbmns * sint[j];
            bth += -mj * bmnc * sint[j] + mj * bmns * cost[j];
        }
    }
    *bmod = bm;
    bder[0] = db / bm;
    bder[1] = 0.0;
    bder[2] = bth / bm;
    free(cost);
    free(sint);

    double sqgbmod2 = SIGN_THETA * field_psi_pr * (field_Bphcov + field_iota * field_Bthcov);
    double sqgbmod = sqgbmod2 / bm;
    *sqrtg = sqgbmod / bm;

    hcovar[0] = 0.0;
    hcovar[1] = field_Bphcov / bm;
    hcovar[2] = field_Bthcov / bm;
    hctrvr[0] = 0.0;
    hctrvr[1] = SIGN_THETA * field_psi_pr / sqgbmod;
    hctrvr[2] = SIGN_THETA * field_iota * field_psi_pr / sqgbmod;
    hcurl[0] = hcurl[1] = hcurl[2] = 0.0;
}

void do_magfie_pert_init(const char *path)
{
    double flux, a, R0;
    p_ncol1 = 5;
    p_ncol2 = (field_inp_swi == 8) ? 4 : 8;
    /* boozer_read uses the do_magfie ncol1/ncol2 macros; set them for the call. */
    int save_nmode = field_nmode;
    ncol1 = p_ncol1;
    ncol2 = p_ncol2;
    boozer_read(path, p_ncol2, &p_nflux, &p_nmode, &p_nfp, &flux, &a, &R0, &p_params,
                &p_modes);
    p_nseg = p_nflux - 1;
    double *dummy1;
    /* build splines into p_spl2 (we only need the mode splines for the amplitude). */
    int seg = p_nseg, pcols = p_ncol1 + 1, mcols = p_ncol2 + 2;
    double *x = (double *)malloc(p_nflux * sizeof(double));
    double *y = (double *)malloc(p_nflux * sizeof(double));
    for (int i = 0; i < p_nflux; i++)
        x[i] = p_params[i * pcols + 0];
    dummy1 = (double *)malloc((size_t)p_ncol1 * seg * 5 * sizeof(double));
    for (int k = 0; k < p_ncol1; k++) {
        for (int i = 0; i < p_nflux; i++)
            y[i] = p_params[i * pcols + (k + 1)];
        spline_coeff(x, y, p_nflux, dummy1 + k * seg * 5);
    }
    p_spl2 = (double *)malloc((size_t)p_ncol2 * p_nmode * seg * 5 * sizeof(double));
    for (int j = 0; j < p_nmode; j++)
        for (int k = 0; k < p_ncol2; k++) {
            for (int i = 0; i < p_nflux; i++)
                y[i] = p_modes[(i * p_nmode + j) * mcols + (k + 2)];
            spline_coeff(x, y, p_nflux, p_spl2 + (k * p_nmode + j) * seg * 5);
        }
    free(x);
    free(y);
    free(dummy1);
    field_nmode = save_nmode; /* restore axisymmetric nmode */
}

void do_magfie_pert_amp(const double x[3], double _Complex *bamp)
{
    double x1, bf = field_bfac, sv[3];
    x1 = PPAR(0, 0) > x[0] ? PPAR(0, 0) : x[0];
    if (PPAR(p_nflux - 1, 0) < x1)
        x1 = PPAR(p_nflux - 1, 0);

    double _Complex sum = 0.0;
    if (field_inp_swi == 8) {
        for (int j = 0; j < p_nmode; j++) {
            spline_val_0(PS2(3, j), p_nseg, x1, &jstart, sv);
            double bmnc = 1.0e4 * sv[0] * bf;
            sum += bmnc * cos(PMOD(0, j, 0) * x[2]);
        }
    } else {
        double m0 = PMOD(0, 0, 0), m1 = PMOD(0, 1, 0), dm = m1 - m0;
        double _Complex ff = cexp(I * m0 * x[2]);
        double _Complex rot = cexp(I * dm * x[2]);
        for (int j = 0; j < p_nmode; j++) {
            spline_val_0(PS2(6, j), p_nseg, x1, &jstart, sv);
            double bmnc = 1.0e4 * sv[0] * bf;
            spline_val_0(PS2(7, j), p_nseg, x1, &jstart, sv);
            double bmns = 1.0e4 * sv[0] * bf;
            sum += (bmnc - I * bmns) * ff;
            ff *= rot;
        }
    }
    *bamp = sum;
}
