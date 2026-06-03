#ifndef NEORT_FIELD_H
#define NEORT_FIELD_H

/* Literal C port of do_magfie_mod / do_magfie_pert_mod (do_magfie_standalone.f90).
 * Single-threaded: the Fortran threadprivate state becomes file-scope statics,
 * which is correct under the OMP_NUM_THREADS=1 golden gate. */

/* Shared axisymmetric-field parameters exposed to callers (mirror module vars). */
extern double field_bfac;     /* B scaling factor (set before reading) */
extern int field_inp_swi;     /* input switch: 8 tok_circ, 9 ASDEX */
extern double field_psi_pr, field_R0, field_a, field_B00;
extern double field_Bthcov, field_Bphcov, field_dBthcovds, field_dBphcovds;
extern double field_q, field_dqds, field_iota, field_eps, field_B0h;
extern int field_nflux, field_nmode;

#define FIELD_SIGN_THETA (-1.0)  /* do_magfie_mod sign_theta parameter */

/* Read axisymmetric Boozer file and build splines (call once). */
void do_magfie_init(const char *path);

/* Evaluate axisymmetric field at x=(s, ph, th).
 * Outputs bmod, sqrtg, and 3-vectors bder, hcovar, hctrvr, hcurl. */
void do_magfie(const double x[3], double *bmod, double *sqrtg, double bder[3],
               double hcovar[3], double hctrvr[3], double hcurl[3]);

/* Perturbation field. */
void do_magfie_pert_init(const char *path);
void do_magfie_pert_amp(const double x[3], double _Complex *bamp);

#endif
