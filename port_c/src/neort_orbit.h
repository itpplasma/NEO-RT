#ifndef NEORT_ORBIT_H
#define NEORT_ORBIT_H

/* Port of neort_orbit (orbit.f90). Bounce integration via SUNDIALS CVODE
 * (CV_ADAMS + fixed-point solver + CVodeRootInit), the faithful match for the
 * Fortran DVODE Adams path with g_fcn=bounceroots (nevents=2). */

#define ORBIT_NVAR 7

extern double orbit_th0;     /* starting poloidal angle */
extern int orbit_noshear;    /* neglect magnetic shear */

/* Set/get the flux-surface coordinate s used by the orbit RHS (do_magfie_mod s). */
void orbit_set_s(double s);
double orbit_get_s(void);

double vpar(double v, double eta, double bmod);
double vperp(double v, double eta, double bmod);

/* Bounce time: integrate poloidal motion (2 vars) to one turn via root-finding. */
double bounce_time(double v, double eta, double taub_estimate, int have_estimate);

/* Fixed-interval bounce average: integrate the nvar system from 0 to taub (no
 * root-finding), return bounceavg = y/taub. Mirrors orbit.f90 bounce_fast. */
void bounce_fast(double v, double eta, double taub, double bounceavg[ORBIT_NVAR]);

/* RHS function type for a custom timestep (v, eta, neq, t, y, ydot). */
typedef void (*orbit_ts_fn)(double v, double eta, int neq, double t,
                            const double *y, double *ydot);

/* Generalized bounce_fast with a custom RHS and CVODE istate out (2 = success).
 * Used by transport with its perturbation-Hamiltonian timestep. */
void bounce_fast_ext(double v, double eta, double taub, double bounceavg[ORBIT_NVAR],
                     orbit_ts_fn ts, int *istate_out);

/* Poloidal drift velocity RHS piece (orbit.f90 poloidal_velocity), reused by
 * the transport timestep. */
void poloidal_velocity(double v, double eta, double bmod, double hthctr,
                       double hderth, double v_par, double ydot[2]);

/* Full bounce: taub + bounce averages of the nvar integrands. */
void bounce(double v, double eta, double *taub, double bounceavg[ORBIT_NVAR],
            double taub_estimate, int have_estimate);

#endif
