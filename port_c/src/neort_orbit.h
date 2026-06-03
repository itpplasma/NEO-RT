#ifndef NEORT_ORBIT_H
#define NEORT_ORBIT_H

/* Port of neort_orbit (orbit.f90). Bounce integration via SUNDIALS CVODE
 * (CV_ADAMS + fixed-point solver + CVodeRootInit), the faithful match for the
 * Fortran DVODE Adams path with g_fcn=bounceroots (nevents=2). */

#define ORBIT_NVAR 7

extern double orbit_th0;     /* starting poloidal angle */
extern int orbit_noshear;    /* neglect magnetic shear */

/* Set the flux-surface coordinate s used by the orbit RHS (do_magfie_mod s). */
void orbit_set_s(double s);

double vpar(double v, double eta, double bmod);
double vperp(double v, double eta, double bmod);

/* Bounce time: integrate poloidal motion (2 vars) to one turn via root-finding. */
double bounce_time(double v, double eta, double taub_estimate, int have_estimate);

/* Full bounce: taub + bounce averages of the nvar integrands. */
void bounce(double v, double eta, double *taub, double bounceavg[ORBIT_NVAR],
            double taub_estimate, int have_estimate);

#endif
