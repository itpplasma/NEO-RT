#ifndef NEORT_RESONANCE_H
#define NEORT_RESONANCE_H

/* Port of neort_resonance (resonance.f90): resonance-line root finding over eta,
 * res(eta) = mth*Om_th + mph*Om_ph, built on the freq layer. */

/* Coarse bracket: fill roots[nroots][2] with sign-change intervals of res(eta)
 * over [eta_min,eta_max] using ninterv subintervals. roots must hold ninterv rows
 * x (at least 2 cols, stride 3 to mirror the Fortran roots(nlev,3)). */
void driftorbit_coarse(double v, double eta_min, double eta_max, double *roots,
                       int ninterv, int *nroots);

int driftorbit_nroot(double v, double eta_min, double eta_max);

/* Bisection root within [eta_min,eta_max]; out[0]=resonant eta, out[1]=dres/deta. */
void driftorbit_root(double v, double tol, double eta_min, double eta_max,
                     double out[2]);

#endif
