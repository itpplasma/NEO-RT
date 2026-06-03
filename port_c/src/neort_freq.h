#ifndef NEORT_FREQ_H
#define NEORT_FREQ_H

/* Port of neort_freq (freq.f90): canonical-frequency splines over eta and their
 * evaluation with near-boundary logarithmic extrapolation. Built on the verified
 * orbit (bounce_fast/bounce_time) and spline layers. */

#define FREQ_NETASPL 100
#define FREQ_NETASPL_PASS 100

void init_canon_freq_trapped_spline(void);
void init_canon_freq_passing_spline(void);

void Om_th(double v, double eta, double *Omth, double *dOmthdv, double *dOmthdeta);
void Om_tB(double v, double eta, double *OmtB, double *dOmtBdv, double *dOmtBdeta);
void Om_ph(double v, double eta, double *Omph, double *dOmphdv, double *dOmphdeta);
void d_Om_ds(double v, double eta, double taub_estimate, double *dOmthds,
             double *dOmphds);

#endif
