#ifndef NEORT_TRANSPORT_H
#define NEORT_TRANSPORT_H

/* Port of neort_transport (transport.f90): velocity-space integration of the
 * resonant transport coefficients and torque density for one harmonic. */

/* Midpoint-rule integral over v in [vmin,vmax] (vsteps steps): fills D[0]=D11,
 * D[1]=D12 (normalized to the plateau coefficient) and *T = torque density. */
void compute_transport_integral(double vmin, double vmax, int vsteps, double D[2],
                                double *T);

#endif
