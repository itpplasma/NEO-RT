#ifndef NEORT_MAGFIE_H
#define NEORT_MAGFIE_H

/* Port of neort_magfie (magfie.f90): flux-surface characterization. Sets the
 * driftorbit globals B0, Bmin, Bmax, etatp, etadt, dVds, the field eps, and the
 * orbit th0 (poloidal angle of minimum B). */

void init_flux_surface_average(double s);

#endif
