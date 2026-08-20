#ifndef NEORT_UTIL_H
#define NEORT_UTIL_H

/* Physical constants and shared mutable plasma scalars, mirroring util.f90. */

#define NEORT_PI 3.14159265358979323846

/* CGS constants (util.f90) */
#define NEORT_QE 4.803204e-10   /* elementary charge */
#define NEORT_MU 1.660538e-24   /* 1u */
#define NEORT_C 2.997925e+10    /* speed of light */
#define NEORT_EV 1.602176e-12   /* 1 electron volt */

/* Flux-surface-dependent charge/mass, set by the profiles module (util qi, mi). */
extern double util_qi;
extern double util_mi;

#endif
