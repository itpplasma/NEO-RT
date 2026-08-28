#ifndef NEORT_DRIFTORBIT_H
#define NEORT_DRIFTORBIT_H

/* State module driftorbit.f90 (declarations only; the resonance/velocity-space
 * logic lives in freq/resonance/transport). File-scope globals mirror the
 * module variables under the single-threaded gate. */

extern double do_efac;       /* radial E-field normalization */
extern double do_epsmn;      /* perturbation amplitude B1/B0 */
extern int do_m0;            /* Boozer poloidal perturbation mode */
extern int do_mth;           /* canonical poloidal mode */
extern int do_mph;           /* toroidal perturbation mode (from pert file) */
extern int do_magdrift;      /* consider magnetic drift */
extern int do_nopassing;     /* neglect passing particles */
extern int do_pertfile;      /* read perturbation from file */
extern int do_comptorque;    /* compute torque */
extern int do_nonlin;        /* nonlinear calculation */

extern double do_dVds, do_etadt, do_etatp, do_etamin, do_etamax;
extern double do_B0, do_Bmin, do_Bmax;
extern double do_sign_vpar, do_sign_vpar_htheta;

/* spline-distance and turning-point tolerances (driftorbit parameters) */
#define DO_EPST_SPL 1.0e-6
#define DO_EPSP_SPL 1.0e-6
#define DO_EPSST_SPL 1.0e-3
#define DO_EPSSP_SPL 1.0e-3
#define DO_EPST 1.0e-8
#define DO_EPSP 1.0e-8
#define DO_NLEV 100

#endif
