#ifndef NEORT_COLLIS_H
#define NEORT_COLLIS_H

/* Literal C port of collis_alp (collis_nbi.f90). Module state efcolf/velrat/enrat
 * is file-scope (single-threaded gate). */

#define COLLIS_NSORTS 3

extern double collis_efcolf[COLLIS_NSORTS];
extern double collis_velrat[COLLIS_NSORTS];
extern double collis_enrat[COLLIS_NSORTS];

void loacol_nbi(double amb, double am1, double am2, double Zb, double Z1, double Z2,
                double densi1, double densi2, double tempi1, double tempi2,
                double tempe, double ebeam, double *v0, double *dchichi,
                double *slowrate, double *dchichi_norm, double *slowrate_norm);

void coleff(double p, double *dpp, double *dhh, double *fpeff);

#endif
