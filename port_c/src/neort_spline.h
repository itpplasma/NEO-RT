#ifndef NEORT_SPLINE_H
#define NEORT_SPLINE_H

/* Cubic spline per Sormann, literal port of itpplasma/spline (spline.f90).
 * Coefficient layout matches the Fortran (n-1, 5) array stored row-major:
 *   coeff[i][0]=x_i, coeff[i][1]=y_i, coeff[i][2..4]=c1,c2,c3 (poly coeffs).
 * The Fortran coefficient solve uses LAPACK dptsv; we call the same routine so
 * the result is bit-identical to the Fortran path. */

/* Number of spline intervals for n knots is n-1. coeff must hold (n-1)*5 doubles. */
void spline_coeff(const double *x, const double *y, int n, double *coeff);

/* Evaluate spline at x. out[0]=value, out[1]=first deriv, out[2]=second deriv.
 * j_start carries the saved interval index between calls (as in the Fortran
 * threadprivate j_start); pass a pointer to a persistent int initialised to 1. */
void spline_val_0(const double *coeff, int nseg, double x, int *j_start, double out[3]);

#endif
