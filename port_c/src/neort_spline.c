#include "neort_spline.h"

/* LAPACK symmetric positive-definite tridiagonal solver (Fortran ABI). */
extern void dptsv_(const int *n, const int *nrhs, double *d, double *e,
                   double *b, const int *ldb, int *info);

/* coeff is row-major (n-1) x 5. Index helper. */
#define C(i, k) coeff[(i) * 5 + (k)]

void spline_coeff(const double *x, const double *y, int n, double *coeff)
{
    /* h(i)=x(i+1)-x(i), r(i)=y(i+1)-y(i), i=0..n-2 (n-1 values). */
    int m = n - 1; /* number of segments */
    double *h = (double *)__builtin_alloca(m * sizeof(double));
    double *r = (double *)__builtin_alloca(m * sizeof(double));
    for (int i = 0; i < m; i++) {
        h[i] = x[i + 1] - x[i];
        r[i] = y[i + 1] - y[i];
    }

    /* Tridiagonal SPD system of size n-2.
     * d (diag, n-2), dl (sub-diag, n-3), c (rhs/solution, n-2). */
    int sys = n - 2;
    double *d = (double *)__builtin_alloca(sys * sizeof(double));
    double *dl = (double *)__builtin_alloca((sys > 1 ? sys - 1 : 1) * sizeof(double));
    double *c = (double *)__builtin_alloca(sys * sizeof(double));

    for (int k = 0; k < sys; k++) {
        d[k] = 2.0 * (h[k] + h[k + 1]);
        c[k] = 3.0 * (r[k + 1] / h[k + 1] - r[k] / h[k]);
    }
    for (int k = 0; k < sys - 1; k++)
        dl[k] = h[k + 1];

    int nrhs = 1, ldb = sys, info = 0;
    dptsv_(&sys, &nrhs, d, dl, c, &ldb, &info);

    /* Assemble coefficients exactly as the Fortran does. */
    for (int i = 0; i < m; i++) {
        C(i, 0) = x[i];
        C(i, 1) = y[i];
        C(i, 2) = 0.0;
        C(i, 3) = 0.0;
        C(i, 4) = 0.0;
    }

    /* column 3 (index 2): first-derivative-like coefficient */
    C(0, 2) = r[0] / h[0] - h[0] / 3.0 * c[0];
    for (int i = 1; i <= n - 3; i++)
        C(i, 2) = r[i] / h[i] - h[i] / 3.0 * (c[i] + 2.0 * c[i - 1]);
    C(m - 1, 2) = r[m - 1] / h[m - 1] - h[m - 1] / 3.0 * (2.0 * c[sys - 1]);

    /* column 4 (index 3): coeff(1,4)=0; coeff(2:,4)=c */
    C(0, 3) = 0.0;
    for (int i = 1; i <= m - 1; i++)
        C(i, 3) = c[i - 1];

    /* column 5 (index 4): second-difference coefficient */
    C(0, 4) = 1.0 / (3.0 * h[0]) * c[0];
    for (int i = 1; i <= n - 3; i++)
        C(i, 4) = 1.0 / (3.0 * h[i]) * (c[i] - c[i - 1]);
    C(m - 1, 4) = 1.0 / (3.0 * h[m - 1]) * (-c[sys - 1]);
}

void spline_val_0(const double *coeff, int nseg, double x, int *j_start, double out[3])
{
    int n = nseg + 1; /* matches Fortran n = size(coeff,1)+1 */
    int j = *j_start;
    double z;

    if (j < 1 || j >= n)
        j = 1;
    if (n <= 1) {
        out[0] = out[1] = out[2] = 0.0;
        return;
    }
    /* Fortran uses 1-based j; convert search to 0-based interval index ji=j-1. */
    int ji = j - 1;
    if (x < C(ji, 0)) {
        while (ji > 0) {
            if (x < C(ji, 0))
                ji--;
            else
                break;
        }
    } else {
        while (ji < n - 2) {
            if (x >= C(ji + 1, 0))
                ji++;
            else
                break;
        }
    }
    *j_start = ji + 1;

    z = x - C(ji, 0);
    out[0] = ((C(ji, 4) * z + C(ji, 3)) * z + C(ji, 2)) * z + C(ji, 1);
    out[1] = (3.0 * C(ji, 4) * z + 2.0 * C(ji, 3)) * z + C(ji, 2);
    out[2] = 6.0 * C(ji, 4) * z + 2.0 * C(ji, 3);
}
