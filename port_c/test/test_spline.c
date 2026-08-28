/* Cross-checks the C spline port against reference values produced by the
 * Fortran itpplasma spline (spline_ref). Reads "COEFF" and "VAL" lines from
 * stdin and asserts the C port reproduces them to rtol 1e-12. */
#include <math.h>
#include <stdio.h>
#include "../src/neort_spline.h"

#define NP 9

static int close_enough(double a, double b)
{
    double rtol = 1e-12, atol = 1e-14;
    return fabs(a - b) <= atol + rtol * fabs(b);
}

int main(void)
{
    double x[NP], y[NP], coeff[(NP - 1) * 5];
    for (int i = 0; i < NP; i++) {
        x[i] = (double)i * 0.37;
        y[i] = sin(1.3 * x[i]) + 0.5 * x[i];
    }
    spline_coeff(x, y, NP, coeff);

    char tag[16];
    int fails = 0, checks = 0;
    /* compare coefficients */
    for (int line = 0; line < NP - 1; line++) {
        int idx;
        double ref[5];
        if (scanf("%15s", tag) != 1) {
            fprintf(stderr, "unexpected EOF reading COEFF\n");
            return 2;
        }
        if (scanf("%d %lf %lf %lf %lf %lf", &idx, &ref[0], &ref[1], &ref[2],
                  &ref[3], &ref[4]) != 6)
            return 2;
        for (int k = 0; k < 5; k++) {
            checks++;
            if (!close_enough(coeff[idx * 5 + k], ref[k])) {
                fprintf(stderr, "COEFF[%d][%d] C=%.16e ref=%.16e\n", idx, k,
                        coeff[idx * 5 + k], ref[k]);
                fails++;
            }
        }
    }

    /* compare evaluations */
    int jstart = 1;
    while (scanf("%15s", tag) == 1) {
        double xe, ref[3], out[3];
        if (scanf("%lf %lf %lf %lf", &xe, &ref[0], &ref[1], &ref[2]) != 4)
            return 2;
        spline_val_0(coeff, NP - 1, xe, &jstart, out);
        for (int k = 0; k < 3; k++) {
            checks++;
            if (!close_enough(out[k], ref[k])) {
                fprintf(stderr, "VAL(x=%.6f)[%d] C=%.16e ref=%.16e\n", xe, k,
                        out[k], ref[k]);
                fails++;
            }
        }
    }

    printf("spline: %d checks, %d failures\n", checks, fails);
    return fails == 0 ? 0 : 1;
}
