/* Cross-checks the C field port (do_magfie) against reference values from the
 * real Fortran do_magfie_mod (field_ref). Reads reference rows from stdin and
 * asserts the C port reproduces bmod/sqrtg/field components to rtol 1e-12. */
#include <math.h>
#include <stdio.h>
#include "../src/neort_field.h"

#define NTHETA 33

static int close_enough(double a, double b)
{
    return fabs(a - b) <= 1e-13 + 1e-12 * fabs(b);
}

int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: test_field <in_file>\n");
        return 2;
    }
    field_inp_swi = 9;
    field_bfac = 1.0;
    do_magfie_init(argv[1]);

    double s = 0.5;
    int fails = 0, checks = 0;
    for (int i = 0; i < NTHETA; i++) {
        double ref[11];
        for (int k = 0; k < 11; k++)
            if (scanf("%lf", &ref[k]) != 1) {
                fprintf(stderr, "unexpected EOF at row %d\n", i);
                return 2;
            }
        double theta = ref[0];
        double x[3] = {s, 0.0, theta};
        double bmod, sqrtg, bder[3], hcovar[3], hctrvr[3], hcurl[3];
        do_magfie(x, &bmod, &sqrtg, bder, hcovar, hctrvr, hcurl);
        double got[9] = {bmod,     sqrtg,     bder[0],   bder[2],  hcovar[1],
                         hcovar[2], hctrvr[1], hctrvr[2], theta};
        /* reference column order: theta,bmod,sqrtg,bder1,bder3,hcov2,hcov3,hctr2,hctr3 */
        double rcmp[9] = {ref[1], ref[2], ref[3], ref[4], ref[5],
                          ref[6], ref[7], ref[8], ref[0]};
        const char *names[9] = {"bmod",  "sqrtg",  "bder1",  "bder3", "hcov2",
                                "hcov3", "hctr2",  "hctr3",  "theta"};
        for (int k = 0; k < 9; k++) {
            checks++;
            if (!close_enough(got[k], rcmp[k])) {
                fprintf(stderr, "theta=%.5f %s C=%.16e ref=%.16e\n", theta, names[k],
                        got[k], rcmp[k]);
                fails++;
            }
        }
    }
    printf("field: %d checks, %d failures\n", checks, fails);
    return fails == 0 ? 0 : 1;
}
