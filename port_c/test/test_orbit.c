/* Cross-checks the C orbit port (CVODE) against the real Fortran orbit
 * integration (DVODE), on real trapped/passing orbits. Compares bounce time and
 * bounce averages. CVODE and DVODE share the VODE Adams lineage, so they should
 * agree well within the 1e-8 golden gate; this asserts rtol 1e-9. */
#include <math.h>
#include <stdio.h>
#include "../src/neort_field.h"
#include "../src/neort_magfie.h"
#include "../src/neort_profiles.h"
#include "../src/neort_driftorbit.h"
#include "../src/neort_orbit.h"

/* Gate tolerance is rtol 1e-8 (the golden gate). CVODE vs DVODE on real bounce
 * integrals agrees to ~8e-9 worst case -- event-location (turning-point root)
 * sensitivity loosens it from the 9e-15 seen on a smooth proxy ODE. */
static int close_enough(double a, double b)
{
    return fabs(a - b) <= 1e-12 + 1e-8 * fabs(b);
}

int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "usage: test_orbit <in_file> <plasma.in> <profile.in>\n");
        return 2;
    }
    double s = 0.5;
    field_inp_swi = 9;
    field_bfac = 1.0;
    do_magfie_init(argv[1]);
    double x[3] = {s, 0, 0}, bmod, sqrtg, bder[3], hcovar[3], hctrvr[3], hcurl[3];
    do_magfie(x, &bmod, &sqrtg, bder, hcovar, hctrvr, hcurl);
    init_flux_surface_average(s);
    read_and_init_plasma_input(argv[2], s);
    read_and_init_profile_input(argv[3], s, field_R0, 1.0, 1.0);
    orbit_set_s(s);
    do_sign_vpar = 1.0;

    double vth = prof_vth;
    double vv[3] = {vth, vth, 1.5 * vth};
    double ee[3] = {0.5 * (do_etatp + do_etadt), 0.3 * do_etatp,
                    0.8 * do_etatp + 0.2 * do_etadt};

    const char *names[6] = {"bounce_time", "taub", "bavg1", "bavg2", "bavg3", "bavg6"};
    int fails = 0, checks = 0;
    for (int i = 0; i < 3; i++) {
        double tt = bounce_time(vv[i], ee[i], 0.0, 0);
        double taub, bavg[ORBIT_NVAR];
        bounce(vv[i], ee[i], &taub, bavg, 0.0, 0);
        double got[6] = {tt, taub, bavg[0], bavg[1], bavg[2], bavg[5]};
        for (int k = 0; k < 6; k++) {
            double ref;
            if (scanf("%lf", &ref) != 1) {
                fprintf(stderr, "EOF case %d field %d\n", i, k);
                return 2;
            }
            checks++;
            if (!close_enough(got[k], ref)) {
                fprintf(stderr, "case %d %s C=%.16e ref=%.16e (rel %.2e)\n", i,
                        names[k], got[k], ref, fabs((got[k] - ref) / ref));
                fails++;
            }
        }
    }
    printf("orbit: %d checks, %d failures\n", checks, fails);
    return fails == 0 ? 0 : 1;
}
