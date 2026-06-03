/* Cross-checks the C profiles/collis port against the real Fortran modules
 * (profiles_ref). Same init sequence; compares 21 scalars to rtol 1e-12. */
#include <math.h>
#include <stdio.h>
#include "../src/neort_field.h"
#include "../src/neort_profiles.h"
#include "../src/neort_util.h"
#include "../src/neort_collis.h"

static int close_enough(double a, double b)
{
    return fabs(a - b) <= 1e-13 + 1e-12 * fabs(b);
}

int main(int argc, char **argv)
{
    if (argc < 4) {
        fprintf(stderr, "usage: test_profiles <in_file> <plasma.in> <profile.in>\n");
        return 2;
    }
    double s = 0.5;
    field_inp_swi = 9;
    field_bfac = 1.0;
    do_magfie_init(argv[1]);
    double x[3] = {s, 0.0, 0.0}, bmod, sqrtg, bder[3], hcovar[3], hctrvr[3], hcurl[3];
    do_magfie(x, &bmod, &sqrtg, bder, hcovar, hctrvr, hcurl);

    read_and_init_plasma_input(argv[2], s);
    read_and_init_profile_input(argv[3], s, field_R0, 1.0, 1.0);
    init_thermodynamic_forces(field_psi_pr, field_q);

    double got[21] = {
        prof_vth, prof_M_t, prof_Om_tE, prof_A1, prof_A2,
        prof_ni1, prof_ni2, prof_Ti1, prof_Ti2, prof_Te,
        util_qi, util_mi,
        collis_efcolf[0], collis_efcolf[1], collis_efcolf[2],
        collis_velrat[0], collis_velrat[1], collis_velrat[2],
        collis_enrat[0], collis_enrat[1], collis_enrat[2]};
    const char *names[21] = {
        "vth", "M_t", "Om_tE", "A1", "A2", "ni1", "ni2", "Ti1", "Ti2", "Te",
        "qi", "mi", "efcolf1", "efcolf2", "efcolf3", "velrat1", "velrat2",
        "velrat3", "enrat1", "enrat2", "enrat3"};

    int fails = 0;
    for (int k = 0; k < 21; k++) {
        double ref;
        if (scanf("%lf", &ref) != 1) {
            fprintf(stderr, "EOF at %d\n", k);
            return 2;
        }
        if (!close_enough(got[k], ref)) {
            fprintf(stderr, "%s C=%.16e ref=%.16e\n", names[k], got[k], ref);
            fails++;
        }
    }
    printf("profiles: 21 checks, %d failures\n", fails);
    return fails == 0 ? 0 : 1;
}
