#ifndef NEORT_PROFILES_H
#define NEORT_PROFILES_H

/* Literal C port of neort_profiles (profiles.f90). File-scope state mirrors the
 * module variables (single-threaded gate). */

extern double prof_vth, prof_dvthds, prof_M_t, prof_dM_tds;
extern double prof_Om_tE, prof_dOm_tEds;
extern double prof_ni1, prof_ni2, prof_Ti1, prof_Ti2, prof_Te;
extern double prof_dni1ds, prof_dni2ds, prof_dTi1ds, prof_dTi2ds, prof_dTeds;
extern double prof_A1, prof_A2;

/* Read plasma.in, build splines, interpolate at s (sets qi, mi, collision state). */
void read_and_init_plasma_input(const char *path, double s);

/* Read profile.in (rotation), build spline, interpolate M_t at s. */
void read_and_init_profile_input(const char *path, double s, double R0, double efac,
                                 double bfac);

/* Thermodynamic forces A1, A2 from psi_pr and q. */
void init_thermodynamic_forces(double psi_pr, double q);

#endif
