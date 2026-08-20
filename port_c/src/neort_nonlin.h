#ifndef NEORT_NONLIN_H
#define NEORT_NONLIN_H

#include "neort_orbit.h" /* ORBIT_NVAR */

/* Port of neort_nonlin (nonlin.f90). The golden gate runs nonlin=false, where
 * nonlinear_attenuation == 1. The nonlin=true path additionally needs
 * attenuation_factor (Lagrange lookup over thetafun_inp.dat) which is outside
 * the golden path; it is not yet ported and errors loudly if invoked. */

double nonlinear_attenuation(double ux, double eta, const double bounceavg[ORBIT_NVAR],
                             double Omth, double dOmthdv, double dOmthdeta, double Hmn2);

double omega_prime(double ux, double eta, const double bounceavg[ORBIT_NVAR],
                   double Omth, double dOmdv, double dOmdeta, double dOmdpph);

#endif
