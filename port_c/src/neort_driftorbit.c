#include "neort_driftorbit.h"

double do_efac = 1.0;
double do_epsmn = 1.0;
int do_m0 = 1;
int do_mth = 1;
int do_mph = 1;
int do_magdrift = 1;
int do_nopassing = 0;
int do_pertfile = 0;
int do_comptorque = 1;
int do_nonlin = 0;

double do_dVds = 0.0, do_etadt = 0.0, do_etatp = 0.0, do_etamin = 0.0, do_etamax = 0.0;
double do_B0 = 0.0, do_Bmin = 0.0, do_Bmax = 0.0;
double do_sign_vpar = 1.0, do_sign_vpar_htheta = 1.0;
