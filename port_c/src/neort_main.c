/* NEO-RT C port driver: parse <case>.in, run the golden-path computation
 * (mirrors neort_lib::neort_compute_at_s + neort::compute_transport), and write
 * the 5 golden output files. Single-threaded (matches the deterministic gate). */
#include "neort_field.h"
#include "neort_profiles.h"
#include "neort_magfie.h"
#include "neort_orbit.h"
#include "neort_freq.h"
#include "neort_transport.h"
#include "neort_driftorbit.h"
#include "neort_util.h"
#include <complex.h>
#include <ctype.h>
#include <math.h>
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* ---- config (namelist /params/) ---- */
struct config {
    double s, M_t, qs, ms, vth, epsmn, bfac, efac;
    int m0, mph, comptorque, magdrift, nopassing, noshear, pertfile, nonlin;
    int inp_swi, vsteps, log_level;
};

static int parse_bool(const char *v)
{
    return (v[0] == 't' || v[0] == 'T' || strstr(v, ".true.") || strstr(v, ".TRUE."));
}

static void read_config(const char *path, struct config *c)
{
    FILE *f = fopen(path, "r");
    if (!f) {
        fprintf(stderr, "cannot open %s\n", path);
        exit(1);
    }
    /* defaults match the Fortran module defaults (driftorbit/orbit); the namelist
     * only overrides keys it contains. comptorque/magdrift default true. */
    c->bfac = c->efac = 1.0;
    c->s = c->M_t = c->qs = c->ms = c->vth = c->epsmn = 0.0;
    c->m0 = c->mph = c->inp_swi = c->vsteps = c->log_level = 0;
    c->comptorque = 1;
    c->magdrift = 1;
    c->nopassing = c->noshear = c->pertfile = c->nonlin = 0;

    char line[512];
    while (fgets(line, sizeof line, f)) {
        char *eq = strchr(line, '=');
        if (!eq)
            continue;
        char key[64];
        int n = 0;
        for (char *p = line; p < eq && n < 63; p++)
            if (!isspace((unsigned char)*p))
                key[n++] = (char)tolower((unsigned char)*p);
        key[n] = 0;
        char *val = eq + 1;
        while (*val && isspace((unsigned char)*val))
            val++;
        /* strip trailing comment */
        char *bang = strchr(val, '!');
        if (bang)
            *bang = 0;

#define DI(name, field) else if (!strcmp(key, name)) c->field = atoi(val)
#define DF(name, field) else if (!strcmp(key, name)) c->field = atof(val)
#define DB(name, field) else if (!strcmp(key, name)) c->field = parse_bool(val)
        if (!strcmp(key, "s")) c->s = atof(val);
        DF("m_t", M_t); DF("qs", qs); DF("ms", ms); DF("vth", vth);
        DF("epsmn", epsmn); DF("bfac", bfac); DF("efac", efac);
        DI("m0", m0); DI("mph", mph); DI("inp_swi", inp_swi);
        DI("vsteps", vsteps); DI("log_level", log_level);
        DB("comptorque", comptorque); DB("magdrift", magdrift);
        DB("nopassing", nopassing); DB("noshear", noshear);
        DB("pertfile", pertfile); DB("nonlin", nonlin);
    }
    fclose(f);
}

/* ---- per-harmonic transport accumulation ---- */
struct harmonic {
    int mth;
    double Dresco[2], Dresctr[2], Drest[2];
    double Tresco, Tresctr, Trest;
    double vminp_vth, vmaxp_vth, vmint_vth, vmaxt_vth;
};

static void set_to_trapped_region(double *emin, double *emax)
{
    *emin = (1.0 + DO_EPST) * do_etatp;
    *emax = (1.0 - DO_EPST) * do_etadt;
}
static void set_to_passing_region(double *emin, double *emax)
{
    *emin = DO_EPSP * do_etatp;
    *emax = (1.0 - DO_EPSP) * do_etatp;
}

/* ---- main computation, mirrors neort_compute_at_s + compute_transport ---- */
int main(int argc, char **argv)
{
    if (argc < 2) {
        fprintf(stderr, "usage: neo_rt_c <runname>\n");
        return 2;
    }
    const char *run = argv[1];
    char cfgpath[1024];
    snprintf(cfgpath, sizeof cfgpath, "%s.in", run);
    struct config c;
    read_config(cfgpath, &c);

    /* set globals from config */
    field_inp_swi = c.inp_swi;
    field_bfac = c.bfac;
    do_epsmn = c.epsmn;
    do_m0 = c.m0;
    do_efac = c.efac;
    do_comptorque = c.comptorque;
    do_magdrift = c.magdrift;
    do_nopassing = c.nopassing;
    do_pertfile = c.pertfile;
    do_nonlin = c.nonlin;
    orbit_noshear = c.noshear;
    int vsteps = c.vsteps;
    double s = c.s;

    /* init field + perturbation */
    do_magfie_init("in_file");
    if (do_pertfile) {
        do_magfie_pert_init("in_file_pert");
        do_mph = field_pert_mph;
    } else {
        do_mph = c.mph;
    }

    /* plasma + rotation profiles (spline mode) */
    read_and_init_plasma_input("plasma.in", s);
    read_and_init_profile_input("profile.in", s, field_R0, c.efac, c.bfac);

    /* geometry + frequency splines */
    init_flux_surface_average(s);
    orbit_set_s(s);
    init_canon_freq_trapped_spline();
    if (!do_nopassing)
        init_canon_freq_passing_spline();
    do_sign_vpar = 1.0;
    set_to_trapped_region(&do_etamin, &do_etamax);
    if (do_comptorque)
        init_thermodynamic_forces(field_psi_pr, field_q);

    /* recompute Om_tE/dOm_tEds (compute_transport head) */
    prof_Om_tE = prof_vth * prof_M_t / field_R0;
    prof_dOm_tEds = prof_vth * prof_dM_tds / field_R0 + prof_M_t * prof_dvthds / field_R0;

    /* harmonic loop */
    int mthmin = -(int)ceil(2.0 * fabs(do_mph * field_q));
    int mthmax = (int)ceil(2.0 * fabs(do_mph * field_q));
    int nh = (mthmax >= mthmin) ? (mthmax - mthmin + 1) : 0;
    struct harmonic *H = (struct harmonic *)calloc(nh > 0 ? nh : 1, sizeof *H);
    double Dco[2] = {0, 0}, Dctr[2] = {0, 0}, Dt[2] = {0, 0};
    double Tco = 0, Tctr = 0, Tt = 0;
    double vth = prof_vth;

    int idx = 0;
    for (int j = mthmin; j <= mthmax; j++, idx++) {
        do_mth = j;
        struct harmonic *h = &H[idx];
        h->mth = j;
        double vminp = 1.0e-6 * vth, vmaxp = 3.0 * vth, vmint = vminp, vmaxt = vmaxp;
        double D[2], T;

        if (!do_nopassing) {
            do_sign_vpar = 1.0;
            set_to_passing_region(&do_etamin, &do_etamax);
            compute_transport_integral(vminp, vmaxp, vsteps, D, &T);
            h->Dresco[0] = D[0]; h->Dresco[1] = D[1]; h->Tresco = T;
            Dco[0] += D[0]; Dco[1] += D[1]; Tco += T;

            do_sign_vpar = -1.0;
            set_to_passing_region(&do_etamin, &do_etamax);
            compute_transport_integral(vminp, vmaxp, vsteps, D, &T);
            h->Dresctr[0] = D[0]; h->Dresctr[1] = D[1]; h->Tresctr = T;
            Dctr[0] += D[0]; Dctr[1] += D[1]; Tctr += T;
        }
        do_sign_vpar = 1.0;
        set_to_trapped_region(&do_etamin, &do_etamax);
        compute_transport_integral(vmint, vmaxt, vsteps, D, &T);
        h->Drest[0] = D[0]; h->Drest[1] = D[1]; h->Trest = T;
        Dt[0] += D[0]; Dt[1] += D[1]; Tt += T;

        h->vminp_vth = vminp / vth; h->vmaxp_vth = vmaxp / vth;
        h->vmint_vth = vmint / vth; h->vmaxt_vth = vmaxt / vth;
    }

    /* ---- write transport output files ---- */
    char path[1100];
    FILE *fp;
    double tot1 = Dco[0] + Dctr[0] + Dt[0], tot2 = Dco[1] + Dctr[1] + Dt[1];

    snprintf(path, sizeof path, "%s.out", run);
    fp = fopen(path, "w");
    fprintf(fp, " # M_t D11co D11ctr D11t D11 D12co D12ctr D12t D12\n");
    fprintf(fp, " %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
            prof_M_t, Dco[0], Dctr[0], Dt[0], tot1, Dco[1], Dctr[1], Dt[1], tot2);
    fclose(fp);

    if (do_comptorque) {
        snprintf(path, sizeof path, "%s_torque.out", run);
        fp = fopen(path, "w");
        fprintf(fp, " # s dVds M_t Tco Tctr Tt\n");
        fprintf(fp, " %.17e %.17e %.17e %.17e %.17e %.17e\n", s, do_dVds, prof_M_t,
                Tco, Tctr, Tt);
        fclose(fp);
    }

    snprintf(path, sizeof path, "%s_integral.out", run);
    fp = fopen(path, "w");
    char path2[1100];
    snprintf(path2, sizeof path2, "%s_torque_integral.out", run);
    FILE *fp2 = fopen(path2, "w");
    for (int k = 0; k < nh; k++) {
        struct harmonic *h = &H[k];
        double t1 = h->Dresco[0] + h->Dresctr[0] + h->Drest[0];
        double t2 = h->Dresco[1] + h->Dresctr[1] + h->Drest[1];
        fprintf(fp, " %.17e %d %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
                prof_M_t, h->mth, h->Dresco[0], h->Dresctr[0], h->Drest[0], t1,
                h->Dresco[1], h->Dresctr[1], h->Drest[1], t2, h->vminp_vth,
                h->vmaxp_vth, h->vmint_vth, h->vmaxt_vth);
        fprintf(fp2, " %d %.17e %.17e %.17e\n", h->mth, h->Tresco, h->Tresctr, h->Trest);
    }
    fclose(fp);
    fclose(fp2);

    /* ---- check_magfie -> _magfie.out (50 theta samples) ---- */
    snprintf(path, sizeof path, "%s_magfie.out", run);
    fp = fopen(path, "w");
    int nth = 50;
    for (int k = 0; k < nth; k++) {
        double th = -NEORT_PI + k * (2.0 * NEORT_PI) / (nth - 1);
        double x[3] = {s, 0.0, th}, bmod, sqrtg, hder[3], hcovar[3], hctrvr[3], hcurl[3];
        do_magfie(x, &bmod, &sqrtg, hder, hcovar, hctrvr, hcurl);
        double _Complex bn;
        if (do_pertfile) {
            do_magfie_pert_amp(x, &bn);
            bn = do_epsmn * bn / bmod;
        } else {
            bn = do_epsmn * cexp(I * do_m0 * th);
        }
        double _Complex eps_exp = do_epsmn * cexp(I * do_m0 * th);
        fprintf(fp,
                " %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e %.17e\n",
                th, bmod, sqrtg, hder[0], hder[1], hder[2], hcovar[0], hcovar[1],
                hcovar[2], hctrvr[0], hctrvr[1], hctrvr[2], hcurl[0], hcurl[1],
                hcurl[2], creal(bn), cimag(bn), creal(eps_exp), cimag(eps_exp));
    }
    fclose(fp);

    free(H);
    return 0;
}
