#include "neort_freq.h"
#include "neort_orbit.h"
#include "neort_field.h"
#include "neort_profiles.h"
#include "neort_driftorbit.h"
#include "neort_util.h"
#include "neort_spline.h"
#include <math.h>

/* trapped-region splines */
static double Omth_spl[(FREQ_NETASPL - 1) * 5];
static double OmtB_spl[(FREQ_NETASPL - 1) * 5];
/* passing-region splines */
static double Omth_pass_spl[(FREQ_NETASPL_PASS - 1) * 5];
static double OmtB_pass_spl[(FREQ_NETASPL_PASS - 1) * 5];

static double k_taub_p = 0, d_taub_p = 0, k_taub_t = 0, d_taub_t = 0;
static double k_OmtB_p = 0, d_OmtB_p = 0, k_OmtB_t = 0, d_OmtB_t = 0;

static int js_omth_t = 1, js_omtb_t = 1, js_omth_p = 1, js_omtb_p = 1;

void init_canon_freq_trapped_spline(void)
{
    const int n = FREQ_NETASPL;
    double etarange[FREQ_NETASPL], Om_tB_v[FREQ_NETASPL], Omth_v[FREQ_NETASPL];
    double taub0 = 0, taub1 = 0, leta0 = 0, leta1 = 0, OmtB0 = 0, OmtB1 = 0;
    double v = prof_vth, taub_est = 0;

    do_etamin = (1.0 + DO_EPST) * do_etatp;
    do_etamax = do_etatp + (do_etadt - do_etatp) * (1.0 - DO_EPSST_SPL);

    double b = log(DO_EPST_SPL);
    double aa = 1.0 / (n - 1.0) * (log(do_etamax / do_etamin - 1.0) - b);

    for (int k = n - 1; k >= 0; k--) {
        double eta = do_etamin * (1.0 + exp(aa * k + b));
        etarange[k] = eta;
        taub_est = (k == n - 1) ? bounce_time(v, eta, 0.0, 0)
                                : bounce_time(v, eta, taub_est, 1);
        double taub = taub_est, bavg[ORBIT_NVAR];
        bounce_fast(v, eta, taub, bavg);
        if (do_magdrift)
            Om_tB_v[k] = bavg[2];
        Omth_v[k] = 2.0 * NEORT_PI / (v * taub);
        if (k == 0) {
            leta0 = log(eta - do_etatp);
            taub0 = v * taub;
            if (do_magdrift)
                OmtB0 = Om_tB_v[k] / Omth_v[k];
        }
        if (k == 1) {
            leta1 = log(eta - do_etatp);
            taub1 = v * taub;
            if (do_magdrift)
                OmtB1 = Om_tB_v[k] / Omth_v[k];
        }
    }
    k_taub_t = (taub1 - taub0) / (leta1 - leta0);
    d_taub_t = taub0 - k_taub_t * leta0;
    spline_coeff(etarange, Omth_v, n, Omth_spl);
    if (do_magdrift) {
        k_OmtB_t = (OmtB1 - OmtB0) / (leta1 - leta0);
        d_OmtB_t = OmtB0 - k_OmtB_t * leta0;
        spline_coeff(etarange, Om_tB_v, n, OmtB_spl);
    }
}

void init_canon_freq_passing_spline(void)
{
    const int n = FREQ_NETASPL_PASS;
    double etarange[FREQ_NETASPL_PASS], Om_tB_v[FREQ_NETASPL_PASS], Omth_v[FREQ_NETASPL_PASS];
    double taub0 = 0, taub1 = 0, leta0 = 0, leta1 = 0, OmtB0 = 0, OmtB1 = 0;
    double v = prof_vth, taub_est = 0;

    do_etamin = do_etatp * DO_EPSSP_SPL;
    do_etamax = do_etatp;

    double b = log((do_etamax - do_etamin) / do_etamax);
    double aa = 1.0 / (n - 1.0) * (log(DO_EPSP_SPL) - b);

    for (int k = n - 1; k >= 0; k--) {
        double eta = do_etamax * (1.0 - exp(aa * k + b));
        etarange[k] = eta;
        taub_est = (k == n - 1) ? bounce_time(v, eta, 0.0, 0)
                                : bounce_time(v, eta, taub_est, 1);
        double taub = taub_est, bavg[ORBIT_NVAR];
        bounce_fast(v, eta, taub, bavg);
        if (do_magdrift)
            Om_tB_v[k] = bavg[2];
        Omth_v[k] = 2.0 * NEORT_PI / (v * taub);
        if (k == n - 2) {
            leta0 = log(do_etatp - eta);
            taub0 = v * taub;
            if (do_magdrift)
                OmtB0 = Om_tB_v[k] / Omth_v[k];
        }
        if (k == n - 1) {
            leta1 = log(do_etatp - eta);
            taub1 = v * taub;
            if (do_magdrift)
                OmtB1 = Om_tB_v[k] / Omth_v[k];
        }
    }
    k_taub_p = (taub1 - taub0) / (leta1 - leta0);
    d_taub_p = taub0 - k_taub_p * leta0;
    spline_coeff(etarange, Omth_v, n, Omth_pass_spl);
    if (do_magdrift) {
        k_OmtB_p = (OmtB1 - OmtB0) / (leta1 - leta0);
        d_OmtB_p = OmtB0 - k_OmtB_p * leta0;
        spline_coeff(etarange, Om_tB_v, n, OmtB_pass_spl);
    }
}

void Om_th(double v, double eta, double *Omth, double *dOmthdv, double *dOmthdeta)
{
    double sv[3];
    if (eta > do_etatp) {
        if (eta > do_etatp * (1.0 + DO_EPST_SPL)) {
            spline_val_0(Omth_spl, FREQ_NETASPL - 1, eta, &js_omth_t, sv);
        } else {
            sv[0] = 2.0 * NEORT_PI / (k_taub_t * log(eta - do_etatp) + d_taub_t);
            sv[1] = -sv[0] * sv[0] / (2.0 * NEORT_PI) * k_taub_t / (eta - do_etatp);
        }
    } else {
        if (eta < do_etatp * (1.0 - DO_EPSP_SPL)) {
            spline_val_0(Omth_pass_spl, FREQ_NETASPL_PASS - 1, eta, &js_omth_p, sv);
        } else {
            sv[0] = 2.0 * NEORT_PI / (k_taub_p * log(do_etatp - eta) + d_taub_p);
            sv[1] = -sv[0] * sv[0] / (2.0 * NEORT_PI) * k_taub_p / (eta - do_etatp);
        }
    }
    *Omth = do_sign_vpar * sv[0] * v;
    *dOmthdv = do_sign_vpar * sv[0];
    *dOmthdeta = do_sign_vpar * sv[1] * v;
}

void Om_tB(double v, double eta, double *OmtB, double *dOmtBdv, double *dOmtBdeta)
{
    double sv[3], Omth, dOmthdv, dOmthdeta;
    if (eta > do_etatp) {
        if (eta > do_etatp * (1.0 + DO_EPST_SPL)) {
            spline_val_0(OmtB_spl, FREQ_NETASPL - 1, eta, &js_omtb_t, sv);
        } else {
            Om_th(v, eta, &Omth, &dOmthdv, &dOmthdeta);
            sv[0] = do_sign_vpar * (k_OmtB_t * log(eta - do_etatp) + d_OmtB_t) * Omth / v;
            sv[1] = do_sign_vpar * (Omth / v * k_OmtB_t / (eta - do_etatp) +
                                    dOmthdeta / v *
                                        (k_OmtB_t * log(eta - do_etatp) + d_OmtB_t));
        }
    } else {
        if (eta < do_etatp * (1.0 - DO_EPSP_SPL)) {
            spline_val_0(OmtB_pass_spl, FREQ_NETASPL_PASS - 1, eta, &js_omtb_p, sv);
        } else {
            Om_th(v, eta, &Omth, &dOmthdv, &dOmthdeta);
            sv[0] = do_sign_vpar * (k_OmtB_p * log(do_etatp - eta) + d_OmtB_p) * Omth / v;
            sv[1] = do_sign_vpar * (Omth / v * k_OmtB_p / (eta - do_etatp) +
                                    dOmthdeta / v *
                                        (k_OmtB_p * log(do_etatp - eta) + d_OmtB_p));
        }
    }
    *OmtB = sv[0] * v * v;
    *dOmtBdv = 2.0 * sv[0] * v;
    *dOmtBdeta = sv[1] * v * v;
}

void Om_ph(double v, double eta, double *Omph, double *dOmphdv, double *dOmphdeta)
{
    double Omth, dOmthdv, dOmthdeta, OmtB, dOmtBdv, dOmtBdeta;
    if (eta > do_etatp) {
        *Omph = prof_Om_tE;
        *dOmphdv = 0.0;
        *dOmphdeta = 0.0;
        if (do_magdrift) {
            Om_tB(v, eta, &OmtB, &dOmtBdv, &dOmtBdeta);
            *Omph += OmtB;
            *dOmphdv += dOmtBdv;
            *dOmphdeta += dOmtBdeta;
        }
    } else {
        Om_th(v, eta, &Omth, &dOmthdv, &dOmthdeta);
        *Omph = prof_Om_tE + Omth / field_iota;
        *dOmphdv = dOmthdv / field_iota;
        *dOmphdeta = dOmthdeta / field_iota;
        if (do_magdrift) {
            Om_tB(v, eta, &OmtB, &dOmtBdv, &dOmtBdeta);
            *Omph += OmtB;
            *dOmphdv += dOmtBdv;
            *dOmphdeta += dOmtBdeta;
        }
    }
}

void d_Om_ds(double v, double eta, double taub_estimate, double *dOmthds,
             double *dOmphds)
{
    /* current s is held by the orbit module; nudge it by +-ds/2. */
    double ds = 2.0e-8;
    double s0 = orbit_get_s();

    orbit_set_s(s0 - ds / 2.0);
    double taub = bounce_time(v, eta, taub_estimate, 1), bavg[ORBIT_NVAR];
    bounce_fast(v, eta, taub, bavg);
    double Omth = do_sign_vpar_htheta * 2.0 * NEORT_PI / taub;
    double Omph_noE;
    if (do_magdrift)
        Omph_noE = (eta > do_etatp) ? bavg[2] * v * v
                                    : bavg[2] * v * v + Omth / field_iota;
    else
        Omph_noE = (eta > do_etatp) ? 0.0 : Omth / field_iota;

    orbit_set_s(s0 + ds / 2.0);
    taub = bounce_time(v, eta, taub_estimate, 1);
    bounce_fast(v, eta, taub, bavg);
    *dOmthds = do_sign_vpar_htheta *
               (2.0 * NEORT_PI / taub - do_sign_vpar_htheta * Omth) / ds;
    if (do_magdrift) {
        if (eta > do_etatp)
            *dOmphds = prof_dOm_tEds + (bavg[2] * v * v - Omph_noE) / ds;
        else
            *dOmphds = prof_dOm_tEds +
                       (bavg[2] * v * v + (2.0 * NEORT_PI / taub) / field_iota - Omph_noE) / ds;
    } else {
        if (eta > do_etatp)
            *dOmphds = prof_dOm_tEds;
        else
            *dOmphds = prof_dOm_tEds +
                       ((2.0 * NEORT_PI / taub) / field_iota - Omph_noE) / ds;
    }
    orbit_set_s(s0);
}
