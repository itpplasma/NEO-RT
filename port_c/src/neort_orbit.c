#include "neort_orbit.h"
#include "neort_field.h"
#include "neort_driftorbit.h"
#include "neort_util.h"
#include <math.h>
#include <stdio.h>
#include <cvode/cvode.h>
#include <nvector/nvector_serial.h>
#include <sunnonlinsol/sunnonlinsol_fixedpoint.h>
#include <sundials/sundials_context.h>

double orbit_th0 = 0.0;
int orbit_noshear = 0;

static double orbit_s = 0.0;
void orbit_set_s(double s) { orbit_s = s; }
double orbit_get_s(void) { return orbit_s; }

double vpar(double v, double eta, double bmod)
{
    double r = v * sqrt(1.0 - eta * bmod);
    return isnan(r) ? 0.0 : r;
}
double vperp(double v, double eta, double bmod)
{
    double r = v * sqrt(eta * bmod);
    return isnan(r) ? 0.0 : r;
}

static void evaluate_bfield_local(double *bmod, double *htheta)
{
    double x[3] = {orbit_s, 0.0, orbit_th0}, sqrtg, hder[3], hcovar[3], hctrvr[3], hcurl[3];
    do_magfie(x, bmod, &sqrtg, hder, hcovar, hctrvr, hcurl);
    *htheta = hctrvr[2];
}

/* ---- ODE right-hand sides (mirror timestep / timestep_poloidal_motion) ---- */
typedef void (*ts_fn)(double v, double eta, int neq, double t, const double *y,
                      double *ydot);

void poloidal_velocity(double v, double eta, double bmod, double hthctr,
                       double hderth, double v_par, double ydot[2])
{
    ydot[0] = v_par * hthctr;
    ydot[1] = -v * v * eta / 2.0 * hthctr * hderth * bmod;
}

static void timestep_poloidal_motion(double v, double eta, int neq, double t,
                                     const double *y, double *ydot)
{
    (void)neq; (void)t;
    double x[3] = {orbit_s, 0.0, y[0]}, bmod, sqrtg, hder[3], hcovar[3], hctrvr[3], hcurl[3];
    do_magfie(x, &bmod, &sqrtg, hder, hcovar, hctrvr, hcurl);
    poloidal_velocity(v, eta, bmod, hctrvr[2], hder[2], y[1], ydot);
}

static void timestep(double v, double eta, int neq, double t, const double *y,
                     double *ydot)
{
    (void)t;
    double x[3] = {orbit_s, 0.0, y[0]}, bmod, sqrtg, hder[3], hcovar[3], hctrvr[3], hcurl[3];
    do_magfie(x, &bmod, &sqrtg, hder, hcovar, hctrvr, hcurl);
    double shearterm = orbit_noshear ? 0.0 : field_Bphcov * field_dqds;
    double Om_tB_v = util_mi * NEORT_C * field_q /
                     (2.0 * util_qi * FIELD_SIGN_THETA * field_psi_pr * bmod) *
                     (-(2.0 - eta * bmod) * bmod * hder[0] +
                      2.0 * (1.0 - eta * bmod) * hctrvr[2] *
                          (field_dBthcovds + field_q * field_dBphcovds + shearterm));
    ydot[0] = y[1] * hctrvr[2];
    ydot[1] = -0.5 * v * v * eta * hctrvr[2] * hder[2] * bmod;
    ydot[2] = Om_tB_v;
    for (int i = 3; i < neq; i++)
        ydot[i] = 0.0;
}

/* ---- CVODE glue ---- */
static double g_sign_vpar_htheta = 1.0;

struct ode_ctx {
    double v, eta;
    int neq;
    ts_fn ts;
};

static int cv_rhs(sunrealtype t, N_Vector y, N_Vector ydot, void *ud)
{
    struct ode_ctx *c = (struct ode_ctx *)ud;
    c->ts(c->v, c->eta, c->neq, t, N_VGetArrayPointer(y), N_VGetArrayPointer(ydot));
    return 0;
}

static int cv_roots(sunrealtype t, N_Vector y, sunrealtype *gout, void *ud)
{
    (void)t; (void)ud;
    double th = N_VGetArrayPointer(y)[0];
    gout[0] = g_sign_vpar_htheta * (th - orbit_th0);
    gout[1] = g_sign_vpar_htheta * (2.0 * NEORT_PI - (th - orbit_th0));
    return 0;
}

/* Mirrors bounce_integral: integrate in dt chunks with root-finding, accept the
 * turning-point event. Returns ti in out[0], final state in out[1..neq]. */
static void bounce_integral(double v, double eta, int neq, const double *y0,
                            double dt, ts_fn ts, double *out)
{
    SUNContext sctx;
    SUNContext_Create(SUN_COMM_NULL, &sctx);
    N_Vector y = N_VNew_Serial(neq, sctx);
    double *yv = N_VGetArrayPointer(y);
    for (int i = 0; i < neq; i++)
        yv[i] = y0[i];

    struct ode_ctx oc = {v, eta, neq, ts};
    void *cv = CVodeCreate(CV_ADAMS, sctx);
    CVodeInit(cv, cv_rhs, 0.0, y);
    CVodeSStolerances(cv, 1.0e-9, 1.0e-10);
    CVodeSetMaxNumSteps(cv, 50000);
    CVodeSetUserData(cv, &oc);
    SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(y, 0, sctx);
    CVodeSetNonlinearSolver(cv, NLS);
    CVodeRootInit(cv, 2, cv_roots);

    int passing = (eta < do_etatp);
    double ti = 0.0;
    double yold[ORBIT_NVAR];
    for (int k = 2; k <= 500; k++) {
        for (int i = 0; i < neq; i++)
            yold[i] = yv[i];
        double tout = ti + dt;
        int flag = CVode(cv, tout, y, &ti, CV_NORMAL);
        if (flag < 0) {
            fprintf(stderr, "CVODE error %d in bounce_integral\n", flag);
            break;
        }
        if (flag == CV_ROOT_RETURN) {
            if (passing || (yold[0] - orbit_th0) < 0.0)
                break;
        }
    }
    out[0] = ti;
    for (int i = 0; i < neq; i++)
        out[1 + i] = yv[i];

    SUNNonlinSolFree(NLS);
    N_VDestroy(y);
    CVodeFree(&cv);
    SUNContext_Free(&sctx);
}

void bounce_fast_ext(double v, double eta, double taub, double bounceavg[ORBIT_NVAR],
                     orbit_ts_fn ts, int *istate_out)
{
    double bmod, htheta;
    evaluate_bfield_local(&bmod, &htheta);
    g_sign_vpar_htheta = (htheta >= 0 ? 1.0 : -1.0) * do_sign_vpar;
    do_sign_vpar_htheta = g_sign_vpar_htheta;

    SUNContext sctx;
    SUNContext_Create(SUN_COMM_NULL, &sctx);
    N_Vector y = N_VNew_Serial(ORBIT_NVAR, sctx);
    double *yv = N_VGetArrayPointer(y);
    for (int i = 0; i < ORBIT_NVAR; i++)
        yv[i] = 1.0e-15;
    yv[0] = orbit_th0;
    yv[1] = g_sign_vpar_htheta * vpar(v, eta, bmod);
    for (int i = 2; i < 6; i++)
        yv[i] = 0.0;

    struct ode_ctx oc = {v, eta, ORBIT_NVAR, ts};
    void *cv = CVodeCreate(CV_ADAMS, sctx);
    CVodeInit(cv, cv_rhs, 0.0, y);
    CVodeSStolerances(cv, 1.0e-9, 1.0e-10);
    CVodeSetMaxNumSteps(cv, 50000);
    CVodeSetUserData(cv, &oc);
    SUNNonlinearSolver NLS = SUNNonlinSol_FixedPoint(y, 0, sctx);
    CVodeSetNonlinearSolver(cv, NLS);

    double tret;
    int flag = CVode(cv, taub, y, &tret, CV_NORMAL);
    for (int i = 0; i < ORBIT_NVAR; i++)
        bounceavg[i] = yv[i] / taub;
    if (istate_out)
        *istate_out = (flag < 0) ? -1 : 2; /* mirror DVODE istate (2=ok, -1=MXSTEP) */

    SUNNonlinSolFree(NLS);
    N_VDestroy(y);
    CVodeFree(&cv);
    SUNContext_Free(&sctx);
}

void bounce_fast(double v, double eta, double taub, double bounceavg[ORBIT_NVAR])
{
    bounce_fast_ext(v, eta, taub, bounceavg, timestep, NULL);
}

double bounce_time(double v, double eta, double taub_estimate, int have_estimate)
{
    double bmod, htheta;
    evaluate_bfield_local(&bmod, &htheta);
    g_sign_vpar_htheta = (htheta >= 0 ? 1.0 : -1.0) * do_sign_vpar;
    do_sign_vpar_htheta = g_sign_vpar_htheta;

    double y0[2] = {orbit_th0, g_sign_vpar_htheta * vpar(v, eta, bmod)};
    double taub = have_estimate
                      ? taub_estimate
                      : 2.0 * NEORT_PI /
                            fabs(vperp(v, eta, bmod) * field_iota / field_R0 *
                                 sqrt(field_eps / 2.0));
    double out[3];
    bounce_integral(v, eta, 2, y0, taub, timestep_poloidal_motion, out);
    return out[0];
}

void bounce(double v, double eta, double *taub, double bounceavg[ORBIT_NVAR],
            double taub_estimate, int have_estimate)
{
    double bmod, htheta;
    evaluate_bfield_local(&bmod, &htheta);
    g_sign_vpar_htheta = (htheta >= 0 ? 1.0 : -1.0) * do_sign_vpar;
    do_sign_vpar_htheta = g_sign_vpar_htheta;

    double y0[ORBIT_NVAR];
    for (int i = 0; i < ORBIT_NVAR; i++)
        y0[i] = 1.0e-15;
    y0[0] = orbit_th0;
    y0[1] = g_sign_vpar_htheta * vpar(v, eta, bmod);
    y0[2] = y0[3] = y0[4] = y0[5] = 0.0;

    double tb = have_estimate
                    ? taub_estimate
                    : 2.0 * NEORT_PI /
                          fabs(vperp(v, eta, bmod) * field_iota / field_R0 *
                               sqrt(field_eps / 2.0));
    double out[ORBIT_NVAR + 1];
    bounce_integral(v, eta, ORBIT_NVAR, y0, tb / 5.0, timestep, out);
    *taub = out[0];
    for (int i = 0; i < ORBIT_NVAR; i++)
        bounceavg[i] = out[1 + i] / out[0];
}
