//! Orbit bounce integration, literal port of neort_orbit (orbit.f90), using
//! SUNDIALS CVODE (CV_ADAMS + fixed-point + root finding) -- the same library
//! and algorithm as the C port, so trajectories are bit-identical to C.

use crate::cvode::*;
use crate::driftorbit as dorb;
use crate::field::{do_magfie, with_field, SIGN_THETA};
use crate::profiles;
use std::cell::Cell;
use std::os::raw::{c_int, c_long, c_void};
use std::ptr;

pub const NVAR: usize = 7;
const PI: f64 = std::f64::consts::PI;

thread_local! {
    static ORBIT_S: Cell<f64> = Cell::new(0.0);
}
pub fn set_s(s: f64) {
    ORBIT_S.with(|c| c.set(s));
}
pub fn get_s() -> f64 {
    ORBIT_S.with(|c| c.get())
}

pub fn vpar(v: f64, eta: f64, bmod: f64) -> f64 {
    let r = v * (1.0 - eta * bmod).sqrt();
    if r.is_nan() {
        0.0
    } else {
        r
    }
}
pub fn vperp(v: f64, eta: f64, bmod: f64) -> f64 {
    let r = v * (eta * bmod).sqrt();
    if r.is_nan() {
        0.0
    } else {
        r
    }
}

fn evaluate_bfield_local() -> (f64, f64) {
    let s = get_s();
    let th0 = dorb::get().th0;
    let (bmod, _sqrtg, _bder, _hcov, hctrvr) = do_magfie(&[s, 0.0, th0]);
    (bmod, hctrvr[2])
}

pub fn poloidal_velocity(v: f64, eta: f64, bmod: f64, hthctr: f64, hderth: f64, v_par: f64, ydot: &mut [f64]) {
    ydot[0] = v_par * hthctr;
    ydot[1] = -v * v * eta / 2.0 * hthctr * hderth * bmod;
}

pub fn timestep_poloidal_motion(v: f64, eta: f64, _neq: usize, _t: f64, y: &[f64], ydot: &mut [f64]) {
    let s = get_s();
    let (bmod, _sqrtg, bder, _hcov, hctrvr) = do_magfie(&[s, 0.0, y[0]]);
    poloidal_velocity(v, eta, bmod, hctrvr[2], bder[2], y[1], ydot);
}

pub fn timestep(v: f64, eta: f64, neq: usize, _t: f64, y: &[f64], ydot: &mut [f64]) {
    let s = get_s();
    let (bmod, _sqrtg, bder, _hcov, hctrvr) = do_magfie(&[s, 0.0, y[0]]);
    let d = dorb::get();
    let (q, dqds, psi_pr, bphcov, dbthcovds, dbphcovds) =
        with_field(|f| (f.q, f.dqds, f.psi_pr, f.bphcov, f.dbthcovds, f.dbphcovds));
    let p = profiles::get();
    let shearterm = if d.noshear { 0.0 } else { bphcov * dqds };
    let om_tb_v = p.mi * profiles::C * q / (2.0 * p.qi * SIGN_THETA * psi_pr * bmod)
        * (-(2.0 - eta * bmod) * bmod * bder[0]
            + 2.0 * (1.0 - eta * bmod) * hctrvr[2] * (dbthcovds + q * dbphcovds + shearterm));
    ydot[0] = y[1] * hctrvr[2];
    ydot[1] = -0.5 * v * v * eta * hctrvr[2] * bder[2] * bmod;
    ydot[2] = om_tb_v;
    for i in 3..neq {
        ydot[i] = 0.0;
    }
}

pub type TsFn = fn(f64, f64, usize, f64, &[f64], &mut [f64]);

struct OdeCtx {
    v: f64,
    eta: f64,
    neq: usize,
    ts: TsFn,
}

extern "C" fn cv_rhs(t: f64, y: NVector, ydot: NVector, ud: *mut c_void) -> c_int {
    let ctx = unsafe { &*(ud as *const OdeCtx) };
    let yp = unsafe { N_VGetArrayPointer(y) };
    let ydp = unsafe { N_VGetArrayPointer(ydot) };
    let ys = unsafe { std::slice::from_raw_parts(yp, ctx.neq) };
    let yds = unsafe { std::slice::from_raw_parts_mut(ydp, ctx.neq) };
    (ctx.ts)(ctx.v, ctx.eta, ctx.neq, t, ys, yds);
    0
}

extern "C" fn cv_roots(_t: f64, y: NVector, gout: *mut f64, _ud: *mut c_void) -> c_int {
    let d = dorb::get();
    let th = unsafe { *N_VGetArrayPointer(y) };
    unsafe {
        *gout = d.sign_vpar_htheta * (th - d.th0);
        *gout.add(1) = d.sign_vpar_htheta * (2.0 * PI - (th - d.th0));
    }
    0
}

/// bounce_integral: integrate in dt chunks with root finding; returns (ti, y).
fn bounce_integral(v: f64, eta: f64, neq: usize, y0: &[f64], dt: f64, ts: TsFn) -> (f64, Vec<f64>) {
    let etatp = dorb::get().etatp;
    let th0 = dorb::get().th0;
    unsafe {
        let mut ctx: SunContext = ptr::null_mut();
        SUNContext_Create(sun_comm_null(), &mut ctx);
        let y = N_VNew_Serial(neq as c_long, ctx);
        let yv = N_VGetArrayPointer(y);
        for i in 0..neq {
            *yv.add(i) = y0[i];
        }
        let oc = Box::new(OdeCtx { v, eta, neq, ts });
        let mem = CVodeCreate(CV_ADAMS, ctx);
        CVodeInit(mem, cv_rhs, 0.0, y);
        CVodeSStolerances(mem, 1.0e-9, 1.0e-10);
        CVodeSetMaxNumSteps(mem, 50000);
        CVodeSetUserData(mem, &*oc as *const OdeCtx as *mut c_void);
        let nls = SUNNonlinSol_FixedPoint(y, 0, ctx);
        CVodeSetNonlinearSolver(mem, nls);
        CVodeRootInit(mem, 2, cv_roots);

        let passing = eta < etatp;
        let mut ti = 0.0f64;
        let mut yold = vec![0.0f64; neq];
        for _k in 2..=500 {
            for i in 0..neq {
                yold[i] = *yv.add(i);
            }
            let tout = ti + dt;
            let flag = CVode(mem, tout, y, &mut ti, CV_NORMAL);
            if flag < 0 {
                eprintln!("CVODE error {} in bounce_integral", flag);
                break;
            }
            if flag == CV_ROOT_RETURN && (passing || (yold[0] - th0) < 0.0) {
                break;
            }
        }
        let out: Vec<f64> = (0..neq).map(|i| *yv.add(i)).collect();
        SUNNonlinSolFree(nls);
        N_VDestroy(y);
        let mut m = mem;
        CVodeFree(&mut m);
        SUNContext_Free(&mut ctx);
        drop(oc);
        (ti, out)
    }
}

fn set_sign_vpar_htheta() -> f64 {
    let (_bmod, htheta) = evaluate_bfield_local();
    let d = dorb::get();
    let svh = (if htheta >= 0.0 { 1.0 } else { -1.0 }) * d.sign_vpar;
    dorb::update(|x| x.sign_vpar_htheta = svh);
    svh
}

pub fn bounce_time(v: f64, eta: f64, taub_estimate: f64, have: bool) -> f64 {
    let (bmod, _) = evaluate_bfield_local();
    let svh = set_sign_vpar_htheta();
    let th0 = dorb::get().th0;
    let y0 = [th0, svh * vpar(v, eta, bmod)];
    let taub = if have {
        taub_estimate
    } else {
        let (iota, r0, eps) = with_field(|f| (f.iota, f.r0, f.eps));
        2.0 * PI / (vperp(v, eta, bmod) * iota / r0 * (eps / 2.0).sqrt()).abs()
    };
    bounce_integral(v, eta, 2, &y0, taub, timestep_poloidal_motion).0
}

pub fn bounce_fast_ext(v: f64, eta: f64, taub: f64, ts: TsFn) -> ([f64; NVAR], i32) {
    let (bmod, _) = evaluate_bfield_local();
    let svh = set_sign_vpar_htheta();
    let th0 = dorb::get().th0;
    let mut y0 = [1.0e-15f64; NVAR];
    y0[0] = th0;
    y0[1] = svh * vpar(v, eta, bmod);
    for i in 2..6 {
        y0[i] = 0.0;
    }
    // integrate to taub without events (large dt, single chunk via CV_NORMAL)
    let (_ti, y) = integrate_fixed(v, eta, NVAR, &y0, taub, ts);
    let mut avg = [0.0f64; NVAR];
    for i in 0..NVAR {
        avg[i] = y[i] / taub;
    }
    (avg, 2)
}

pub fn bounce_fast(v: f64, eta: f64, taub: f64) -> [f64; NVAR] {
    bounce_fast_ext(v, eta, taub, timestep).0
}

fn integrate_fixed(v: f64, eta: f64, neq: usize, y0: &[f64], taub: f64, ts: TsFn) -> (f64, Vec<f64>) {
    unsafe {
        let mut ctx: SunContext = ptr::null_mut();
        SUNContext_Create(sun_comm_null(), &mut ctx);
        let y = N_VNew_Serial(neq as c_long, ctx);
        let yv = N_VGetArrayPointer(y);
        for i in 0..neq {
            *yv.add(i) = y0[i];
        }
        let oc = Box::new(OdeCtx { v, eta, neq, ts });
        let mem = CVodeCreate(CV_ADAMS, ctx);
        CVodeInit(mem, cv_rhs, 0.0, y);
        CVodeSStolerances(mem, 1.0e-9, 1.0e-10);
        CVodeSetMaxNumSteps(mem, 50000);
        CVodeSetUserData(mem, &*oc as *const OdeCtx as *mut c_void);
        let nls = SUNNonlinSol_FixedPoint(y, 0, ctx);
        CVodeSetNonlinearSolver(mem, nls);
        let mut tret = 0.0f64;
        CVode(mem, taub, y, &mut tret, CV_NORMAL);
        let out: Vec<f64> = (0..neq).map(|i| *yv.add(i)).collect();
        SUNNonlinSolFree(nls);
        N_VDestroy(y);
        let mut m = mem;
        CVodeFree(&mut m);
        SUNContext_Free(&mut ctx);
        drop(oc);
        (tret, out)
    }
}

pub fn bounce(v: f64, eta: f64, taub_estimate: f64, have: bool) -> (f64, [f64; NVAR]) {
    let (bmod, _) = evaluate_bfield_local();
    let svh = set_sign_vpar_htheta();
    let th0 = dorb::get().th0;
    let mut y0 = [1.0e-15f64; NVAR];
    y0[0] = th0;
    y0[1] = svh * vpar(v, eta, bmod);
    for i in 2..6 {
        y0[i] = 0.0;
    }
    let tb = if have {
        taub_estimate
    } else {
        let (iota, r0, eps) = with_field(|f| (f.iota, f.r0, f.eps));
        2.0 * PI / (vperp(v, eta, bmod) * iota / r0 * (eps / 2.0).sqrt()).abs()
    };
    let (taub, y) = bounce_integral(v, eta, NVAR, &y0, tb / 5.0, timestep);
    let mut avg = [0.0f64; NVAR];
    for i in 0..NVAR {
        avg[i] = y[i] / taub;
    }
    (taub, avg)
}
