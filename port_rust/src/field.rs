//! Axisymmetric + perturbation magnetic field, literal port of do_magfie_mod /
//! do_magfie_pert_mod (do_magfie_standalone.f90), mirroring the C port. Module
//! state lives in a thread_local (single-threaded gate), with accessors for the
//! orbit/transport layers that read the shared field quantities.

use crate::spline::{spline_coeff, spline_val_0};
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub const SIGN_THETA: f64 = -1.0;
const ITOB: f64 = 2.0e-1 * SIGN_THETA;
const PI: f64 = std::f64::consts::PI;

#[derive(Default)]
pub struct Field {
    pub bfac: f64,
    pub inp_swi: i32,
    pub psi_pr: f64,
    pub r0: f64,
    pub a: f64,
    pub b00: f64,
    pub bthcov: f64,
    pub bphcov: f64,
    pub dbthcovds: f64,
    pub dbphcovds: f64,
    pub q: f64,
    pub dqds: f64,
    pub iota: f64,
    pub eps: f64,
    pub b0h: f64,
    pub nflux: usize,
    pub nmode: usize,
    ncol1: usize,
    ncol2: usize,
    params0: Vec<f64>, // nflux x (ncol1+1)
    modes0: Vec<f64>,  // nflux x nmode x (ncol2+2)
    spl1: Vec<f64>,    // ncol1 splines of nseg*5
    spl2: Vec<f64>,    // ncol2*nmode splines of nseg*5
    nseg: usize,
    jstart: i32,
    // perturbation
    pub pert_mph: i32,
    p_ncol2: usize,
    p_nflux: usize,
    p_nmode: usize,
    p_nfp: i32,
    p_nseg: usize,
    p_params: Vec<f64>,
    p_modes: Vec<f64>,
    p_spl2: Vec<f64>,
}

thread_local! {
    static F: RefCell<Field> = RefCell::new(Field { bfac: 1.0, jstart: 1, ..Default::default() });
}

pub fn with_field<R>(f: impl FnOnce(&Field) -> R) -> R {
    F.with(|c| f(&c.borrow()))
}
pub fn set_bfac(v: f64) {
    F.with(|c| c.borrow_mut().bfac = v);
}
pub fn set_inp_swi(v: i32) {
    F.with(|c| c.borrow_mut().inp_swi = v);
}
pub fn set_eps(v: f64) {
    F.with(|c| c.borrow_mut().eps = v);
}

/// Read a Boozer file, mirroring boozer_read (header skip 5, per-surface skip 2
/// before params and 1 before modes). Returns (nflux, nmode, nfp, flux, a, R0,
/// params, modes) with params nflux x (ncol1+1), modes nflux x nmode x (nc2+2).
fn boozer_read(path: &str, ncol1: usize, nc2: usize) -> (usize, usize, i32, f64, f64, f64, Vec<f64>, Vec<f64>) {
    let file = File::open(path).unwrap_or_else(|_| panic!("cannot open {}", path));
    let mut lines = BufReader::new(file).lines().map(|l| l.unwrap());
    for _ in 0..5 {
        lines.next();
    }
    let hdr: Vec<f64> = lines
        .next()
        .unwrap()
        .split_whitespace()
        .map(|t| t.parse().unwrap())
        .collect();
    let (m0b, n0b, nflux, nfp) = (hdr[0] as i32, hdr[1] as i32, hdr[2] as usize, hdr[3] as i32);
    let (flux, a, r0) = (hdr[4], hdr[5], hdr[6]);
    let nmode = ((m0b + 1) * (n0b + 1)) as usize;
    let pcols = ncol1 + 1;
    let mcols = nc2 + 2;
    let mut params = vec![0.0f64; nflux * pcols];
    let mut modes = vec![0.0f64; nflux * nmode * mcols];
    for ks in 0..nflux {
        lines.next();
        lines.next(); // two label lines
        let pv: Vec<f64> = lines
            .next()
            .unwrap()
            .split_whitespace()
            .map(|t| t.parse().unwrap())
            .collect();
        for k in 0..pcols {
            params[ks * pcols + k] = pv[k];
        }
        lines.next(); // mode-label line
        for j in 0..nmode {
            let mv: Vec<f64> = lines
                .next()
                .unwrap()
                .split_whitespace()
                .map(|t| t.parse().unwrap())
                .collect();
            for k in 0..mcols {
                modes[(ks * nmode + j) * mcols + k] = mv[k];
            }
        }
    }
    (nflux, nmode, nfp, flux, a, r0, params, modes)
}

fn build_splines(params: &[f64], modes: &[f64], nflux: usize, nmode: usize, nc2: usize, pcols: usize, mcols: usize) -> (Vec<f64>, Vec<f64>) {
    let ncol1 = pcols - 1;
    let seg = nflux - 1;
    let x: Vec<f64> = (0..nflux).map(|i| params[i * pcols]).collect();
    let mut s1 = vec![0.0f64; ncol1 * seg * 5];
    for k in 0..ncol1 {
        let y: Vec<f64> = (0..nflux).map(|i| params[i * pcols + k + 1]).collect();
        let c = spline_coeff(&x, &y);
        s1[k * seg * 5..(k + 1) * seg * 5].copy_from_slice(&c);
    }
    let mut s2 = vec![0.0f64; nc2 * nmode * seg * 5];
    for j in 0..nmode {
        for k in 0..nc2 {
            let y: Vec<f64> = (0..nflux).map(|i| modes[(i * nmode + j) * mcols + k + 2]).collect();
            let c = spline_coeff(&x, &y);
            let off = (k * nmode + j) * seg * 5;
            s2[off..off + seg * 5].copy_from_slice(&c);
        }
    }
    (s1, s2)
}

pub fn do_magfie_init(path: &str) {
    F.with(|c| {
        let mut f = c.borrow_mut();
        f.ncol1 = 5;
        f.ncol2 = if f.inp_swi == 8 { 4 } else { 8 };
        let (nflux, nmode, _nfp, flux, a, r0, params, modes) = boozer_read(path, f.ncol1, f.ncol2);
        f.nflux = nflux;
        f.nmode = nmode;
        f.a = 100.0 * a;
        f.r0 = 100.0 * r0;
        f.psi_pr = 1.0e8 * flux / (2.0 * PI) * f.bfac;
        f.nseg = nflux - 1;
        let (s1, s2) = build_splines(&params, &modes, nflux, nmode, f.ncol2, f.ncol1 + 1, f.ncol2 + 2);
        f.params0 = params;
        f.modes0 = modes;
        f.spl1 = s1;
        f.spl2 = s2;
        let mcols = f.ncol2 + 2;
        f.r0 = f.modes0[0 * nmode * mcols + 0 * mcols + 2] * 100.0;
        f.b00 = 1.0e4 * f.modes0[0 * nmode * mcols + 0 * mcols + 5] * f.bfac;
    });
}

/// do_magfie at x=(s,ph,th): returns (bmod, sqrtg, bder[3], hcovar[3], hctrvr[3]).
/// Updates the shared bthcov/bphcov/iota/q/dqds etc.
pub fn do_magfie(x: &[f64; 3]) -> (f64, f64, [f64; 3], [f64; 3], [f64; 3]) {
    F.with(|c| {
        let mut f = c.borrow_mut();
        let pcols = f.ncol1 + 1;
        let mcols = f.ncol2 + 2;
        let nseg = f.nseg;
        let nm = f.nmode;
        let bf = f.bfac;
        let mut x1 = if f.params0[0] > x[0] { f.params0[0] } else { x[0] };
        let last = f.params0[(f.nflux - 1) * pcols];
        if last < x1 {
            x1 = last;
        }
        let mut js = f.jstart;
        let sv = spline_val_0(&f.spl1[2 * nseg * 5..3 * nseg * 5], nseg, x1, &mut js);
        f.bthcov = ITOB * sv[0] * bf;
        f.dbthcovds = ITOB * sv[1] * bf;
        let sv = spline_val_0(&f.spl1[1 * nseg * 5..2 * nseg * 5], nseg, x1, &mut js);
        f.bphcov = ITOB * sv[0] * bf;
        f.dbphcovds = ITOB * sv[1] * bf;
        let sv = spline_val_0(&f.spl1[0..nseg * 5], nseg, x1, &mut js);
        f.iota = sv[0];
        f.q = 1.0 / f.iota;
        f.dqds = -sv[1] / (f.iota * f.iota);

        let m0 = f.modes0[(0 * nm + 0) * mcols + 0];
        let m1 = f.modes0[(0 * nm + 1) * mcols + 0];
        let dm = m1 - m0;
        let mut cost = vec![0.0f64; nm];
        let mut sint = vec![0.0f64; nm];
        // fast_sin_cos
        let (mut ffr, mut ffi) = ((m0 * x[2]).cos(), (m0 * x[2]).sin());
        let (rotr, roti) = ((dm * x[2]).cos(), (dm * x[2]).sin());
        for j in 0..nm {
            cost[j] = ffr;
            sint[j] = ffi;
            let nr = ffr * rotr - ffi * roti;
            let ni = ffr * roti + ffi * rotr;
            ffr = nr;
            ffi = ni;
        }

        let (mut bm, mut db, mut bth) = (0.0, 0.0, 0.0);
        if f.inp_swi == 8 {
            for j in 0..nm {
                let s2 = &f.spl2[(3 * nm + j) * nseg * 5..(3 * nm + j + 1) * nseg * 5];
                let sv = spline_val_0(s2, nseg, x1, &mut js);
                let bmnc = 1.0e4 * sv[0] * bf;
                let dbmnc = 1.0e4 * sv[1] * bf;
                if j == 0 {
                    f.b0h = bmnc;
                }
                bm += bmnc * cost[j];
                db += dbmnc * cost[j];
                bth += -f.modes0[(0 * nm + j) * mcols] * bmnc * sint[j];
            }
        } else {
            for j in 0..nm {
                let s2c = &f.spl2[(6 * nm + j) * nseg * 5..(6 * nm + j + 1) * nseg * 5];
                let sv = spline_val_0(s2c, nseg, x1, &mut js);
                let bmnc = 1.0e4 * sv[0] * bf;
                let dbmnc = 1.0e4 * sv[1] * bf;
                let s2s = &f.spl2[(7 * nm + j) * nseg * 5..(7 * nm + j + 1) * nseg * 5];
                let sv = spline_val_0(s2s, nseg, x1, &mut js);
                let bmns = 1.0e4 * sv[0] * bf;
                let dbmns = 1.0e4 * sv[1] * bf;
                if j == 0 {
                    f.b0h = bmnc;
                }
                let mj = f.modes0[(0 * nm + j) * mcols];
                bm += bmnc * cost[j] + bmns * sint[j];
                db += dbmnc * cost[j] + dbmns * sint[j];
                bth += -mj * bmnc * sint[j] + mj * bmns * cost[j];
            }
        }
        f.jstart = js;
        let bder = [db / bm, 0.0, bth / bm];
        let sqgbmod2 = SIGN_THETA * f.psi_pr * (f.bphcov + f.iota * f.bthcov);
        let sqgbmod = sqgbmod2 / bm;
        let sqrtg = sqgbmod / bm;
        let hcovar = [0.0, f.bphcov / bm, f.bthcov / bm];
        let hctrvr = [
            0.0,
            SIGN_THETA * f.psi_pr / sqgbmod,
            SIGN_THETA * f.iota * f.psi_pr / sqgbmod,
        ];
        (bm, sqrtg, bder, hcovar, hctrvr)
    })
}

pub fn do_magfie_pert_init(path: &str) {
    F.with(|c| {
        let mut f = c.borrow_mut();
        let ncol1 = 5usize;
        let nc2 = if f.inp_swi == 8 { 4 } else { 8 };
        let (nflux, nmode, nfp, _flux, _a, _r0, params, modes) = boozer_read(path, ncol1, nc2);
        f.p_ncol2 = nc2;
        f.p_nflux = nflux;
        f.p_nmode = nmode;
        f.p_nfp = nfp;
        f.p_nseg = nflux - 1;
        let mcols = nc2 + 2;
        let nv = modes[(0 * nmode + 0) * mcols + 1];
        f.pert_mph = (nfp as f64 * nv).round() as i32;
        let (_s1, s2) = build_splines(&params, &modes, nflux, nmode, nc2, ncol1 + 1, mcols);
        f.p_params = params;
        f.p_modes = modes;
        f.p_spl2 = s2;
    });
}

/// Perturbation amplitude bamp at x=(s,ph,th) as (re, im).
pub fn do_magfie_pert_amp(x: &[f64; 3]) -> (f64, f64) {
    F.with(|c| {
        let mut f = c.borrow_mut();
        let pcols = 6usize;
        let mcols = f.p_ncol2 + 2;
        let nseg = f.p_nseg;
        let nm = f.p_nmode;
        let bf = f.bfac;
        let mut x1 = if f.p_params[0] > x[0] { f.p_params[0] } else { x[0] };
        let last = f.p_params[(f.p_nflux - 1) * pcols];
        if last < x1 {
            x1 = last;
        }
        let mut js = f.jstart;
        let (mut sr, mut si) = (0.0, 0.0);
        if f.inp_swi == 8 {
            for j in 0..nm {
                let s2 = &f.p_spl2[(3 * nm + j) * nseg * 5..(3 * nm + j + 1) * nseg * 5];
                let sv = spline_val_0(s2, nseg, x1, &mut js);
                let bmnc = 1.0e4 * sv[0] * bf;
                let m = f.p_modes[(0 * nm + j) * mcols];
                sr += bmnc * (m * x[2]).cos();
            }
        } else {
            let m0 = f.p_modes[(0 * nm + 0) * mcols];
            let m1 = f.p_modes[(0 * nm + 1) * mcols];
            let dm = m1 - m0;
            let (mut ffr, mut ffi) = ((m0 * x[2]).cos(), (m0 * x[2]).sin());
            let (rotr, roti) = ((dm * x[2]).cos(), (dm * x[2]).sin());
            for j in 0..nm {
                let s2c = &f.p_spl2[(6 * nm + j) * nseg * 5..(6 * nm + j + 1) * nseg * 5];
                let sv = spline_val_0(s2c, nseg, x1, &mut js);
                let bmnc = 1.0e4 * sv[0] * bf;
                let s2s = &f.p_spl2[(7 * nm + j) * nseg * 5..(7 * nm + j + 1) * nseg * 5];
                let sv = spline_val_0(s2s, nseg, x1, &mut js);
                let bmns = 1.0e4 * sv[0] * bf;
                // (bmnc - i bmns) * (ffr + i ffi)
                sr += bmnc * ffr + bmns * ffi;
                si += bmnc * ffi - bmns * ffr;
                let nr = ffr * rotr - ffi * roti;
                let ni = ffr * roti + ffi * rotr;
                ffr = nr;
                ffi = ni;
            }
        }
        f.jstart = js;
        (sr, si)
    })
}
