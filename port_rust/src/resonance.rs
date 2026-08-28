//! Resonance-line root finding, literal port of neort_resonance (resonance.f90).

use crate::driftorbit as dorb;
use crate::freq::{om_ph, om_th};

/// Fill roots (flat ninterv x 3, stride 3) with sign-change intervals of
/// res(eta)=mth*Om_th+mph*Om_ph; returns nroots.
pub fn driftorbit_coarse(v: f64, eta_min: f64, eta_max: f64, roots: &mut [f64], ninterv: usize) -> usize {
    let deta = (eta_max - eta_min) / ninterv as f64;
    let d = dorb::get();
    let mut resold = 0.0;
    let mut nroots = 0usize;
    for k in 0..=ninterv {
        let eta = eta_min + k as f64 * deta;
        let (omth, ..) = om_th(v, eta);
        let (omph, ..) = om_ph(v, eta);
        let res = d.mth as f64 * omth + d.mph as f64 * omph;
        if k > 0 {
            let snew = if res < 0.0 { -1.0 } else { 1.0 };
            let sold = if resold < 0.0 { -1.0 } else { 1.0 };
            if snew != sold {
                roots[nroots * 3] = eta - deta;
                roots[nroots * 3 + 1] = eta;
                nroots += 1;
            }
        }
        resold = res;
    }
    nroots
}

pub fn driftorbit_nroot(v: f64, eta_min: f64, eta_max: f64) -> usize {
    let mut roots = vec![0.0f64; dorb::NLEV * 3];
    driftorbit_coarse(v, eta_min, eta_max, &mut roots, dorb::NLEV)
}

/// Bisection root in [eta_min,eta_max]; returns (eta, dres/deta).
pub fn driftorbit_root(v: f64, tol: f64, eta_min: f64, eta_max: f64) -> (f64, f64) {
    let d = dorb::get();
    let (mth, mph) = (d.mth as f64, d.mph as f64);
    let maxit = 100;
    let mut etamin2 = eta_min;
    let mut etamax2 = eta_max;
    let mut eta = etamin2;

    let (omth, _, _) = om_th(v, eta);
    let (omph, _, _) = om_ph(v, eta);
    let resmin = mph * omph + mth * omth;
    eta = etamax2;
    let (omth, _, _) = om_th(v, eta);
    let (omph, _, _) = om_ph(v, eta);
    let resmax = mph * omph + mth * omth;
    let slope_pos = resmax - resmin > 0.0;

    if driftorbit_nroot(v, etamin2, etamax2) == 0 {
        let (_, _, deth) = om_th(v, eta);
        let (_, _, deph) = om_ph(v, eta);
        return (eta, mph * deph + mth * deth);
    }

    let mut out_eta = eta;
    let mut out_d = 0.0;
    let mut state = -2;
    for _ in 1..=maxit {
        let (omth, _, deth) = om_th(v, eta);
        let (omph, _, deph) = om_ph(v, eta);
        let res = mph * omph + mth * omth;
        out_eta = eta;
        if res.abs() < tol {
            state = 1;
            out_d = mph * deph + mth * deth;
            break;
        } else if (slope_pos && res > 0.0) || (!slope_pos && res < 0.0) {
            etamax2 = eta;
            eta = (eta + etamin2) / 2.0;
        } else {
            etamin2 = eta;
            eta = (eta + etamax2) / 2.0;
        }
    }
    if state < 0 {
        let (_, _, deth) = om_th(v, out_eta);
        let (_, _, deph) = om_ph(v, out_eta);
        out_d = mph * deph + mth * deth;
    }
    (out_eta, out_d)
}
