//! Canonical-frequency splines + evaluation, literal port of neort_freq
//! (freq.f90), built on the verified orbit + spline layers. thread_local state.

use crate::driftorbit as dorb;
use crate::field::with_field;
use crate::orbit::{bounce_fast, bounce_time, set_s, NVAR};
use crate::profiles;
use crate::spline::{spline_coeff, spline_val_0};
use std::cell::RefCell;

pub const NETASPL: usize = 100;
const PI: f64 = std::f64::consts::PI;

#[derive(Default)]
struct FreqState {
    omth_spl: Vec<f64>,
    omtb_spl: Vec<f64>,
    omth_pass_spl: Vec<f64>,
    omtb_pass_spl: Vec<f64>,
    k_taub_p: f64,
    d_taub_p: f64,
    k_taub_t: f64,
    d_taub_t: f64,
    k_omtb_p: f64,
    d_omtb_p: f64,
    k_omtb_t: f64,
    d_omtb_t: f64,
    js: [i32; 4],
}

thread_local! {
    static FS: RefCell<FreqState> = RefCell::new(FreqState { js: [1; 4], ..Default::default() });
}

pub fn init_canon_freq_trapped_spline() {
    let n = NETASPL;
    let v = profiles::get().vth;
    let d = dorb::get();
    let etatp = d.etatp;
    let etamin = (1.0 + dorb::EPST) * etatp;
    let etamax = etatp + (d.etadt - etatp) * (1.0 - dorb::EPSST_SPL);
    dorb::update(|x| {
        x.etamin = etamin;
        x.etamax = etamax;
    });
    let magdrift = d.magdrift;

    let b = (dorb::EPST_SPL).ln();
    let aa = 1.0 / (n as f64 - 1.0) * ((etamax / etamin - 1.0).ln() - b);
    let mut etarange = vec![0.0f64; n];
    let mut omtb_v = vec![0.0f64; n];
    let mut omth_v = vec![0.0f64; n];
    let (mut taub0, mut taub1, mut leta0, mut leta1, mut otb0, mut otb1) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut est = 0.0f64;
    for k in (0..n).rev() {
        let eta = etamin * (1.0 + (aa * k as f64 + b).exp());
        etarange[k] = eta;
        est = if k == n - 1 { bounce_time(v, eta, 0.0, false) } else { bounce_time(v, eta, est, true) };
        let taub = est;
        let bavg = bounce_fast(v, eta, taub);
        if magdrift {
            omtb_v[k] = bavg[2];
        }
        omth_v[k] = 2.0 * PI / (v * taub);
        if k == 0 {
            leta0 = (eta - etatp).ln();
            taub0 = v * taub;
            if magdrift {
                otb0 = omtb_v[k] / omth_v[k];
            }
        }
        if k == 1 {
            leta1 = (eta - etatp).ln();
            taub1 = v * taub;
            if magdrift {
                otb1 = omtb_v[k] / omth_v[k];
            }
        }
    }
    let omth_spl = spline_coeff(&etarange, &omth_v);
    let omtb_spl = if magdrift { spline_coeff(&etarange, &omtb_v) } else { Vec::new() };
    FS.with(|c| {
        let mut f = c.borrow_mut();
        f.k_taub_t = (taub1 - taub0) / (leta1 - leta0);
        f.d_taub_t = taub0 - f.k_taub_t * leta0;
        f.omth_spl = omth_spl;
        if magdrift {
            f.k_omtb_t = (otb1 - otb0) / (leta1 - leta0);
            f.d_omtb_t = otb0 - f.k_omtb_t * leta0;
            f.omtb_spl = omtb_spl;
        }
    });
}

pub fn init_canon_freq_passing_spline() {
    let n = NETASPL;
    let v = profiles::get().vth;
    let d = dorb::get();
    let etatp = d.etatp;
    let etamin = etatp * dorb::EPSSP_SPL;
    let etamax = etatp;
    dorb::update(|x| {
        x.etamin = etamin;
        x.etamax = etamax;
    });
    let magdrift = d.magdrift;

    let b = ((etamax - etamin) / etamax).ln();
    let aa = 1.0 / (n as f64 - 1.0) * ((dorb::EPSP_SPL).ln() - b);
    let mut etarange = vec![0.0f64; n];
    let mut omtb_v = vec![0.0f64; n];
    let mut omth_v = vec![0.0f64; n];
    let (mut taub0, mut taub1, mut leta0, mut leta1, mut otb0, mut otb1) = (0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
    let mut est = 0.0f64;
    for k in (0..n).rev() {
        let eta = etamax * (1.0 - (aa * k as f64 + b).exp());
        etarange[k] = eta;
        est = if k == n - 1 { bounce_time(v, eta, 0.0, false) } else { bounce_time(v, eta, est, true) };
        let taub = est;
        let bavg = bounce_fast(v, eta, taub);
        if magdrift {
            omtb_v[k] = bavg[2];
        }
        omth_v[k] = 2.0 * PI / (v * taub);
        if k == n - 2 {
            leta0 = (etatp - eta).ln();
            taub0 = v * taub;
            if magdrift {
                otb0 = omtb_v[k] / omth_v[k];
            }
        }
        if k == n - 1 {
            leta1 = (etatp - eta).ln();
            taub1 = v * taub;
            if magdrift {
                otb1 = omtb_v[k] / omth_v[k];
            }
        }
    }
    let omth_pass_spl = spline_coeff(&etarange, &omth_v);
    let omtb_pass_spl = if magdrift { spline_coeff(&etarange, &omtb_v) } else { Vec::new() };
    FS.with(|c| {
        let mut f = c.borrow_mut();
        f.k_taub_p = (taub1 - taub0) / (leta1 - leta0);
        f.d_taub_p = taub0 - f.k_taub_p * leta0;
        f.omth_pass_spl = omth_pass_spl;
        if magdrift {
            f.k_omtb_p = (otb1 - otb0) / (leta1 - leta0);
            f.d_omtb_p = otb0 - f.k_omtb_p * leta0;
            f.omtb_pass_spl = omtb_pass_spl;
        }
    });
}

pub fn om_th(v: f64, eta: f64) -> (f64, f64, f64) {
    let d = dorb::get();
    let etatp = d.etatp;
    let sv = FS.with(|c| {
        let mut f = c.borrow_mut();
        if eta > etatp {
            if eta > etatp * (1.0 + dorb::EPST_SPL) {
                let mut js = f.js[0];
                let r = spline_val_0(&f.omth_spl, NETASPL - 1, eta, &mut js);
                f.js[0] = js;
                r
            } else {
                let v0 = 2.0 * PI / (f.k_taub_t * (eta - etatp).ln() + f.d_taub_t);
                [v0, -v0 * v0 / (2.0 * PI) * f.k_taub_t / (eta - etatp), 0.0]
            }
        } else if eta < etatp * (1.0 - dorb::EPSP_SPL) {
            let mut js = f.js[2];
            let r = spline_val_0(&f.omth_pass_spl, NETASPL - 1, eta, &mut js);
            f.js[2] = js;
            r
        } else {
            let v0 = 2.0 * PI / (f.k_taub_p * (etatp - eta).ln() + f.d_taub_p);
            [v0, -v0 * v0 / (2.0 * PI) * f.k_taub_p / (eta - etatp), 0.0]
        }
    });
    (d.sign_vpar * sv[0] * v, d.sign_vpar * sv[0], d.sign_vpar * sv[1] * v)
}

pub fn om_tb(v: f64, eta: f64) -> (f64, f64, f64) {
    let d = dorb::get();
    let etatp = d.etatp;
    let sv: [f64; 2] = if eta > etatp {
        if eta > etatp * (1.0 + dorb::EPST_SPL) {
            FS.with(|c| {
                let mut f = c.borrow_mut();
                let mut js = f.js[1];
                let r = spline_val_0(&f.omtb_spl, NETASPL - 1, eta, &mut js);
                f.js[1] = js;
                [r[0], r[1]]
            })
        } else {
            let (omth, _dv, domthdeta) = om_th(v, eta);
            let (k, dd) = FS.with(|c| { let f = c.borrow(); (f.k_omtb_t, f.d_omtb_t) });
            [
                d.sign_vpar * (k * (eta - etatp).ln() + dd) * omth / v,
                d.sign_vpar * (omth / v * k / (eta - etatp) + domthdeta / v * (k * (eta - etatp).ln() + dd)),
            ]
        }
    } else if eta < etatp * (1.0 - dorb::EPSP_SPL) {
        FS.with(|c| {
            let mut f = c.borrow_mut();
            let mut js = f.js[3];
            let r = spline_val_0(&f.omtb_pass_spl, NETASPL - 1, eta, &mut js);
            f.js[3] = js;
            [r[0], r[1]]
        })
    } else {
        let (omth, _dv, domthdeta) = om_th(v, eta);
        let (k, dd) = FS.with(|c| { let f = c.borrow(); (f.k_omtb_p, f.d_omtb_p) });
        [
            d.sign_vpar * (k * (etatp - eta).ln() + dd) * omth / v,
            d.sign_vpar * (omth / v * k / (eta - etatp) + domthdeta / v * (k * (etatp - eta).ln() + dd)),
        ]
    };
    (sv[0] * v * v, 2.0 * sv[0] * v, sv[1] * v * v)
}

pub fn om_ph(v: f64, eta: f64) -> (f64, f64, f64) {
    let d = dorb::get();
    let om_te = profiles::get().om_te;
    let iota = with_field(|f| f.iota);
    if eta > d.etatp {
        let (mut omph, mut dv, mut de) = (om_te, 0.0, 0.0);
        if d.magdrift {
            let (otb, dvtb, detb) = om_tb(v, eta);
            omph += otb;
            dv += dvtb;
            de += detb;
        }
        (omph, dv, de)
    } else {
        let (omth, dvth, deth) = om_th(v, eta);
        let (mut omph, mut dv, mut de) = (om_te + omth / iota, dvth / iota, deth / iota);
        if d.magdrift {
            let (otb, dvtb, detb) = om_tb(v, eta);
            omph += otb;
            dv += dvtb;
            de += detb;
        }
        (omph, dv, de)
    }
}

pub fn d_om_ds(v: f64, eta: f64, taub_estimate: f64) -> (f64, f64) {
    let ds = 2.0e-8;
    let s0 = crate::orbit::get_s();
    let d = dorb::get();
    let iota = with_field(|f| f.iota);
    let dom_teds = profiles::get().dom_teds;

    set_s(s0 - ds / 2.0);
    let taub = bounce_time(v, eta, taub_estimate, true);
    let bavg: [f64; NVAR] = bounce_fast(v, eta, taub);
    let omth = d.sign_vpar_htheta * 2.0 * PI / taub;
    let omph_noe = if d.magdrift {
        if eta > d.etatp { bavg[2] * v * v } else { bavg[2] * v * v + omth / iota }
    } else if eta > d.etatp {
        0.0
    } else {
        omth / iota
    };

    set_s(s0 + ds / 2.0);
    let taub2 = bounce_time(v, eta, taub_estimate, true);
    let bavg2 = bounce_fast(v, eta, taub2);
    let domthds = d.sign_vpar_htheta * (2.0 * PI / taub2 - d.sign_vpar_htheta * omth) / ds;
    let domphds = if d.magdrift {
        if eta > d.etatp {
            dom_teds + (bavg2[2] * v * v - omph_noe) / ds
        } else {
            dom_teds + (bavg2[2] * v * v + (2.0 * PI / taub2) / iota - omph_noe) / ds
        }
    } else if eta > d.etatp {
        dom_teds
    } else {
        dom_teds + ((2.0 * PI / taub2) / iota - omph_noe) / ds
    };
    set_s(s0);
    (domthds, domphds)
}
