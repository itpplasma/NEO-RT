//! Velocity-space transport integration, literal port of neort_transport
//! (transport.f90). Produces D11/D12 and torque density for one harmonic.

use crate::driftorbit as dorb;
use crate::field::{do_magfie, do_magfie_pert_amp, with_field, SIGN_THETA};
use crate::orbit::{bounce_fast_ext, get_s, poloidal_velocity, NVAR};
use crate::profiles;
use crate::freq::om_th;
use crate::resonance::{driftorbit_coarse, driftorbit_root};
use std::cell::Cell;

const PI: f64 = std::f64::consts::PI;

thread_local! {
    static OMTH: Cell<f64> = Cell::new(0.0);
}

fn d11int(ux: f64, taub: f64, hmn2: f64) -> f64 {
    let d = dorb::get();
    let q = with_field(|f| f.q);
    let psi_pr = with_field(|f| f.psi_pr);
    let p = profiles::get();
    PI.powf(1.5) * (d.mph * d.mph) as f64 * profiles::C * profiles::C * q * p.vth
        / (p.qi * p.qi * d.dvds * psi_pr.abs())
        * ux * ux * ux * (-ux * ux).exp() * taub * hmn2
}
fn d12int(ux: f64, taub: f64, hmn2: f64) -> f64 {
    d11int(ux, taub, hmn2) * ux * ux
}
fn tphi_int(ux: f64, taub: f64, hmn2: f64) -> f64 {
    let d = dorb::get();
    let (q, psi_pr) = with_field(|f| (f.q, f.psi_pr));
    let p = profiles::get();
    let sgn = if psi_pr * q * SIGN_THETA >= 0.0 { 1.0 } else { -1.0 };
    sgn * PI.powf(1.5) * (d.mph * d.mph) as f64 * p.ni1 * profiles::C * p.vth / p.qi
        * ux * ux * ux * (-ux * ux).exp() * taub * hmn2 * (p.a1 + p.a2 * ux * ux)
}

pub fn timestep_transport(v: f64, eta: f64, neq: usize, t: f64, y: &[f64], ydot: &mut [f64]) {
    let _ = neq;
    let s = get_s();
    let (bmod, _sqrtg, bder, _hcov, hctrvr) = do_magfie(&[s, 0.0, y[0]]);
    poloidal_velocity(v, eta, bmod, hctrvr[2], bder[2], y[1], ydot);
    let d = dorb::get();
    let q = with_field(|f| f.q);
    let omth = OMTH.with(|c| c.get());

    // epsn (complex)
    let (er, ei) = if d.pertfile {
        let (ar, ai) = do_magfie_pert_amp(&[s, 0.0, y[0]]);
        (d.epsmn * ar / bmod, d.epsmn * ai / bmod)
    } else {
        (d.epsmn * (d.m0 as f64 * y[0]).cos(), d.epsmn * (d.m0 as f64 * y[0]).sin())
    };
    // phase
    let ph = if eta > d.etatp {
        q * d.mph as f64 * y[0] - d.mth as f64 * t * omth
    } else {
        q * d.mph as f64 * y[0] - (d.mth as f64 + q * d.mph as f64) * t * omth
    };
    let (cr, ci) = (ph.cos(), ph.sin());
    let amp = 2.0 - eta * bmod;
    // Hn = amp * epsn * exp(i ph)
    let hr = amp * (er * cr - ei * ci);
    let hi = amp * (er * ci + ei * cr);
    ydot[2] = hr;
    ydot[3] = hi;
    if d.nonlin {
        ydot[4] = 1.0 / bmod;
        ydot[5] = bmod;
    } else {
        ydot[4] = 0.0;
        ydot[5] = 0.0;
    }
    ydot[6] = 0.0;
}

/// Midpoint integral over v: returns (D[2], T).
pub fn compute_transport_integral(vmin: f64, vmax: f64, vsteps: usize) -> ([f64; 2], f64) {
    let p = profiles::get();
    let vth = p.vth;
    let d0 = dorb::get();
    let mut roots = vec![0.0f64; dorb::NLEV * 3];
    let mut d = [0.0f64; 2];
    let mut t = 0.0f64;
    let du = (vmax - vmin) / (vsteps as f64 * vth);
    let mut ux = vmin / vth + du / 2.0;
    let om_te = p.om_te;
    let mi = p.mi;

    for _ku in 1..=vsteps {
        let v = ux * vth;
        let nroots = driftorbit_coarse(v, d0.etamin, d0.etamax, &mut roots, dorb::NLEV);
        for kr in 0..nroots {
            let (eta, dres) = driftorbit_root(v, 1.0e-8 * om_te.abs(), roots[kr * 3], roots[kr * 3 + 1]);
            let (omth, _dv, _de) = om_th(v, eta);
            OMTH.with(|c| c.set(omth));
            let taub = 2.0 * PI / omth.abs();
            let (bavg, _ist): ([f64; NVAR], i32) = bounce_fast_ext(v, eta, taub, timestep_transport);
            let hmn2 = (bavg[2] * bavg[2] + bavg[3] * bavg[3]) * (mi * (ux * vth) * (ux * vth) / 2.0).powi(2);
            let atten = if d0.nonlin {
                panic!("nonlinear attenuation path not ported (outside golden gate)");
            } else {
                1.0
            };
            let dd11 = du * d11int(ux, taub, hmn2) / dres.abs();
            let dd12 = du * d12int(ux, taub, hmn2) / dres.abs();
            d[0] += dd11 * atten;
            d[1] += dd12 * atten;
            if d0.comptorque {
                let dt = du * tphi_int(ux, taub, hmn2) / dres.abs();
                t += dt * atten;
            }
        }
        ux += du;
    }

    let (r0, iota) = with_field(|f| (f.r0, f.iota));
    let a = with_field(|f| f.a);
    let d_plateau = PI * vth.powi(3) / (16.0 * r0 * iota * (p.qi * d0.b0 / (mi * profiles::C)).powi(2));
    let dsdreff = 2.0 / a * get_s().sqrt();
    d[0] = dsdreff.powi(-2) * d[0] / d_plateau;
    d[1] = dsdreff.powi(-2) * d[1] / d_plateau;
    (d, t)
}
