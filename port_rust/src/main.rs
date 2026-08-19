//! NEO-RT Rust port driver: parse <case>.in, run the golden-path computation
//! (mirrors neort_compute_at_s + compute_transport), write the 5 output files.
//! Single-threaded (deterministic gate). Uses the same CVODE/LAPACK as C.

use neo_rt_rust::driftorbit as dorb;
use neo_rt_rust::field::{do_magfie, do_magfie_init, do_magfie_pert_amp, do_magfie_pert_init, set_bfac, set_inp_swi, with_field};
use neo_rt_rust::magfie::init_flux_surface_average;
use neo_rt_rust::orbit::set_s;
use neo_rt_rust::profiles::{self, init_thermodynamic_forces, read_and_init_plasma_input, read_and_init_profile_input};
use neo_rt_rust::freq::{init_canon_freq_passing_spline, init_canon_freq_trapped_spline};
use neo_rt_rust::transport::compute_transport_integral;
use std::fs::File;
use std::io::{BufRead, BufReader, Write};

const PI: f64 = std::f64::consts::PI;

#[derive(Default)]
struct Config {
    s: f64, epsmn: f64, bfac: f64, efac: f64,
    m0: i32, mph: i32, inp_swi: i32, vsteps: i32,
    comptorque: bool, magdrift: bool, nopassing: bool, noshear: bool, pertfile: bool, nonlin: bool,
}

fn read_config(path: &str) -> Config {
    let mut c = Config { bfac: 1.0, efac: 1.0, comptorque: true, magdrift: true, ..Default::default() };
    let f = File::open(path).unwrap_or_else(|_| panic!("cannot open {}", path));
    for line in BufReader::new(f).lines().map(|l| l.unwrap()) {
        let line = line.split('!').next().unwrap();
        if let Some((k, v)) = line.split_once('=') {
            let key = k.trim().to_lowercase();
            let val = v.trim();
            let b = || val.starts_with('t') || val.starts_with('T') || val.contains(".true.") || val.contains(".TRUE.");
            match key.as_str() {
                "s" => c.s = val.parse().unwrap_or(0.0),
                "epsmn" => c.epsmn = val.parse().unwrap_or(0.0),
                "bfac" => c.bfac = val.parse().unwrap_or(1.0),
                "efac" => c.efac = val.parse().unwrap_or(1.0),
                "m0" => c.m0 = val.parse().unwrap_or(0),
                "mph" => c.mph = val.parse().unwrap_or(0),
                "inp_swi" => c.inp_swi = val.parse().unwrap_or(0),
                "vsteps" => c.vsteps = val.parse().unwrap_or(0),
                "comptorque" => c.comptorque = b(),
                "magdrift" => c.magdrift = b(),
                "nopassing" => c.nopassing = b(),
                "noshear" => c.noshear = b(),
                "pertfile" => c.pertfile = b(),
                "nonlin" => c.nonlin = b(),
                _ => {}
            }
        }
    }
    c
}

#[derive(Clone, Copy, Default)]
struct Harmonic {
    mth: i32,
    dresco: [f64; 2], dresctr: [f64; 2], drest: [f64; 2],
    tresco: f64, tresctr: f64, trest: f64,
    vminp_vth: f64, vmaxp_vth: f64, vmint_vth: f64, vmaxt_vth: f64,
}

fn set_trapped() {
    let d = dorb::get();
    dorb::update(|x| {
        x.etamin = (1.0 + dorb::EPST) * d.etatp;
        x.etamax = (1.0 - dorb::EPST) * d.etadt;
    });
}
fn set_passing() {
    let d = dorb::get();
    dorb::update(|x| {
        x.etamin = dorb::EPSP * d.etatp;
        x.etamax = (1.0 - dorb::EPSP) * d.etatp;
    });
}

fn main() {
    let run = std::env::args().nth(1).expect("usage: neo_rt_rust <runname>");
    let cfg = read_config(&format!("{}.in", run));
    set_inp_swi(cfg.inp_swi);
    set_bfac(cfg.bfac);
    dorb::update(|d| {
        d.epsmn = cfg.epsmn;
        d.m0 = cfg.m0;
        d.efac = cfg.efac;
        d.comptorque = cfg.comptorque;
        d.magdrift = cfg.magdrift;
        d.nopassing = cfg.nopassing;
        d.pertfile = cfg.pertfile;
        d.nonlin = cfg.nonlin;
        d.noshear = cfg.noshear;
    });
    let s = cfg.s;

    do_magfie_init("in_file");
    if cfg.pertfile {
        do_magfie_pert_init("in_file_pert");
        let mph = with_field(|f| f.pert_mph);
        dorb::update(|d| d.mph = mph);
    } else {
        dorb::update(|d| d.mph = cfg.mph);
    }
    read_and_init_plasma_input("plasma.in", s);
    let r0 = with_field(|f| f.r0);
    read_and_init_profile_input("profile.in", s, r0, cfg.efac, cfg.bfac);
    init_flux_surface_average(s);
    set_s(s);
    init_canon_freq_trapped_spline();
    if !cfg.nopassing {
        init_canon_freq_passing_spline();
    }
    dorb::update(|d| d.sign_vpar = 1.0);
    set_trapped();
    let (psi_pr, q) = with_field(|f| (f.psi_pr, f.q));
    if cfg.comptorque {
        init_thermodynamic_forces(psi_pr, q);
    }

    let mthmax = (2.0 * (dorb::get().mph as f64 * q).abs()).ceil() as i32;
    let mthmin = -mthmax;
    let mph = dorb::get().mph;
    let _ = mph;
    let mut harmonics: Vec<Harmonic> = Vec::new();
    let mut dco = [0.0; 2];
    let mut dctr = [0.0; 2];
    let mut dt = [0.0; 2];
    let (mut tco, mut tctr, mut tt) = (0.0, 0.0, 0.0);
    let vth = profiles::get().vth;

    for j in mthmin..=mthmax {
        dorb::update(|d| d.mth = j);
        let mut h = Harmonic { mth: j, ..Default::default() };
        let (vminp, vmaxp) = (1.0e-6 * vth, 3.0 * vth);
        if !cfg.nopassing {
            dorb::update(|d| d.sign_vpar = 1.0);
            set_passing();
            let (dd, t) = compute_transport_integral(vminp, vmaxp, cfg.vsteps as usize);
            h.dresco = dd; h.tresco = t;
            dco[0] += dd[0]; dco[1] += dd[1]; tco += t;
            dorb::update(|d| d.sign_vpar = -1.0);
            set_passing();
            let (dd, t) = compute_transport_integral(vminp, vmaxp, cfg.vsteps as usize);
            h.dresctr = dd; h.tresctr = t;
            dctr[0] += dd[0]; dctr[1] += dd[1]; tctr += t;
        }
        dorb::update(|d| d.sign_vpar = 1.0);
        set_trapped();
        let (dd, t) = compute_transport_integral(vminp, vmaxp, cfg.vsteps as usize);
        h.drest = dd; h.trest = t;
        dt[0] += dd[0]; dt[1] += dd[1]; tt += t;
        h.vminp_vth = vminp / vth; h.vmaxp_vth = vmaxp / vth;
        h.vmint_vth = vminp / vth; h.vmaxt_vth = vmaxp / vth;
        harmonics.push(h);
    }

    let p = profiles::get();
    let m_t = p.m_t;
    let dvds = dorb::get().dvds;
    let tot1 = dco[0] + dctr[0] + dt[0];
    let tot2 = dco[1] + dctr[1] + dt[1];

    let mut f = File::create(format!("{}.out", run)).unwrap();
    writeln!(f, " # M_t D11co D11ctr D11t D11 D12co D12ctr D12t D12").unwrap();
    writeln!(f, " {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}",
        m_t, dco[0], dctr[0], dt[0], tot1, dco[1], dctr[1], dt[1], tot2).unwrap();

    if cfg.comptorque {
        let mut f = File::create(format!("{}_torque.out", run)).unwrap();
        writeln!(f, " # s dVds M_t Tco Tctr Tt").unwrap();
        writeln!(f, " {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}", s, dvds, m_t, tco, tctr, tt).unwrap();
    }

    let mut fi = File::create(format!("{}_integral.out", run)).unwrap();
    let mut ft = File::create(format!("{}_torque_integral.out", run)).unwrap();
    for h in &harmonics {
        let t1 = h.dresco[0] + h.dresctr[0] + h.drest[0];
        let t2 = h.dresco[1] + h.dresctr[1] + h.drest[1];
        writeln!(fi, " {:.17e} {} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}",
            m_t, h.mth, h.dresco[0], h.dresctr[0], h.drest[0], t1, h.dresco[1], h.dresctr[1], h.drest[1], t2,
            h.vminp_vth, h.vmaxp_vth, h.vmint_vth, h.vmaxt_vth).unwrap();
        writeln!(ft, " {} {:.17e} {:.17e} {:.17e}", h.mth, h.tresco, h.tresctr, h.trest).unwrap();
    }

    // check_magfie -> _magfie.out
    let mut fm = File::create(format!("{}_magfie.out", run)).unwrap();
    let nth = 50;
    let d = dorb::get();
    for k in 0..nth {
        let th = -PI + k as f64 * (2.0 * PI) / (nth - 1) as f64;
        let (bmod, sqrtg, bder, hcovar, hctrvr) = do_magfie(&[s, 0.0, th]);
        let (bnr, bni) = if cfg.pertfile {
            let (ar, ai) = do_magfie_pert_amp(&[s, 0.0, th]);
            (d.epsmn * ar / bmod, d.epsmn * ai / bmod)
        } else {
            (d.epsmn * (d.m0 as f64 * th).cos(), d.epsmn * (d.m0 as f64 * th).sin())
        };
        let (er, ei) = (d.epsmn * (d.m0 as f64 * th).cos(), d.epsmn * (d.m0 as f64 * th).sin());
        writeln!(fm, " {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e} {:.17e}",
            th, bmod, sqrtg, bder[0], bder[1], bder[2], hcovar[0], hcovar[1], hcovar[2],
            hctrvr[0], hctrvr[1], hctrvr[2], 0.0, 0.0, 0.0, bnr, bni, er, ei).unwrap();
    }
}
