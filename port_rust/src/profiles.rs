//! Plasma + rotation profiles, literal port of neort_profiles (profiles.f90).
//! thread_local state; uses collis::loacol_nbi and the spline layer.

use crate::collis::loacol_nbi;
use crate::spline::{spline_coeff, spline_val_0};
use std::cell::RefCell;
use std::fs::File;
use std::io::{BufRead, BufReader};

pub const QE: f64 = 4.803204e-10;
pub const MU: f64 = 1.660538e-24;
pub const C: f64 = 2.997925e+10;
pub const EV: f64 = 1.602176e-12;
const SIGN_THETA: f64 = -1.0;

#[derive(Default, Clone, Copy)]
pub struct Profiles {
    pub vth: f64,
    pub dvthds: f64,
    pub m_t: f64,
    pub dm_tds: f64,
    pub om_te: f64,
    pub dom_teds: f64,
    pub ni1: f64,
    pub ni2: f64,
    pub ti1: f64,
    pub ti2: f64,
    pub te: f64,
    pub dni1ds: f64,
    pub dti1ds: f64,
    pub a1: f64,
    pub a2: f64,
    pub qi: f64,
    pub mi: f64,
}

thread_local! {
    static P: RefCell<Profiles> = RefCell::new(Profiles { qi: QE, mi: 2.014 * MU, ..Default::default() });
}

pub fn get() -> Profiles {
    P.with(|c| *c.borrow())
}

fn read_rows(path: &str, skip: usize, ncol: usize) -> Vec<Vec<f64>> {
    let f = File::open(path).unwrap_or_else(|_| panic!("cannot open {}", path));
    let mut out = Vec::new();
    for (i, line) in BufReader::new(f).lines().enumerate() {
        let line = line.unwrap();
        if i < skip {
            continue;
        }
        let v: Vec<f64> = line.split_whitespace().filter_map(|t| t.parse().ok()).collect();
        if v.len() >= ncol {
            out.push(v[..ncol].to_vec());
        }
    }
    out
}

pub fn read_and_init_plasma_input(path: &str, s: f64) {
    // header: line0 "% N am1...", line1 "N am1 am2 Z1 Z2", line2 "% s ni...", rows
    let f = File::open(path).unwrap_or_else(|_| panic!("cannot open {}", path));
    let lines: Vec<String> = BufReader::new(f).lines().map(|l| l.unwrap()).collect();
    let hdr: Vec<f64> = lines[1].split_whitespace().map(|t| t.parse().unwrap()).collect();
    let nplasma = hdr[0] as usize;
    let (am1, am2, z1, z2) = (hdr[1], hdr[2], hdr[3], hdr[4]);
    let mut plasma = vec![vec![0.0f64; 6]; nplasma];
    for k in 0..nplasma {
        let v: Vec<f64> = lines[3 + k].split_whitespace().map(|t| t.parse().unwrap()).collect();
        plasma[k] = v[..6].to_vec();
    }
    let x: Vec<f64> = (0..nplasma).map(|i| plasma[i][0]).collect();
    let mut spl = Vec::with_capacity(5);
    for kk in 0..5 {
        let y: Vec<f64> = (0..nplasma).map(|i| plasma[i][kk + 1]).collect();
        spl.push(spline_coeff(&x, &y));
    }
    let nseg = nplasma - 1;
    let mut js = 1i32;
    let ev0 = spline_val_0(&spl[0], nseg, s, &mut js);
    let ev1 = spline_val_0(&spl[1], nseg, s, &mut js);
    let ev2 = spline_val_0(&spl[2], nseg, s, &mut js);
    let ev3 = spline_val_0(&spl[3], nseg, s, &mut js);
    let ev4 = spline_val_0(&spl[4], nseg, s, &mut js);

    let qi = z1 * QE;
    let mi = am1 * MU;
    let vth = (2.0 * ev2[0] * EV / mi).sqrt();
    let dvthds = 0.5 * (2.0 * EV / (mi * ev2[0])).sqrt() * ev2[1];

    let pmass = 1.6726e-24;
    let v0 = vth;
    let ebeam = 2.0 * pmass * v0 * v0 / (2.0 * EV);
    loacol_nbi(2.0, am1, am2, 1.0, z1, z2, ev0[0], ev1[0], ev2[0], ev3[0], ev4[0], ebeam);

    P.with(|c| {
        let mut p = c.borrow_mut();
        p.ni1 = ev0[0]; p.dni1ds = ev0[1];
        p.ni2 = ev1[0];
        p.ti1 = ev2[0]; p.dti1ds = ev2[1];
        p.ti2 = ev3[0];
        p.te = ev4[0];
        p.qi = qi; p.mi = mi; p.vth = vth; p.dvthds = dvthds;
    });
}

pub fn read_and_init_profile_input(path: &str, s: f64, r0: f64, efac: f64, bfac: f64) {
    let rows = read_rows(path, 0, 3);
    let n = rows.len();
    let sx: Vec<f64> = rows.iter().map(|r| r[0]).collect();
    let mt: Vec<f64> = rows.iter().map(|r| r[1]).collect();
    let spl = spline_coeff(&sx, &mt);
    let mut js = 1i32;
    let sv = spline_val_0(&spl, n - 1, s, &mut js);
    P.with(|c| {
        let mut p = c.borrow_mut();
        p.m_t = sv[0] * efac / bfac;
        p.dm_tds = sv[1] * efac / bfac;
        p.om_te = p.vth * p.m_t / r0;
        p.dom_teds = p.vth * p.dm_tds / r0 + p.m_t * p.dvthds / r0;
    });
}

pub fn init_thermodynamic_forces(psi_pr: f64, q: f64) {
    P.with(|c| {
        let mut p = c.borrow_mut();
        p.a1 = p.dni1ds / p.ni1
            - p.qi / (p.ti1 * EV) * SIGN_THETA * psi_pr / (q * C) * p.om_te
            - 1.5 * p.dti1ds / p.ti1;
        p.a2 = p.dti1ds / p.ti1;
    });
}
