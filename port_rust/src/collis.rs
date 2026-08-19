//! Collision coefficients, literal port of collis_alp (collis_nbi.f90). Module
//! state efcolf/velrat/enrat in a thread_local (single-threaded gate). collis
//! keeps its local ev=1.6022e-12 distinct from util ev, faithful to Fortran.

use std::cell::RefCell;

pub const NSORTS: usize = 3;

thread_local! {
    static ST: RefCell<[[f64; NSORTS]; 3]> = RefCell::new([[0.0; NSORTS]; 3]);
    // [0]=efcolf, [1]=velrat, [2]=enrat
}

pub fn efcolf() -> [f64; NSORTS] {
    ST.with(|s| s.borrow()[0])
}
pub fn velrat() -> [f64; NSORTS] {
    ST.with(|s| s.borrow()[1])
}
pub fn enrat() -> [f64; NSORTS] {
    ST.with(|s| s.borrow()[2])
}

fn onseff(v: f64) -> (f64, f64, f64) {
    let sqp = 1.7724538;
    let cons = 0.75225278;
    let v2 = v * v;
    let v3 = v2 * v;
    if v < 0.01 {
        (cons * (1.0 - 0.6 * v2), cons * (1.0 - 0.2 * v2), 2.0 * cons * (1.0 - 1.2 * v2))
    } else if v > 6.0 {
        (1.0 / v3, (1.0 - 0.5 / v2) / v, -1.0 / v3)
    } else {
        let ex = (-v2).exp() / sqp;
        let er = erf(v);
        let d_p = er / v3 - 2.0 * ex / v2;
        (d_p, er * (1.0 - 0.5 / v2) / v + ex / v2, 4.0 * ex - d_p)
    }
}

/// dpp, dhh, fpeff (coleff).
pub fn coleff(p: f64) -> (f64, f64, f64) {
    let plim = if p > 1.0e-8 { p } else { 1.0e-8 };
    let (ef, vr, en) = (efcolf(), velrat(), enrat());
    let (mut dpp, mut dhh, mut fpeff) = (0.0, 0.0, 0.0);
    for i in 0..NSORTS {
        let (d_p, dh, dpd) = onseff(p * vr[i]);
        dpp += d_p * ef[i];
        dhh += dh * ef[i];
        fpeff += (dpd / plim - 2.0 * d_p * p * en[i]) * ef[i];
    }
    dhh /= plim * plim;
    (dpp, dhh, fpeff)
}

#[allow(clippy::too_many_arguments)]
pub fn loacol_nbi(amb: f64, am1: f64, am2: f64, zb: f64, z1: f64, z2: f64, densi1: f64, densi2: f64, tempi1: f64, tempi2: f64, tempe: f64, ebeam: f64) -> f64 {
    let pi = 3.14159265358979;
    let (pmass, emass, e, ev): (f64, f64, f64, f64) = (1.6726e-24, 9.1094e-28, 4.8032e-10, 1.6022e-12);
    let eps = f64::EPSILON;
    let v0 = (2.0 * ebeam * ev / (amb * pmass)).sqrt();
    let vti1 = (2.0 * tempi1 * ev / (pmass * am1)).sqrt();
    let vti2 = (2.0 * tempi2 * ev / (pmass * am2)).sqrt();
    let vte = (2.0 * tempe * ev / emass).sqrt();
    let dense = densi1 * z1 + densi2 * z2;
    let a1 = (densi1 * z1 * z1 / tempi1).sqrt() * zb * z1 * (amb + am1) / (amb * tempi1 + am1 * ebeam);
    let a2 = (densi2 * z2 * z2 / tempi2).sqrt() * zb * z2 * (amb + am2) / (amb * tempi2 + am2 * ebeam);
    let alami1 = 23.0 - (if a1 > eps { a1 } else { eps }).ln();
    let alami2 = 23.0 - (if a2 > eps { a2 } else { eps }).ln();
    let alame = 24.0 - (dense.sqrt() / tempe).ln();
    let mut frecol = 2.0 * pi * dense * e.powi(4) * zb * zb / ((amb * pmass).powi(2) * v0.powi(3));
    frecol /= v0;
    ST.with(|s| {
        let mut st = s.borrow_mut();
        st[2] = [ebeam / tempi1, ebeam / tempi2, ebeam / tempe];
        st[1] = [v0 / vti1, v0 / vti2, v0 / vte];
        st[0] = [
            frecol * z1 * z1 * alami1 * densi1 / dense,
            frecol * z2 * z2 * alami2 * densi2 / dense,
            frecol * alame,
        ];
        for i in 0..NSORTS {
            st[0][i] *= st[1][i];
        }
    });
    v0
}

/// erf via libm (matches Fortran intrinsic erf).
fn erf(x: f64) -> f64 {
    extern "C" {
        fn erf(x: f64) -> f64;
    }
    unsafe { erf(x) }
}
