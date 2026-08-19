//! Shared driftorbit state (driftorbit.f90) + orbit th0/noshear, as a
//! thread_local struct (single-threaded gate). Mirrors the C globals.

use std::cell::RefCell;

pub const EPST_SPL: f64 = 1.0e-6;
pub const EPSP_SPL: f64 = 1.0e-6;
pub const EPSST_SPL: f64 = 1.0e-3;
pub const EPSSP_SPL: f64 = 1.0e-3;
pub const EPST: f64 = 1.0e-8;
pub const EPSP: f64 = 1.0e-8;
pub const NLEV: usize = 100;

#[derive(Clone, Copy)]
pub struct DriftOrbit {
    pub efac: f64,
    pub epsmn: f64,
    pub m0: i32,
    pub mth: i32,
    pub mph: i32,
    pub magdrift: bool,
    pub nopassing: bool,
    pub pertfile: bool,
    pub comptorque: bool,
    pub nonlin: bool,
    pub dvds: f64,
    pub etadt: f64,
    pub etatp: f64,
    pub etamin: f64,
    pub etamax: f64,
    pub b0: f64,
    pub bmin: f64,
    pub bmax: f64,
    pub sign_vpar: f64,
    pub sign_vpar_htheta: f64,
    pub noshear: bool,
    pub th0: f64,
}

impl Default for DriftOrbit {
    fn default() -> Self {
        DriftOrbit {
            efac: 1.0, epsmn: 1.0, m0: 1, mth: 1, mph: 1,
            magdrift: true, nopassing: false, pertfile: false, comptorque: true,
            nonlin: false, dvds: 0.0, etadt: 0.0, etatp: 0.0, etamin: 0.0, etamax: 0.0,
            b0: 0.0, bmin: 0.0, bmax: 0.0, sign_vpar: 1.0, sign_vpar_htheta: 1.0,
            noshear: false, th0: 0.0,
        }
    }
}

thread_local! {
    static D: RefCell<DriftOrbit> = RefCell::new(DriftOrbit::default());
}

pub fn get() -> DriftOrbit {
    D.with(|c| *c.borrow())
}
pub fn update(f: impl FnOnce(&mut DriftOrbit)) {
    D.with(|c| f(&mut c.borrow_mut()))
}
