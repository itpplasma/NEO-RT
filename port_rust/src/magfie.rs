//! Flux-surface characterization, literal port of neort_magfie (magfie.f90).
//! Sets driftorbit B0/Bmin/Bmax/etatp/etadt/dVds/th0 and field eps.

use crate::driftorbit::update;
use crate::field::{do_magfie, set_eps};

const PI: f64 = std::f64::consts::PI;

pub fn init_flux_surface_average(s: f64) {
    let nth = 1000;
    let dth = 2.0 * PI / nth as f64;
    let (mut dvds, mut b0, mut eps) = (0.0, 0.0, 0.0);
    let (mut bmin, mut bmax, mut th0) = (-1.0f64, 0.0f64, 0.0f64);
    for k in 1..=nth {
        let th = -PI + k as f64 * 2.0 * PI / nth as f64;
        let (bmod, sqrtg, _bder, _hcov, _hctr) = do_magfie(&[s, 0.0, th]);
        dvds += sqrtg.abs() * dth;
        b0 += bmod * dth;
        eps -= th.cos() * bmod * dth;
        if bmin < 0.0 || bmod < bmin {
            bmin = bmod;
            th0 = th;
        }
        if bmod > bmax {
            bmax = bmod;
        }
    }
    dvds *= 2.0 * PI;
    b0 /= 2.0 * PI;
    eps /= b0 * PI;
    set_eps(eps);
    update(|d| {
        d.dvds = dvds;
        d.b0 = b0;
        d.bmin = bmin;
        d.bmax = bmax;
        d.etatp = 1.0 / bmax;
        d.etadt = 1.0 / bmin;
        d.th0 = th0;
    });
}
