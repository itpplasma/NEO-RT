//! Cross-checks the Rust orbit port (CVODE) against the Fortran DVODE orbit
//! reference (orbit_ref dump on stdin: 3 rows x 6). CVODE and DVODE share the
//! VODE Adams lineage, so within the 1e-8 gate.
use neo_rt_rust::driftorbit as dorb;
use neo_rt_rust::field::{do_magfie, do_magfie_init, set_bfac, set_inp_swi, with_field};
use neo_rt_rust::magfie::init_flux_surface_average;
use neo_rt_rust::orbit::{bounce, bounce_time, set_s};
use neo_rt_rust::profiles::{get, read_and_init_plasma_input, read_and_init_profile_input};
use std::io::{self, Read};

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-12 + 1e-8 * b.abs()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let s = 0.5;
    set_inp_swi(9);
    set_bfac(1.0);
    do_magfie_init(&a[1]);
    let _ = do_magfie(&[s, 0.0, 0.0]);
    let r0 = with_field(|f| f.r0);
    init_flux_surface_average(s);
    read_and_init_plasma_input(&a[2], s);
    read_and_init_profile_input(&a[3], s, r0, 1.0, 1.0);
    set_s(s);
    dorb::update(|d| d.sign_vpar = 1.0);

    let d = dorb::get();
    let vth = get().vth;
    let vv = [vth, vth, 1.5 * vth];
    let ee = [0.5 * (d.etatp + d.etadt), 0.3 * d.etatp, 0.8 * d.etatp + 0.2 * d.etadt];

    let mut buf = String::new();
    io::stdin().read_to_string(&mut buf).unwrap();
    let refv: Vec<f64> = buf.split_whitespace().map(|t| t.parse().unwrap()).collect();
    let names = ["bounce_time", "taub", "bavg1", "bavg2", "bavg3", "bavg6"];
    let mut fails = 0;
    let mut checks = 0;
    for i in 0..3 {
        let tt = bounce_time(vv[i], ee[i], 0.0, false);
        let (taub, bavg) = bounce(vv[i], ee[i], 0.0, false);
        let got = [tt, taub, bavg[0], bavg[1], bavg[2], bavg[5]];
        for k in 0..6 {
            let r = refv[i * 6 + k];
            checks += 1;
            if !close(got[k], r) {
                eprintln!("case {} {} R={:.16e} ref={:.16e} rel={:.2e}", i, names[k], got[k], r, ((got[k] - r) / r).abs());
                fails += 1;
            }
        }
    }
    println!("orbit: {} checks, {} failures", checks, fails);
    std::process::exit(if fails == 0 { 0 } else { 1 });
}
