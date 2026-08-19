//! Cross-checks the Rust profiles/collis port against the Fortran modules
//! (profiles_ref dump on stdin: 21 scalars), to rtol 1e-12.
use neo_rt_rust::collis::{efcolf, enrat, velrat};
use neo_rt_rust::field::{do_magfie, do_magfie_init, set_bfac, set_inp_swi, with_field};
use neo_rt_rust::profiles::{get, init_thermodynamic_forces, read_and_init_plasma_input, read_and_init_profile_input};
use std::io::{self, Read};

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-13 + 1e-12 * b.abs()
}

fn main() {
    let a: Vec<String> = std::env::args().collect();
    let s = 0.5;
    set_inp_swi(9);
    set_bfac(1.0);
    do_magfie_init(&a[1]);
    let _ = do_magfie(&[s, 0.0, 0.0]);
    let (r0, psi_pr, q) = with_field(|f| (f.r0, f.psi_pr, f.q));
    read_and_init_plasma_input(&a[2], s);
    read_and_init_profile_input(&a[3], s, r0, 1.0, 1.0);
    init_thermodynamic_forces(psi_pr, q);

    let p = get();
    let ef = efcolf();
    let vr = velrat();
    let en = enrat();
    let got = [
        p.vth, p.m_t, p.om_te, p.a1, p.a2, p.ni1, p.ni2, p.ti1, p.ti2, p.te, p.qi, p.mi,
        ef[0], ef[1], ef[2], vr[0], vr[1], vr[2], en[0], en[1], en[2],
    ];
    let names = [
        "vth", "M_t", "Om_tE", "A1", "A2", "ni1", "ni2", "Ti1", "Ti2", "Te", "qi", "mi",
        "efcolf1", "efcolf2", "efcolf3", "velrat1", "velrat2", "velrat3", "enrat1", "enrat2", "enrat3",
    ];
    let mut buf = String::new();
    io::stdin().read_to_string(&mut buf).unwrap();
    let refv: Vec<f64> = buf.split_whitespace().map(|t| t.parse().unwrap()).collect();
    let mut fails = 0;
    for k in 0..21 {
        if !close(got[k], refv[k]) {
            eprintln!("{} R={:.16e} ref={:.16e}", names[k], got[k], refv[k]);
            fails += 1;
        }
    }
    println!("profiles: 21 checks, {} failures", fails);
    std::process::exit(if fails == 0 { 0 } else { 1 });
}
