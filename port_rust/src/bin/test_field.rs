//! Cross-checks the Rust field port (do_magfie) against the Fortran do_magfie_mod
//! reference (field_ref dump on stdin), to rtol 1e-12.
use neo_rt_rust::field::{do_magfie, do_magfie_init, set_bfac, set_inp_swi};
use std::io::{self, Read};

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-13 + 1e-12 * b.abs()
}

fn main() {
    let in_file = std::env::args().nth(1).expect("usage: test_field <in_file>");
    set_inp_swi(9);
    set_bfac(1.0);
    do_magfie_init(&in_file);

    let mut buf = String::new();
    io::stdin().read_to_string(&mut buf).unwrap();
    let nums: Vec<f64> = buf.split_whitespace().map(|t| t.parse().unwrap()).collect();
    let s = 0.5;
    let names = ["bmod", "sqrtg", "bder1", "bder3", "hcov2", "hcov3", "hctr2", "hctr3"];
    let mut checks = 0;
    let mut fails = 0;
    let nrow = nums.len() / 11;
    for i in 0..nrow {
        let r = &nums[i * 11..i * 11 + 11];
        let theta = r[0];
        let (bmod, sqrtg, bder, hcovar, hctrvr) = do_magfie(&[s, 0.0, theta]);
        let got = [bmod, sqrtg, bder[0], bder[2], hcovar[1], hcovar[2], hctrvr[1], hctrvr[2]];
        let refv = [r[1], r[2], r[3], r[4], r[5], r[6], r[7], r[8]];
        for k in 0..8 {
            checks += 1;
            if !close(got[k], refv[k]) {
                eprintln!("theta={:.5} {} R={:.16e} ref={:.16e}", theta, names[k], got[k], refv[k]);
                fails += 1;
            }
        }
    }
    println!("field: {} checks, {} failures", checks, fails);
    std::process::exit(if fails == 0 { 0 } else { 1 });
}
