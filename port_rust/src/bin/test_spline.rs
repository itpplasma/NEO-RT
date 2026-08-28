//! Cross-checks the Rust spline against the Fortran itpplasma spline reference
//! (spline_ref dump on stdin: COEFF and VAL lines), to rtol 1e-12.
use neo_rt_rust::spline::{spline_coeff, spline_val_0};
use std::io::{self, Read};

const NP: usize = 9;

fn close(a: f64, b: f64) -> bool {
    (a - b).abs() <= 1e-14 + 1e-12 * b.abs()
}

fn main() {
    let x: Vec<f64> = (0..NP).map(|i| i as f64 * 0.37).collect();
    let y: Vec<f64> = x.iter().map(|&xi| (1.3 * xi).sin() + 0.5 * xi).collect();
    let coeff = spline_coeff(&x, &y);

    let mut input = String::new();
    io::stdin().read_to_string(&mut input).unwrap();
    let mut tok = input.split_whitespace();

    let mut checks = 0;
    let mut fails = 0;
    // coefficients
    for _ in 0..(NP - 1) {
        assert_eq!(tok.next().unwrap(), "COEFF");
        let idx: usize = tok.next().unwrap().parse().unwrap();
        for k in 0..5 {
            let refv: f64 = tok.next().unwrap().parse().unwrap();
            checks += 1;
            if !close(coeff[idx * 5 + k], refv) {
                eprintln!("COEFF[{}][{}] R={:.16e} ref={:.16e}", idx, k, coeff[idx * 5 + k], refv);
                fails += 1;
            }
        }
    }
    // evaluations
    let mut jstart = 1i32;
    while let Some(tag) = tok.next() {
        assert_eq!(tag, "VAL");
        let xe: f64 = tok.next().unwrap().parse().unwrap();
        let refv: [f64; 3] = [
            tok.next().unwrap().parse().unwrap(),
            tok.next().unwrap().parse().unwrap(),
            tok.next().unwrap().parse().unwrap(),
        ];
        let out = spline_val_0(&coeff, NP - 1, xe, &mut jstart);
        for k in 0..3 {
            checks += 1;
            if !close(out[k], refv[k]) {
                eprintln!("VAL(x={:.6})[{}] R={:.16e} ref={:.16e}", xe, k, out[k], refv[k]);
                fails += 1;
            }
        }
    }
    println!("spline: {} checks, {} failures", checks, fails);
    std::process::exit(if fails == 0 { 0 } else { 1 });
}
