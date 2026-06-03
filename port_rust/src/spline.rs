//! Cubic spline per Sormann, literal port of itpplasma/spline (spline.f90),
//! matching the C port. Coefficient layout: (n-1) rows x 5, row-major:
//! [x_i, y_i, c1, c2, c3]. Uses LAPACK dptsv for the tridiagonal solve so the
//! coefficients are bit-identical to the Fortran path.

extern "C" {
    fn dptsv_(
        n: *const i32,
        nrhs: *const i32,
        d: *mut f64,
        e: *mut f64,
        b: *mut f64,
        ldb: *const i32,
        info: *mut i32,
    );
}

/// Build spline coefficients for knots (x, y). Returns a flat (n-1)*5 vector.
pub fn spline_coeff(x: &[f64], y: &[f64]) -> Vec<f64> {
    let n = x.len();
    let m = n - 1; // segments
    let mut h = vec![0.0f64; m];
    let mut r = vec![0.0f64; m];
    for i in 0..m {
        h[i] = x[i + 1] - x[i];
        r[i] = y[i + 1] - y[i];
    }

    let sys = n - 2;
    let mut d = vec![0.0f64; sys];
    let mut dl = vec![0.0f64; if sys > 1 { sys - 1 } else { 1 }];
    let mut c = vec![0.0f64; sys];
    for k in 0..sys {
        d[k] = 2.0 * (h[k] + h[k + 1]);
        c[k] = 3.0 * (r[k + 1] / h[k + 1] - r[k] / h[k]);
    }
    for k in 0..sys.saturating_sub(1) {
        dl[k] = h[k + 1];
    }
    let (nn, nrhs, ldb) = (sys as i32, 1i32, sys as i32);
    let mut info = 0i32;
    unsafe {
        dptsv_(&nn, &nrhs, d.as_mut_ptr(), dl.as_mut_ptr(), c.as_mut_ptr(), &ldb, &mut info);
    }

    let mut coeff = vec![0.0f64; m * 5];
    let cset = |co: &mut [f64], i: usize, k: usize, v: f64| co[i * 5 + k] = v;
    for i in 0..m {
        cset(&mut coeff, i, 0, x[i]);
        cset(&mut coeff, i, 1, y[i]);
    }
    // column 3 (index 2)
    coeff[0 * 5 + 2] = r[0] / h[0] - h[0] / 3.0 * c[0];
    for i in 1..=(n - 3) {
        coeff[i * 5 + 2] = r[i] / h[i] - h[i] / 3.0 * (c[i] + 2.0 * c[i - 1]);
    }
    coeff[(m - 1) * 5 + 2] = r[m - 1] / h[m - 1] - h[m - 1] / 3.0 * (2.0 * c[sys - 1]);
    // column 4 (index 3)
    coeff[0 * 5 + 3] = 0.0;
    for i in 1..=(m - 1) {
        coeff[i * 5 + 3] = c[i - 1];
    }
    // column 5 (index 4)
    coeff[0 * 5 + 4] = 1.0 / (3.0 * h[0]) * c[0];
    for i in 1..=(n - 3) {
        coeff[i * 5 + 4] = 1.0 / (3.0 * h[i]) * (c[i] - c[i - 1]);
    }
    coeff[(m - 1) * 5 + 4] = 1.0 / (3.0 * h[m - 1]) * (-c[sys - 1]);
    coeff
}

/// Evaluate the spline at `x`. Returns [value, 1st deriv, 2nd deriv]. `j_start`
/// carries the saved interval guess between calls (value-independent; x clamped).
pub fn spline_val_0(coeff: &[f64], nseg: usize, x: f64, j_start: &mut i32) -> [f64; 3] {
    let n = nseg + 1;
    let mut j = *j_start;
    if j < 1 || j >= n as i32 {
        j = 1;
    }
    if n <= 1 {
        return [0.0, 0.0, 0.0];
    }
    let c = |i: usize, k: usize| coeff[i * 5 + k];
    let mut ji = (j - 1) as usize;
    if x < c(ji, 0) {
        while ji > 0 {
            if x < c(ji, 0) {
                ji -= 1;
            } else {
                break;
            }
        }
    } else {
        while ji < n - 2 {
            if x >= c(ji + 1, 0) {
                ji += 1;
            } else {
                break;
            }
        }
    }
    *j_start = ji as i32 + 1;
    let z = x - c(ji, 0);
    [
        ((c(ji, 4) * z + c(ji, 3)) * z + c(ji, 2)) * z + c(ji, 1),
        (3.0 * c(ji, 4) * z + 2.0 * c(ji, 3)) * z + c(ji, 2),
        6.0 * c(ji, 4) * z + 2.0 * c(ji, 3),
    ]
}
