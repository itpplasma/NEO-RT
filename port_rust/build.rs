// Link the same numerical packages as the C port: SUNDIALS CVODE (Adams + root
// finding) and LAPACK/BLAS (dptsv for the spline solve). SUNDIALS here is
// MPI-enabled, so libmpi is required for SUN_COMM_NULL.
fn main() {
    println!("cargo:rustc-link-lib=dylib=sundials_cvode");
    println!("cargo:rustc-link-lib=dylib=sundials_core");
    println!("cargo:rustc-link-lib=dylib=mpi");
    println!("cargo:rustc-link-lib=dylib=lapack");
    println!("cargo:rustc-link-lib=dylib=blas");
    println!("cargo:rustc-link-lib=dylib=m");
}
