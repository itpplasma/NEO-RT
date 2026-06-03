//! Rust port of NEO-RT (golden path). Modules mirror the Fortran/C structure;
//! numerics use SUNDIALS CVODE + LAPACK via FFI for parity with the C port.

pub mod spline;
