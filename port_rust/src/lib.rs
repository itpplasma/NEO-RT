//! Rust port of NEO-RT (golden path). Modules mirror the Fortran/C structure;
//! numerics use SUNDIALS CVODE + LAPACK via FFI for parity with the C port.

pub mod spline;
pub mod field;
pub mod collis;
pub mod profiles;
pub mod cvode;
pub mod driftorbit;
pub mod magfie;
pub mod orbit;
pub mod freq;
pub mod resonance;
pub mod transport;
