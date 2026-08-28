//! Minimal FFI to SUNDIALS CVODE (Adams) + fixed-point nonlinear solver +
//! root finding -- the same library the C port uses, so the Rust orbit
//! integration is bit-identical to C. Single-threaded use only.
//!
//! sunrealtype is f64 in the default SUNDIALS build. Opaque handles are raw
//! pointers. SUN_COMM_NULL is a null comm (SUNDIALS here is MPI-enabled, hence
//! libmpi is linked in build.rs); passing null works for the serial context.

use std::os::raw::{c_int, c_long, c_void};

pub const CV_ADAMS: c_int = 1;
pub const CV_NORMAL: c_int = 1;
pub const CV_ROOT_RETURN: c_int = 2;

pub type SunContext = *mut c_void;
pub type NVector = *mut c_void;
pub type CvodeMem = *mut c_void;
pub type SunNonlinSol = *mut c_void;

pub type CvRhsFn = extern "C" fn(t: f64, y: NVector, ydot: NVector, user_data: *mut c_void) -> c_int;
pub type CvRootFn = extern "C" fn(t: f64, y: NVector, gout: *mut f64, user_data: *mut c_void) -> c_int;

extern "C" {
    // OpenMPI's MPI_COMM_NULL is &ompi_mpi_comm_null; SUN_COMM_NULL maps to it.
    // Passing a plain null pointer makes SUNDIALS call MPI_Comm_dup -> abort.
    pub static ompi_mpi_comm_null: u8;
}

/// The SUNComm value equal to C's SUN_COMM_NULL (MPI_COMM_NULL) on this build.
pub fn sun_comm_null() -> *mut c_void {
    unsafe { &ompi_mpi_comm_null as *const u8 as *mut c_void }
}

extern "C" {
    pub fn SUNContext_Create(comm: *mut c_void, ctx: *mut SunContext) -> c_int;
    pub fn SUNContext_Free(ctx: *mut SunContext) -> c_int;
    pub fn N_VNew_Serial(len: c_long, ctx: SunContext) -> NVector;
    pub fn N_VGetArrayPointer(v: NVector) -> *mut f64;
    pub fn N_VDestroy(v: NVector);
    pub fn CVodeCreate(lmm: c_int, ctx: SunContext) -> CvodeMem;
    pub fn CVodeInit(mem: CvodeMem, f: CvRhsFn, t0: f64, y0: NVector) -> c_int;
    pub fn CVodeSStolerances(mem: CvodeMem, reltol: f64, abstol: f64) -> c_int;
    pub fn CVodeSetMaxNumSteps(mem: CvodeMem, mxsteps: c_long) -> c_int;
    pub fn CVodeSetUserData(mem: CvodeMem, user_data: *mut c_void) -> c_int;
    pub fn SUNNonlinSol_FixedPoint(y: NVector, m: c_int, ctx: SunContext) -> SunNonlinSol;
    pub fn CVodeSetNonlinearSolver(mem: CvodeMem, nls: SunNonlinSol) -> c_int;
    pub fn SUNNonlinSolFree(nls: SunNonlinSol) -> c_int;
    pub fn CVodeRootInit(mem: CvodeMem, nrtfn: c_int, g: CvRootFn) -> c_int;
    pub fn CVode(mem: CvodeMem, tout: f64, yout: NVector, tret: *mut f64, itask: c_int) -> c_int;
    pub fn CVodeFree(mem: *mut CvodeMem);
}
