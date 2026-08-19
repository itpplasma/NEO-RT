"""Minimal ccall bindings to SUNDIALS CVODE (Adams) + fixed-point solver + root
finding -- the same library as the C/Rust ports, so trajectories are identical.
Single-threaded. SUN_COMM_NULL = MPI_COMM_NULL = &ompi_mpi_comm_null."""
module Cvode

const LIBCV = "libsundials_cvode"
const LIBCORE = "libsundials_core"
const LIBNVEC = "libsundials_nvecserial"
const LIBFP = "libsundials_sunnonlinsolfixedpoint"
const CV_ADAMS = Cint(1)
const CV_NORMAL = Cint(1)
const CV_ROOT_RETURN = Cint(2)

# MPI_COMM_NULL handle (OpenMPI): address of the global ompi_mpi_comm_null.
sun_comm_null() = cglobal((:ompi_mpi_comm_null, "libmpi"))

ctx_create() = begin
    ctx = Ref{Ptr{Cvoid}}(C_NULL)
    ccall((:SUNContext_Create, LIBCORE), Cint, (Ptr{Cvoid}, Ref{Ptr{Cvoid}}), sun_comm_null(), ctx)
    ctx[]
end
ctx_free(ctx) = (r = Ref(ctx); ccall((:SUNContext_Free, LIBCORE), Cint, (Ref{Ptr{Cvoid}},), r))

nvnew(n, ctx) = ccall((:N_VNew_Serial, LIBNVEC), Ptr{Cvoid}, (Clong, Ptr{Cvoid}), n, ctx)
nvdata(v) = ccall((:N_VGetArrayPointer, LIBNVEC), Ptr{Float64}, (Ptr{Cvoid},), v)
nvdestroy(v) = ccall((:N_VDestroy, LIBNVEC), Cvoid, (Ptr{Cvoid},), v)

cvcreate(lmm, ctx) = ccall((:CVodeCreate, LIBCV), Ptr{Cvoid}, (Cint, Ptr{Cvoid}), lmm, ctx)
cvinit(mem, f, t0, y0) = ccall((:CVodeInit, LIBCV), Cint, (Ptr{Cvoid}, Ptr{Cvoid}, Float64, Ptr{Cvoid}), mem, f, t0, y0)
cvtol(mem, rt, at) = ccall((:CVodeSStolerances, LIBCV), Cint, (Ptr{Cvoid}, Float64, Float64), mem, rt, at)
cvmaxsteps(mem, n) = ccall((:CVodeSetMaxNumSteps, LIBCV), Cint, (Ptr{Cvoid}, Clong), mem, n)
nls_fixedpoint(y, m, ctx) = ccall((:SUNNonlinSol_FixedPoint, LIBFP), Ptr{Cvoid}, (Ptr{Cvoid}, Cint, Ptr{Cvoid}), y, m, ctx)
cvsetnls(mem, nls) = ccall((:CVodeSetNonlinearSolver, LIBCV), Cint, (Ptr{Cvoid}, Ptr{Cvoid}), mem, nls)
nlsfree(nls) = ccall((:SUNNonlinSolFree, LIBCORE), Cint, (Ptr{Cvoid},), nls)
cvrootinit(mem, n, g) = ccall((:CVodeRootInit, LIBCV), Cint, (Ptr{Cvoid}, Cint, Ptr{Cvoid}), mem, n, g)
cvode(mem, tout, y, tret, itask) = ccall((:CVode, LIBCV), Cint, (Ptr{Cvoid}, Float64, Ptr{Cvoid}, Ref{Float64}, Cint), mem, tout, y, tret, itask)
cvfree(mem) = (r = Ref(mem); ccall((:CVodeFree, LIBCV), Cvoid, (Ref{Ptr{Cvoid}},), r))

end # module
