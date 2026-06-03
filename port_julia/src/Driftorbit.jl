"""Shared driftorbit state (driftorbit.f90) + orbit th0/noshear, as a mutable
module struct (single-threaded gate). Mirrors the C/Rust globals."""
module Driftorbit

const EPST_SPL = 1.0e-6
const EPSP_SPL = 1.0e-6
const EPSST_SPL = 1.0e-3
const EPSSP_SPL = 1.0e-3
const EPST = 1.0e-8
const EPSP = 1.0e-8
const NLEV = 100

mutable struct DriftOrbitState
    efac::Float64; epsmn::Float64
    m0::Int; mth::Int; mph::Int
    magdrift::Bool; nopassing::Bool; pertfile::Bool; comptorque::Bool; nonlin::Bool
    dvds::Float64; etadt::Float64; etatp::Float64; etamin::Float64; etamax::Float64
    b0::Float64; bmin::Float64; bmax::Float64
    sign_vpar::Float64; sign_vpar_htheta::Float64; noshear::Bool; th0::Float64
end
DriftOrbitState() = DriftOrbitState(1.0, 1.0, 1, 1, 1, true, false, false, true, false,
    0,0,0,0,0, 0,0,0, 1.0, 1.0, false, 0.0)
const DS = DriftOrbitState()

end # module
