"""Plasma + rotation profiles, port of neort_profiles (profiles.f90). Mutable
module state; uses Collis.loacol_nbi and the spline layer."""
module Profiles

using ..Spline
using ..Collis

export read_and_init_plasma_input, read_and_init_profile_input, init_thermodynamic_forces, PS, ProfState

const QE = 4.803204e-10
const MU = 1.660538e-24
const C = 2.997925e+10
const EV = 1.602176e-12
const SIGN_THETA = -1.0

mutable struct ProfState
    vth::Float64; dvthds::Float64; m_t::Float64; dm_tds::Float64
    om_te::Float64; dom_teds::Float64
    ni1::Float64; ni2::Float64; ti1::Float64; ti2::Float64; te::Float64
    dni1ds::Float64; dti1ds::Float64; a1::Float64; a2::Float64
    qi::Float64; mi::Float64
end
ProfState() = ProfState(0,0,0,0, 0,0, 0,0,0,0,0, 0,0,0,0, QE, 2.014*MU)
const PS = ProfState()

function read_and_init_plasma_input(path::String, s::Float64)
    lines = readlines(path)
    hdr = parse.(Float64, split(lines[2]))
    nplasma = Int(hdr[1]); am1 = hdr[2]; am2 = hdr[3]; z1 = hdr[4]; z2 = hdr[5]
    plasma = zeros(nplasma, 6)
    for k in 1:nplasma
        plasma[k, :] .= parse.(Float64, split(lines[3+k]))[1:6]
    end
    x = plasma[:, 1]
    spl = [spline_coeff(x, plasma[:, kk+1]) for kk in 1:5]
    js = Ref(1)
    e0 = spline_val_0(spl[1], s, js); e1 = spline_val_0(spl[2], s, js)
    e2 = spline_val_0(spl[3], s, js); e3 = spline_val_0(spl[4], s, js)
    e4 = spline_val_0(spl[5], s, js)
    qi = z1*QE; mi = am1*MU
    vth = sqrt(2.0*e2[1]*EV/mi)
    dvthds = 0.5*sqrt(2.0*EV/(mi*e2[1]))*e2[2]
    pmass = 1.6726e-24; v0 = vth
    ebeam = 2.0*pmass*v0^2/(2.0*EV)
    loacol_nbi(2.0, am1, am2, 1.0, z1, z2, e0[1], e1[1], e2[1], e3[1], e4[1], ebeam)
    PS.ni1 = e0[1]; PS.dni1ds = e0[2]; PS.ni2 = e1[1]
    PS.ti1 = e2[1]; PS.dti1ds = e2[2]; PS.ti2 = e3[1]; PS.te = e4[1]
    PS.qi = qi; PS.mi = mi; PS.vth = vth; PS.dvthds = dvthds
end

function read_and_init_profile_input(path::String, s::Float64, r0::Float64, efac::Float64, bfac::Float64)
    rows = [parse.(Float64, split(l)) for l in readlines(path) if length(split(l)) >= 3]
    n = length(rows)
    sx = [r[1] for r in rows]; mt = [r[2] for r in rows]
    spl = spline_coeff(sx, mt)
    js = Ref(1)
    sv = spline_val_0(spl, s, js)
    PS.m_t = sv[1]*efac/bfac
    PS.dm_tds = sv[2]*efac/bfac
    PS.om_te = PS.vth*PS.m_t/r0
    PS.dom_teds = PS.vth*PS.dm_tds/r0 + PS.m_t*PS.dvthds/r0
end

function init_thermodynamic_forces(psi_pr::Float64, q::Float64)
    PS.a1 = PS.dni1ds/PS.ni1 - PS.qi/(PS.ti1*EV)*SIGN_THETA*psi_pr/(q*C)*PS.om_te - 1.5*PS.dti1ds/PS.ti1
    PS.a2 = PS.dti1ds/PS.ti1
end

end # module
