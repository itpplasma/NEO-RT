"""
Axisymmetric + perturbation magnetic field, port of do_magfie_mod /
do_magfie_pert_mod (do_magfie_standalone.f90). Mutable module state (single
threaded gate). Mirrors the verified C/Rust field.
"""
module Field

include("Spline.jl")
using .Spline

export do_magfie_init, do_magfie, do_magfie_pert_init, do_magfie_pert_amp, FS, FieldState

const SIGN_THETA = -1.0
const ITOB = 2.0e-1 * SIGN_THETA
const PI = pi

mutable struct FieldState
    bfac::Float64
    inp_swi::Int
    psi_pr::Float64; r0::Float64; a::Float64; b00::Float64
    bthcov::Float64; bphcov::Float64; dbthcovds::Float64; dbphcovds::Float64
    q::Float64; dqds::Float64; iota::Float64; eps::Float64; b0h::Float64
    nflux::Int; nmode::Int; ncol1::Int; ncol2::Int
    params0::Matrix{Float64}; modes0::Array{Float64,3}
    spl1::Vector{Matrix{Float64}}; spl2::Vector{Matrix{Float64}}
    nseg::Int; jstart::Base.RefValue{Int}
    pert_mph::Int; p_ncol2::Int; p_nflux::Int; p_nmode::Int; p_nfp::Int; p_nseg::Int
    p_params::Matrix{Float64}; p_modes::Array{Float64,3}; p_spl2::Vector{Matrix{Float64}}
end

FieldState() = FieldState(1.0, 0, 0,0,0,0, 0,0,0,0, 0,0,0,0,0, 0,0,5,0,
    zeros(0,0), zeros(0,0,0), Matrix{Float64}[], Matrix{Float64}[], 0, Ref(1),
    0,0,0,0,0,0, zeros(0,0), zeros(0,0,0), Matrix{Float64}[])

const FS = FieldState()

function boozer_read(path::String, ncol1::Int, nc2::Int)
    lines = readlines(path)
    idx = 6  # skip 5 header lines
    hdr = parse.(Float64, split(lines[idx])); idx += 1
    m0b = Int(hdr[1]); n0b = Int(hdr[2]); nflux = Int(hdr[3]); nfp = Int(hdr[4])
    flux = hdr[5]; a = hdr[6]; r0 = hdr[7]
    nmode = (m0b + 1) * (n0b + 1)
    pcols = ncol1 + 1; mcols = nc2 + 2
    params = zeros(Float64, nflux, pcols)
    modes = zeros(Float64, nflux, nmode, mcols)
    for ks in 1:nflux
        idx += 2  # two label lines
        pv = parse.(Float64, split(lines[idx])); idx += 1
        params[ks, :] .= pv[1:pcols]
        idx += 1  # mode-label line
        for j in 1:nmode
            mv = parse.(Float64, split(lines[idx])); idx += 1
            modes[ks, j, :] .= mv[1:mcols]
        end
    end
    return nflux, nmode, nfp, flux, a, r0, params, modes
end

function build_splines(params, modes, nflux, nmode, ncol1, nc2)
    x = params[:, 1]
    s1 = Matrix{Float64}[]
    for k in 1:ncol1
        push!(s1, spline_coeff(x, params[:, k+1]))
    end
    s2 = Vector{Matrix{Float64}}(undef, nc2 * nmode)
    for j in 1:nmode, k in 1:nc2
        s2[(k-1)*nmode + j] = spline_coeff(x, modes[:, j, k+2])
    end
    return s1, s2
end

function do_magfie_init(path::String)
    FS.ncol1 = 5
    FS.ncol2 = FS.inp_swi == 8 ? 4 : 8
    nflux, nmode, _nfp, flux, a, r0, params, modes = boozer_read(path, FS.ncol1, FS.ncol2)
    FS.nflux = nflux; FS.nmode = nmode
    FS.a = 100.0 * a; FS.r0 = 100.0 * r0
    FS.psi_pr = 1.0e8 * flux / (2.0 * PI) * FS.bfac
    FS.nseg = nflux - 1
    FS.params0 = params; FS.modes0 = modes
    FS.spl1, FS.spl2 = build_splines(params, modes, nflux, nmode, FS.ncol1, FS.ncol2)
    FS.r0 = modes[1, 1, 3] * 100.0
    FS.b00 = 1.0e4 * modes[1, 1, 6] * FS.bfac
end

function do_magfie(x::NTuple{3,Float64})
    nm = FS.nmode; bf = FS.bfac
    x1 = max(FS.params0[1, 1], x[1]); x1 = min(FS.params0[FS.nflux, 1], x1)
    js = FS.jstart
    sv = spline_val_0(FS.spl1[3], x1, js)
    FS.bthcov = ITOB * sv[1] * bf; FS.dbthcovds = ITOB * sv[2] * bf
    sv = spline_val_0(FS.spl1[2], x1, js)
    FS.bphcov = ITOB * sv[1] * bf; FS.dbphcovds = ITOB * sv[2] * bf
    sv = spline_val_0(FS.spl1[1], x1, js)
    FS.iota = sv[1]; FS.q = 1.0 / FS.iota; FS.dqds = -sv[2] / FS.iota^2

    m0 = FS.modes0[1, 1, 1]; m1 = FS.modes0[1, 2, 1]; dm = m1 - m0
    cost = zeros(nm); sint = zeros(nm)
    ffr = cos(m0 * x[3]); ffi = sin(m0 * x[3])
    rotr = cos(dm * x[3]); roti = sin(dm * x[3])
    for j in 1:nm
        cost[j] = ffr; sint[j] = ffi
        ffr, ffi = ffr*rotr - ffi*roti, ffr*roti + ffi*rotr
    end

    bm = 0.0; db = 0.0; bth = 0.0
    if FS.inp_swi == 8
        for j in 1:nm
            sv = spline_val_0(FS.spl2[(4-1)*nm + j], x1, js)
            bmnc = 1.0e4*sv[1]*bf; dbmnc = 1.0e4*sv[2]*bf
            j == 1 && (FS.b0h = bmnc)
            bm += bmnc*cost[j]; db += dbmnc*cost[j]
            bth += -FS.modes0[1, j, 1]*bmnc*sint[j]
        end
    else
        for j in 1:nm
            sv = spline_val_0(FS.spl2[(7-1)*nm + j], x1, js)
            bmnc = 1.0e4*sv[1]*bf; dbmnc = 1.0e4*sv[2]*bf
            sv = spline_val_0(FS.spl2[(8-1)*nm + j], x1, js)
            bmns = 1.0e4*sv[1]*bf; dbmns = 1.0e4*sv[2]*bf
            j == 1 && (FS.b0h = bmnc)
            mj = FS.modes0[1, j, 1]
            bm += bmnc*cost[j] + bmns*sint[j]
            db += dbmnc*cost[j] + dbmns*sint[j]
            bth += -mj*bmnc*sint[j] + mj*bmns*cost[j]
        end
    end
    bder = (db/bm, 0.0, bth/bm)
    sqgbmod2 = SIGN_THETA * FS.psi_pr * (FS.bphcov + FS.iota * FS.bthcov)
    sqgbmod = sqgbmod2 / bm
    sqrtg = sqgbmod / bm
    hcovar = (0.0, FS.bphcov/bm, FS.bthcov/bm)
    hctrvr = (0.0, SIGN_THETA*FS.psi_pr/sqgbmod, SIGN_THETA*FS.iota*FS.psi_pr/sqgbmod)
    return bm, sqrtg, bder, hcovar, hctrvr
end

function do_magfie_pert_init(path::String)
    nc2 = FS.inp_swi == 8 ? 4 : 8
    nflux, nmode, nfp, _flux, _a, _r0, params, modes = boozer_read(path, 5, nc2)
    FS.p_ncol2 = nc2; FS.p_nflux = nflux; FS.p_nmode = nmode; FS.p_nfp = nfp
    FS.p_nseg = nflux - 1
    FS.pert_mph = round(Int, nfp * modes[1, 1, 2])
    _s1, s2 = build_splines(params, modes, nflux, nmode, 5, nc2)
    FS.p_params = params; FS.p_modes = modes; FS.p_spl2 = s2
end

function do_magfie_pert_amp(x::NTuple{3,Float64})
    nm = FS.p_nmode; bf = FS.bfac
    x1 = max(FS.p_params[1, 1], x[1]); x1 = min(FS.p_params[FS.p_nflux, 1], x1)
    js = FS.jstart
    sr = 0.0; si = 0.0
    if FS.inp_swi == 8
        for j in 1:nm
            sv = spline_val_0(FS.p_spl2[(4-1)*nm + j], x1, js)
            bmnc = 1.0e4*sv[1]*bf
            sr += bmnc * cos(FS.p_modes[1, j, 1] * x[3])
        end
    else
        m0 = FS.p_modes[1, 1, 1]; m1 = FS.p_modes[1, 2, 1]; dm = m1 - m0
        ffr = cos(m0*x[3]); ffi = sin(m0*x[3]); rotr = cos(dm*x[3]); roti = sin(dm*x[3])
        for j in 1:nm
            sv = spline_val_0(FS.p_spl2[(7-1)*nm + j], x1, js); bmnc = 1.0e4*sv[1]*bf
            sv = spline_val_0(FS.p_spl2[(8-1)*nm + j], x1, js); bmns = 1.0e4*sv[1]*bf
            sr += bmnc*ffr + bmns*ffi
            si += bmnc*ffi - bmns*ffr
            ffr, ffi = ffr*rotr - ffi*roti, ffr*roti + ffi*rotr
        end
    end
    return sr, si
end

end # module
