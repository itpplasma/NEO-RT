"""Canonical-frequency splines + evaluation, port of neort_freq (freq.f90).
Built on the verified orbit + spline layers."""
module Freq

using ..Spline
using ..Field
using ..Driftorbit
using ..Driftorbit: DS
using ..Profiles
using ..Profiles: PS
using ..Orbit

export init_canon_freq_trapped_spline, init_canon_freq_passing_spline, om_th, om_tb, om_ph, d_om_ds

const NETASPL = 100
const PI = pi

mutable struct FreqState
    omth_spl::Matrix{Float64}; omtb_spl::Matrix{Float64}
    omth_pass_spl::Matrix{Float64}; omtb_pass_spl::Matrix{Float64}
    k_taub_p::Float64; d_taub_p::Float64; k_taub_t::Float64; d_taub_t::Float64
    k_omtb_p::Float64; d_omtb_p::Float64; k_omtb_t::Float64; d_omtb_t::Float64
    js::Vector{Base.RefValue{Int}}
end
FreqState() = FreqState(zeros(0,0), zeros(0,0), zeros(0,0), zeros(0,0),
    0,0,0,0, 0,0,0,0, [Ref(1), Ref(1), Ref(1), Ref(1)])
const FQ = FreqState()

function init_canon_freq_trapped_spline()
    n = NETASPL; v = PS.vth; etatp = DS.etatp
    etamin = (1.0 + Driftorbit.EPST) * etatp
    etamax = etatp + (DS.etadt - etatp) * (1.0 - Driftorbit.EPSST_SPL)
    DS.etamin = etamin; DS.etamax = etamax
    magdrift = DS.magdrift
    b = log(Driftorbit.EPST_SPL)
    aa = 1.0 / (n - 1.0) * (log(etamax / etamin - 1.0) - b)
    etarange = zeros(n); omtb_v = zeros(n); omth_v = zeros(n)
    taub0 = 0.0; taub1 = 0.0; leta0 = 0.0; leta1 = 0.0; otb0 = 0.0; otb1 = 0.0
    est = 0.0
    for k in (n-1):-1:0
        eta = etamin * (1.0 + exp(aa * k + b))
        etarange[k+1] = eta
        est = k == n - 1 ? bounce_time(v, eta, 0.0, false) : bounce_time(v, eta, est, true)
        taub = est
        bavg = bounce_fast(v, eta, taub)
        magdrift && (omtb_v[k+1] = bavg[3])
        omth_v[k+1] = 2.0 * PI / (v * taub)
        if k == 0
            leta0 = log(eta - etatp); taub0 = v * taub
            magdrift && (otb0 = omtb_v[k+1] / omth_v[k+1])
        end
        if k == 1
            leta1 = log(eta - etatp); taub1 = v * taub
            magdrift && (otb1 = omtb_v[k+1] / omth_v[k+1])
        end
    end
    FQ.k_taub_t = (taub1 - taub0) / (leta1 - leta0)
    FQ.d_taub_t = taub0 - FQ.k_taub_t * leta0
    FQ.omth_spl = spline_coeff(etarange, omth_v)
    if magdrift
        FQ.k_omtb_t = (otb1 - otb0) / (leta1 - leta0)
        FQ.d_omtb_t = otb0 - FQ.k_omtb_t * leta0
        FQ.omtb_spl = spline_coeff(etarange, omtb_v)
    end
end

function init_canon_freq_passing_spline()
    n = NETASPL; v = PS.vth; etatp = DS.etatp
    etamin = etatp * Driftorbit.EPSSP_SPL
    etamax = etatp
    DS.etamin = etamin; DS.etamax = etamax
    magdrift = DS.magdrift
    b = log((etamax - etamin) / etamax)
    aa = 1.0 / (n - 1.0) * (log(Driftorbit.EPSP_SPL) - b)
    etarange = zeros(n); omtb_v = zeros(n); omth_v = zeros(n)
    taub0 = 0.0; taub1 = 0.0; leta0 = 0.0; leta1 = 0.0; otb0 = 0.0; otb1 = 0.0
    est = 0.0
    for k in (n-1):-1:0
        eta = etamax * (1.0 - exp(aa * k + b))
        etarange[k+1] = eta
        est = k == n - 1 ? bounce_time(v, eta, 0.0, false) : bounce_time(v, eta, est, true)
        taub = est
        bavg = bounce_fast(v, eta, taub)
        magdrift && (omtb_v[k+1] = bavg[3])
        omth_v[k+1] = 2.0 * PI / (v * taub)
        if k == n - 2
            leta0 = log(etatp - eta); taub0 = v * taub
            magdrift && (otb0 = omtb_v[k+1] / omth_v[k+1])
        end
        if k == n - 1
            leta1 = log(etatp - eta); taub1 = v * taub
            magdrift && (otb1 = omtb_v[k+1] / omth_v[k+1])
        end
    end
    FQ.k_taub_p = (taub1 - taub0) / (leta1 - leta0)
    FQ.d_taub_p = taub0 - FQ.k_taub_p * leta0
    FQ.omth_pass_spl = spline_coeff(etarange, omth_v)
    if magdrift
        FQ.k_omtb_p = (otb1 - otb0) / (leta1 - leta0)
        FQ.d_omtb_p = otb0 - FQ.k_omtb_p * leta0
        FQ.omtb_pass_spl = spline_coeff(etarange, omtb_v)
    end
end

function om_th(v, eta)
    etatp = DS.etatp
    if eta > etatp
        if eta > etatp * (1.0 + Driftorbit.EPST_SPL)
            sv = spline_val_0(FQ.omth_spl, eta, FQ.js[1])
        else
            v0 = 2.0 * PI / (FQ.k_taub_t * log(eta - etatp) + FQ.d_taub_t)
            sv = (v0, -v0^2 / (2.0 * PI) * FQ.k_taub_t / (eta - etatp), 0.0)
        end
    elseif eta < etatp * (1.0 - Driftorbit.EPSP_SPL)
        sv = spline_val_0(FQ.omth_pass_spl, eta, FQ.js[3])
    else
        v0 = 2.0 * PI / (FQ.k_taub_p * log(etatp - eta) + FQ.d_taub_p)
        sv = (v0, -v0^2 / (2.0 * PI) * FQ.k_taub_p / (eta - etatp), 0.0)
    end
    return DS.sign_vpar * sv[1] * v, DS.sign_vpar * sv[1], DS.sign_vpar * sv[2] * v
end

function om_tb(v, eta)
    etatp = DS.etatp
    if eta > etatp
        if eta > etatp * (1.0 + Driftorbit.EPST_SPL)
            sv = spline_val_0(FQ.omtb_spl, eta, FQ.js[2])
            s1, s2 = sv[1], sv[2]
        else
            omth, _, deth = om_th(v, eta)
            k = FQ.k_omtb_t; dd = FQ.d_omtb_t
            s1 = DS.sign_vpar * (k * log(eta - etatp) + dd) * omth / v
            s2 = DS.sign_vpar * (omth / v * k / (eta - etatp) + deth / v * (k * log(eta - etatp) + dd))
        end
    elseif eta < etatp * (1.0 - Driftorbit.EPSP_SPL)
        sv = spline_val_0(FQ.omtb_pass_spl, eta, FQ.js[4])
        s1, s2 = sv[1], sv[2]
    else
        omth, _, deth = om_th(v, eta)
        k = FQ.k_omtb_p; dd = FQ.d_omtb_p
        s1 = DS.sign_vpar * (k * log(etatp - eta) + dd) * omth / v
        s2 = DS.sign_vpar * (omth / v * k / (eta - etatp) + deth / v * (k * log(etatp - eta) + dd))
    end
    return s1 * v * v, 2.0 * s1 * v, s2 * v * v
end

function om_ph(v, eta)
    iota = Field.FS.iota
    if eta > DS.etatp
        omph = PS.om_te; dv = 0.0; de = 0.0
        if DS.magdrift
            otb, dvtb, detb = om_tb(v, eta)
            omph += otb; dv += dvtb; de += detb
        end
        return omph, dv, de
    else
        omth, dvth, deth = om_th(v, eta)
        omph = PS.om_te + omth / iota; dv = dvth / iota; de = deth / iota
        if DS.magdrift
            otb, dvtb, detb = om_tb(v, eta)
            omph += otb; dv += dvtb; de += detb
        end
        return omph, dv, de
    end
end

function d_om_ds(v, eta, taub_estimate)
    ds = 2.0e-8; s0 = get_s(); iota = Field.FS.iota
    set_s(s0 - ds / 2.0)
    taub = bounce_time(v, eta, taub_estimate, true)
    bavg = bounce_fast(v, eta, taub)
    omth = DS.sign_vpar_htheta * 2.0 * PI / taub
    omph_noe = DS.magdrift ?
        (eta > DS.etatp ? bavg[3] * v * v : bavg[3] * v * v + omth / iota) :
        (eta > DS.etatp ? 0.0 : omth / iota)
    set_s(s0 + ds / 2.0)
    taub2 = bounce_time(v, eta, taub_estimate, true)
    bavg2 = bounce_fast(v, eta, taub2)
    domthds = DS.sign_vpar_htheta * (2.0 * PI / taub2 - DS.sign_vpar_htheta * omth) / ds
    domphds = DS.magdrift ?
        (eta > DS.etatp ? PS.dom_teds + (bavg2[3] * v * v - omph_noe) / ds :
            PS.dom_teds + (bavg2[3] * v * v + (2.0 * PI / taub2) / iota - omph_noe) / ds) :
        (eta > DS.etatp ? PS.dom_teds : PS.dom_teds + ((2.0 * PI / taub2) / iota - omph_noe) / ds)
    set_s(s0)
    return domthds, domphds
end

end # module
