"""Resonance-line root finding, port of neort_resonance (resonance.f90)."""
module Resonance

using ..Driftorbit
using ..Driftorbit: DS
using ..Freq

export driftorbit_coarse, driftorbit_nroot, driftorbit_root

function driftorbit_coarse!(roots, v, eta_min, eta_max, ninterv)
    deta = (eta_max - eta_min) / ninterv
    resold = 0.0; nroots = 0
    for k in 0:ninterv
        eta = eta_min + k * deta
        omth, = om_th(v, eta)
        omph, = om_ph(v, eta)
        res = DS.mth * omth + DS.mph * omph
        if k > 0
            snew = res < 0.0 ? -1.0 : 1.0
            sold = resold < 0.0 ? -1.0 : 1.0
            if snew != sold
                roots[nroots*3 + 1] = eta - deta
                roots[nroots*3 + 2] = eta
                nroots += 1
            end
        end
        resold = res
    end
    return nroots
end

driftorbit_coarse(v, eta_min, eta_max, roots, ninterv) = driftorbit_coarse!(roots, v, eta_min, eta_max, ninterv)

function driftorbit_nroot(v, eta_min, eta_max)
    roots = zeros(Driftorbit.NLEV * 3)
    driftorbit_coarse!(roots, v, eta_min, eta_max, Driftorbit.NLEV)
end

function driftorbit_root(v, tol, eta_min, eta_max)
    mth = Float64(DS.mth); mph = Float64(DS.mph)
    maxit = 100
    etamin2 = eta_min; etamax2 = eta_max; eta = etamin2
    omth, = om_th(v, eta); omph, = om_ph(v, eta)
    resmin = mph * omph + mth * omth
    eta = etamax2
    omth, = om_th(v, eta); omph, = om_ph(v, eta)
    resmax = mph * omph + mth * omth
    slope_pos = resmax - resmin > 0.0
    if driftorbit_nroot(v, etamin2, etamax2) == 0
        _, _, deth = om_th(v, eta); _, _, deph = om_ph(v, eta)
        return eta, mph * deph + mth * deth
    end
    out_eta = eta; out_d = 0.0; state = -2
    for _ in 1:maxit
        omth, _, deth = om_th(v, eta); omph, _, deph = om_ph(v, eta)
        res = mph * omph + mth * omth
        out_eta = eta
        if abs(res) < tol
            state = 1; out_d = mph * deph + mth * deth; break
        elseif (slope_pos && res > 0.0) || (!slope_pos && res < 0.0)
            etamax2 = eta; eta = (eta + etamin2) / 2.0
        else
            etamin2 = eta; eta = (eta + etamax2) / 2.0
        end
    end
    if state < 0
        _, _, deth = om_th(v, out_eta); _, _, deph = om_ph(v, out_eta)
        out_d = mph * deph + mth * deth
    end
    return out_eta, out_d
end

end # module
