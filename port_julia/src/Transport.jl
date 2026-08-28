"""Velocity-space transport integration, port of neort_transport (transport.f90)."""
module Transport

using ..Field
using ..Field: FS
using ..Driftorbit
using ..Driftorbit: DS
using ..Profiles
using ..Profiles: PS
using ..Orbit
using ..Freq
using ..Resonance

export compute_transport_integral

const PI = pi
const SIGN_THETA = -1.0
const OMTH = Ref(0.0)

function d11int(ux, taub, hmn2)
    PI^1.5 * Float64(DS.mph)^2 * Profiles.C^2 * FS.q * PS.vth /
        (PS.qi^2 * DS.dvds * abs(FS.psi_pr)) * ux^3 * exp(-ux^2) * taub * hmn2
end
d12int(ux, taub, hmn2) = d11int(ux, taub, hmn2) * ux^2
function tphi_int(ux, taub, hmn2)
    sgn = FS.psi_pr * FS.q * SIGN_THETA >= 0 ? 1.0 : -1.0
    sgn * PI^1.5 * Float64(DS.mph)^2 * PS.ni1 * Profiles.C * PS.vth / PS.qi *
        ux^3 * exp(-ux^2) * taub * hmn2 * (PS.a1 + PS.a2 * ux^2)
end

function timestep_transport(v, eta, neq, t, y, ydot)
    bmod, _sqrtg, bder, _hcov, hctrvr = do_magfie((get_s(), 0.0, y[1]))
    poloidal_velocity(v, eta, bmod, hctrvr[3], bder[3], y[2], ydot)
    q = FS.q; omth = OMTH[]
    if DS.pertfile
        ar, ai = do_magfie_pert_amp((get_s(), 0.0, y[1]))
        er = DS.epsmn * ar / bmod; ei = DS.epsmn * ai / bmod
    else
        er = DS.epsmn * cos(DS.m0 * y[1]); ei = DS.epsmn * sin(DS.m0 * y[1])
    end
    ph = eta > DS.etatp ? q * DS.mph * y[1] - DS.mth * t * omth :
                          q * DS.mph * y[1] - (DS.mth + q * DS.mph) * t * omth
    cr = cos(ph); ci = sin(ph)
    amp = 2.0 - eta * bmod
    ydot[3] = amp * (er * cr - ei * ci)
    ydot[4] = amp * (er * ci + ei * cr)
    if DS.nonlin
        ydot[5] = 1.0 / bmod; ydot[6] = bmod
    else
        ydot[5] = 0.0; ydot[6] = 0.0
    end
    ydot[7] = 0.0
end

function compute_transport_integral(vmin, vmax, vsteps)
    vth = PS.vth
    roots = zeros(Driftorbit.NLEV * 3)
    d = [0.0, 0.0]; t = 0.0
    du = (vmax - vmin) / (vsteps * vth)
    ux = vmin / vth + du / 2.0
    om_te = PS.om_te; mi = PS.mi
    for _ku in 1:vsteps
        v = ux * vth
        nroots = driftorbit_coarse(v, DS.etamin, DS.etamax, roots, Driftorbit.NLEV)
        for kr in 0:nroots-1
            eta, dres = driftorbit_root(v, 1.0e-8 * abs(om_te), roots[kr*3+1], roots[kr*3+2])
            omth, = om_th(v, eta)
            OMTH[] = omth
            taub = 2.0 * PI / abs(omth)
            bavg, _ = bounce_fast_ext(v, eta, taub, timestep_transport)
            hmn2 = (bavg[3]^2 + bavg[4]^2) * (mi * (ux * vth)^2 / 2.0)^2
            DS.nonlin && error("nonlinear path not ported (outside golden gate)")
            atten = 1.0
            d[1] += du * d11int(ux, taub, hmn2) / abs(dres) * atten
            d[2] += du * d12int(ux, taub, hmn2) / abs(dres) * atten
            if DS.comptorque
                t += du * tphi_int(ux, taub, hmn2) / abs(dres) * atten
            end
        end
        ux += du
    end
    r0 = FS.r0; iota = FS.iota; a = FS.a
    d_plateau = PI * vth^3 / (16.0 * r0 * iota * (PS.qi * DS.b0 / (mi * Profiles.C))^2)
    dsdreff = 2.0 / a * sqrt(get_s())
    d[1] = dsdreff^(-2) * d[1] / d_plateau
    d[2] = dsdreff^(-2) * d[2] / d_plateau
    return d, t
end

end # module
