"""Collision coefficients, port of collis_alp (collis_nbi.f90). Module state in a
mutable struct. Keeps local ev=1.6022e-12 distinct from util ev, per Fortran."""
module Collis

export loacol_nbi, coleff, CS

const NSORTS = 3

mutable struct CollisState
    efcolf::Vector{Float64}
    velrat::Vector{Float64}
    enrat::Vector{Float64}
end
CollisState() = CollisState(zeros(NSORTS), zeros(NSORTS), zeros(NSORTS))
const CS = CollisState()

function onseff(v::Float64)
    sqp = 1.7724538; cons = 0.75225278
    v2 = v*v; v3 = v2*v
    if v < 0.01
        return cons*(1.0-0.6*v2), cons*(1.0-0.2*v2), 2.0*cons*(1.0-1.2*v2)
    elseif v > 6.0
        return 1.0/v3, (1.0-0.5/v2)/v, -1.0/v3
    else
        ex = exp(-v2)/sqp; er = erf_(v)
        d_p = er/v3 - 2.0*ex/v2
        return d_p, er*(1.0-0.5/v2)/v + ex/v2, 4.0*ex - d_p
    end
end

erf_(x::Float64) = ccall((:erf, "libm"), Float64, (Float64,), x)

function coleff(p::Float64)
    plim = max(p, 1.0e-8)
    dpp = 0.0; dhh = 0.0; fpeff = 0.0
    for i in 1:NSORTS
        d_p, dh, dpd = onseff(p * CS.velrat[i])
        dpp += d_p * CS.efcolf[i]
        dhh += dh * CS.efcolf[i]
        fpeff += (dpd/plim - 2.0*d_p*p*CS.enrat[i]) * CS.efcolf[i]
    end
    dhh /= plim^2
    return dpp, dhh, fpeff
end

function loacol_nbi(amb, am1, am2, zb, z1, z2, densi1, densi2, tempi1, tempi2, tempe, ebeam)
    pi_ = 3.14159265358979
    pmass = 1.6726e-24; emass = 9.1094e-28; e = 4.8032e-10; ev = 1.6022e-12
    eps = eps_f64()
    v0 = sqrt(2.0*ebeam*ev/(amb*pmass))
    vti1 = sqrt(2.0*tempi1*ev/(pmass*am1))
    vti2 = sqrt(2.0*tempi2*ev/(pmass*am2))
    vte = sqrt(2.0*tempe*ev/emass)
    dense = densi1*z1 + densi2*z2
    a1 = sqrt(densi1*z1^2/tempi1)*zb*z1*(amb+am1)/(amb*tempi1+am1*ebeam)
    a2 = sqrt(densi2*z2^2/tempi2)*zb*z2*(amb+am2)/(amb*tempi2+am2*ebeam)
    alami1 = 23.0 - log(max(a1, eps))
    alami2 = 23.0 - log(max(a2, eps))
    alame = 24.0 - log(sqrt(dense)/tempe)
    frecol = 2.0*pi_*dense*e^4*zb^2/((amb*pmass)^2*v0^3)
    frecol /= v0
    CS.enrat = [ebeam/tempi1, ebeam/tempi2, ebeam/tempe]
    CS.velrat = [v0/vti1, v0/vti2, v0/vte]
    CS.efcolf = [frecol*z1^2*alami1*densi1/dense, frecol*z2^2*alami2*densi2/dense, frecol*alame]
    for i in 1:NSORTS
        CS.efcolf[i] *= CS.velrat[i]
    end
    return v0
end

eps_f64() = eps(1.0)

end # module
