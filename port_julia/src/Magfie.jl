"""Flux-surface characterization, port of neort_magfie (magfie.f90). Sets
driftorbit B0/Bmin/Bmax/etatp/etadt/dVds/th0 and field eps."""
module Magfie

using ..Field
using ..Driftorbit: DS

export init_flux_surface_average

const PI = pi

function init_flux_surface_average(s::Float64)
    nth = 1000
    dth = 2.0 * PI / nth
    dvds = 0.0; b0 = 0.0; eps = 0.0
    bmin = -1.0; bmax = 0.0; th0 = 0.0
    for k in 1:nth
        th = -PI + k * 2.0 * PI / nth
        bmod, sqrtg, _bder, _hcov, _hctr = do_magfie((s, 0.0, th))
        dvds += abs(sqrtg) * dth
        b0 += bmod * dth
        eps -= cos(th) * bmod * dth
        if bmin < 0.0 || bmod < bmin
            bmin = bmod; th0 = th
        end
        bmod > bmax && (bmax = bmod)
    end
    dvds *= 2.0 * PI
    b0 /= 2.0 * PI
    eps /= b0 * PI
    Field.FS.eps = eps
    DS.dvds = dvds; DS.b0 = b0; DS.bmin = bmin; DS.bmax = bmax
    DS.etatp = 1.0 / bmax; DS.etadt = 1.0 / bmin; DS.th0 = th0
end

end # module
