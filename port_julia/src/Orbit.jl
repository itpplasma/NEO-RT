"""Orbit bounce integration, port of neort_orbit (orbit.f90), using SUNDIALS
CVODE via ccall -- same library/algorithm as C/Rust. The C callbacks read
module globals for the current (v, eta, ts); valid because the gate is
single-threaded and CVode is synchronous."""
module Orbit

using ..Field
using ..Driftorbit
using ..Driftorbit: DS
using ..Profiles
using ..Profiles: PS
using ..Cvode
using ..Cvode: CV_ADAMS, CV_NORMAL, CV_ROOT_RETURN

export bounce, bounce_time, bounce_fast, bounce_fast_ext, set_s, get_s, NVAR,
    poloidal_velocity, timestep, set_sign_vpar_htheta

const NVAR = 7
const PI = pi
const SIGN_THETA = -1.0

const ORBIT_S = Ref(0.0)
set_s(s::Float64) = (ORBIT_S[] = s)
get_s() = ORBIT_S[]

vpar(v, eta, bmod) = (r = v * sqrt(1.0 - eta * bmod); isnan(r) ? 0.0 : r)
vperp(v, eta, bmod) = (r = v * sqrt(eta * bmod); isnan(r) ? 0.0 : r)

function evaluate_bfield_local()
    bmod, _sqrtg, _bder, _hcov, hctrvr = do_magfie((get_s(), 0.0, DS.th0))
    return bmod, hctrvr[3]
end

function poloidal_velocity(v, eta, bmod, hthctr, hderth, v_par, ydot)
    ydot[1] = v_par * hthctr
    ydot[2] = -v * v * eta / 2.0 * hthctr * hderth * bmod
end

function timestep_poloidal_motion(v, eta, neq, t, y, ydot)
    bmod, _sqrtg, bder, _hcov, hctrvr = do_magfie((get_s(), 0.0, y[1]))
    poloidal_velocity(v, eta, bmod, hctrvr[3], bder[3], y[2], ydot)
end

function timestep(v, eta, neq, t, y, ydot)
    bmod, _sqrtg, bder, _hcov, hctrvr = do_magfie((get_s(), 0.0, y[1]))
    shearterm = DS.noshear ? 0.0 : Field.FS.bphcov * Field.FS.dqds
    q = Field.FS.q
    om_tb_v = PS.mi * Profiles.C * q / (2.0 * PS.qi * SIGN_THETA * Field.FS.psi_pr * bmod) *
              (-(2.0 - eta * bmod) * bmod * bder[1] +
               2.0 * (1.0 - eta * bmod) * hctrvr[3] * (Field.FS.dbthcovds + q * Field.FS.dbphcovds + shearterm))
    ydot[1] = y[2] * hctrvr[3]
    ydot[2] = -0.5 * v * v * eta * hctrvr[3] * bder[3] * bmod
    ydot[3] = om_tb_v
    for i in 4:neq
        ydot[i] = 0.0
    end
end

# --- CVODE callback plumbing (module globals for the current integration) ---
mutable struct CurOde
    v::Float64; eta::Float64; neq::Int; ts::Function
end
const CUR = CurOde(0.0, 0.0, 0, timestep)

function rhs_jl(t::Float64, y::Ptr{Cvoid}, ydot::Ptr{Cvoid}, ud::Ptr{Cvoid})::Cint
    yv = unsafe_wrap(Array, Cvode.nvdata(y), CUR.neq)
    ydv = unsafe_wrap(Array, Cvode.nvdata(ydot), CUR.neq)
    CUR.ts(CUR.v, CUR.eta, CUR.neq, t, yv, ydv)
    return Cint(0)
end

function roots_jl(t::Float64, y::Ptr{Cvoid}, gout::Ptr{Float64}, ud::Ptr{Cvoid})::Cint
    th = unsafe_load(Cvode.nvdata(y), 1)
    svh = DS.sign_vpar_htheta
    unsafe_store!(gout, svh * (th - DS.th0), 1)
    unsafe_store!(gout, svh * (2.0 * PI - (th - DS.th0)), 2)
    return Cint(0)
end

const RHS_C = @cfunction(rhs_jl, Cint, (Float64, Ptr{Cvoid}, Ptr{Cvoid}, Ptr{Cvoid}))
const ROOTS_C = @cfunction(roots_jl, Cint, (Float64, Ptr{Cvoid}, Ptr{Float64}, Ptr{Cvoid}))

function bounce_integral(v, eta, neq, y0, dt, ts)
    CUR.v = v; CUR.eta = eta; CUR.neq = neq; CUR.ts = ts
    etatp = DS.etatp; th0 = DS.th0
    ctx = Cvode.ctx_create()
    y = Cvode.nvnew(neq, ctx)
    yp = Cvode.nvdata(y)
    for i in 1:neq
        unsafe_store!(yp, y0[i], i)
    end
    mem = Cvode.cvcreate(CV_ADAMS, ctx)
    Cvode.cvinit(mem, RHS_C, 0.0, y)
    Cvode.cvtol(mem, 1.0e-9, 1.0e-10)
    Cvode.cvmaxsteps(mem, 50000)
    nls = Cvode.nls_fixedpoint(y, 0, ctx)
    Cvode.cvsetnls(mem, nls)
    Cvode.cvrootinit(mem, 2, ROOTS_C)

    passing = eta < etatp
    ti = Ref(0.0)
    yold = zeros(neq)
    for _k in 2:500
        for i in 1:neq
            yold[i] = unsafe_load(yp, i)
        end
        tout = ti[] + dt
        flag = Cvode.cvode(mem, tout, y, ti, CV_NORMAL)
        if flag < 0
            println(stderr, "CVODE error $flag")
            break
        end
        if flag == CV_ROOT_RETURN && (passing || (yold[1] - th0) < 0.0)
            break
        end
    end
    out = [unsafe_load(yp, i) for i in 1:neq]
    Cvode.nlsfree(nls)
    Cvode.nvdestroy(y)
    Cvode.cvfree(mem)
    Cvode.ctx_free(ctx)
    return ti[], out
end

function integrate_fixed(v, eta, neq, y0, taub, ts)
    CUR.v = v; CUR.eta = eta; CUR.neq = neq; CUR.ts = ts
    ctx = Cvode.ctx_create()
    y = Cvode.nvnew(neq, ctx)
    yp = Cvode.nvdata(y)
    for i in 1:neq
        unsafe_store!(yp, y0[i], i)
    end
    mem = Cvode.cvcreate(CV_ADAMS, ctx)
    Cvode.cvinit(mem, RHS_C, 0.0, y)
    Cvode.cvtol(mem, 1.0e-9, 1.0e-10)
    Cvode.cvmaxsteps(mem, 50000)
    nls = Cvode.nls_fixedpoint(y, 0, ctx)
    Cvode.cvsetnls(mem, nls)
    tret = Ref(0.0)
    Cvode.cvode(mem, taub, y, tret, CV_NORMAL)
    out = [unsafe_load(yp, i) for i in 1:neq]
    Cvode.nlsfree(nls)
    Cvode.nvdestroy(y)
    Cvode.cvfree(mem)
    Cvode.ctx_free(ctx)
    return out
end

function set_sign_vpar_htheta()
    _bmod, htheta = evaluate_bfield_local()
    svh = (htheta >= 0.0 ? 1.0 : -1.0) * DS.sign_vpar
    DS.sign_vpar_htheta = svh
    return svh
end

function bounce_time(v, eta, taub_estimate, have)
    bmod, _ = evaluate_bfield_local()
    svh = set_sign_vpar_htheta()
    y0 = [DS.th0, svh * vpar(v, eta, bmod)]
    taub = have ? taub_estimate :
           2.0 * PI / abs(vperp(v, eta, bmod) * Field.FS.iota / Field.FS.r0 * sqrt(Field.FS.eps / 2.0))
    ti, _ = bounce_integral(v, eta, 2, y0, taub, timestep_poloidal_motion)
    return ti
end

function bounce_fast_ext(v, eta, taub, ts)
    bmod, _ = evaluate_bfield_local()
    svh = set_sign_vpar_htheta()
    y0 = fill(1.0e-15, NVAR)
    y0[1] = DS.th0
    y0[2] = svh * vpar(v, eta, bmod)
    for i in 3:6
        y0[i] = 0.0
    end
    y = integrate_fixed(v, eta, NVAR, y0, taub, ts)
    return [y[i] / taub for i in 1:NVAR], 2
end

bounce_fast(v, eta, taub) = bounce_fast_ext(v, eta, taub, timestep)[1]

function bounce(v, eta, taub_estimate, have)
    bmod, _ = evaluate_bfield_local()
    svh = set_sign_vpar_htheta()
    y0 = fill(1.0e-15, NVAR)
    y0[1] = DS.th0
    y0[2] = svh * vpar(v, eta, bmod)
    for i in 3:6
        y0[i] = 0.0
    end
    tb = have ? taub_estimate :
         2.0 * PI / abs(vperp(v, eta, bmod) * Field.FS.iota / Field.FS.r0 * sqrt(Field.FS.eps / 2.0))
    taub, y = bounce_integral(v, eta, NVAR, y0, tb / 5.0, timestep)
    return taub, [y[i] / taub for i in 1:NVAR]
end

end # module
