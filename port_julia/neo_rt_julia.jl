#!/usr/bin/env julia
# NEO-RT Julia port driver: parse <case>.in, run the golden-path computation
# (mirrors neort_compute_at_s + compute_transport), write the 5 output files.
# Single-threaded (deterministic gate). Same CVODE/LAPACK as C/Rust.
include(joinpath(@__DIR__, "src", "NeoRT.jl"))
using .NeoRT.Field, .NeoRT.Magfie, .NeoRT.Profiles, .NeoRT.Orbit
using .NeoRT.Freq, .NeoRT.Transport
using .NeoRT.Driftorbit: DS
using Printf

const PI = pi

function read_config(path)
    cfg = Dict{String,Any}("bfac"=>1.0, "efac"=>1.0, "comptorque"=>true, "magdrift"=>true,
        "nopassing"=>false, "noshear"=>false, "pertfile"=>false, "nonlin"=>false,
        "s"=>0.0, "epsmn"=>0.0, "m0"=>0, "mph"=>0, "inp_swi"=>0, "vsteps"=>0)
    for line in readlines(path)
        line = split(line, '!')[1]
        if occursin('=', line)
            k, v = split(line, '=', limit=2)
            key = strip(lowercase(k)); val = strip(v)
            if key in ("comptorque","magdrift","nopassing","noshear","pertfile","nonlin")
                cfg[key] = startswith(val, "t") || startswith(val, "T") || occursin(".true.", val) || occursin(".TRUE.", val)
            elseif key in ("m0","mph","inp_swi","vsteps")
                cfg[key] = parse(Int, split(val)[1])
            elseif key in ("s","epsmn","bfac","efac")
                cfg[key] = parse(Float64, split(val)[1])
            end
        end
    end
    cfg
end

function main()
    run = ARGS[1]
    cfg = read_config("$run.in")
    NeoRT.Field.FS.inp_swi = cfg["inp_swi"]
    NeoRT.Field.FS.bfac = cfg["bfac"]
    DS.epsmn = cfg["epsmn"]; DS.m0 = cfg["m0"]; DS.efac = cfg["efac"]
    DS.comptorque = cfg["comptorque"]; DS.magdrift = cfg["magdrift"]
    DS.nopassing = cfg["nopassing"]; DS.pertfile = cfg["pertfile"]
    DS.nonlin = cfg["nonlin"]; DS.noshear = cfg["noshear"]
    s = cfg["s"]; vsteps = cfg["vsteps"]

    do_magfie_init("in_file")
    if cfg["pertfile"]
        do_magfie_pert_init("in_file_pert")
        DS.mph = NeoRT.Field.FS.pert_mph
    else
        DS.mph = cfg["mph"]
    end
    read_and_init_plasma_input("plasma.in", s)
    r0 = NeoRT.Field.FS.r0
    read_and_init_profile_input("profile.in", s, r0, cfg["efac"], cfg["bfac"])
    init_flux_surface_average(s)
    set_s(s)
    init_canon_freq_trapped_spline()
    cfg["nopassing"] || init_canon_freq_passing_spline()
    DS.sign_vpar = 1.0
    DS.etamin = (1.0 + NeoRT.Driftorbit.EPST) * DS.etatp
    DS.etamax = (1.0 - NeoRT.Driftorbit.EPST) * DS.etadt
    psi_pr = NeoRT.Field.FS.psi_pr; q = NeoRT.Field.FS.q
    cfg["comptorque"] && init_thermodynamic_forces(psi_pr, q)

    PS = NeoRT.Profiles.PS
    PS.om_te = PS.vth * PS.m_t / r0
    PS.dom_teds = PS.vth * PS.dm_tds / r0 + PS.m_t * PS.dvthds / r0

    set_trapped() = (DS.etamin = (1.0 + NeoRT.Driftorbit.EPST) * DS.etatp; DS.etamax = (1.0 - NeoRT.Driftorbit.EPST) * DS.etadt)
    set_passing() = (DS.etamin = NeoRT.Driftorbit.EPSP * DS.etatp; DS.etamax = (1.0 - NeoRT.Driftorbit.EPSP) * DS.etatp)

    mthmax = ceil(Int, 2.0 * abs(DS.mph * q)); mthmin = -mthmax
    vth = PS.vth
    dco = [0.0, 0.0]; dctr = [0.0, 0.0]; dt = [0.0, 0.0]
    tco = 0.0; tctr = 0.0; tt = 0.0
    harmonics = Vector{Any}()
    vminp = 1.0e-6 * vth; vmaxp = 3.0 * vth
    for j in mthmin:mthmax
        DS.mth = j
        h = Dict{String,Any}("mth"=>j, "dresco"=>[0.0,0.0], "dresctr"=>[0.0,0.0], "drest"=>[0.0,0.0],
            "tresco"=>0.0, "tresctr"=>0.0, "trest"=>0.0)
        if !cfg["nopassing"]
            DS.sign_vpar = 1.0; set_passing()
            dd, t = compute_transport_integral(vminp, vmaxp, vsteps)
            h["dresco"] = dd; h["tresco"] = t; dco .+= dd; tco += t
            DS.sign_vpar = -1.0; set_passing()
            dd, t = compute_transport_integral(vminp, vmaxp, vsteps)
            h["dresctr"] = dd; h["tresctr"] = t; dctr .+= dd; tctr += t
        end
        DS.sign_vpar = 1.0; set_trapped()
        dd, t = compute_transport_integral(vminp, vmaxp, vsteps)
        h["drest"] = dd; h["trest"] = t; dt .+= dd; tt += t
        push!(harmonics, h)
    end

    m_t = PS.m_t; dvds = DS.dvds
    tot1 = dco[1] + dctr[1] + dt[1]; tot2 = dco[2] + dctr[2] + dt[2]
    f(x) = @sprintf("%.17e", x)
    open("$run.out", "w") do io
        println(io, " # M_t D11co D11ctr D11t D11 D12co D12ctr D12t D12")
        println(io, " ", join(f.([m_t, dco[1], dctr[1], dt[1], tot1, dco[2], dctr[2], dt[2], tot2]), " "))
    end
    if cfg["comptorque"]
        open("$(run)_torque.out", "w") do io
            println(io, " # s dVds M_t Tco Tctr Tt")
            println(io, " ", join(f.([s, dvds, m_t, tco, tctr, tt]), " "))
        end
    end
    open("$(run)_integral.out", "w") do io
        open("$(run)_torque_integral.out", "w") do io2
            for h in harmonics
                t1 = h["dresco"][1] + h["dresctr"][1] + h["drest"][1]
                t2 = h["dresco"][2] + h["dresctr"][2] + h["drest"][2]
                println(io, " ", f(m_t), " ", h["mth"], " ",
                    join(f.([h["dresco"][1], h["dresctr"][1], h["drest"][1], t1,
                             h["dresco"][2], h["dresctr"][2], h["drest"][2], t2,
                             vminp/vth, vmaxp/vth, vminp/vth, vmaxp/vth]), " "))
                println(io2, " ", h["mth"], " ", join(f.([h["tresco"], h["tresctr"], h["trest"]]), " "))
            end
        end
    end
    open("$(run)_magfie.out", "w") do io
        nth = 50
        for k in 0:nth-1
            th = -PI + k * (2.0 * PI) / (nth - 1)
            bmod, sqrtg, bder, hcovar, hctrvr = do_magfie((s, 0.0, th))
            if cfg["pertfile"]
                ar, ai = do_magfie_pert_amp((s, 0.0, th)); bnr = DS.epsmn*ar/bmod; bni = DS.epsmn*ai/bmod
            else
                bnr = DS.epsmn*cos(DS.m0*th); bni = DS.epsmn*sin(DS.m0*th)
            end
            er = DS.epsmn*cos(DS.m0*th); ei = DS.epsmn*sin(DS.m0*th)
            println(io, " ", join(f.([th, bmod, sqrtg, bder[1], bder[2], bder[3],
                hcovar[1], hcovar[2], hcovar[3], hctrvr[1], hctrvr[2], hctrvr[3],
                0.0, 0.0, 0.0, bnr, bni, er, ei]), " "))
        end
    end
end

main()
