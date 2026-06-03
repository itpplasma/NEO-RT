include("../src/NeoRT.jl")
using .NeoRT.Field, .NeoRT.Magfie, .NeoRT.Profiles, .NeoRT.Orbit
using .NeoRT.Driftorbit: DS

close_enough(a,b) = abs(a-b) <= 1e-12 + 1e-8*abs(b)

function main()
    s = 0.5
    NeoRT.Field.FS.inp_swi = 9; NeoRT.Field.FS.bfac = 1.0
    do_magfie_init(ARGS[1])
    _ = do_magfie((s,0.0,0.0))
    r0 = NeoRT.Field.FS.r0
    init_flux_surface_average(s)
    read_and_init_plasma_input(ARGS[2], s)
    read_and_init_profile_input(ARGS[3], s, r0, 1.0, 1.0)
    set_s(s)
    DS.sign_vpar = 1.0
    vth = NeoRT.Profiles.PS.vth
    vv = [vth, vth, 1.5*vth]
    ee = [0.5*(DS.etatp+DS.etadt), 0.3*DS.etatp, 0.8*DS.etatp+0.2*DS.etadt]
    refv = parse.(Float64, split(read(stdin, String)))
    names = ["bounce_time","taub","bavg1","bavg2","bavg3","bavg6"]
    fails = 0; checks = 0
    for i in 0:2
        tt = bounce_time(vv[i+1], ee[i+1], 0.0, false)
        taub, bavg = bounce(vv[i+1], ee[i+1], 0.0, false)
        got = [tt, taub, bavg[1], bavg[2], bavg[3], bavg[6]]
        for k in 1:6
            r = refv[i*6+k]; checks += 1
            if !close_enough(got[k], r)
                println(stderr, "case $i $(names[k]) J=$(got[k]) ref=$r rel=$(abs((got[k]-r)/r))"); fails += 1
            end
        end
    end
    println("orbit: $checks checks, $fails failures")
    exit(fails == 0 ? 0 : 1)
end
main()
