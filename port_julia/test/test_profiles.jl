include("../src/NeoRT.jl")
using .NeoRT.Field
using .NeoRT.Profiles
using .NeoRT.Collis

close_enough(a,b) = abs(a-b) <= 1e-13 + 1e-12*abs(b)

function main()
    s = 0.5
    NeoRT.Field.FS.inp_swi = 9
    NeoRT.Field.FS.bfac = 1.0
    do_magfie_init(ARGS[1])
    _ = do_magfie((s, 0.0, 0.0))
    r0 = NeoRT.Field.FS.r0; psi_pr = NeoRT.Field.FS.psi_pr; q = NeoRT.Field.FS.q
    read_and_init_plasma_input(ARGS[2], s)
    read_and_init_profile_input(ARGS[3], s, r0, 1.0, 1.0)
    init_thermodynamic_forces(psi_pr, q)
    p = NeoRT.Profiles.PS; cs = NeoRT.Collis.CS
    got = [p.vth, p.m_t, p.om_te, p.a1, p.a2, p.ni1, p.ni2, p.ti1, p.ti2, p.te, p.qi, p.mi,
           cs.efcolf[1], cs.efcolf[2], cs.efcolf[3], cs.velrat[1], cs.velrat[2], cs.velrat[3],
           cs.enrat[1], cs.enrat[2], cs.enrat[3]]
    names = ["vth","M_t","Om_tE","A1","A2","ni1","ni2","Ti1","Ti2","Te","qi","mi",
             "efcolf1","efcolf2","efcolf3","velrat1","velrat2","velrat3","enrat1","enrat2","enrat3"]
    refv = parse.(Float64, split(read(stdin, String)))
    fails = 0
    for k in 1:21
        if !close_enough(got[k], refv[k])
            println(stderr, "$(names[k]) J=$(got[k]) ref=$(refv[k])"); fails += 1
        end
    end
    println("profiles: 21 checks, $fails failures")
    exit(fails == 0 ? 0 : 1)
end
main()
