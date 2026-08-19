# Cross-checks the Julia field port (do_magfie) against the Fortran do_magfie_mod
# reference (field_ref dump on stdin), to rtol 1e-12.
include("../src/NeoRT.jl")
using .NeoRT.Field

close_enough(a, b) = abs(a - b) <= 1e-13 + 1e-12 * abs(b)

function main()
    in_file = ARGS[1]
    NeoRT.Field.FS.inp_swi = 9
    NeoRT.Field.FS.bfac = 1.0
    do_magfie_init(in_file)

    nums = parse.(Float64, split(read(stdin, String)))
    s = 0.5
    names = ["bmod","sqrtg","bder1","bder3","hcov2","hcov3","hctr2","hctr3"]
    nrow = div(length(nums), 11)
    checks = 0; fails = 0
    for i in 0:nrow-1
        r = nums[i*11+1 : i*11+11]
        theta = r[1]
        bmod, sqrtg, bder, hcovar, hctrvr = do_magfie((s, 0.0, theta))
        got = [bmod, sqrtg, bder[1], bder[3], hcovar[2], hcovar[3], hctrvr[2], hctrvr[3]]
        refv = [r[2], r[3], r[4], r[5], r[6], r[7], r[8], r[9]]
        for k in 1:8
            checks += 1
            if !close_enough(got[k], refv[k])
                println(stderr, "theta=$theta $(names[k]) J=$(got[k]) ref=$(refv[k])")
                fails += 1
            end
        end
    end
    println("field: $checks checks, $fails failures")
    exit(fails == 0 ? 0 : 1)
end

main()
