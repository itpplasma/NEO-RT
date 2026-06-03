# Cross-checks the Julia spline against the Fortran itpplasma spline reference
# (spline_ref dump on stdin: COEFF and VAL lines), to rtol 1e-12.
include("../src/NeoRT.jl")
using .NeoRT.Spline

const NP = 9

close_enough(a, b) = abs(a - b) <= 1e-14 + 1e-12 * abs(b)

function main()
    x = [Float64(i) * 0.37 for i in 0:NP-1]
    y = [sin(1.3 * xi) + 0.5 * xi for xi in x]
    coeff = spline_coeff(x, y)

    toks = split(read(stdin, String))
    idx = 1
    checks = 0; fails = 0
    for _ in 1:(NP-1)
        @assert toks[idx] == "COEFF"; idx += 1
        row = parse(Int, toks[idx]) + 1; idx += 1
        for k in 1:5
            refv = parse(Float64, toks[idx]); idx += 1
            checks += 1
            if !close_enough(coeff[row, k], refv)
                println(stderr, "COEFF[$row][$k] J=$(coeff[row,k]) ref=$refv")
                fails += 1
            end
        end
    end
    jstart = Ref(1)
    while idx <= length(toks)
        @assert toks[idx] == "VAL"; idx += 1
        xe = parse(Float64, toks[idx]); idx += 1
        ref = (parse(Float64, toks[idx]), parse(Float64, toks[idx+1]), parse(Float64, toks[idx+2])); idx += 3
        out = spline_val_0(coeff, xe, jstart)
        for k in 1:3
            checks += 1
            if !close_enough(out[k], ref[k])
                println(stderr, "VAL(x=$xe)[$k] J=$(out[k]) ref=$(ref[k])")
                fails += 1
            end
        end
    end
    println("spline: $checks checks, $fails failures")
    exit(fails == 0 ? 0 : 1)
end

main()
