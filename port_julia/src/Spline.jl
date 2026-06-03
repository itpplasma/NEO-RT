"""
Cubic spline per Sormann, port of itpplasma/spline (spline.f90). Julia is
1-based like Fortran, so this mirrors the Fortran closely. Uses LAPACK dptsv
(same routine as the C/Rust ports) for the tridiagonal solve.
"""
module Spline

export spline_coeff, spline_val_0

function spline_coeff(x::Vector{Float64}, y::Vector{Float64})
    n = length(x)
    h = x[2:end] .- x[1:n-1]
    r = y[2:end] .- y[1:n-1]
    dl = h[2:n-2]                       # length n-3
    d = 2.0 .* (h[1:n-2] .+ h[2:end])   # length n-2
    c = 3.0 .* (r[2:end] ./ h[2:end] .- r[1:n-2] ./ h[1:n-2])  # length n-2

    nn = Int32(n - 2); nrhs = Int32(1); ldb = Int32(n - 2); info = Ref{Int32}(0)
    ccall((:dptsv_, "liblapack"), Cvoid,
          (Ref{Int32}, Ref{Int32}, Ptr{Float64}, Ptr{Float64}, Ptr{Float64}, Ref{Int32}, Ref{Int32}),
          nn, nrhs, d, dl, c, ldb, info)

    coeff = zeros(Float64, n - 1, 5)
    coeff[:, 1] .= x[1:n-1]
    coeff[:, 2] .= y[1:n-1]
    coeff[1, 3] = r[1] / h[1] - h[1] / 3.0 * c[1]
    coeff[2:n-2, 3] .= r[2:n-2] ./ h[2:n-2] .- h[2:n-2] ./ 3.0 .* (c[2:n-2] .+ 2.0 .* c[1:n-3])
    coeff[n-1, 3] = r[n-1] / h[n-1] - h[n-1] / 3.0 * (2.0 * c[n-2])
    coeff[1, 4] = 0.0
    coeff[2:end, 4] .= c
    coeff[1, 5] = 1.0 / (3.0 * h[1]) * c[1]
    coeff[2:n-2, 5] .= 1.0 ./ (3.0 .* h[2:n-2]) .* (c[2:n-2] .- c[1:n-3])
    coeff[n-1, 5] = 1.0 / (3.0 * h[n-1]) * (-c[n-2])
    return coeff
end

"""Evaluate spline at x; returns (value, 1st deriv, 2nd deriv). jstart is a Ref
carrying the saved interval index (value-independent; x clamped in range)."""
function spline_val_0(coeff::Matrix{Float64}, x::Float64, jstart::Ref{Int})
    n = size(coeff, 1) + 1
    j = jstart[]
    if j < 1 || j >= n
        j = 1
    end
    if n <= 1
        return (0.0, 0.0, 0.0)
    end
    if x < coeff[j, 1]
        while j > 1
            if x < coeff[j, 1]
                j -= 1
            else
                break
            end
        end
    else
        while j < n - 1
            if x >= coeff[j+1, 1]
                j += 1
            else
                break
            end
        end
    end
    jstart[] = j
    z = x - coeff[j, 1]
    v0 = ((coeff[j, 5] * z + coeff[j, 4]) * z + coeff[j, 3]) * z + coeff[j, 2]
    v1 = (3.0 * coeff[j, 5] * z + 2.0 * coeff[j, 4]) * z + coeff[j, 3]
    v2 = 6.0 * coeff[j, 5] * z + 2.0 * coeff[j, 4]
    return (v0, v1, v2)
end

end # module
