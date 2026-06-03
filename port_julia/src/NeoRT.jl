"""NEO-RT Julia port (golden path). Submodules mirror the Fortran/C/Rust
structure; numerics use LAPACK dptsv and SUNDIALS CVODE via ccall for parity
with the C/Rust ports."""
module NeoRT

include("Spline.jl")
include("Collis.jl")
include("Field.jl")
include("Profiles.jl")

end # module
