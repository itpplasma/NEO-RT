The S3 diagnostic uses the same physical ideal-helical perturbation for both
drive routes on `examples/base`: `m0=2`, `epsmn=1e-3`, no gauge shift, and the
Boozer-collapse displacement enabled. The background and toroidal mode come
from `examples/base/driftorbit.in`. All drive values are divided by the kinetic
energy `E = mi*vth**2/2`; they are dimensionless.

The generator writes two complete sampled orbits (257 equally spaced time
points) and a 36-point pitch scan. Passing particles use `m_b=-4`; trapped
particles use `m_b=0`. At each pitch the E x B rotation is chosen to satisfy
the selected resonance. This is a harmonic-equivalence diagnostic, rather than
a torque profile at fixed rotation. The chosen rotation is retained in the CSV.

The pitch fraction is `eta/(2*etatp)` for passing particles and
`0.5 + (eta-etatp)/(2*(etadt-etatp))` for trapped particles. The scan excludes
the zero-pitch endpoint, the trapped-passing boundary at 0.5, and the deeply
trapped endpoint. The orbit CSVs retain time in seconds, angle in radians,
parallel speed in cm/s, and field strength in gauss. Their plotted complex
integrands include the harmonic phase used by `line_bounce`.

`diagnose_line_drive.x` must be registered as an executable linked to `neo_rt`
in `test/CMakeLists.txt`. Build and execution belong to `fo`; the wrapper
requires a Slurm or PBS allocation and limits shared-memory libraries to one
thread. It needs a Python environment with NumPy and Matplotlib; version 3.6
is sufficient for the plotting script.
Inside an allocation, run from the repository root:

```sh
bash test/line_drive/run_diagnostics.sh "$PWD/artifacts/line-drive"
```

Set `PYTHON` to an existing suitable interpreter if the system Python lacks
the plotting dependencies. The output includes `harmonics.csv`,
`integrand_passing.csv`, `integrand_trapped.csv`, `metadata.txt`, and PNG/PDF
versions of the integrand and harmonic figures. Additional grayscale figures
support visual inspection. No data are smoothed or cropped. Harmonic power is
shown on logarithmic axes; complex harmonic errors use a symmetric-log scale
that is linear below `1e-12`. The Okabe-Ito blue/orange palette and solid/dashed
lines with circle markers distinguish the routes in color and grayscale.

This generator is evidence production, not a substitute for the independent
physical oracles in `test_line_drive`, `test_line_prestudy`, and the end-to-end
line/Boozer torque benchmark. Inspect the actual PNG/PDF outputs and the error
table before publishing or recording stage completion.

The circular pre-study test uses a thin orbit. Its period is checked against
independent conserved-energy quadrature, and its electrostatic toggle against
Cartesian circular geometry and time quadrature. The historical 9–11% potential
misalignment range came from full orbits with finite radial width; it does not
define a pass band for the thin path. Retain those original inputs and results
when comparing with the general raw-field/FOW extension. The thin tests and
these ideal-helical figures do not certify that full-width reproduction.
