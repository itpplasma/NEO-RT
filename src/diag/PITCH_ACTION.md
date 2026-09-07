# Off-root pitch and action export

The drift-kinetic code NEO-RT normally integrates over velocity while locating
resonant pitch roots. This diagnostic exports the off-root data needed to
check the same linear, zero-orbit-width integral in the opposite order. It
does not change the production transport path or compute an alternative
torque integral itself.

From a prepared input directory, run through the build driver:

```sh
fo exec --cwd /path/to/inputs neo_rt_diag.x pitch_action case points.in
```

The input directory must contain `case.in`, `in_file`, `in_file_pert`,
`plasma.in` and `profile.in`, with the same identities as the reference lane.
A campaign must hash the actual executable and linked libraries before
execution and retain all input hashes; a source revision alone is insufficient.

## Point file and model guard

Each nonblank, non-comment line has exactly five fields:

```text
s_tor  branch  mth  eta  ux
```

`branch` is 1 for co-passing, 2 for counter-passing, and 3 for trapped.
`eta` is inverse magnetic field, in inverse gauss for the standalone input
path; `ux=v/vth`. Points may repeat a pitch at different speeds, and may
contain an adaptive pitch mesh. Lines beginning with `#` are comments.
Harmonic labels retain their native signs. Points must lie within the
executed branch's pitch and velocity support.

The diagnostic requires `nonlin=false`, `supban=false` and `comptorque=true`.
Passing points are refused when `nopassing=true`. It uses
`neort_setup_at_s` with the ordinary positive-sign frequency initialization,
then changes only the branch's evaluation sign. It verifies `q*iota=1`
within 1e-10 and records the magnetic-drift and shear switches.

## Output and independent checks

`case_pitch_action.dat` names every column in its schema header. It records
the complex normalized bounce action, bounce time, native `Hmn2` and
`Tphi_int`, drive coefficients, toroidal covector signs, branch endpoints,
direct frequencies and derivatives, and quadratic resonance coefficients.
The production callback reads a thread-local orbit frequency: the diagnostic
sets that exact module variable before calling `bounce_fast`. A local copy
would not provide an equivalent action. Solver status must be 2 and every
written value finite.

For normalized speed `u`, the ordinary resonance is
`g=a*u^2+b*u+c`. The coefficients come from the native unit-speed magnetic
drift and poloidal frequency, with the passing `mph/iota` transit term.
Both direct `g,dg/du` and their polynomial evaluations are written so a
consumer can reject a mismatch. No agreement is inferred from the export
itself. Repeated-pitch measurements at two speeds must independently check
the complex bounce amplitude, `u*taub`, `Hmn2/u^4`, and native weight before
using continuum speed scaling in a quadrature.

Tangencies, root births, branch endpoints and the remaining pitch integral
need a separate event-partitioned numerical method. This point exporter
establishes none of those gates and does not authorize a sign, phase or gain
fit to MARS or another code.
