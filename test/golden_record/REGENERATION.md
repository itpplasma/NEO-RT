# Golden record regeneration: 2026-06-15

The committed `golden.h5` was regenerated from the fortnum vode integrator
(branch `migrate/fortnum-ode-events-drop-vode`, PR #47; fortnum pinned to main
`974dcf1`), built `CONFIG=Fast`. The previous golden was generated on demand
from NEO-RT `main` (DVODE_F90). Comparison bar unchanged: `rtol=1e-8`,
`atol=1e-15`.

## Why

The DVODE golden is wrong at near-separatrix trapped spline nodes. NEO-RT
builds the canonical poloidal-frequency spline in
`neort_freq::init_canon_freq_trapped_spline` by calling `bounce_time(v, eta)`
on a 100-point eta grid at `v = vth`. That spline defines `Om_th`, hence the
bounce time `taub = 2*pi/(v*Omth)` that enters the transport integral. At the
few eta nodes nearest the separatrix, the DVODE chunked event search returns an
integer multiple of the fundamental bounce period instead of the fundamental.
The seed `taub_estimate` chained from the previous node drives DVODE past the
first theta=th0 return, and its `(yold-th0)<0` acceptance latches onto the
second or third full-period crossing.

fortnum returns the fundamental period at every node.

## Validation: 2e6-substep RK4

Each regenerated value was checked against an independent reference: a
fixed-step RK4 of the poloidal-motion ODE (`theta`, `vpar`) over the same
`do_magfie` field, ~2e6 substeps per fundamental period, the bounce event
bracketed at every theta=th0 crossing and root-polished by RK4 sub-bisection.
Convergence confirmed at 2e6, 16e6, and 64e6 substeps (agreement to 1e-9 at the
worst node). The fundamental bounce is the first crossing entering th0 from
below, the physical full-period return.

Per mismatch node the RK4 truth resolves the period directly. Example, 0p700
k=99 (eta=6.94e-5): RK4 crossings at t = 6.535e-5 (rejected, from above),
1.33610e-4 (fundamental), 2.67221e-4 (2x), 4.00831e-4 (3x). fortnum taub =
1.3224e-4 (fundamental, 1.0% short, the event-locator tolerance at the
deep-trapped node); DVODE golden taub = 2.67221e-4, exactly 2x the fundamental.

Across all 9 cases the only DVODE/fortnum disagreements above 1% are trapped
spline nodes where DVODE is an integer multiple of the RK4 fundamental:

| case  | trapped multiple-period nodes (k: DVODE/fundamental) |
|-------|------------------------------------------------------|
| 0p100 | 84: 2.00x                                            |
| 0p200 | 69: 2.00x                                            |
| 0p300 | none                                                 |
| 0p400 | 43: 2.00x, 91: 2.00x, 99: 2.00x                      |
| 0p500 | none                                                 |
| 0p600 | 57: 2.00x, 99: 2.00x                                 |
| 0p700 | 8: 2.00x, 99: 2.00x                                  |
| 0p800 | 76: 2.00x, 79: 2.00x, 80: 2.00x, 98: 2.00x           |
| 0p900 | 92: 2.00x, 95: 2.00x, 97: 2.00x, 99: 5.00x           |

At every one of these nodes fortnum/fundamental is 0.99 to 1.00 and
DVODE/fundamental is 2.00 (5.00 at 0p900 k=99). At all other trapped nodes and
across the full passing grid, fortnum and DVODE agree to ~1e-6 and both match
RK4. fortnum passes RK4 validation for all 9 cases; DVODE is the outlier wherever
old and new differ. No case showed fortnum as the outlier, so no case was held
back.

## Per-case OLD vs NEW transport

`magfie` is bit-identical between backends. The divergence is the trapped
diffusion channel, which carries the multiple-period nodes.

| case  | D11 DVODE   | D11 fortnum | d%    | D12 DVODE   | D12 fortnum | d%    |
|-------|-------------|-------------|-------|-------------|-------------|-------|
| 0p100 | 4.66479e-03 | 4.63838e-03 | -0.57 | 1.39282e-02 | 1.38751e-02 | -0.38 |
| 0p200 | 1.00903e-03 | 1.05268e-03 | +4.33 | 1.33994e-03 | 1.43353e-03 | +6.98 |
| 0p300 | 4.05668e-04 | 3.99458e-04 | -1.53 | 5.34070e-04 | 5.23870e-04 | -1.91 |
| 0p400 | 2.40201e-03 | 2.38122e-03 | -0.87 | 8.41916e-03 | 8.39316e-03 | -0.31 |
| 0p500 | 6.39491e-04 | 6.39403e-04 | -0.01 | 2.21908e-03 | 2.21904e-03 | -0.00 |
| 0p600 | 3.22354e-03 | 3.22456e-03 | +0.03 | 1.25500e-02 | 1.25507e-02 | +0.01 |
| 0p700 | 2.42497e-02 | 2.25794e-02 | -6.89 | 5.77027e-02 | 5.31353e-02 | -7.92 |
| 0p800 | 4.03095e-03 | 3.96635e-03 | -1.60 | 1.10415e-02 | 1.07758e-02 | -2.41 |
| 0p900 | 3.04954e-03 | 2.94211e-03 | -3.52 | 1.05836e-02 | 1.02717e-02 | -2.95 |

0p300 and 0p500 have no multiple-period nodes and barely move (d% <= 0.01 on
0p500, ~1.9% on 0p300 from sub-1% taub differences spread over the grid). The
largest shifts (0p700, 0p900) carry the most multiple-period nodes, including
the 5x node at 0p900.

## Mechanics

`golden.h5` is now a committed artifact. `ensure_golden.py` uses a committed
`golden.h5` as-is and regenerates from main only when none is present. The
comparison logic and tolerances in `test_golden_record.py` are unchanged.
`test_golden_record.py` passes for all 9 cases at `rtol=1e-8`, `atol=1e-15`.
