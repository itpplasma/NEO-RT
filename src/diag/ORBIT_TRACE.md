# Common orbit trace

`neo_rt_diag.x orbit_trace <runname> <ux> <eta> [nsteps] [mth]` emits
`<runname>_orbit_trace.dat` for one fixed surface and one resonant orbit.
The diagnostic is default-off and does not alter a production torque run.

The packet is an ordered, same-callback trace.  Its position columns are
Boozer coordinates `(s_tor, phi, theta)`, with `rho_tor=sqrt(s_tor)`.  The
phase gauge is `t=0` at the configured `th0` and `phi=0`.

Schema v2 makes the orientation map explicit.  `vpar_state` is the signed
second orbit state and `hctrvr_theta` is the native `hctrvr(3)` factor in
`ydot(theta)=vpar_state*hctrvr_theta`.  Therefore
`orientation_vpar=sign(vpar_state*hctrvr_theta)` is the physical parallel
direction, while `orientation_state=sign(vpar_state)` is retained as a
coordinate-state diagnostic.  Either orientation is `0` at an exactly zero
factor.  A v1 packet is legacy evidence: its `orientation` column was only the
state sign despite the old header label and must not be treated as a physical
parallel-direction field.

`H_inst_re/im` is evaluated by the same `evaluate_hamiltonian` routine used by
`timestep_transport`, while `H_action_re/im` is the running integral stored in
the transport state.

The packet deliberately records the native NEO-RT operator only.  It does not
claim a MARS-equivalent coarea measure, toroidal covector, finite-width
resonance prescription, or physical cylindrical embedding.  Those quantities
must be supplied by the corresponding producer before a cross-code action
gate can close.
