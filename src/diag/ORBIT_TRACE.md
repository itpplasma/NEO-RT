# Common orbit trace

`neo_rt_diag.x orbit_trace <runname> <ux> <eta> [nsteps] [mth]` emits
`<runname>_orbit_trace.dat` for one fixed surface and one resonant orbit.
The diagnostic is default-off and does not alter a production torque run.

The packet is an ordered, same-callback trace.  Its position columns are
Boozer coordinates `(s_tor, phi, theta)`, with `rho_tor=sqrt(s_tor)`.  The
phase gauge is `t=0` at the configured `th0` and `phi=0`; `orientation` is the
sign of the parallel velocity.  `H_inst_re/im` is evaluated by the same
`evaluate_hamiltonian` routine used by `timestep_transport`, while
`H_action_re/im` is the running integral stored in the transport state.

The packet deliberately records the native NEO-RT operator only.  It does not
claim a MARS-equivalent coarea measure, toroidal covector, finite-width
resonance prescription, or physical cylindrical embedding.  Those quantities
must be supplied by the corresponding producer before a cross-code action
gate can close.
