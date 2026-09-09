# Common orbit trace

`neo_rt_diag.x orbit_trace <runname> <ux> <eta> [nsteps] [mth]` emits
`<runname>_orbit_trace.dat` for one fixed surface and one resonant orbit.
The diagnostic is default-off and does not alter a production torque run.

The packet is an ordered, same-callback trace.  Its radial and poloidal
columns are evaluated in the Boozer chart with `rho_tor=sqrt(s_tor)`.  The
`phi` column is named `phi_trace` in the metadata because it is the
drift-free field-line phase

```text
phi_trace(t) = q * (theta_orb(t) - th0).
```

It is suitable for the native axisymmetric field and Fourier-amplitude
callbacks, but it is **not** a complete physical guiding-centre toroidal
position.  The trace does not carry the canonical toroidal angle, periodic
`Delta phi` correction, or the secular `Omega_t` precession.  Consequently,
the old shorthand `position_coordinates = Boozer(s_tor,phi,theta)` must not
be used as evidence of a physical MARS/NEO-RT embedding.

The distinction follows the canonical-angle construction in Albert et al.
(2016), Eqs. (angles), (amn), and (resonance):

```text
theta^1 = phi - Delta_phi(theta^2,J)
theta^2 = Omega_theta * tau
theta^3 = phi_H + q*theta^2*delta_tp - q*theta_orb(theta0,tau)
Omega_phi = q*Omega_theta*delta_tp + Omega_t.
```

The `evaluate_hamiltonian` callback intentionally uses the canonical Fourier
coefficient `exp(i*q*mph*theta - i*(mth+q*mph*delta_tp)*Omega_theta*t)`;
adding `Omega_t*t` to that coefficient would change the native operator rather
than repair this diagnostic label.  A physical position/action join requires
an additional producer for `phi_H`, `Delta_phi`, the secular precession and
the associated metric/covector map.

The phase gauge is `t=0` at the configured `th0` and `phi_trace=0`.

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
