# Common orbit trace

`neo_rt_diag.x orbit_trace <runname> <ux> <eta> [nsteps] [mth]` emits
`<runname>_orbit_trace.dat` for one fixed surface and one resonant orbit.
The diagnostic is default-off and does not alter a production torque run.

`neo_rt_diag.x orbit_trace_physical <runname> <ux> <eta> [nsteps] [mth]` emits
the same trace with schema v3 and an additional, explicitly scoped candidate
toroidal trajectory.  The v2 command remains byte-compatible with its prior
field-line-only output; v3 is a separate diagnostic lane and is not consumed by
the production transport path.

`neo_rt_diag.x orbit_trace_hcovar <runname> <ux> <eta> [nsteps] [mth]` emits
schema v4 with the source-bound signed `h_s=B_s/|B|` column in the native cgs
length unit (cm).  It calls the standalone `.bc` geometry hook without feeding
the value back into `do_magfie`; the native operator and torque therefore remain
unchanged.  The command is restricted to `inp_swi=9` and is a diagnostic source
export, not a canonical toroidal embedding or a cross-code acceptance lane.

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

## Thin-orbit toroidal candidate (schema v3)

The optional physical trace integrates the source-faithful thin-orbit rate

```text
d phi_gc / d t = v_parallel hctrvr_phi + v^2 (Omega_tB / v^2) + Omega_tE.
```

Here `v_parallel = sign(hctrvr_theta)*vpar_state`, `hctrvr_phi` is the native
contravariant toroidal field-line component, and the magnetic term is the exact
expression used by `neort_orbit:timestep` for its bounce-averaged
`Omega_tB/v^2`.  `phi_gc(0)=0`; the packet also emits the three addends,
`phi_gc`, and the running `phi_gc/t - Omph` residual so the period identity can
be checked without reconstructing a drift term from a plot.  The diagnostic
rejects the analytic superbanana (`supban`) path because its native `Om_ph` is
not the local thin-orbit drift expression integrated here.

This is a candidate physical guiding-centre coordinate in the thin-orbit model,
not the canonical angle itself.  It still omits the periodic `Delta_phi` chart
term, the canonical origin `phi_H`, and finite-orbit/non-axisymmetric
corrections.  A v3 period identity therefore verifies only the local source
decomposition; it does not close the MARS reflection, phase, coarea, covector,
or work map and cannot authorize a sign, gain, harmonic relabel, or torque
patch.

## Radial covariant export (schema v4)

The v4 trace retains the field-line phase and orientation columns from schema v2
and adds the source-bound finite-orbit inputs after `bmod`.  The radial source
expression is

```text
B_s = (psi'/R)*(R_s*R_theta + Z_s*Z_theta)/J_pol
      + F*(2*pi/nper)*d(v_shift)/ds,
J_pol = R_s*Z_theta - R_theta*Z_s,
F = Bphcov = ItoB*(Jpol/nper).
```

`Jpol/nper` is the toroidal covariant field function (`Bphcov=F`); `Itor` is
the poloidal covariant component (`Bthcov`) and is not substituted into the
shift term.  The geometry and shift derivatives come from the same Sormann
natural radial splines and Fourier convention as the standalone field reader.
The exported component is a signed source quantity with length units.  Its
potential derivative uses the same gauge identity as NEO-RT's chartmap path,

```text
A_phi(s)-A_phi(s_ref) = psi_pr*integral(iota(s),s_ref,s) ds,
A_phi'(s) = psi_pr*iota(s),
delta_phi_H = -c*mi*v_parallel*h_s/(qi*A_phi').
```

The trace emits `vpar_physical`, `A_phi'`, the signed `delta_phi_H`, and the
toroidal-harmonic phase increment `mph*delta_phi_H`.  It also emits the
mechanical canonical-momentum component `mi*v_parallel*h_phi` and factors the
exact signed mode multiplier in native `Tphi_int` into
`toroidal_covector_native=sign(psi_pr*q*sign_theta)*mph` and
`drive_mode_native=mph`.  Their product is the native signed `mph**2` factor.

These are producer values, not a completed canonical map.  The local zero-FOW
state does not contain the periodic action-angle function `Delta_phi` or the
radial derivatives at fixed canonical actions needed for a complete covector.
The trace also lacks the reflected MARS half-bounce and common phase gauge.
Those are exact upstream obstructions to a cross-code phase, coarea, covector,
or work admission; none of the v4 fields is fed back into the native operator.
