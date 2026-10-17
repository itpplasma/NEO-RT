# Equal-temperature test-particle Fokker–Planck response

Set `tp_fp_collisions=.true.` in `&params` to solve pitch and energy diffusion
against every Maxwellian background in an explicit species mixture. D, T and
electrons use the same `species_file` interface and respond independently in
one invocation. The background temperatures must be exactly equal at the
evaluated surface. Affine profiles with equal reference temperatures can
therefore still be rejected away from their reference surface.

```fortran
&params
 s=.5, species_file='species.in', inp_swi=8
 epsmn=.001, m0=-6, mph=3
 comptorque=.true., magdrift=.false., noshear=.true.
 tp_fp_collisions=.true.
 fp_rate_scale=1, fp_speed_steps=16, fp_pitch_steps=16, fp_umax=4
 output_format='both'
/
```

Use the explicit D/T/e mixture format in [species.md](species.md), with equal
temperatures. `fp_rate_scale` multiplies both physical diffusion coefficients
and must be finite and strictly positive. `fp_speed_steps` and `fp_pitch_steps`
are Q1 cell counts in speed and on each pitch branch; both must be at least 2.
`fp_umax` is the finite positive upper speed in units of the responding
species' thermal speed. These controls default to the values shown. The mode
defaults to disabled and cannot be combined with Krook or Lorentz activation.

This implementation retains the [Lorentz mode](lorentz.md) geometry contract:
standalone format 8 with exactly the `(m,n)=(0,0),(1,0)` field modes, one circular
cosine magnetic well, constant safety factor and integer `n*q`. Linear thin
orbits, both passing signs and the trapped population are required. Magnetic
shear, shaped wells, finite radial orbit width, oscillating drift phase,
nonlinear attenuation and trapped-only operation are rejected. A common magnetic
drift model must apply to all populations. Harmonic output always has the two
common-phase records: zero for trapped particles and `-n*q` for both passing
signs. `mth_max_abs` does not alter this model.

Native orbit sampling retains the existing integration tolerances. The focused
circular oracle bounds magnetic-drift samples against converged independent
angle quadrature to `1e-6` relative; its largest measured error is `4.17e-7`.
This sampling bound is separate from response mesh and tail convergence.

## Physical operator and velocity response

For test species `a` and Maxwellian background `b`, the local operator is

```text
C_ab[f] = nuD_ab/2 * d_xi[(1-xi^2)*d_xi f]
        + v^-2*d_v{v^2*Dvv_ab*[d_v f + ma*v*f/Tb]}.
```

Temperature is an energy in this formula. The independently evaluated pair
coefficients are the Maxwellian deflection frequency `nuD_ab(v)` and radial
velocity diffusion `Dvv_ab(v)`; energy diffusion is not inferred from the
deflection frequency. Every explicit background contributes, including self,
ion–ion, ion–electron and electron–electron pairs, with the supplied common
Coulomb logarithm. `Dvv(0)` has its finite physical limit; the singular
deflection frequency is sampled only at positive interior quadrature speeds.
There is no collision-frequency floor.

Write the nonadiabatic perturbation as `h=F0a*g`, where
`g=(delta_f-delta_H*dF0/dH)/F0a`. At equal temperature the Maxwellian drag is
contained in the conservative energy flux `F0a*Dvv*d_v g`; it has not been
omitted. For `u=v/vth_a`, `eta=(1-xi^2)/B` and a positive native orbit period
`tau`, the weak form uses

```text
W = u^3*exp(-u^2)*tau,
P = 2*eta*W*sum_b(nuD_ab)*<1/B-eta>,
E = W*sum_b(Dvv_ab)/vth_a^2,
integral [P*conj(a_eta)*g_eta + E*conj(a_u)*g_u
          + i*W*R*conj(a)*g] = integral W*conj(a)*S.
```

The orbit sampler supplies the actual complex Hamiltonian `H`, drift frequency,
period and inverse-field average. In this restricted common phase, `u*tau`,
`H/u^2` and the magnetic drift divided by `u^2` are independent of speed. The
solver uses their native orbit values, with the physical electric rotation in
`R=n*OmegaE+u^2*R_magnetic(u=1)`.

The source is `S=-i*beta*H*(A1+A2*u^2)`, with
`beta=n*q/[(qi/c)*sign_theta*psi_pr]` and native thermodynamic forces `A1,A2`.
Energy diffusion couples speeds. The two responses solve sources `H` and
`u^2*H` separately, followed by `g=-i*beta*(A1*z1+A2*z2)`. Multiplying one
response by `A1+A2*u^2` would give a different operator. The collision pullback
`delta C*F0=-C0[H*d_H F0]` is retained by this nonadiabatic source convention.

The tensor Q1 space joins the three pitch intervals at every positive speed.
It imposes continuity and summed collisional flux balance at the trapped/passing
junction, with natural zero flux at the outer pitch endpoints. All pitch and
branch values at `u=0` collapse to one physical degree of freedom. The upper
boundary at `fp_umax` is reflecting; it is a declared finite-domain
approximation, so speed, pitch and upper-tail convergence must be established
for a physical comparison.

With zero detuning, only one global density null mode remains for each test
species. The physical combined source must have zero global `W` mean, and the
response uses one global zero-mean gauge. The auxiliary columns may each be
incompatible even when their physical combination is compatible. In that case
the solver uses the physical combined source directly; an incompatible physical
source fails. Exact zero physical forcing gives exactly zero response and torque.

## Torque and output

The reported electromagnetic torque is direct perturbation work,

```text
dTphi_int/ds = -n*sqrt(pi)*ns*vth_a*abs(psi_pr/q)
              * sum_branches integral W*Im(conj(H)*g) du d_eta.
```

This expression uses native cgs `ns`, `vth_a`, `H` and flux, and returns erg.
It is not a momentum-conservation diagnostic. The `D11` and `D12` output fields
are zero placeholders with `transport/susceptibility_valid=false`; no direct
particle-flux or legacy diffusion-coefficient interpretation is supplied.
Text output marks these placeholders in `<runname>_<species>_collision_model.out`.

For each `/species/<name>`, HDF5 records the five FP controls in `config` and
the consumed `fp_thermal_energy_J`, `fp_vth_m_s` and `fp_charge_over_c` in
`magfie`. This mode uses exact SI elementary charge and the rate provider's
thermal speed consistently. In native formulas, `qi` is chosen so that
`qi/c=Z*e_SI/10`; `qi` is the native electromagnetic coefficient, not a literal
physical statcoulomb value. Species charge metadata remains the physical `Z`.
Other modes retain their historical constants and default behavior.

`transport/fp_source_work_real` and `fp_source_work_imag` contain
`integral W*conj(H)*response_j` for each branch and source column. In Fortran
their shape is `(2,3)`, ordered `(H,u^2*H)` and `(trapped,co,counter)`; h5py sees
shape `(3,2)`. `fp_auxiliary_valid=true` identifies the two auxiliary responses.
If `fp_physical_response=true`, column one instead contains physical `g` work,
column two is zero and `fp_auxiliary_valid=false`. The receiving observable
remains the original Hamiltonian in both cases. Mixture torque sums physical
species torques; species diffusion placeholders are never summed as fluxes.

## Model limits

This is a test-particle operator with maintained Maxwellian backgrounds. It
conserves the test-particle number and has nonnegative diffusion dissipation;
it does not conserve test-species momentum or energy. Background perturbations,
field-particle terms and reciprocal interspecies momentum and energy exchange
are absent. Summing background rates therefore does not make this a coupled,
conserving multispecies Fokker–Planck calculation.

Orbit averaging still requires collisions slow compared with bounce/transit
motion. This ordering fails near the trapped/passing boundary and near zero
speed; finite numerical solutions do not qualify arbitrary collisionalities in
those layers. Comparisons with an unaveraged solver must resolve the layer and
tail limitations. The current equal-temperature restriction is implementation
scope; unequal-temperature Maxwellians need additional energy-drift terms and
a declared maintained-background ordering.

The local test-particle coefficients follow
[Abel et al. (2008), section II](https://arxiv.org/html/0808.1300), with the
conservative transformed operator in
[Albert (2020), section 4.9](https://doi.org/10.3217/978-3-85125-746-5).
The joined weak projection is the model implemented here. The unequal-temperature
term and the distinction from a conserving collision operator are discussed in
[Sugama et al. (2019), appendix B](https://arxiv.org/html/1906.07427).
