# Resonant electromagnetic (`j×b`) torque

NEO-RT keeps this response-field torque separate from neoclassical toroidal
viscosity (NTV). The existing Boozer perturbation file contains the scalar
Lagrangian field-strength perturbation needed by the particle Hamiltonian. It
does not contain the vector magnetic perturbation or induced current required
to reconstruct electromagnetic torque.

## Implemented kernels

For physical cylindrical components, `cylindrical_toroidal_torque` evaluates

\[
 \left\langle R(\mathbf j_1\times\mathbf b_1)_\phi\right\rangle
 = {1\over2}\sum_k w_k R_k\,\Re\left[
 j_{Z,k}b_{R,k}^*-j_{R,k}b_{Z,k}^*\right].
\]

The factor \(1/2\) is the time or toroidal average of complex single-mode
amplitudes. The weights must include the desired normalized surface measure.
With current in A/m², magnetic field in T, radius in m, and dimensionless
normalized weights, the result is in N/m².

`mars_surface_torque` evaluates the collocated MARS contravariant harmonic
contraction

\[
 T_{j\times b}={1\over2}\sum_m\Re\left(j_1^*b_2-j_2^*b_1\right).
\]

Its inputs and result retain MARS normalization. `mars_half_mesh_torque`
implements the active MARS `KCTORQ=2` discretization: the \(j_2^*b_1\) product
is averaged from the two adjacent full-mesh surfaces, while \(j_1^*b_2\) is
evaluated on the half mesh. No coordinate or SI conversion is guessed.

`integrate_mars_profile` reproduces MARS's midpoint radial integration,

\[
 T_\mathrm{net}=4\pi^2 C\sum_i T_i(s_{i+1}-s_i),
\]

where \(C=1\) gives native normalization and a caller-supplied `scale` performs
an explicitly documented dimensional conversion.

## Command-line profile integration

The companion executable reads the first column of a MARS `PROFEQ.OUT` as full
mesh edges and the second column of `TORQUEJXB.OUT` as half-mesh density:

```sh
fo exec neo_rt_jxb_profile.x PROFEQ.OUT TORQUEJXB.OUT [SI_SCALE]
```

It rejects profiles unless the full mesh has exactly one more row than the
half mesh.

## Verification and benchmark

The unit tests have independent behavioral oracles:

- direct cylindrical vector cross products;
- the published MARS harmonic contraction;
- invariance under a common complex Fourier phase;
- exact integration of a linear manufactured profile on a nonuniform mesh;
- command-line parsing of MARS-shaped multi-column files.

An archived MARS validation case provides a profile-integration compatibility
benchmark:

```sh
fo exec neo_rt_jxb_profile.x \
  /mnt/storage/codex-mars/mastu-native-reload-validation-20260716/workspace/pair_a/native/PROFEQ.OUT \
  /mnt/storage/codex-mars/mastu-native-reload-validation-20260716/workspace/pair_a/native/TORQUEJXB.OUT
```

NEO-RT returns `-5.2440799808745e-5`; the MARS log prints `-5.24408e-5`.
This checks mesh interpretation, normalization, and file compatibility. It is
not an independent shielding-response physics benchmark because both values
come from the MARS-computed torque profile.

The complete radial profile and cumulative-integral comparison are reproduced
with:

```sh
python python/plot_jxb_mars_comparison.py \
  PROFEQ.OUT TORQUEJXB.OUT reconstructed.out mars_neort_jxb_profile.png
```

The figure uses MARS native `rho_pol=sqrt(psi_pol)` and torque normalization.
It applies no smoothing, fit, sign change, normalization, or radial remapping.

When raw MARS harmonics are available, NEO-RT independently reconstructs the
profile from `BPLASMA.OUT` and `JPLASMA.OUT`:

```sh
fo exec neo_rt_jxb_from_mars.x \
  PROFEQ.OUT BPLASMA.OUT JPLASMA.OUT reconstructed.out \
  5 5 0.999
```

The last three values reproduce the executed MARS `NTORQ`, smoothing-pass
count, and `CTEDGE`. The output columns are radial coordinate, raw harmonic
contraction, and postprocessed torque. For the archived MAST-U Pair A run, the
postprocessed reconstruction matches all 350 printed `TORQUEJXB.OUT` values
with maximum absolute error `4.46e-11` and relative L2 error `1.69e-8`.

The ITER TC24 MARS log reports two native evaluations,
`-2.30520e-7` and `-6.28797e-7`, corresponding to `-124.072 N m` and
`-338.436 N m`. The supplied TC24 archive does not include `TORQUEJXB.OUT` or
the complex current harmonics, so those totals cannot be independently
reintegrated in NEO-RT. They remain reference results, not claimed benchmark
agreement.

## Code and literature survey

- MARS-F/MARS-Q computes the flux-surface-averaged toroidal electromagnetic
  torque from perturbed current and field and exposes it separately as
  `TORQUEJXB.OUT`. Liu et al., *Physics of Plasmas* **19**, 102507 (2012),
  derive the contraction above and show continuum-resonance splitting near
  rational surfaces:
  <https://scientific-publications.ukaea.uk/wp-content/uploads/Published/POPVOL19P102507.pdf>.
- GPEC/DCON represents ideal singular shielding currents and can expose
  rational-surface fields and currents. PENTRC torque is a kinetic
  nonambipolar/NTV channel and must not be treated as `j×b`. GPEC also offers
  resistive DCON/RMATCH inner-layer matching; its output documentation lists
  the distinct singular-field and kinetic-torque products:
  <https://princetonuniversity.github.io/GPEC/outputs.html>.
- KiLCA resolves a kinetic-Maxwell rational layer in cylindrical geometry.
  QL-Balance obtains species electromagnetic forcing from the absorbed power
  and evolves slow profiles. See Heyn et al., *Nuclear Fusion* **48**, 024005
  (2008), <https://doi.org/10.1088/0029-5515/48/2/024005>, and Heyn et al.,
  *Nuclear Fusion* **54**, 064005 (2014),
  <https://doi.org/10.1088/0029-5515/54/6/064005>.
- MEPHIT iterates a toroidal Ampère solve with localized ideal shielding
  currents and is useful for GPEC/KiLCA geometry comparisons, but an ideal
  delta-like shielding construction is not a finite resistive or kinetic
  layer: <https://arxiv.org/abs/2211.12167>.
- Fitzpatrick's response-matrix and nonlinear layer theories express the
  torque through the complex rational-surface response. They are appropriate
  when only layer fluxes and matching indices—not volume-resolved
  \(\mathbf j_1,\mathbf b_1\)—are available:
  <https://arxiv.org/abs/1808.04482>.

Other MHD codes can evaluate the same coordinate-free volume integral if they
export the complex vector current, magnetic field, geometry, and surface
measure. A code name or a scalar resonant-field spectrum alone is not enough.

## Applicability limits

This implementation calculates torque from a supplied response; it does not
solve for shielding currents. Consequently:

- scalar `Delta|B|` input is insufficient;
- vacuum fields without plasma current give zero plasma `j×b` torque;
- ideal singular currents require a specified distributional or layer
  integration convention;
- smoothing, rational-surface masks, and sign flips are never applied
  implicitly;
- all Fourier, coordinate, Jacobian, radial, and unit conventions remain the
  responsibility of the response-field adapter.
