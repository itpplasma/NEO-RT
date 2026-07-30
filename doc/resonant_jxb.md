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

## Boozer Ampère kernel

`boozer_harmonic_current` moves the independent calculation one step upstream:
given the three complex covariant components of a perturbation harmonic and
the radial derivatives of its angular components, it evaluates
`curl(B)/mu0` and returns all three contravariant current components.  The
Fourier convention is

```text
B_hat(s) exp(i * (n*phi + m*theta))
```

for coordinates ordered `(s, phi, theta)`.  The signed coordinate Jacobian is
an explicit input, so left- and right-handed charts are not silently folded
together.  The source kernel is generated from the fortsym chart curl; see
`derivations/fortsym/README.md`.  Its independent behavioral test uses the
analytic curl of a general complex cylindrical vector field.

The current chartmap NetCDF contract contains geometry, `A_phi`, `B_theta`,
`B_phi`, and scalar `Bmod`.  The current perturbation `.bc` contract likewise
contains harmonics of scalar `delta|B|`.  Neither contains the three vector
components of `delta B`, and scalar field magnitude alone does not uniquely
determine `curl(delta B)`.  Consequently this kernel is ready for a vector
perturbation input, but it does not pretend to infer shielding current from
the existing scalar file.  The shielding-response layer must either supply a
documented vector reconstruction or extend the NetCDF contract with vector
components or a vector potential.

The new vector perturbation contract takes that explicit extension route.  A
NetCDF file contains:

- dimensions `s` and `mode`;
- coordinates `s(s)` and signed poloidal integers `m(mode)`;
- global integer `toroidal_mode`;
- complex covariant components split into
  `B_s_real/imag(mode,s)`, `B_phi_real/imag(mode,s)`, and
  `B_theta_real/imag(mode,s)`;
- required convention attributes `coordinate_order="s,phi,theta"`,
  `component_variance="covariant"`, `radial_coordinate="s_tor"`,
  `magnetic_component_units="T m"`, and the Fourier convention above.

NEO-RT differentiates the angular covariant components on the nonuniform
radial grid and emits the Jacobian-weighted contravariant current
`J * J^i`.  Keeping this weighted quantity avoids artificial poloidal-mode
coupling by the angle-dependent Boozer Jacobian:

```sh
fo exec neo_rt_boozer_current.x perturbation_vector.nc current.out
```

The optional third argument overrides SI `mu0`.  Output rows contain `s`, `m`,
and the real and imaginary parts of all three weighted current components.
`test_boozer_vector_io` verifies the complete reader-and-current path against
the analytic `curl(curl(A))` of two manufactured vector-potential harmonics on
a nonuniform radial grid.

### Paired axisymmetric/vector NetCDF torque profile

The full SI path takes the existing axisymmetric Boozer chartmap as its
geometry input and the vector perturbation contract above as its response-field
input:

```sh
fo exec neo_rt_boozer_jxb.x \
  axisymmetric_chartmap.nc perturbation_vector.nc jxb_profile.out
```

The optional fourth command-line value overrides SI `mu0`. NEO-RT reads the
chartmap's `x`, `y`, and `z` arrays in metres, differentiates both periodic
angles spectrally, and interpolates geometry in `rho_tor=sqrt(s_tor)`. For
coordinate order `(s,phi,theta)`, it uses the signed metric determinant to
raise the covariant perturbation:

\[
C^i=\mathcal J B^i=\mathcal Jg^{ij}B_j.
\]

With \(K^i=\mathcal JJ^i\) from the Ampère kernel, the volume-weighted
covariant toroidal force is

\[
\mathcal J(\mathbf J\times\mathbf B)_\phi
 =K^\theta C^s-K^sC^\theta.
\]

The executable reconstructs every delivered poloidal harmonic on the
chartmap theta grid, evaluates the complex phase average, performs the full
angular quadrature, and writes `s_tor`, `rho_tor`, and
`dT_phi/ds [N m]` on every perturbation surface. The metric raise and
contraction are generated by fortsym. Independent tests compare the local
kernel with explicit Cartesian phase sampling in an oblique basis and compare
the complete radial profile with analytic Cartesian `J cross B` sampling in a
circular torus.

This path reads no MARS `JPLASMA`, `TORQUEJXB`, `PROFEQ`, or `RMZM` output.
If its vector perturbation was converted from `BPLASMA`, MARS still supplied
the response field itself. A vector perturbation produced by another
shielding-response model uses exactly the same NEO-RT current and torque path.

The production ITER chartmap stores two toroidal planes. Because it is
axisymmetric, NEO-RT uses the exact rotational derivative
`e_phi=(-y,x,0)` rather than trying to differentiate those two samples.

### ITER paired-NetCDF versus native-MARS benchmark

The real TC24 chartmap and all three converted perturbations—total, applied
vacuum, and `PLS-VAC` plasma response—were run through the paired-NetCDF
executable. `python/plot_iter_boozer_mars_jxb_comparison.py` compares every
resulting surface against a native-staggered calculation made from the same
MARS `BPLASMA` outputs:

```sh
python python/plot_iter_boozer_mars_jxb_comparison.py \
  MARSF_results/conversion_phiI000.toml OUTRMAR run/iter_tc24_phiI000 \
  iter_boozer_mars_jxb_comparison.png
```

`OUTRMAR` supplies exactly the `ASPCT`, `R0EXP`, and `B0EXP` values used by
MARS's own
`TORQFAC=R0EXP**3*B0EXP**2/(mu0*ASPCT**2)`. The script changes the native
abscissa from `sqrt(psi_pol)` to `s_tor` by a cellwise Jacobian, so the
integrated torque is invariant under the radial remapping.

The benchmark deliberately exposes a failed approximation. Integrated values
for native staggering versus the collocated Boozer contract are:

- total: `+7.100e3 N m` versus `-5.838e5 N m`;
- applied vacuum: `+5.49e-5 N m` versus `-3.570e3 N m`;
- plasma response: `+1.292e3 N m` versus `-6.410e5 N m`.

The whole-radius curves, cumulative integrals, CSV, and provenance JSON are in
`doc/figures/iter_boozer_mars_jxb_comparison.*`. In particular, the native
vacuum result is numerically zero while the collocated path invents a finite
vacuum torque. This is decisive evidence that a MARS-derived perturbation must
retain the full/half radial staggering through the curl. These collocated
ITER values are a diagnostic of the lossy adapter, not physical torque
predictions. The native curves reconstruct current from exported `BPLASMA`;
they are not claimed to reproduce unavailable combined-coil `JPLASMA` or
`TORQUEJXB` files.

### Full-vector MARS adapter

`python/mars_vector_to_boozer_netcdf.py` reads the three covariant MARS
`BPLASMA` components rather than projecting them immediately to
`delta|B|`. It reconstructs the cylindrical vector, maps equal-arc MARS
surfaces to the axisymmetric Boozer chart, transforms all three components to
the covariant Boozer basis, and writes the vector NetCDF contract above. The
component can be selected explicitly as MARS total, applied vacuum, or their
difference (plasma response).

For the ITER TC24 zero-phasing case:

```sh
python python/mars_vector_to_boozer_netcdf.py \
  MARSF_results/conversion_phiI000.toml iter_phiI000_total_vector.nc \
  --component total \
  --scalar-reference MARSF_results/validation_phiI000/components.npz
fo exec neo_rt_boozer_current.x \
  iter_phiI000_total_vector.nc iter_phiI000_total_current.out
```

The conversion retains 641 poloidal modes on 237 complete surfaces. The
largest discarded vector power fraction is `1.17e-10`, and the largest
relative geometry mismatch is `9.37e-3`. Projecting the transformed vector
back onto the equilibrium field reproduces the separately converted
Eulerian `delta|B|` spectrum with relative L2 difference `1.53e-2`. That
projection is a round-trip validation of the coordinate and vector-basis
transformation; it is not used when NEO-RT calculates the curl.

The resulting NEO-RT current uses no `JPLASMA` data. If the selected
`BPLASMA` is the MARS plasma solution, however, MARS still supplied the
physical shielding response in that magnetic field. Moreover, this first
vector contract is collocated: it is appropriate for smooth vector NetCDF
inputs and for validating the MARS-to-Boozer field transformation, but not
for differentiating exported MARS shielding layers. The native-staggered
path below is the authoritative MARS current reconstruction. A future
staggered vector contract can move that exact derivative through the
Boozer adapter without damping it.

### Independent MARS `B -> J` benchmark

The native-coordinate benchmark uses `RMZM_F.OUT` and `BPLASMA.OUT` to lower
MARS's exported Jacobian-weighted contravariant magnetic components,

\[
\begin{aligned}
B_s&=(g_{ss}C^s+g_{s\chi}C^\chi)/\mathcal J,\\
B_\chi&=(g_{s\chi}C^s+g_{\chi\chi}C^\chi)/\mathcal J,\\
B_\phi&=R^2C^\phi/\mathcal J,
\end{aligned}
\]

then evaluates the same fortsym-derived coordinate curl on the staggered MARS
mesh. `JPLASMA.OUT` is read only after that calculation and serves as the
independent reference:

```sh
python python/plot_mars_ampere_comparison.py \
  RMZM_F.OUT BPLASMA.OUT JPLASMA.OUT mars_ampere_current_profile.png \
  --case-name "MAST-U Pair A"
```

The figure and its machine-readable curves cover every surface on which the
exported staggering determines a current: all 350 half-mesh `J1` points and
349 interior full-mesh `J2/J3` points. No smoothing, fitted scale, phase
rotation, radial remapping, or mode selection is applied. On archived MAST-U
Pair A, the interior relative L2 differences are `2.67e-4`, `9.15e-2`, and
`7.02e-2` for `J1`, `J2`, and `J3`. The angular-only `J1` reconstruction is
nearly exact. `J2/J3` require radial derivatives of exported half-mesh fields;
their remaining difference measures this cubic reconstruction against MARS's
internal finite-element projection.

The TC24 postprocessing archives contain `JPLASMA.OUT` files whose plasma
rows are identically zero, so they cannot serve as a current oracle. The
ITER vector conversion remains valid as a field-transformation benchmark,
but current agreement is benchmarked against the nonzero MAST-U reference
rather than reported as a meaningless division by zero.

### ITER shielding-current profiles

For MARS output, the curl must precede radial collocation. `B1` is exported
on the full mesh while `B2/B3` are exported on the half mesh; centering all
three fields first preserves a scalar field projection but suppresses the
sharp radial derivatives that constitute shielding current. The ITER profile
utility therefore combines coil rows and phases first, keeps the native
staggering, lowers the vector with the MARS metric, and differentiates:

```sh
python python/plot_iter_mars_shielding_current.py \
  MARSF_results/conversion_phiI000.toml \
  iter_mars_shielding_current_profile.png \
  iter_mars_resonant_current_modes.png
```

The first figure contains the complete resolved profiles of all three
Jacobian-weighted current components for the applied vacuum field, pure
plasma response (`PLS-VAC`), and total field. The three independently
calculated pieces close linearly to `6.03e-16` relative L2. The applied
vacuum curl is `1.45e-3` of the total current norm, providing a numerical
zero-current check without `JPLASMA`.

The second figure contains every rational harmonic in the delivered ITER
range: `m=4,...,10` for `n=-3`, with the corresponding `q=m/3` surfaces
marked at `rho_tor=sqrt(s_tor)` values
`0.5803, 0.6989, 0.7815, 0.8441, 0.8921, 0.9287, 0.9847`.
It plots the coordinate-weighted tangential amplitude
`sqrt(|sqrt(g) J^chi_m|^2 + |sqrt(g) J^phi_m|^2)`. This is a transparent
native-MARS diagnostic, not a coordinate-invariant physical current norm.
Both figures cover the full resolved radial domain, apply no smoothing or
mode filtering, and have companion CSV and provenance JSON files. The CSV
retains both `s_tor` and the published `rho_tor` abscissa.

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

This implementation calculates current and torque from a supplied vector
magnetic response; it does not yet solve the plasma shielding response that
produces that vector field. Consequently:

- scalar `Delta|B|` input is insufficient;
- vacuum fields without plasma current give zero plasma `j×b` torque;
- ideal singular currents require a specified distributional or layer
  integration convention;
- smoothing, rational-surface masks, and sign flips are never applied
  implicitly;
- all Fourier, coordinate, Jacobian, radial, and unit conventions remain the
  responsibility of the response-field adapter.
