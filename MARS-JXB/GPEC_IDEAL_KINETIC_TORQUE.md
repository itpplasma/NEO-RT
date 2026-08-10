# Ideal shielding current and GPEC kinetic torque

## Result in one sentence

An ideal GPEC response can contain large local shielding currents and local
Lorentz forces, but its flux-surface-averaged and volume-integrated toroidal
torque is zero. Finite GPEC torque appears when the kinetic pressure operator
is enabled, and that torque is neoclassical toroidal viscosity rather than a
resonant-layer \(J\times B\) torque.

## Ideal torque null

For a static ideal-MHD perturbation, the linear force balance is

\[
\delta\mathbf{J}\times\mathbf{B}_0
+\mathbf{J}_0\times\delta\mathbf{B}
=\nabla\delta p .
\]

Neither term on the left must vanish locally. Their toroidal moment can have
large positive and negative structure. On a closed axisymmetric equilibrium,
however, the volume integral of the toroidal moment of the pressure gradient
vanishes. The ideal force operator is self-adjoint, so its quadratic response
energy is real. There is no anti-Hermitian part that can absorb power or
produce an average torque.

GPEC exposes the same identity in its control-surface output. The source prints

\[
T_\phi=-2n\,\operatorname{Im}(p_y).
\]

For an ideal response, \(p_y\) should be real. A small printed value is an
anti-Hermitian numerical residual. It is not evidence for physical braking.
The old TC24 ideal calculation at source commit
`a430ad2ac27d0acc8cc09b7c261b161b09a3ef53` printed a torque of order
\(1\ \mathrm{N\,m}\) beside a plasma energy of order
\(1.25\times10^4\ \mathrm{J}\). The relative imaginary part was about
\(10^{-5}\), consistent with a numerical residual.

This also explains why summing torque from ideal shielding-current sheets can
reproduce the same small global number. Such a sum reconstructs the
anti-Hermitian residue of that discrete response. It does not turn the residue
into a physical torque.

The cancellation can be lost if a calculation uses only one part of the
linearized Lorentz force, crosses a response current with the wrong magnetic
field, changes the native staggering, or mixes independently interpolated
fields and currents. A physical production claim therefore needs either the
native torque diagnostic of the source code or a reconstruction that preserves
all terms and reproduces that diagnostic.

Markl et al. state the result directly in Appendix B. The ideal response is
omitted from their torque integral because it produces no average torque. They
also note that the ideal Lorentz force equals a pressure gradient and cannot
produce integral torque [Markl et al., Nuclear Fusion 63, 126007
(2023)](https://doi.org/10.1088/1741-4326/acf587).

## What changes in kinetic GPEC

Kinetic GPEC sets `kin_flag = .true.` and continues through the resonant
surfaces with `con_flag = .true.`. PENTRC supplies a drift-kinetic anisotropic
pressure response. The resulting force operator is not self-adjoint when
neoclassical toroidal viscosity is finite. Its anti-Hermitian energy gives a
finite torque.

Park and Logan derive this construction. The first-order anisotropic-pressure
force changes the local perturbed equilibrium, while its second-order
flux-surface average gives NTV. Their energy and torque integral is the basis
of the GPEC torque response matrix [Park and Logan, Physics of Plasmas 24,
032505 (2017)](https://doi.org/10.1063/1.4977898).

The public GPEC documentation confirms the executable contract:

- `kin_flag` selects the kinetic Euler-Lagrange equation.
- `con_flag` continues the integration through layers.
- `kinfac2` scales the torque contribution.
- `dw_flag` writes self-consistent energy and torque profiles only when the
  DCON result is kinetic.

See the [GPEC input and output
documentation](https://princetonuniversity.github.io/GPEC/outputs.html).

## Three torque channels that must stay separate

| Channel | Nonideal mechanism | What the current calculation supplies |
| --- | --- | --- |
| GPEC and PENTRC NTV | Nonambipolar drift-kinetic pressure response in a three-dimensional field | A self-consistent NTV profile from `gpec_dw_profile_n3.out` |
| KiLCA and KAMEL resonant torque | Finite-width kinetic resonant layer, with a phase-shifted parallel current and absorbed power | KiLCA supplies the torque. GPEC ideal shielding current supplies only the toroidal amplitude correction |
| MARS-F electromagnetic torque | Resistive and rotating MHD response, evaluated from the native perturbed current and magnetic field | Native MARS \(J\times B\) output is still required for the TC24 electromagnetic benchmark |

The first and third rows are distinct momentum channels. They must not be
overlaid as if one reproduced the other.

## KAMEL current rescaling

KAMEL does not assign torque to the ideal GPEC current by itself. The local
QL-Balance implementation computes

\[
A_{\mathrm{ant}}=
\left(\frac{I_{\parallel,\mathrm{GPEC}}}
{I_{\parallel,\mathrm{KiLCA}}}\right)^2
\]

after converting the two current conventions. It multiplies KiLCA
quasilinear transport coefficients by this factor. The square is required
because quasilinear power, transport, and torque are quadratic in perturbation
amplitude. The nonideal phase and the torque remain those of the KiLCA layer
solution.

This follows the coupling used by Markl et al. Their equation for the toroidal
correction is a ratio of the ideal toroidal GPEC shielding current to the
cylindrical current. The correction then enters quadratically.

The underlying KiLCA model computes a finite electromagnetic force in a
kinetic resonant layer. Heyn et al. relate the mode torque to absorbed power,
with toroidal torque proportional to \(nP/\omega\). The power absorption is
the required nonideal ingredient [Heyn et al., Physics of Plasmas 15, 052509
(2008)](https://doi.org/10.1063/1.2913264).

## Kasilov NTV distinction

Kasilov et al. call the applied three-dimensional magnetic perturbation
nonresonant and ideal. Here, ideal describes the imposed magnetic geometry
used for guiding-center motion. NEO-2 then solves a drift-kinetic equation.
The torque follows from nonambipolar radial particle flux through the
force-flux relation. It is NTV, not the volume integral of the ideal-MHD
shielding-current Lorentz force [Kasilov et al., Physics of Plasmas 21, 092506
(2014)](https://doi.org/10.1063/1.4894479).

NEO-RT evaluates this same nonresonant kinetic torque family in low
collisionality resonant transport regimes. Calling the magnetic perturbation
ideal does not make the particle response self-adjoint or ambipolar.

## Interpretation rules for the TC24 evidence

1. The kinetic GPEC profile is accepted only from the same executable and
   unchanged physical inputs as its ideal twin.
2. No radial interpolation, smoothing, phase fitting, or amplitude fitting is
   used in the comparison.
3. The profile endpoint must close against the independent global torque
   printed by GPEC.
4. The plotted sign is the native GPEC sign. A machine-coordinate CCW sign is
   not claimed until the handedness conversion has an independent
   manufactured-mode oracle.
5. A successful kinetic GPEC run validates the NTV path and the ideal torque
   null. It does not replace the missing native MARS current needed for the
   electromagnetic \(J\times B\) benchmark.

## Executed TC24 kinetic check

The same-source twin campaign
`tc24-gpec-kinetic-2e031ae6-20260730-r3` completed on 2026-07-30. It used GPEC
`v1.5.7-581-g2e031ae6` and the accepted supplier-profile correction. The ideal
and kinetic decks differ only in `kin_flag` and `con_flag`.

The ideal run produced no kinetic torque profile and printed
`-1.24265549 N m`. This reproduces the small residue of the archived ideal
run. The kinetic run printed `255.783088 N m`. Its native profile endpoint is
`255.782098 N m`, a relative closure difference of `3.87e-6`.

This is accepted evidence for a finite self-consistent GPEC NTV torque. The
unsmoothed local density has very large sign-changing spikes near ideal
rational surfaces. Their amplitudes are not promoted as converged local
physics without a radial and layer resolution study. The finite global
integral is independently checked by the control-surface scalar.

The full evidence and warning inventory is in the [run
report](results/tc24-gpec-kinetic-2e031ae6-20260730-r3/REPORT.md).

The run also demonstrates an output-only source defect. At the executed
commit, the NetCDF writer packs a `0:mstep` spline array into a shape of
`mstep`, shifting values between the real and imaginary planes. The accepted
profile therefore comes from `gpec_dw_profile_n3.out`, whose endpoint closes
the global scalar. The NetCDF `T` and `dTdpsi` arrays are rejected.
