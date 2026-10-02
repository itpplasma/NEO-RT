# Periodic-cylinder perturbation adapter

`python/periodic_cylinder.py` maps complex KiLCA/KAMEL harmonics through an
explicit chart into physical orthonormal `(R, phi, Z)` components. Its sample
NPZ output is an intermediate adapter format; it is not the libneo NetCDF
perturbation container or a rectangular R-Z grid.

## Native contract

`EB.dat` has 13 columns: radius, then real/imaginary pairs for
`Er, Etheta, Ez, Br, Btheta, Bz`. Lengths are cm, B is G, A is G cm,
E is statV/cm, and Phi is statV. Native modes use
`exp(i*(m*theta+n*z/R0-omega_lab*t))`, with signed integer m/n and machine
R0 independent of the magnetic-axis radius. The laboratory frequency in
`mode_data.dat` is in Hz; `laboratory_omega` multiplies its laboratory entry
by 2 pi and does not substitute the Doppler-shifted moving-frame frequency.
The adapter reconstructs potentials for nonzero n.

An explicit continuous radial interval is selected before checking duplicate
radii. Exact repeated rows are accepted. Conflicting interface or antenna
rows are rejected, including when the full archived domain is selected;
the adapter never averages a discontinuity. For AUG 33353 at 2900 ms,
`--r-max-cm 61.8` selects the closed-flux interval before the conflicting
r=67/70 cm interfaces. Interpolation and chart evaluation reject extrapolation.

## Geometry, orientation, and curl

The supplied embedding is
`R=R(r,theta), Z=Z(r,theta), phi=z/R0+nu(r,theta)`.
The physical deformation, with columns in native `(er, etheta, ez)` and rows
in target `(eR, ephi, eZ)`, is

```text
M = [[R_r,       R_theta/r,       0],
     [R*nu_r,    R*nu_theta/r,    R/R0],
     [Z_r,       Z_theta/r,       0]]
J = R/(R0*r) * (R_theta*Z_r - R_r*Z_theta)
A_target = inverse(transpose(M)) A_native
E_target = inverse(transpose(M)) E_native
B_target = M B_native / J
```

These are the one-form and magnetic-flux transformations. Every exported
single-n amplitude also contains `exp(i*(m*theta-n*nu))`; its toroidal/time
factor is `exp(i*n*phi-i*omega_lab*t)`. The adapter requires J>0. The analytic
clockwise circle `R=R0+r*cos(theta), Z=Z0-r*sin(theta), nu=0` has J=R/R0.
Reversing the poloidal orientation without updating the native convention
is rejected. Radial labels are supplied explicitly and are never interpreted
as geometric minor radius. `SurfaceChart` interpolates coordinate functions
and differentiates those same functions to provide an integrable chart.

## Magnetic and paired electromagnetic exports

The magnetic reconstruction fixes `Ar=0`, `Az(r_min)=0`, and
`Atheta(r_min)=i*Br(r_min)/(n/R0)`. It integrates Btheta and r times the
interpolated Bz exactly, then independently checks the remaining Br equation
at radial nodes and midpoints. A magnetic export includes diagnostic E but
does not provide a scalar potential or assert that the supplied E is static.

The paired laboratory reconstruction fixes Az=0 and uses
`Ar=-i*Btheta/kz, Atheta=i*Br/kz, Phi=i*Ez/kz`, where kz=n/R0.
All supplied components must pass explicit curl-A, paired-E, Faraday, and
div-B limits before export. Its Gaussian convention is
`E=-grad(Phi)+i*omega_lab*A/c`. Independent tests reject a reversed time sign,
a missing 2 pi conversion, and substitution of the moving-frame frequency.
Low-frequency Faraday diagnostics retain absolute residuals and the sum-of-
terms conditioning factor; cancellation does not justify relaxing a limit.

```sh
python3 python/kilca_to_toroidal.py EB.dat --m 5 --n 2 --R0-cm 165 \
  --r-max-cm 61.8 --geometry chart.npz --output magnetic.npz \
  --max-relative-br-residual 1e-6
```

For `--gauge axial-electromagnetic`, supply all four explicit limits:
`--max-relative-B-residual`, `--max-relative-E-residual`,
`--max-relative-Faraday-residual`, and `--max-relative-div-B`.
`--diagnostics-only` reports source residuals without exporting fields.

## Reconstructed AUG chart and remaining limits

`python/geqdsk_cylinder_chart.py` uses the recovered equilibrium table's
outboard anchor and radial label to regenerate a clockwise straight-field-line
angle from a specified GEQDSK. It requires an explicit sign relating native q
to the stored table, validates signed q and positive J, and records all three
input hashes. Independent circular finite-aspect formulas test both signed q
and enclosed toroidal flux.

The archived label is `r_eff=sqrt(2*abs(Phi_label/Btf))`. KAMEL's legacy
`field_line_rhs` accumulator uses the contour value F(surface) over the whole
enclosed weighted area. At finite beta this differs from true enclosed flux,
which integrates F(psi)/R over the interior. The builder preserves the native
label exactly and reports separate legacy-label and true-flux discrepancies;
neither quantity is silently substituted. For the matched AUG data the
difference is about 0.53%, while the legacy-label reconstruction agrees much
more closely. Mapped cylinder B0 can consequently differ from GEQDSK B0.

The historical angle table and byte-identical original GEQDSK were not
recovered. The reconstructed chart has declared physical matching checks;
it does not certify the missing historical chart. Archived FLRE lab E has
large low-frequency derivative cancellation and remains unverified for a
full paired electromagnetic drive. Its magnetic channel and the vacuum
paired channel have separate acceptance gates. These limits must accompany
downstream orbit-drive comparisons.
