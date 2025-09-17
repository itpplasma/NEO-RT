---
title: Output file formats
---

# Output file formats

This reference describes the structure of the plain-text files written by `neo_rt.x`. Each file is written in free-format (values separated by spaces). Unless noted otherwise, numeric quantities use cgs units consistent with the solver internals.

## `<runname>.out`

Summary of transport coefficients evaluated at a single flux surface. The first line is a header. The data line contains:

1. `M_t` – Toroidal Mach number used in the run.
2. `D11co` – Contribution of co-passing particles to the radial particle flux coefficient `D11`.
3. `D11ctr` – Contribution of counter-passing particles to `D11`.
4. `D11t` – Contribution of trapped particles to `D11`.
5. `D11` – Sum of the three `D11` contributions.
6. `D12co` – Contribution of co-passing particles to the momentum transport coefficient `D12`.
7. `D12ctr` – Contribution of counter-passing particles to `D12`.
8. `D12t` – Contribution of trapped particles to `D12`.
9. `D12` – Sum of the three `D12` contributions.

`D11` corresponds to the particle flux coefficient and `D12` to the momentum flux coefficient in the effective radius approximation implemented in `compute_transport`.

## `<runname>_integral.out`

Velocity-space integrals for each evaluated poloidal harmonic `mth`. The file appends one line per harmonic after the solver finishes:

1. `M_t` – Toroidal Mach number for the current evaluation.
2. `mth` – Poloidal mode number associated with the resonance under consideration.
3. `D11co` – Co-passing contribution to `D11` for the current harmonic.
4. `D11ctr` – Counter-passing contribution to `D11`.
5. `D11t` – Trapped-particle contribution to `D11`.
6. `D11` – Sum of the three `D11` contributions for the current harmonic.
7. `D12co` – Co-passing contribution to `D12` for the current harmonic.
8. `D12ctr` – Counter-passing contribution to `D12`.
9. `D12t` – Trapped-particle contribution to `D12`.
10. `D12` – Sum of the three `D12` contributions for the current harmonic.
11. `vminp/vth` – Lower bound of the passing velocity grid normalised by the thermal speed.
12. `vmaxp/vth` – Upper bound of the passing velocity grid normalised by the thermal speed.
13. `vmint/vth` – Lower bound of the trapped velocity grid normalised by the thermal speed.
14. `vmaxt/vth` – Upper bound of the trapped velocity grid normalised by the thermal speed.

## `<runname>_magfie_param.out`

Magnetic-field parameters and derived quantities collected at the start of the run. The file is human-readable and contains key-value pairs such as `s`, `R0`, `psi_pr`, `B0`, `q`, `iota`, `dVds`, and geometric scaling factors. When `nonlin=.true.` additional diagnostics (`dpp`, `dhh`, `dfpeff`) are printed. All quantities correspond to the initial flux surface and match the variables written in `neort:check_magfie`.

## `<runname>_magfie.out`

Sampling of the magnetic field along the Boozer poloidal angle for inspection and plotting. Each line contains 19 columns:

1. `theta` – Boozer poloidal angle [rad].
2. `|B|` – Magnetic-field magnitude.
3. `sqrt(g)` – Jacobian determinant.
4–6. `hder(1:3)` – Contravariant basis vectors.
7–9. `hcovar(1:3)` – Covariant basis vectors.
10–12. `hctrvr(1:3)` – Contravariant components of the field-line curvature.
13–15. `hcurl(1:3)` – Components of the curl of the basis vectors.
16. `Re(bn)` – Real part of the perturbation amplitude `bn` normalised by `|B|`.
17. `Im(bn)` – Imaginary part of `bn`.
18. `Re(epsmn*exp(i*m0*theta))` – Reference analytic perturbation.
19. `Im(epsmn*exp(i*m0*theta))` – Imaginary part of the reference analytic perturbation.

When `pertfile=.true.` the `bn` columns correspond to the perturbation read from `in_file_pert`; otherwise the analytic perturbation is used.

## `<runname>_torque.out`

Written when `comptorque=.true.`. Contains one line with:

1. `s` – Normalised toroidal flux.
2. `dV/ds` – Derivative of the flux-surface volume.
3. `M_t` – Toroidal Mach number.
4. `Tco` – Co-passing contribution to the torque density.
5. `Tctr` – Counter-passing contribution to the torque density.
6. `Tt` – Trapped-particle contribution to the torque density.

## `<runname>_torque_integral.out`

Also written when `comptorque=.true.`. Appends one line per harmonic with:

1. `mth` – Poloidal mode number.
2. `Tco` – Co-passing contribution to the torque integral.
3. `Tctr` – Counter-passing contribution to the torque integral.
4. `Tt` – Trapped-particle contribution to the torque integral.

The sum of the three torque components matches the values written in `<runname>_torque.out` once all harmonics have been processed.
