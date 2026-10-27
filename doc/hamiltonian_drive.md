# Perturbation Hamiltonian: Boozer form and general ideal drive

NEO-RT evolves the guiding centre in the canonical action-angle variables of the
unperturbed axisymmetric field. A perturbation enters only through the
first-order Hamiltonian `H_n`. Per unit energy `E`, with `eta*B = mu*B/E`:

| Input | `H_n / E` |
| --- | --- |
| `in_file_pert` only (default) | `(2 - eta*B) * dB_L/B` |
| plus `pertfile_curvature = .true.` and `in_file_pert_kappa` | `eta*B * dB_L/B + 2*(1 - eta*B) * kappa.xi_perp` |

Here `dB_L` is the Lagrangian change of `|B|`, `kappa = (b.grad) b` and
`xi_perp` is the perpendicular displacement.

## Derivation

The guiding-centre one-form is `gamma = (e/c) A + m v_par b + (mc/e) mu dzeta` and
`H = m v_par^2/2 + mu |B| + e Phi`. A perturbation changes `gamma` by `delta gamma`
and `H` by `delta H`. In unperturbed canonical variables the first-order
Hamiltonian is

    H1 = delta H - delta gamma(zdot_0) + dg/dt.

The near-identity Darboux shift that removes the non-exact part of
`delta gamma` gives this result, and `g` is a free gauge function.

1. **Boozer form.** In Boozer variables `rho_par = m c v_par/(e B)`,
   `H = (e B rho_par)^2/(2 m c^2) + mu B`. At fixed canonical momenta, and
   therefore fixed `rho_par`, a change `dB` of `|B|` gives
   `dH = (m v_par^2/B + mu) dB = E (2 - eta B) dB/B`. This form is exact when
   the perturbed field, pulled back along the displacement, has an unchanged
   covariant parallel component, `beta_b = 0` below. The field is then
   Boozer-representable with the same `I` and `G`.
2. **Pullback along an ideal displacement `xi`.** For `delta A = xi x B0` the
   pullback of `A1` is `A0 + grad(xi.A0)`, which is an exact gauge. The pullback
   of `|B1|` is `|B0| + dB_L`. The pullback of the covariant field has the
   parallel component `|B0| beta_b`, with `beta_b = dB_L/B + b.(b.grad)xi`.
   Hence `H1 = E (2 - eta B) dB_L/B - m v_par^2 beta_b`.
3. **Streaming identity.** With `u = xi.b`,
   `b.(b.grad)xi = b.grad(u) - kappa.xi_perp`. Along the unperturbed orbit
   (`xdot = v_par b`, `m dv_par/dt = -mu b.grad B`) the Hamiltonian of step 2 is

       H1 = mu dB_Lperp + m v_par^2 kappa.xi_perp - dS/dt,   S = m v_par u.

   This holds exactly. `dB_Lperp = b.deltaB + xi_perp.grad|B|`. The Eulerian
   route (`delta A = xi_perp x B0` at fixed position) gives the same
   `mu dB_Lperp + m v_par^2 kappa.xi_perp` directly. That route uses the
   magnetic drift, an `O(rho*)` relative remainder, and a potential frozen to the
   displaced surfaces: the `E x B` term cancels `delta Phi = -xi.grad Phi0`.
4. **Torque.** The orbit harmonic of `dS/dt` is `i w_k S_k`, where `w_k` is the
   full detuning. The resonant torque (`nu -> 0`, `w_k = 0`) is therefore
   independent of `S`. The Krook torque has the weight
   `|H_k|^2 nu/(w_k^2 + nu^2)` and changes by
   `nu (w_k^2 |S_k|^2 + 2 w_k Im(conj(H_k) S_k)) / (w_k^2 + nu^2)`.
   A finite Krook rate therefore depends on the representation.

The default drive equals the general one exactly when
`kappa.xi_perp = dB_L/B`. At low beta and in static force balance the
difference is `2(1 - eta B)(mu0 xi.grad p/B^2 - b.deltaB/B)`. Otherwise the
difference `m v_par^2 (dB_L/B - kappa.xi_perp)` is a total derivative `dS/dt`
only when `b.grad(u B) = B (kappa.xi_perp - dB_L/B)` has a single-valued
solution `u`, and only up to the `u b.grad|B|` part of `dB_L`. That part is
`O(eps)` relative, but a small parallel wavenumber enhances it. A strongly
compressive test field has an O(1) difference. There the two forms agree at
`nu -> 0` only to the order of this remainder. At finite Krook rate they differ
by the off-resonant `i w_k S_k` term, which reaches deep-passing phase space.

Krook relaxes the `f1` that is canonically conjugate to the chosen `H1`.
`pertfile_curvature` selects the representation with no parallel relabel,
`u = 0`. In it, `H1` contains only the physical perpendicular displacement.
This is also the form used by PENTRC and MARS-K. The Krook model does not fix
the parallel relabel. Only a conserving collision operator transformed together
with `H1` would make the finite-rate torque independent of it.

## Input

Set `pertfile = .true.` and `pertfile_curvature = .true.`. Provide
`in_file_pert_kappa` in the layout of `in_file_pert`, with the same harmonic
rows. The radial grid may differ. Its `bmnc`/`bmns` columns hold
`|B0| kappa.xi_perp` in tesla on the unperturbed Boozer angles, in the complex
convention of `in_file_pert`. The file is scaled like `in_file_pert`
(`1e4`, `bfac`, `pertfile_scale`) and divided by the local `|B|`. Both files
must use the same perpendicular displacement. If a displacement with a parallel
component `u` was used for `dB_L`, remove `u b.grad|B|` first. Passing
`in_file_pert` itself as the curvature table reproduces the default drive to
rounding, as `test_curvature_drive` checks.

Curvature tables require `USE_STANDALONE=ON`; the native NEO-2 adapter rejects
this option explicitly.

The drive weights are generated by fortsym (`gen/gen_drive.f90`,
`src/generated/drive_leaf.inc`).
