---
title: Circular Kasilov offset-rotation showcase
---

# Circular Kasilov offset rotation

NEO-RT evolves neither a fluid toroidal velocity nor an axisymmetric poloidal
flow. It evaluates non-axisymmetric transport and torque for a prescribed
electric precession frequency. The circular test demonstrates how a consumer
converts that response to the intrinsic contravariant toroidal rotation of
Kasilov et al. (2014), without changing solver output.

## Quantities and units

The code input is the electric field-line precession frequency

```text
Omega_tE = -c Phi'(s) / chi'(s) = c E_s / chi'(s),
chi'     = sqrt(g) B^theta = sign_theta psi_tor' / q.
```

The signs follow from `E_s = -Phi'`. A Cartesian `E cross B` construction in
the test independently checks the field-line-label velocity
`V^phi - q V^theta` and the corresponding electric term in `A1`.

Kasilov's fluid quantities are

```text
k       = 5/2 - D32/D31,
k_NA    = D12_NA/D11_NA - 5/2
          + <B_phi^2> k/(<B^2><g_phi_phi>),
Vphi_in = c k_NA/(e_i chi') dTi/ds.
```

These are the single-ion, negligible-electron-flux forms of Eqs. (33)--(35)
in [Kasilov et al., *Physics of Plasmas* 21, 092506
(2014)](https://doi.org/10.1063/1.4894479).

`Vphi_in` has units `1/s`: it is the contravariant velocity over a
dimensionless toroidal angle. The test uses the requested circular
large-aspect-ratio closure `k = 1.17` and

```text
<B_phi^2>/(<B^2><g_phi_phi>) ~= (Bphcov/(B0 R0))^2.
```

This value of `k` is a test-only circular closure. A general-equilibrium
consumer must obtain `D31`, `D32`, and the geometry factor independently, for
example from an axisymmetric NEO-2 calculation.

The fluid rotation corresponding to a NEO-RT evaluation is reconstructed from
radial force balance,

```text
Vphi = Omega_tE
       - c/(e_i chi') (Ti n_i'/n_i + Ti')
       + c/(e_i chi') geometry_factor k Ti'.
```

Temperatures and derivatives are converted from eV to cgs energy. Magnetic
precession is not added to `Vphi`: `Omega_tB(v,eta)` is orbit-dependent and
already enters the resonant `D11_NA` and `D12_NA` when magnetic drift is
enabled.

## Behavioral test

For each of four circular flux surfaces, the test repeatedly runs NEO-RT while
bisection solves the implicit condition

```text
Vphi(Omega_tE) - Vphi_in(D11(Omega_tE), D12(Omega_tE)) = 0.
```

The transport coefficients are recomputed at every trial frequency. The test
then evaluates independent below, at, and above cases in a four-thread OpenMP
loop. It requires positive restoring torque below the offset, vanishing torque
at the offset, and negative restoring torque above it in the native coordinate
convention.

The plotted points are the unsmoothed solver evaluations at `s = 0.075`.
The test generates the source data, and the plotting script in
`examples/circular_offset_rotation` renders it. Generated CSV, PNG, and PDF
files are intentionally not committed.

## Reproduce the showcase

From the repository root:

```bash
NEORT_OFFSET_SHOWCASE="$PWD/examples/circular_offset_rotation/circular_offset_rotation.csv" \
    fo test test_circular_offset_rotation
python examples/circular_offset_rotation/plot_showcase.py
```

Normal test execution creates no data or plot files. Setting
`NEORT_OFFSET_SHOWCASE` only enables the reproducible showcase export from the
test program; production output files and schemas remain unchanged.
