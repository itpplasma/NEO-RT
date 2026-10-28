---
title: Diagnostics and probes
---

# Diagnostics and probes

Besides the solver `neo_rt.x`, the build produces three auxiliary executables.

## neo_rt_diag.x

`neo_rt_diag.x` links the NEO-RT library and the
[fortplot](https://github.com/lazy-fortran/fortplot) plotting backend. The
drivers in `src/diag/` provide post-processing for debugging, harmonic
inspection and attenuation maps.

## Field probe

```bash
neo_rt_field_probe.x BCF INPUT_SWITCH POINTS OUTPUT
```

Evaluates the Boozer field from file `BCF` at the points listed in `POINTS`,
one whitespace-separated triple `(s_tor, phi_B, theta_B)` per line, angles in
radians. Profiles and transport are not initialised. Input switches 8 and 9 are
supported. The CSV output contains the field magnitude, covariant components,
rotational transform, Jacobian and flux derivative in the CGS units stated in
its header. Points outside the radial range of the file are rejected.

`tools/run_field_probe.py` wraps the executable, also accepts a CSV point list,
and records hashes of input, executable and output next to the result.

## Profile probe

```bash
neo_rt_profile_probe.x BCF INPUT_SWITCH PLASMA MACH EFAC BFAC POINTS OUTPUT
```

Runs only the plasma, Mach-number and collision initialisation of the solver,
without orbits or transport, and reports the resulting profile data at the
radii `s_tor` listed one per line in `POINTS`. `PLASMA` is a plasma profile file
in the solver format (species header and six columns), `MACH` a two-column
Mach-number table, and `EFAC` and `BFAC` are the electric- and magnetic-field
scaling factors of the solver namelist (see [running.md](running.md)). All points must lie within the radial range covered by the
field, plasma and Mach data.

The CSV output contains profile values and their derivatives with respect to
`s_tor`, Mach number, thermal speed and rotation, the major radius `R0` of the
field file, the orbit masses and charge, and the collision arrays in their
native normalisation (`efcolf` inverse free path, `velrat` velocity ratios,
`enrat` energy ratios). Electron density, `Zeff` and thermal pressure are
derived from the loaded ion background. Units are eV, cm^-3, cm/s, g and
erg/cm^3. The header lists the mass and energy constants used separately by the
orbit and collision routines.

Positive background ion densities are required; a zero density of the second
ion species is allowed. The probe runs with `nonlin=.false.`.

`tools/run_profile_probe.py` accepts the same CSV point lists as the field
wrapper and records hashes of all inputs, the executable and the output.
