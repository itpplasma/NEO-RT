# NEO-RT full-FOW structural contracts

This directory contains small, proof-carrying structural contracts for the
full real-space FOW path.  They are deliberately separate from the Fortran
production sources: these files describe ownership and bookkeeping interfaces,
not a formalization of floating-point field evaluation, interval arithmetic,
or the executable's numerical results.

## Toolchain and exact dependency status

The package is pinned to Lean `4.14.0` and uses only the installed `Std`
library.  `lean` and `lake` are available.  A direct `import Mathlib` probe
fails with `unknown module prefix 'Mathlib'`; the local `.ltar` cache is not a
usable Lake dependency.  Therefore mathlib-dependent contracts are **not
claimed** and are kept behind the explicit `check-mathlib.sh` gate.

The available proof target is built with:

```text
lake build
lake env lean -o /tmp/neort_full_fow_contracts.olean NeortFullFow.lean
```

The mathlib gate is expected to fail with exit status `2` in the current
environment:

```text
./check-mathlib.sh
```

`MathlibGate.lean` intentionally imports the real `Mathlib` module.  It is not
part of the default target, and there are no `sorry`, `axiom`, or substitute
stubs in the compiled contracts.

## Contract inventory

| Lean module | Contract | Full-FOW source anchor |
| --- | --- | --- |
| `ClassPartition` | Three non-collapsed classes and shared seam equations | `src/gc_eqdsk_composite_r_ownership.f90:4-12, 152-200` |
| `BoundaryOwnership` | First seam is axis-owned; second seam and final endpoint are outboard-owned | `src/gc_eqdsk_composite_r_ownership.f90:7-17, 176-200` |
| `PartnerConsistency` | Partner boundaries have equal canonical momentum and distinct boundary IDs | `src/gc_operational_partner_crossings.f90:3-11, 145-199` |
| `CanonicalVariation` | A variation certificate carries an enclosure and derivative-coverage bit | `src/gc_certified_canonical_variation.f90:90-120, 230-244` |
| `OrientationSign` | Reversing the declared orientation negates the signed quantity | `src/gc_full_fow_canonical_behavior.f90` tests and `src/transport.f90:182` |
| `StateCompleteness` | Runtime, backend, wall, class/measure, provenance, quadrature, harmonic, and refinement fields are all required | `src/gc_full_fow_runtime_metadata.f90:66-115` |

The source anchors are provenance only.  The Lean contracts do not prove that
the Fortran implementation satisfies them; the existing Fortran tests and
runtime metadata remain the executable evidence for that refinement.
