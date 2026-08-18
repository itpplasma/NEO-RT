# Circular offset-rotation showcase

This example visualizes the circular `k = 1.17` behavioral test documented in
[`doc/offset_rotation.md`](../../doc/offset_rotation.md). It does not add or
change any production NEO-RT output.

Generate the data and local plots from the repository root:

```bash
NEORT_OFFSET_SHOWCASE="$PWD/examples/circular_offset_rotation/circular_offset_rotation.csv" \
    fo test test_circular_offset_rotation
python examples/circular_offset_rotation/plot_showcase.py
```

The generated CSV contains unsmoothed NEO-RT evaluations at normalized toroidal flux
`s = 0.075`. Frequencies are in `s^-1`; torque is shown in the solver's native
output units. Its horizontal scan coordinate is the imposed electric-frequency
shift from the self-consistent root; at fixed pressure and poloidal-flow terms,
this is also the imposed shift in fluid `Vphi`.

Generated CSV, PNG, and PDF files are not committed.
