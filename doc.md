project: NEO-RT
summary: Neoclassical toroidal viscosity solver for resonant transport regimes
author: NEO-RT Developers
output_dir: build/doc
src_dir:
  - src
  - src/diag
exclude:
  - doc/driftorbit.lyx
  - **/POTATO/**
  - build
  - build/**
exclude_dir:
  - POTATO
  - build
  - test
  - examples
  - python
md_extensions:
  - markdown.extensions.toc
  - markdown.extensions.tables
graph: true
warn: false
pages:
  doc/index.md: index
  doc/running.md: Running NEO-RT
  doc/file_formats.md: Output file formats
