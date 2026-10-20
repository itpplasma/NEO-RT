project: NEO-RT
summary: Neoclassical toroidal viscosity solver for resonant transport regimes
author: NEO-RT Developers
output_dir: build/doc
src_dir: src
exclude: src/do_magfie_neo.f90
page_dir: doc
md_extensions: markdown.extensions.toc
    markdown.extensions.tables
graph: true
parallel: 8
warn: false

API documentation uses the default standalone field backend, matching
`USE_STANDALONE=ON` in CMake. The alternative NEO-2 backend defines the same
module names and is excluded from this source set. FORD searches `src`
recursively, including diagnostics, without collecting tests or build dependencies.
