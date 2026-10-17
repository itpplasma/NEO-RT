# Common-temperature multispecies Landau response

`conserving_fp_collisions=T` selects a joint response of every species in the
explicit SI mixture. It is exclusive with test-particle FP, Lorentz and Krook
collisions. Use `neort_compute_mixture_at_s`; a selected-species computation
cannot represent this operator. The CLI already dispatches a species input
file to the mixture entry point.

The collision form contains every species pair and its field-particle cross
blocks. The trial response combines the joined trapped/co/counter pitch graph
and speed elements with exact parallel-flow and energy functions. Pairwise
Landau gradient differences conserve each species number and the summed
parallel momentum and energy admitted by this common-phase space. This is a
finite projection of the linearized operator about Maxwellians at one common
local temperature. It does not assert a fully resolved arbitrary-geometry
Fokker–Planck solution or unequal-temperature equilibrium.

The retained native geometry is the same strict cosine well as the Lorentz
and test-particle implementations: constant safety factor with integer `n*q`,
linear thin orbits, a common phase, every orbit branch, and the same magnetic
drift model on trapped and passing particles. Oscillating drift phase,
finite orbit width, nonlinear transport, superbanana-only and trapped-only
switches are rejected. The native geometry guard remains authoritative.

Active physical profiles use exact SI mixture masses, charges and temperature
energies. Native velocity, electromagnetic ratios, source and work use that
same Maxwellian. Other modes retain their historical constants and output
when this switch is disabled.

The physical inner product is `J_SI * F0a * dtheta * d3v`, separately for each
species. Assembly retains this mass Gram instead of dividing rows by species
densities or thermal speeds. The joint equation is
`(C + D + i P) g = -i beta_a * H_a * (A1_a + A2_a*u_a**2)` projected in this
measure, where `D` is the retained native mirror/parallel generator and `P`
contains the native electric and local magnetic precession. The native
Hamiltonian includes the common spatial phase exactly once.

`C` is the Landau Dirichlet form: for each pair, integrate the transverse
projection of `grad(g_a)/m_a - grad(g_b)/m_b` against the two Maxwellians.
Center-of-mass and relative-velocity coordinates cancel the integrable
relative-speed singularity analytically; positive Gaussian quadrature
weights retain the Gram form. Self collisions retain their one-half pair
factor. No velocity softening, collision floor or deleted coincident pairs
are used. Both gyrophases are retained through their common and relative
azimuths. The parallel-flow and energy additions are exact analytic trial
functions, so their admitted common-temperature collision nulls survive
finite pair quadrature. General graph entries still require quadrature
refinement. This is the common-temperature linearized Landau structure
described by [Sugama et al. (2019)](https://arxiv.org/abs/1906.07427); the
native Hamiltonian ordering follows [Albert et al. (2016)](https://arxiv.org/abs/1607.04665).

The starter controls are deliberately small and require convergence before a
comparison scan:

| Input | Default | Meaning |
| --- | ---: | --- |
| `conserving_speed_steps` | 2 | speed elements |
| `conserving_pitch_steps` | 2 | elements on each joined pitch arm |
| `conserving_umax` | 4 | outer speed knot in thermal-speed units |
| `conserving_rate_scale` | 1 | common multiplier of the complete collision matrix |
| `conserving_theta_order` | 8 | native angle quadrature |
| `conserving_velocity_order` | 4 | native velocity projection quadrature |
| `conserving_cm_order` | 4 | pair center-of-mass quadrature |
| `conserving_rel_order` | 4 | relative-speed quadrature |
| `conserving_polar_order` | 4 | pair polar quadrature |
| `conserving_azimuth_order` | 8 | pair azimuth quadrature |
| `conserving_rank_rtol` | 1e-11 | relative singular-value threshold |
| `conserving_compatibility_rtol` | 1e-9 | null-space source tolerance |

The positive finite `conserving_rate_scale` multiplies every pair and species
block uniformly for collisionality scans, preserving the admitted invariants.
Species-dependent multipliers are not supported.

Mesh counts and Gaussian quadrature orders must be at least two; the pair
azimuth order must be at least four. The outer knot must be finite and positive,
and both tolerances must be finite and strictly
between zero and one. Beyond the outer speed knot graph values continue
constantly in speed; the analytic flow and energy functions retain their
physical tails. This is a basis definition, not a reflecting speed boundary.

The solver whitens the physical mass Gram and computes a singular-value
decomposition of the complete collision, streaming and precession matrix.
The reported relative rank threshold defines numerical resolution, with no
replacement of small singular values. A physical source along a left null
mode is rejected. A compatible source uses the minimum physical-mass norm
gauge. Increasing collisionality can require a tighter threshold to resolve
weakly detuned modes; inspect rank, compatibility and residual together with
mesh and quadrature convergence.

The result reports direct electromagnetic perturbation work. For each
receiving species, `conserving_work(3)` is the complex physical work integral
in trapped/co/counter order, with torque `-n*pi*Im(work)` in erg per normalized
radial coordinate. It uses the original receiving Hamiltonian. Optional
`conserving_cross_work(3,Nspecies)` columns resolve physical source species;
`conserving_cross_work_qualified` identifies columns whose individual source
is compatible with the operator null space. A compatible total source can
exist even when its separate driver columns are incompatible. Undefined
columns are placeholders and must not be interpreted as zero coupling.

HDF5 reverses Fortran array dimensions: cross-work datasets are driver by
receiving branch. Real and imaginary parts are stored separately. Rank,
nullity, singular cutoff, residual, compatibility, mass condition and
collision entropy accompany a successful native result. Only a successful
joint native path sets `conserving_active` and advertises
`source_feature=all-pair-linearized-Landau-commonT-v1`. The legacy D columns
are undefined placeholders and `susceptibility_valid` is false; use direct
work and torque. Disabled modes add no conserving metadata.

Rank, nullity, cutoff, mass condition, compatibility, residual and the
number-weighted collision entropy form describe the joint all-species
operator and are repeated for each receiving species. They are not separate
species conservation diagnostics. The singular cutoff has units of inverse
seconds; the entropy form counts the Maxwellian particle weight per second
without a Boltzmann-constant factor.

The old harmonic velocity-range columns are also undefined zero placeholders:
the physical velocity domain extends to infinity, and `conserving_umax` is a
mesh knot. Active metadata states this tail policy explicitly.
