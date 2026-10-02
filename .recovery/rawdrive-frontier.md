# Frozen raw-drive frontier

Controller review input, 2026-10-02. Mode: REPAIR.

Base: `8c429a5fb2f5182295843fb9e31b6b909ab89e92`.
Inherited analytical-S3 patch:
`2a8920ccf9d3b0255e58d14b586032433decdae30192a1df7318bf1b24d2e1df`.
Library override: `/Users/ert/proj/libneo-pertfield`, base
`0c2ea0ea4ee4433902e3d0618804118d1a4f8992`.

## Terminal claim and scope

For declared finite orbit samples in a smooth axisymmetric background, compare
the retained-order coordinate-free one-form drive with an independent oracle.
Use a complete unperturbed Littlejohn guiding-center trajectory, including its
radial excursion. Its velocity in the drive is exactly the velocity advancing
that trajectory. No rational or signed toroidal mode is removed.

Physical cylindrical coordinates are `(R, phi, Z)`, with orthonormal vector
components `(R, phi, Z)`. The source is a complex vector potential and physical
electrostatic potential, with convention `exp(+i*n*phi-i*omega*t)`. In Gaussian
units their dimensions are cm, G cm, statV; background magnetic field is G,
time s, magnetic moment erg/G, charge statC, mass g. SI uses m, T m, V, T,
J/T, C, kg and replaces every `q/c` coupling by `q`.

The retained-order PLAN contract, with NEO-RT's Hamiltonian sign, is
`H1 = mu*delta|B|_E + q*deltaPhi - (q/c)*deltaA.dot(Xdot)`.
Here `delta|B|_E = b0.dot(curl(deltaA))`. No Boozer weight, displaced-field
strength, or additional parallel-energy term is inserted.

This is a first-order guiding-center perturbation contract. Using a complete
background guiding-center orbit does not prove completeness at every order in
rho*. In particular, this packet does not silently infer that omitted changes
of the magnetic-direction symplectic term have vanished at all orders.

For a poloidal period `T`, `Omega_b=2*pi/T` and
`Omega_phi=(phi(T)-phi(0))/T`. Define `Delta=mb*Omega_b+n*Omega_phi-omega`.
The explicitly timed physical `H1` is projected by `exp(-i*Delta*t)`.
The field's explicit time phase therefore cancels in its canonical coefficient;
its physical frequency still enters resonance detuning and the paired scalar
potential. This is one convention, used consistently in every control.

## Established facts and first gap

Analytical S3 supplies useful nonresonant displacement and Boozer-collapse
checks. Its path streams at fixed s while its drive uses a radial drift, and
its gauge construction cannot eliminate rational parallel content. That
path/velocity inconsistency is the first bridge to repair for raw fields.

Complete background dynamics in Gaussian units are
`Bstar=B+(mass*c/q)*vpar*curl(b)`, `Bstarpar=b.dot(Bstar)`,
`F=mu*grad|B|+q*gradPhi0`,
`Xdot=(vpar*Bstar+(c/q)*b.cross(F))/Bstarpar`,
`vpar_dot=-Bstar.dot(F)/(mass*Bstarpar)`.
The domain excludes zero field, zero charge, nonpositive mass, singular
`Bstarpar`, or an orbit leaving the field's validated interpolation region.

## Falsifier and independent controls

For a prescribed nonzero pure gauge, `deltaA=grad(chi)` and
`deltaPhi=-(1/c)*partial_t(chi)`, the raw drive is
`-(q/c)*dchi/dt` along the same independently specified orbit. Its harmonic obeys
`DeltaH=-(q/c)*[chi(T)*exp(-i*Delta*T)-chi(0)]/T
         -i*(q/c)*Delta*chi_mb`.
Measure and retain the endpoint term; a declared closed periodic projected
chi must make it vanish. On resonance the remaining gauge harmonic is zero.
Off resonance, the detuning law is the required result.

Negative controls must detect: dropping radial velocity, mixing time/potential
conventions, suppressing rational parallel modes, and adding a parallel double
count. The independent closed-orbit gauge construction comes first. The S5
worker's independently derived circular fields and full-GC solver supply a
second oracle with actual radial width and finite resonant magnetic content.

## Minimal implementation boundary and nonclaims

Add a separate raw cylindrical one-form interface backed by the local
`perturbation_field_t` override, a complete GC integrator, and focused diagnostics.
Preserve analytical S3. Source evaluators must retain all signed modes and use
explicit physical potentials and frequencies; no scalar potential is invented
from an incomplete electric channel.

Do not insert raw harmonics into existing thin-orbit resonance roots and call
the resulting torque physical FOW transport. FOW frequencies, roots, canonical
actions, and the action-space measure must match before that product claim.
Fixed-unperturbed-orbit continuity in epsilon (often proportional to w squared)
does not establish an O(w) island-population sliver or the universal S5 limit.
Those remain separate obligations requiring perturbed-orbit/action-domain
evidence. Finite diagnostic contractions are labeled accordingly.

Heavy verification uses the controller's bounded scheduler queue. Workers do
not update authoritative status, commit, push, or promote these artifacts.
