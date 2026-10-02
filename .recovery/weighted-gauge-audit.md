# Weighted gauge correction

Review scope: the weighted exact-one-form argument in
`neort-proofs/monograph/chapters/body_coordinates.tex` only.
Exact base: `c4615729cb7f29586f652306c6a22be8f2d1acdc`.
The separate patch digest is regenerated after the controller's final
endpoint-condition clarification; it is recorded in the handoff message.

The first invalid implication replaced the harmonic-weighted integral of a
total derivative by its unweighted endpoint. The repair retains the finite
endpoint and detuning term. The monograph uses the one-form sign
`+(q/c) A.Xdot - q Phi - mu deltaB`; NEO-RT's raw Hamiltonian sign is its
negative, so its gauge change also has the opposite sign.

For `chi(t)=chi0 exp(i m Omega t)`, `Omega*T=2*pi`, integer `m != 0`, and
`omega=0`, the endpoints coincide but the weighted derivative coefficient is
`i*m*Omega*chi0`, nonzero. This counterexample is retained in the chapter.

Existing formal statements already supply the required identity:

- `NeortRealspace.harmonic` in `lean/NeortRealspace/Harmonic.lean` defines the
  exact weighted finite-window coefficient.
- `NeortRealspace.harmonic_deriv` in `lean/NeortRealspace/Gauge.lean` gives
  the derivative coefficient as endpoint plus `i*kappa*chi_hat`.
- `NeortRealspace.drive_gauge_shift` supplies the same change multiplied by
  the electromagnetic coupling.
- `NeortRealspace.harmonic_deriv_periodic` removes the endpoint under its
  declared periodicity and commensurability hypotheses.
- `NeortRealspace.drive_gauge_invariant_on_resonance` preserves the resonant
  suffix. No new Lean identity is required and no Lean source is changed.

Validation: source audit against these theorem statements and the exact
closed-orbit counterexample; whitespace check passes. The prose checker found
three pre-existing adjacent headings/phrases listed by the checker;
they are outside this narrow correction. No theorem build or document export
has been run locally. Controller review and cluster verification own promotion.
