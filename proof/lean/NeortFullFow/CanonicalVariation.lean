import NeortFullFow.Basic

namespace NeortFullFow

def stepDistance (left right : Int) : Nat :=
  Int.natAbs (right - left)

def canonicalVariation : List Int → Nat
  | [] => 0
  | [_] => 0
  | left :: right :: rest => stepDistance left right + canonicalVariation (right :: rest)

structure VariationCertificate where
  lower : Nat
  upper : Nat
  estimate : Nat
  lower_le_estimate : lower ≤ estimate
  estimate_le_upper : estimate ≤ upper
  derivative_coverage : Bool

def VariationCertificate.valid (c : VariationCertificate) : Prop :=
  c.derivative_coverage = true

theorem canonical_variation_is_nonnegative (samples : List Int) :
    0 ≤ canonicalVariation samples :=
  Nat.zero_le _

theorem certified_variation_is_enclosed (c : VariationCertificate)
    (h : c.valid) :
    c.valid ∧ c.lower ≤ c.estimate ∧ c.estimate ≤ c.upper := by
  exact ⟨h, c.lower_le_estimate, c.estimate_le_upper⟩

theorem canonical_variation_singleton (sample : Int) :
    canonicalVariation [sample] = 0 := by
  rfl

end NeortFullFow
