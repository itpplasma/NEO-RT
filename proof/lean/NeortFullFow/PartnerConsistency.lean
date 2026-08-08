import NeortFullFow.Basic

namespace NeortFullFow

structure PartnerBoundary where
  id : Nat
  canonical : Int
  orientation : Sign

structure PartnerPair where
  first : PartnerBoundary
  second : PartnerBoundary

def PartnerPair.consistent (p : PartnerPair) : Prop :=
  p.first.canonical = p.second.canonical ∧ p.first.id ≠ p.second.id

def makePartner (firstId secondId : Nat) (canonical : Int)
    (firstOrientation secondOrientation : Sign) : PartnerPair :=
  { first := { id := firstId, canonical := canonical, orientation := firstOrientation }
    second := { id := secondId, canonical := canonical, orientation := secondOrientation } }

theorem makePartner_consistent {firstId secondId : Nat} {canonical : Int}
    {firstOrientation secondOrientation : Sign} (different : firstId ≠ secondId) :
    (makePartner firstId secondId canonical firstOrientation secondOrientation).consistent := by
  exact ⟨rfl, different⟩

theorem consistent_implies_same_canonical (p : PartnerPair)
    (h : p.consistent) : p.first.canonical = p.second.canonical :=
  h.1

theorem consistent_implies_distinct_boundaries (p : PartnerPair)
    (h : p.consistent) : p.first.id ≠ p.second.id :=
  h.2

end NeortFullFow
