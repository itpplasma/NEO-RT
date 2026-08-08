import NeortFullFow.ClassPartition

namespace NeortFullFow

def ownsLeftClosed (i : Interval) (x : Nat) : Prop :=
  i.lo ≤ x ∧ x < i.hi

def ownsClosedRight (i : Interval) (x : Nat) : Prop :=
  i.lo ≤ x ∧ x ≤ i.hi

theorem first_seam_is_owned_by_axis (p : ClassPartition) :
    ownsLeftClosed p.axis.span p.seamInboardAxis := by
  constructor
  · simp [p.axis_left_seam]
  · simpa [p.axis_left_seam] using p.axis.span.valid

theorem first_seam_is_not_owned_by_inboard (p : ClassPartition) :
    ¬ ownsLeftClosed p.inboard.span p.seamInboardAxis := by
  intro h
  exact Nat.lt_irrefl p.seamInboardAxis (p.inboard_seam ▸ h.2)

theorem second_seam_is_owned_by_outboard (p : ClassPartition) :
    ownsClosedRight p.outboard.span p.seamAxisOutboard := by
  constructor
  · simp [p.outboard_left_seam]
  · simpa [p.outboard_left_seam] using Nat.le_of_lt p.outboard.span.valid

theorem second_seam_is_not_owned_by_axis (p : ClassPartition) :
    ¬ ownsLeftClosed p.axis.span p.seamAxisOutboard := by
  intro h
  exact Nat.lt_irrefl p.seamAxisOutboard (p.axis_seam ▸ h.2)

theorem boundary_policy_keeps_final_endpoint (p : ClassPartition) :
    ownsClosedRight p.outboard.span p.outboard.span.hi := by
  exact ⟨Nat.le_of_lt p.outboard.span.valid, Nat.le_refl _⟩

end NeortFullFow
