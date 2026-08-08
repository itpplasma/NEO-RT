import NeortFullFow.Basic

namespace NeortFullFow

structure ClassPartition where
  inboard : OrbitClass
  axis : OrbitClass
  outboard : OrbitClass
  seamInboardAxis : Nat
  seamAxisOutboard : Nat
  inboard_seam : inboard.span.hi = seamInboardAxis
  axis_left_seam : axis.span.lo = seamInboardAxis
  axis_seam : axis.span.hi = seamAxisOutboard
  outboard_left_seam : outboard.span.lo = seamAxisOutboard
  distinct_ids :
    inboard.id ≠ axis.id ∧ axis.id ≠ outboard.id ∧ inboard.id ≠ outboard.id

def ClassPartition.classes (p : ClassPartition) : List OrbitClass :=
  [p.inboard, p.axis, p.outboard]

theorem class_partition_has_three_classes (p : ClassPartition) :
    p.classes.length = 3 := by
  rfl

theorem class_partition_seams_are_shared (p : ClassPartition) :
    p.inboard.span.hi = p.axis.span.lo ∧
    p.axis.span.hi = p.outboard.span.lo := by
  exact ⟨p.inboard_seam.trans p.axis_left_seam.symm,
    p.axis_seam.trans p.outboard_left_seam.symm⟩

theorem class_partition_ids_are_not_collapsed (p : ClassPartition) :
    p.inboard.id ≠ p.axis.id ∧ p.axis.id ≠ p.outboard.id ∧
    p.inboard.id ≠ p.outboard.id :=
  p.distinct_ids

end NeortFullFow
