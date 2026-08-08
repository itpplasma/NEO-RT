import Std

namespace NeortFullFow

inductive Sign where
  | positive
  | negative
  deriving DecidableEq, Repr

def Sign.flip : Sign → Sign
  | .positive => .negative
  | .negative => .positive

theorem Sign.flip_flip (s : Sign) : s.flip.flip = s := by
  cases s <;> rfl

structure Interval where
  lo : Nat
  hi : Nat
  valid : lo < hi

structure OrbitClass where
  id : Nat
  span : Interval
  canonical : Int

inductive Boundary where
  | leftEndpoint
  | seam (id : Nat)
  | rightEndpoint
  deriving DecidableEq, Repr

end NeortFullFow
