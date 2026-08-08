import NeortFullFow.Basic

namespace NeortFullFow

def signedQuantity (orientation : Sign) (magnitude : Int) : Int :=
  match orientation with
  | .positive => magnitude
  | .negative => -magnitude

structure OrientationLedger where
  coordinate : Sign
  mode : Int
  torque : Int

def OrientationLedger.valid (ledger : OrientationLedger) : Prop :=
  ledger.torque = signedQuantity ledger.coordinate ledger.mode

theorem signed_quantity_under_reversal (orientation : Sign) (magnitude : Int) :
    signedQuantity orientation.flip magnitude = -signedQuantity orientation magnitude := by
  cases orientation <;> simp [Sign.flip, signedQuantity]

theorem orientation_ledger_preserves_declared_sign (ledger : OrientationLedger)
    (h : ledger.valid) :
    ledger.torque = signedQuantity ledger.coordinate ledger.mode :=
  h

end NeortFullFow
