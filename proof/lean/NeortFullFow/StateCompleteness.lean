import NeortFullFow.Basic

namespace NeortFullFow

structure RuntimeState where
  backend : String
  coordinates : String
  model : Nat
  frequencyModel : Nat
  wallHash : String
  runtimeExecutionComplete : Bool
  orbitBackendCertified : Bool
  wallCertified : Bool
  canonicalMeasureCertified : Bool
  componentIdentityCertified : Bool
  perturbationProvenanceCertified : Bool
  quadratureConvergenceCertified : Bool
  harmonicBatchCertified : Bool
  classReconstructionCertified : Bool
  orbitStepRefinementCertified : Bool

def RuntimeState.complete (s : RuntimeState) : Prop :=
  s.backend = "cylindrical_full_fow" ∧
  s.coordinates = "R,Z,phi" ∧
  s.model = 2 ∧
  s.frequencyModel = 2 ∧
  s.wallHash.length = 64 ∧
  s.runtimeExecutionComplete = true ∧
  s.orbitBackendCertified = true ∧
  s.wallCertified = true ∧
  s.canonicalMeasureCertified = true ∧
  s.componentIdentityCertified = true ∧
  s.perturbationProvenanceCertified = true ∧
  s.quadratureConvergenceCertified = true ∧
  s.harmonicBatchCertified = true ∧
  s.classReconstructionCertified = true ∧
  s.orbitStepRefinementCertified = true

theorem complete_state_uses_full_fow_backend (s : RuntimeState)
    (h : s.complete) : s.backend = "cylindrical_full_fow" :=
  h.1

theorem complete_state_has_runtime_execution (s : RuntimeState)
    (h : s.complete) : s.runtimeExecutionComplete = true := by
  rcases h with ⟨_, _, _, _, _, runtime, _, _, _, _, _, _, _, _, _⟩
  exact runtime

theorem complete_state_has_class_and_measure_certificates (s : RuntimeState)
    (h : s.complete) :
    s.canonicalMeasureCertified = true ∧ s.classReconstructionCertified = true := by
  rcases h with ⟨_, _, _, _, _, _, _, _, measure, _, _, _, _, classes, _⟩
  exact ⟨measure, classes⟩

theorem incomplete_launcher_state_is_not_complete (s : RuntimeState)
    (h : s.runtimeExecutionComplete = false) : ¬ s.complete := by
  intro complete
  rcases complete with ⟨_, _, _, _, _, runtime, _, _, _, _, _, _, _, _, _⟩
  exact Bool.noConfusion (runtime.symm.trans h)

end NeortFullFow
