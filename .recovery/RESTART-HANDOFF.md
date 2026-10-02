# Raw-drive and weighted-gauge restart handoff

User requested stop, publish WIP through the domain controller, and conclude.
No worker commit/push, new science, or new submission after that request.
Last submitted job 21805574 completed independently; no watcher is needed.

## Publication boundary and reconstruction

Worktree: `/Users/ert/proj/NEO-RT-rawdrive`, detached exact base
`8c429a5fb2f5182295843fb9e31b6b909ab89e92`.

Inherited original S3 candidate patch SHA256:
`2a8920ccf9d3b0255e58d14b586032433decdae30192a1df7318bf1b24d2e1df`.
It was applied before raw work. This worktree does not contain S3's later
period/prestudy fixes. Do not promote its old analytical S3 copies over the
S3 controller's separately reviewed latest candidate.

The owned-only patch is `.recovery/raw-owned.patch`, SHA256
`643f7cb0bddf3298f7861cbbd97644f794642064cbbf4aa51a02f9f8e8dd31d5`.
It applies against the exact base using a disposable index check. It includes
only the raw feature additions in the two CMake files, excluding inherited
S3 diagnostic/test additions. Publish it in a clean controller checkout, or
stage precisely these hunks; staging the entire current test/CMakeLists.txt
would also stage inherited S3 changes whose sources belong to another worker.

Owned source paths are the two CMake files, src/raw_drive.f90,
src/raw_orbit.f90, and all five files under test/raw_drive/. Exact source
hashes, the 21 complete reconstruction paths and complete-patch digest are in
`.recovery/rawdrive-manifest.json`. `.recovery/worktree-complete.patch`
reproduces this whole candidate including its older inherited S3 subset;
use that for forensic reconstruction, not selective promotion.

`.recovery/warning-cleanup.patch` SHA256:
`e934aa5649bc4a4f02fa154b39f879b63adcd0eba9ec5f011deb0b32c41e9b7c`.
It records the delta from job 21805565 to final job 21805574, including an
explicit full-orbit resonant pure-gauge zero assertion. No tolerance changed.

## Implemented scope and exact tested sources

Experimental NEORT_ENABLE_RAW_DRIVE defaults OFF. A controller-provided local
libneo container override is needed when enabled; no dependency was pushed.
raw_source_t loads/evaluates physical complex cylindrical A and paired Phi,
keeps arbitrary signed n including zero, and uses source exp(+inphi-iomega t).
SI H is mu b.curl(A)+qPhi-q A.Xdot in joules; Gaussian replaces q by q/c in
the vector coupling and yields erg. Projection uses
exp[-i(mb Omega_b+n Omega_phi-omega)t], with source time phase canceling.
Absent Phi is an explicitly declared zero channel, not inferred from E.

raw_orbit integrates the complete radial Littlejohn GC velocity through
fortnum VODE. The same velocity feeds the Hamiltonian. It computes E,
canonical Pphi, reduced meridional action and local symplectic density.
Actions fail closed without matching background A0 or endpoint closure.
Primitive period and passing orientation remain caller obligations; the
oracle supplies an independently solved primitive trapped period.

Final tested production SHA256:
- src/raw_drive.f90:
  `8ccfac8fec28f6bf21f3ed44dfd52aa99f203a1250be1bd9ae97322bd154a52a`
- src/raw_orbit.f90:
  `1884538f62515227645186fd2298ae16e5bc373057da5eab072cb56845a01e67`
- test/raw_drive/test_raw_oneform.f90:
  `2460e8005cc758b21fd6473332b465cd38148fa65ce9be619451e27727fcd464`
- test/raw_drive/test_raw_gc.f90:
  `3f5c157e1705fa2fdf072daf43c37c234c12e354843e7150a3579e1f814b02e8`

Reference provenance is test/raw_drive/raw_gc_reference.json. It transforms
the existing independent S5 1025-point CSV into a text fixture, with no new
orbit computation. The physical SI unit scales are explicit in raw_fixture.

## Validation, resources and durable results

Job 21805574: acluster / node26, COMPLETED 0:0, 4 seconds, one reserved and
actual CPU, hard 1 GiB, Slurm 5 minutes, immutable v2 guard 180 seconds.
Guard SHA256: 49d12f9c6c954d5e48f6f1b088f3f97e5d71a41f59db3d64bbd5fe6185c3c262.
Kernel memory.max=1073741824, memory.swap.max=0, OOM group=1,
cpu.max=100000/100000; no OOM events. Build log is empty: no warnings.
Compilation/testing used fo, not direct binary execution.

Remote stage:
`/home/ert/runs/claude-recovery-20261002/rawdrive/tiny-fullgc-v2`.
Input manifest SHA256:
`a89a6d778f26c1b4a9100a4643e2f42b4cd56813d16e91a69273968f6d6a051e`.
Runner packet SHA256:
`ca1b95ba506d153086eec4a0d09b1a7576655cc450a9854aaf47da1cdfa17fe0`.
Exact accepted 0c2ea0e container consumer archive SHA256:
`7d5d9c58705aa75a46406f3ab534fe24e4dd50ec16e36d934d188c4dd7000511`.
Its producer job was 21805562; source scratch was already retrieved/deleted
by its controller. The raw verifier uses only copied read-only archives.

Private scratch `/tmp/rawdrive-21805574.56mfk6` was copied and checksum-verified
by its independent exit trap, then deleted. Failure preserves scratch.
Durable remote results: stage/results/21805574. Local results:
`.recovery/results-21805574`; all 20 retrieval checksums verified locally.
Slurm log and runner-input manifest also copied locally.

Independent full-GC errors: velocity 2.22e-16; orbit 2.60e-10;
E 7.84e-12; Pphi 4.59e-10; reduced action 1.16e-10;
retained rational magnetic amplitude 1.53e-11; raw physical energy 1.42e-11.
Actual-orbit gauge endpoint+detuning errors <=5.2e-15; resonant pure-gauge
zero 3.39e-11. Orbit radius spans 0.235516 to 0.243559 m.
Nonzero mu*curl(A), SI/Gaussian energy and signed-n phase controls pass.
Full project CMake integration/full fo pipeline have not been run.

Earlier jobs: 21805538 reader argument-order compile failure; 21805543 finite
gauge PASS; 21805565 first complete GC PASS with temporary warnings;
21805573 staging failure before compilation, safely retrieved; 21805574
final warning-clean PASS. No work ran on protected faepop/faepcr hosts.

## Reviewed and unfinished frontiers

Root reviewed .recovery/rawdrive-frontier.md SHA256
`57ac65bff4c1143b273c111578b741554d09bd2e26ff3e8f06b449516e729cb8`:
PASS within retained-order finite-orbit scope. S5 physics_review independently
confirmed source units/phases and controls, then requested the added nonzero
magnetic-channel oracle. It independently confirmed canonical action/volume
with the passing-orientation correction now encoded in raw_gc_actions.
.recovery/canonical-bridge.md records the exact measure and forbidden claims.

Next production gap: discover a primitive GC return period rather than accept
one from the caller. The reviewed proposal uses certified departure from a
non-tangent Z section, first same-direction return with R/Z/vpar closure,
finite time/event caps, and independently integrated geometric poloidal
winding around an explicit magnetic axis. Accept winding 0 or +/-1 with
physical tolerances; reject axis crossings, tangency and unsupported winding.
Do not demand phi closure. This proposal is reviewed but not implemented.

S5 worker owns unexecuted realspace/s5/action_oracle.py: canonical fixed-E/P
shooting, action derivatives and 6x6 symplectic determinant. Requires a new
controller queue release. Full physical transport needs actual FOW roots,
canonical volume, physical distribution normalization and a consistent
collision operator. A bare-H Lorentzian is a fixed-gauge diagnostic only.
The thin-island population/mask and widths vs radial/resonance/collision
scales remain unproved; a uniform limit may be false.

S4 worker owns the real matched VMEC ideal vector/displacement bridge;
scalar collapse alone is insufficient. S6 real KiLCA magnetic vacuum exports
exist, but certified paired A/Phi export is still withheld after the strict
electric-field gate failed. Do not invent Phi or claim physical torque/Dij.
Arbitrary-field higher-order delta-b symplectic completeness remains separate.
Historical full-FOW prestudy reproduction is also unfinished; S3's corrected
thin prestudy result cannot close it.

## Narrow monograph repair

Only owned prose file: neort-proofs/monograph/chapters/body_coordinates.tex.
Base c4615729cb7f29586f652306c6a22be8f2d1acdc. Patch
.recovery/weighted-gauge-doc.patch SHA256
`2038813afaf9fe0cca775cc89082ad540fcdecacd8fa5eaf7937ec8855bb8b99`.
Current exact file diff still matches that digest. Root independently reviewed
the endpoint+detuning repair and required the explicit endpoint condition for
resonant invariance; that wording is present. Durable off-resonance Fourier
counterexample and mapping to existing Lean Gauge/Harmonic theorems are in
.recovery/weighted-gauge-audit.md. No Lean edits/build claimed.

## Exact next commands

From a clean controller checkout at the recorded base:
`git apply --check /Users/ert/proj/NEO-RT-rawdrive/.recovery/raw-owned.patch`.
Publish only after reviewing the owned paths and WIP scope. Do not overwrite
the independently advanced S3 branch or authoritative PLAN/status files.

Read-only status:
`ssh acluster 'sacct -j 21805574 --format=JobID,State,ExitCode,Elapsed,NodeList -P'`.
After a new controller slot release, replay the frozen minimal verifier with
`ssh acluster 'set -e; cd /home/ert/runs/claude-recovery-20261002/rawdrive/tiny-fullgc-v2; sha256sum --check runner-inputs.sha256; sha256sum --check inputs.sha256; sbatch --parsable job.sbatch'`.
For full integration, stage exact local libneo/fortio overrides and use fo with
NEORT_ENABLE_RAW_DRIVE=ON in a new guarded isolated build. Do not build locally
or mutate the cached S3/container source. No resume is automatic.
