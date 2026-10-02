# S1 restart checkpoint, 2026-10-02

Stopped because the user requested commit/push and conclude for restart. No new jobs were submitted after that request. Worker did not commit or push; controller owns isolated WIP publication.

## Source publication

Worktree: `/Users/ert/proj/NEO-RT-eqdsk`
Branch: `feat/eqdsk-thin-orbit`
Base: `0b9472e50e0a960b3a83e432796a1242c34b3d52`
Remote: `git@github.com:itpplasma/NEO-RT.git`
Current complete publication patch: `publication-current.patch`
SHA256: `db570556ab8718a4a88a2207e980f599a59387f480c949d2525ec14cbb418527`
Exact owned paths and per-file hashes: `publication-current.json` (14 paths).
Force-stage only the four explicit `.bc` paths in that manifest. Those inputs are ignored by Git and are required for CI. No CMake/CTest fixture-generation step exists.

Original scientific candidate remains immutable: `candidate.patch`, SHA256 `b88c2b748121804b505a9372098e33f10c4a0c088fa16b70c35e9184adbb56ee`.
Separate complete fixture followup: `fixture-followup.patch`, SHA256 `0b8dbbaec2d076b9abdc865c4537377913faf4e645eb959b7368c668ba7bda92`.
The followup makes normal Boozer/perturbation fixtures explicit tracked CI inputs and expands normal/reversed Solovev GEQDSK vertical coverage from ±0.8 m to ±0.9 m. The generator changes only Solovev coverage. Normal expanded GEQDSK came from guarded job 21805567; reversed raw data use exact toroidal-field reflection.

The complete current patch contains both original science and fixture followup. The fixture followup has not yet had clean focused/full CI validation. Do not call the entire current patch validated based on the original candidate results alone. LIBNEO_REF has not been changed; original source still pins a4620f8. Published new-library checks used an explicit module override. Any required source dependency-pin update remains a controller decision after resumed clean CI.

## Completed independent validation

Remote root: `acluster:/home/ert/runs/claude-recovery-20261002/eqdsk-s1/`
Core job 21805532: `jobs/21805532/evidence/`; local copy `job-21805532-evidence/`.
- Original candidate focused six tests: exact named set PASS.
- Full 21 CTest names PASS, ten golden pytest cases PASS.
- Original implementation RED: actual fresh 11→9 SIGSEGV, reversed-angle-map winding failure, and absent stream-convention rejection.
- Four deliberate phase/sign mutations caught the actual transport behavioral assertion; baseline and restored baseline passed.
- Six default example `.out` files byte-for-byte identical to exact base 0b9472e.
- Bare `fo` rc1 only from unchanged baseline unused imports. All behavioral tests were separately completed. No claim of a clean bare-fo pipeline.
- Bounded kernel limits: actual4 CPU, hard8 GiB, swap0,2580 s. OOM/max events0. Slurm RSS measures the guard process and cannot establish total scope peak.
Exact evidence/named sets: `verification-summary.json`, core `result.json`, mutation `summary.json`, and raw LastTest logs.

Corrected plot-only job21805545: remote `jobs/21805545/`; local `job-21805545-plots/`.
Corrected plotting driver SHA256 `1f08060cb92ff824b247713b1482ae96af35cb3dc0aec97014d062d272ff7ad0`.
Four PNG/PDF figures, grayscale previews, frequency/harmonic/resonance CSVs, metadata and source are under `evidence/figures/`.
Mass/charge correction applied to constructed resonance loci and both positive quadratic branches retained. Color/grayscale rendering reviewed. The figures still describe original pre-coverage fixture data; they are selected samples, not a complete production resonance scan or converged transport profile.
At256 velocity steps: max difference/max Boozer magnitude is9.57e-5(co torque),6.03e-4(counter torque),2.72e-3(trapped torque),2.38e-4(totalD11),2.52e-4(totalD12). Weak trapped D12 at s=.35 differs about9%; total normalization cannot certify that component.

## Resolved newflux frequency failure

Published new geo source commit `5d760dff9aaf981f3c578afa51fa46f2f722e394`, exported module SHA256 `9aaccba64fd3144dc445a387d76cfcf3e6ef07d419d02f833fc9a641abf2d563`.
Independent analytic scalar oracle SHA256 `19d497143eb6a00c0491d3c4d2cffffab7333b2cfbc828af20f50453db3d75ad`.
Job21805556, original clipped fixture: local `flux-consumer/job-21805556/evidence/`, remote `flux-consumer/jobs/21805556/evidence/`.
The unchanged passing frequency bound failed at65/129/257 grids (1.394e-5/2.906e-5/2.220e-5). Actual-eta analytic period matched the returned poloidal label within6.8e-6, but the normalized toroidal label was wrong.
Independent diagnosis: declared LCFS reaches |Z|=.8080033927 m, beyond grid |Z|≤.8 m. Edge contour clamping corrupts total toroidal flux. Cache refinement did not fix this.
Job21805567, only vertical extent expanded to±.9 m: local `flux-consumer-expanded/job-21805567/evidence/`, remote `flux-consumer-expanded/jobs/21805567/evidence/`.
The unchanged complete backend test now PASSES all65/129/257 cases; passing errors5.244e-6/5.381e-6/5.422e-6, below unchanged1e-5.
Returned poloidal flux at s=.35 converges1.520352123e7→1.520351423e7→1.520351421e7 versus exact1.520351421e7. Independent period residual5.2–5.6e-6; independent quadrature converged5.9e-15.
These are cached-object consumer/link checks with exact module/API/compiler/cache provenance. They do not replace the separately completed whole libneo108-test build.
No library code or tolerance change was needed for the clipping repair.

## Staged unfinished work

`clean-fixture-check/`: NOT SUBMITTED. Exact clean archive + original science patch + fixture followup, without private fixture tar; four fresh test programs and exact six named behavioral tests using cached dependency objects plus published geo override. Frozen manifest SHA256 `680b81362daddccca1042310a7704d4bd6799b98e4e1ac51c77f4397916d2eee`. Proposed profile1 CPU/hard1 GiB/180s. Full clean dependency build and compatibility with unchanged old LIBNEO_REF remain unverified.
`transport-refinement/`: NOT SUBMITTED. Same original cached S1 executable/physics as256 comparison,512/1024 sweeps at all three surfaces, all class-resolved torque/D11/D12 residuals and consecutive-resolution differences. Frozen manifest SHA256 `a066d6957273ab22d0af6196d721290557a0db61b5c7bcc589f8f362b04aa158`. Proposed profile1 CPU/hard1 GiB/240s. It remains a pre-coverage control; refreshed wider-grid final artifacts are a subsequent task.

Resume commands only after aggregate slot release:
```sh
lane=clean-fixture-check  # then transport-refinement, sequentially
ssh acluster "mkdir -p /home/ert/runs/claude-recovery-20261002/eqdsk-s1/$lane"
rsync -ac /Users/ert/proj/recovery-20261002/eqdsk-s1/$lane/ acluster:/home/ert/runs/claude-recovery-20261002/eqdsk-s1/$lane/
ssh acluster "cd /home/ert/runs/claude-recovery-20261002/eqdsk-s1/$lane && sbatch job.sbatch"
```
Both staged payloads use bounded-run-v2, node-local scratch, checksum copy-out to an incoming directory and atomic rename on every exit. They do not delete unrelated data. Existing S1 jobs are all terminal; no S1 running payload depends on this agent staying alive.

Infrastructure attempts21805521 (literal-selector fix),21805524 (missing ignored fixtures),21805528 (BeeGFS Remote I/O) were preserved separately; they are not scientific failures. Consumer preflight attempts21805552/21805554 were preserved and ran no numerical test. Original candidate/reviewer logs remain referenced in the recovery task.

Adjacent unrepaired limits: reversed-Bt main Boozer frequency/q convention issue (separate from angle-map repair), trapped turning-event accuracy, counter-passing orbit oracle coverage, outer interpolation clamping and perturbation header validation. Source promotion and PLAN/status updates are controller-owned.
