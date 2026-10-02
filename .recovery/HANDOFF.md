S3 restart checkpoint, 2026-10-02. No worker commit/push or shared status update.
No S3 jobs are still running; successful final job 21805571 has copied results
independently to shared storage before allocation exit.

Source worktree: /Users/ert/proj/NEO-RT-line
Branch: feat/line-integral-drive
Exact base: 8c429a5fb2f5182295843fb9e31b6b909ab89e92
Patch SHA256: e77d75a29e5bd47ba0395c99dd543f6b117511d2a5ab28d2f299dc1d4e68b4b2
Source digest: 7621a3970998c6fd15b2b5176aa11a401e04201db8f1b3ad0046d1c716ec6282
Patch: .recovery/candidate.patch; exact owned paths in candidate-manifest.json.
Base includes the interrupted Claude's existing S3 local commits; root alone
owns promotion of the source and proof checkpoint. The .bc fixture is ignored
and needs explicit force-staging if root commits it. Preserve unrelated edits.

Owned source paths:
- src/line_drive.f90
- doc/running.md
- test/CMakeLists.txt
- test/test_line_drive.f90
- test/test_line_prestudy.f90
- test/gen_prestudy_circ.py
- test/diagnose_line_drive.f90
- test/line_drive/check_mutations.py
- test/line_drive/plot_diagnostics.py
- test/line_drive/run_diagnostics.sh
- test/line_drive/DIAGNOSTICS.md
- test/fixtures/prestudy/prestudy.geqdsk
- test/fixtures/prestudy/prestudy_boozer.bc
- test/line_drive/period_oracle.f90
- test/line_drive/prestudy_potential_oracle.f90

Owned proof artifacts:
/Users/ert/proj/neort-proofs/realspace/s3/historical-prestudy/
/Users/ert/proj/neort-proofs/realspace/s3/validation-e77d75a/
Final evidence payload SHA256:
df6a910a19e6198db6684c73163d82dd86578cb452150fbeb7b53d44809f859e
The 45-file payload is atomically checkpointed, with exact sources, raw CSVs,
PNG/PDF/grayscale figures, reviews, logs and per-file evidence-manifest hashes.
No PLAN, authoritative table, PR or shared integration state was modified.

Scientific checks:
20 registered CTests PASS (full evidence from job 21805560), then 11 golden/torque
pytest checks PASS, 3 independent physical mutations detected, 6 default scalar
outputs bitwise identical to independent exact S1 baseline 0b9472e.
Passing/trapped diagnostic harmonic errors: 1.470e-5 / 8.902e-10.
Independent algebra/interface review and final figure/data review PASS.
Bare fo attempted and returned nonzero for 3 unchanged unused imports:
POTATO/SRC/box_counting:get_stats; test_bounce_integrand:etadt;
test_chartmap_input:pi. Build and scientific tests passed; full log is retained.

Domain and gaps:
Positive-hθ Boozer backgrounds, analytical perturbation, fixed-s thin path and
canonical gauge projection were tested. Arbitrary raw dA/dPhi, rational modes,
radial finite-width transport, and negative-hθ passing frequency/target support
remain with separate lanes. The historical full-FOW prestudy is retained as
reported evidence, not claimed reproduced by the thin-path .078145 potential
ratio; corrected full-FOW reproduction remains with raw/fullGC extension.

Cluster provenance:
Host acluster.tugraz.at. Final jobs requested 8/used 4 CPU, exclusive allocated 64,
hard 8 GiB, swap 0, 45 min Slurm /2580 s v2 guard. All kernel OOM event counters zero.
Runner SHA256 49d12f9c6c954d5e48f6f1b088f3f97e5d71a41f59db3d64bbd5fe6185c3c262.
21805560 node31 /tmp/neort-s3-21805560: 20 CTests PASS, later A/B harness error.
21805566 node33 /tmp/neort-s3-21805566: 2 focused oracles PASS, A/B run-name error.
21805571 node33 same scratch: remaining checks/figures PASS, COMPLETED exit0.
21805564 cancelled before allocation while node31 was unavailable.
Earlier preserved failures: 21805525/31 node26 shared-home build IO; 21805533
literal-selector no-tests; 21805534 old Newton closure; 21805541 exact trace and
old empirical potential heuristic. These are not successful scientific checks.

Durable remote source+driver+dependency snapshot:
/home/ert/work/recovery-20261002/neort-s3
Durable remote job evidence:
/home/ert/work/recovery-20261002/neort-s3/.recovery/results/job-21805560
/home/ert/work/recovery-20261002/neort-s3/.recovery/results/job-21805571
Durable local private evidence: .recovery/results/job-21805560 and job-21805571.
Six dependency refs in .recovery/dependencies.json; compiler GCC 12.2,
CMake 3.25, RelWithDebInfo, FETCHCONTENT_FULLY_DISCONNECTED=ON.
Separate validation-only portability overlay SHA256
dd540918c5aa632bf871f08cfe34dc3677d403f3f640805de61e0a3ca6834f39
is applied only to the remote test snapshot, outside the S3 science patch.

Restart only, after the root releases the aggregate compute budget:
ssh acluster.tugraz.at
cd /home/ert/work/recovery-20261002/neort-s3
sbatch .recovery/s3-exclusive.sbatch
That driver creates fresh node-local source/dependencies/build, verifies the
frozen science digest and separate overlay, and runs the full bounded pipeline.
All process execution uses fo; neo_rt's literal run target needs argument
'driftorbit'. Required figures are in validation-e77d75a/diagnostics.

S5 child checkpoint (all realspace/s5 owned by the physics reviewer):
/Users/ert/proj/neort-proofs/realspace/s5/restart_handoff.md
SHA256 1ebe1fb210971d76535d3b79045f0fb1d67c910264a88bb10395d221312cc11f.
Final action preparation patch adcfd9f2acfcc4334ea96cb0c452074f213429c1d2310a287cb6779c6dab3ca6;
action source manifest 89bf6af0847c5fdc4cdd653561b3081c2d37cf9e29d3b7273452a5e05677694f.
Reduced S5 evidence-final SHA c5013eb0734dd7686504a9495a0d26b8f126289948e213bcb777841f7ffc7bce
remains valid for its narrower prerequisite. New action numerics remain unrun;
physical transport/action-sliver/island-population bounds remain open.
All children received stop/checkpoint; no new submissions/science/commits/pushes.
