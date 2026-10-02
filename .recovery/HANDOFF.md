# S2 restart checkpoint

The user requested commit/push and conclusion for restart. No new jobs or
scientific calculations are authorized. The root controller owns WIP branch
commits and pushes; this worker has not committed, pushed, or promoted results.

`restart-checkpoint.json` records current repositories, branches, remotes,
owned paths, and exact source overlay. `artifact-manifest.json` freezes all
prepared deliverables. Source base 0b9472e plus full overlay dc72f39d remains
unchanged by the final artifact-only native-helicity repair.

The latest repair adapts native BC phase m*theta_B−nper*ncol*phi_B to each code's
plus-n representation. The +2 scalar conjugates the native coefficients and
reverses m. Generated +2 NEO-RT input reverses m and bmns while retaining stored
ncol; generated -2 input retains m/bmns and negates stored ncol. Both preserve
the original real field and the ν toroidal shift. The two new nonzero-angle
native cosine/sine oracles are prepared and syntax checked, **not executed**.

Static energy-domain and convention reviews remain in `static-review.json`.
They do not certify the newly repaired native helicity numerically. The full
physical #107 conjugation prerequisite, composite dependency build, and all
energy/grid/orbit/radial-box/mass convergence checks remain pending. The six-case
torque/resonance scan has produced no data or figures.

The prepared Slurm wrappers use bounded-run-v2, actual-node proportional RAM
reservation, core binding, and a 120-second cleanup margin. Sources/dependencies
and builds are copied to node-local /tmp. Hashed build bundles and results go to
shared storage. Native POTATO binary reuse checks the CPU feature fingerprint.
No wrapper has been submitted by this worker.

For controller checkpointing, create owned WIP branches from the heads recorded
in `restart-checkpoint.json`. Stage only the listed source paths and the files
listed in `artifact-manifest.json`, plus that manifest itself. Preserve other
workers' edits. The two remotes are origin on GitHub for NEO-RT and origin on
TU Graz GitLab for neort-proofs. Existing source/fixture archives and source
overlays are part of the reproducible checkpoint.

After restart and a fresh controller resource release, set the exact paths in
a copy of `job-config.template.json` and run, in order:

```sh
sbatch tiny.sbatch "$PWD/job-config.json"
sbatch --time=00:18:00 job.sbatch prepare "$PWD/job-config.json"
```

The tiny job must execute the two literal Fortran oracles and four independent
Python phase/native-field/binary-layout oracles successfully. Preparation then
builds the composite pinned sources and constructs the native flux and scalar
fixtures. Establish independent nonzero physical representation-invariance
evidence before filling the prerequisite PASS gate and releasing torque runs.
The ordinary benchmark gate deliberately rejects pending evidence.

The separate perturbation-container evidence audit is complete in
`recovery-20261002/pertfield/review-tools/independent-flux-review-21805562.json`:
focused 1/1, suite 108/108, and bare-fo 108/108 named tests passed; thirty durable
file checksums verified. Its sole lint failure is unchanged baseline int64,
with sixteen new array-temporary warning emissions at ten source/test locations
for the source owner to repair. Frozen container evidence was not modified.
